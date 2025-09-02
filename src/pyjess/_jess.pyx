# coding: utf-8
# cython: language_level=3, linetrace=True, binding=True
"""Bindings to Jess, a 3D template matching software.

References:
    - Barker, J. A., & Thornton, J. M. (2003). *An algorithm for
      constraint-based structural template matching: application to
      3D templates with statistical analysis*. Bioinformatics (Oxford,
      England), 19(13), 1644–1649. :doi:`10.1093/bioinformatics/btg226`.
    - Wallace, A. C., Borkakoti, N., & Thornton, J. M. (1997).
      *TESS: a geometric hashing algorithm for deriving 3D coordinate
      templates for searching structural databases. Application to enzyme
      active sites*. Protein science : a publication of the Protein
      Society, 6(11), 2308–2323. :doi:`10.1002/pro.5560061104`.

"""

# --- C imports --------------------------------------------------------------

cimport cython
from cpython.unicode cimport PyUnicode_FromStringAndSize

from libc.math cimport isnan, exp, INFINITY, NAN
from libc.stdio cimport FILE, fclose, fdopen, printf
from libc.stdint cimport uintptr_t
from libc.stdlib cimport calloc, realloc, free, malloc
from libc.string cimport memcpy, memset, strncpy, strdup

cimport jess.atom
cimport jess.jess
cimport jess.molecule
cimport jess.super
cimport jess.tess_template
cimport jess.tess_atom
cimport jess.res_index
from jess.atom cimport Atom as _Atom
from jess.jess cimport Jess as _Jess
from jess.jess cimport JessQuery as _JessQuery
from jess.molecule cimport Molecule as _Molecule
from jess.super cimport Superposition as _Superposition
from jess.template cimport Template as _Template, IgnoreType as _IgnoreType
from jess.tess_template cimport TessTemplate as _TessTemplate
from jess.tess_atom cimport TessAtom as _TessAtom

# --- Python imports ---------------------------------------------------------

import contextlib
import functools
import io
import itertools
import os
import warnings

__version__ = PROJECT_VERSION

# --- Utils ------------------------------------------------------------------

cdef inline void copy_token(char* dst, const char* src, size_t n) noexcept nogil:
    cdef size_t i
    for i in range(n):
        if src[i] == ord(' ') or src[i] == 0:
            dst[i] = ord('_')
        else:
            dst[i] = src[i]
    dst[n] = 0

@contextlib.contextmanager
def nullcontext(return_value=None):
    yield return_value

# --- Classes ----------------------------------------------------------------

cdef class Molecule:
    """A molecule structure, as a sequence of `Atom` objects.

    .. versionadded:: 0.2.2
       Support identifiers of arbitrary length.

    .. versionadded:: 0.4.0
       Equality, hashing and pickle protocol support.

    """
    cdef _Molecule* _mol
    cdef str        _id

    @classmethod
    def loads(cls, text, str id = None, bint ignore_endmdl = False):
        """Load a molecule from a PDB string.

        Arguments:
            file (`str`, `os.PathLike`, or file-like object): Either the path
                to a file, or a file-like object opened in **text mode**
                containing a PDB molecule.
            id (`str`, optional): The identifier of the molecule. If `None`
                given, the parser will attempt to extract it from the
                ``HEADER`` line.
            ignore_endmdl (`bool`): Pass `True` to make the parser read all
                the atoms from the PDB file. By default, the parser only
                reads the atoms of the first model, and stops at the first
                ``ENDMDL`` line.

        Returns:
            `~pyjess.Molecule`: The molecule parsed from the PDB file.

        See Also:
            `Molecule.load` to load a PDB molecule from a file-like
            object or from a path.

        """
        return cls.load(io.StringIO(text), id=id, ignore_endmdl=ignore_endmdl)

    @classmethod
    def load(cls, file, str id = None, bint ignore_endmdl = False):
        """Load a molecule from a PDB file.

        Arguments:
            file (`str`, `os.PathLike`, or file-like object): Either the path
                to a file, or a file-like object opened in **text mode**
                containing a PDB molecule.
            id (`str`, optional): The identifier of the molecule. If `None`
                given, the parser will attempt to extract it from the
                ``HEADER`` line.
            ignore_endmdl (`bool`): Pass `True` to make the parser read all
                the atoms from the PDB file. By default, the parser only
                reads the atoms of the first model, and stops at the first
                ``ENDMDL`` line.

        Returns:
            `~pyjess.Molecule`: The molecule parsed from the PDB file.

        """
        try:
            handle = open(file)
        except TypeError:
            handle = nullcontext(file)
        with handle as f:
            atoms = []
            for line in f:
                if line.startswith("HEADER"):
                    if id is None:
                        id = line[62:66].strip() or None
                elif line.startswith(("ATOM", "HETATM")):
                    atoms.append(Atom.loads(line))
                elif line.startswith("ENDMDL"):
                    if not ignore_endmdl:
                        break
        return cls(atoms, id=id)

    def __cinit__(self):
        self._mol = NULL

    def __dealloc__(self):
        jess.molecule.Molecule_free(self._mol)

    def __init__(self, object atoms = (), str id = None):
        """__init__(self, atoms=(), id=None)\n--\n

        Create a new molecule.

        Arguments:
            atoms (sequence of `~pyjess.Atom`): The atoms of the molecule.
            id (`str`, optional): The identifier of the molecule.

        Raises:
            `MemoryError`: When the system allocator fails to allocate
                enough memory for the molecule storage.

        """
        cdef Atom atom
        cdef int i
        cdef int count = len(atoms)

        self._mol = <_Molecule*> malloc(sizeof(_Molecule) + count * sizeof(_Atom*))
        if self._mol is NULL:
            raise MemoryError("Failed to allocate molecule")

        self._mol.index = NULL
        self._mol.count = count
        for i in range(count):
            self._mol.atom[i] = NULL
        memset(self._mol.id, b' ', 5)
        self._id = id

        for i, atom in enumerate(atoms):
            self._mol.atom[i] = <_Atom*> malloc(sizeof(_Atom))
            if self._mol.atom[i] is NULL:
                raise MemoryError("Failed to allocate atom")
            memcpy(self._mol.atom[i], atom._atom, sizeof(_Atom))

        self._mol.index = jess.res_index.ResIndex_create(self._mol.atom, count)
        if self._mol.index is NULL:
            raise MemoryError("Failed to allocate residue index")

    def __len__(self):
        assert self._mol is not NULL
        return self._mol.count

    def __getitem__(self, object index):
        assert self._mol is not NULL

        cdef Atom    atom
        cdef ssize_t index_
        cdef ssize_t length = self._mol.count

        if isinstance(index, slice):
            indices = range(*index.indices(length))
            return type(self)(atoms=[self[i] for i in indices], id=self.id)
        else:
            index_ = index
            if index_ < 0:
                index_ += length
            if index_ < 0 or index_ >= length:
                raise IndexError(index)
            atom = Atom.__new__(Atom)
            atom.owner = self
            atom.owned = True
            atom._atom = <_Atom*> jess.molecule.Molecule_atom(self._mol, index_)
            return atom

    def __copy__(self):
        return self.copy()

    def __eq__(self, object other):
        cdef Molecule other_
        if not isinstance(other, Molecule):
            return NotImplemented
        other_ = other
        if self._id != other_._id:
            return False
        if self._mol.count != other_._mol.count:
            return False
        return all(x == y for x,y in zip(self, other_))

    def __hash__(self):
        return hash((self._id, *(hash(x) for x in self)))

    def __reduce__(self):
        return type(self), (list(self), self.id)

    def __sizeof__(self):
        assert self._mol is not NULL
        return (
            sizeof(self)
            + sizeof(_Molecule)
            + self._mol.count*(sizeof(_Atom*) + sizeof(_Atom))
        )

    @property
    def id(self):
        return self._id

    cpdef Molecule conserved(self, double cutoff = 0.0):
        assert self._mol is not NULL
        cdef Atom atom
        return type(self)(
            id=self.id,
            atoms=[
                atom
                for atom in self
                if cutoff <= 0.0
                or atom._atom.tempFactor >= cutoff
            ]
        )

    cpdef Molecule copy(self):
        """Create a copy of this molecule and its atoms.

        Returns:
            `~pyjess.Molecule`: A newly allocated molecule with the same
            identifier and atoms.

        .. versionadded:: 0.4.0

        """
        cdef Molecule copy = Molecule.__new__(Molecule)
        cdef size_t   size = sizeof(_Molecule) + self._mol.count * sizeof(_Atom*)

        with nogil:
            # allocate molecule storage
            copy._mol = <_Molecule*> malloc(size)
            if copy._mol is NULL:
                raise MemoryError("Failed to allocate molecule")
            # copy molecule attributes
            copy._mol.index = NULL
            copy._mol.count = self._mol.count
            memset(copy._mol.id, b' ', 5)
            # copy molecule atoms
            for i in range(self._mol.count):
                copy._mol.atom[i] = <_Atom*> malloc(sizeof(_Atom))
                if copy._mol.atom[i] is NULL:
                    raise MemoryError("Failed to allocate atom")
                memcpy(copy._mol.atom[i], self._mol.atom[i], sizeof(_Atom))
            # regenerate index
            copy._mol.index = jess.res_index.ResIndex_create(copy._mol.atom, copy._mol.count)
            if copy._mol.index is NULL:
                raise MemoryError("Failed to allocate residue index")

        copy._id = self._id
        return copy


cdef class Atom:
    """A single atom in a molecule.

    .. versionadded:: 0.4.0
       Equality, hashing and pickle protocol support.

    """
    cdef object owner
    cdef bint   owned
    cdef _Atom* _atom

    @classmethod
    def load(cls, file):
        """Load an atom from the given file.

        Arguments:
            file (file-like object): A file-like object opened in text
                mode to read the atom from.

        """
        return cls.loads(file.read())

    @classmethod
    def loads(cls, text):
        """Load an atom from the given string.

        Arguments:
            text (`str`, `bytes` or `bytearray`): The atom line to read the
                atom metadata from.

        """
        cdef const unsigned char* s
        cdef bytearray            b
        cdef Atom                 atom

        if isinstance(text, str):
            b = bytearray(text, 'utf-8')
        else:
            b = bytearray(text)
        if not b.endswith(b'\n'):
            b.append(b'\n')
        b.append(b'\0')
        s = b

        atom = cls.__new__(cls)
        with nogil:
            atom._atom = <_Atom*> malloc(sizeof(_Atom))
            if atom._atom == NULL:
                raise MemoryError("Failed to allocate atom")
            if not jess.atom.Atom_parse(atom._atom, <const char*> s):
                raise ValueError(f"Failed to parse atom: {text!r}")

        return atom

    def __cinit__(self):
        self._atom = NULL
        self.owner = None
        self.owned = False

    def __dealloc__(self):
        if not self.owned:
            free(self._atom)

    def __init__(
        self,
        *,
        int serial,
        str name,
        str altloc,
        str residue_name,
        str chain_id,
        int residue_number,
        str insertion_code,
        double x,
        double y,
        double z,
        double occupancy = 0.0,
        double temperature_factor = 0.0,
        str segment = '',
        str element = '',
        int charge = 0,
    ):
        """__init__(self, *, serial, name, altloc, residue_name, chain_id, residue_number, insertion_code, x, y, z, occupancy=0.0, temperature_factor=0.0, segment='', element='', charge=0)\n--\n

        Create a new atom.

        Raises:
            `MemoryError`: When the system allocator fails to allocate
                enough memory for the atom storage.
            `ValueError`: When either of the ``name``, ``residue_name``,
                ``segment``, ``element`` or ``chain_id`` strings is too
                long.

        """
        if len(name) > 4:
            raise ValueError(f"Invalid atom name: {name!r}")
        if len(residue_name) > 3:
            raise ValueError(f"Invalid residue name: {residue_name!r}")
        if len(segment) > 3:
            raise ValueError(f"Invalid segment: {segment!r}")
        if len(element) > 2:
            raise ValueError(f"Invalid element: {element!r}")
        if len(chain_id) > 2:
            raise ValueError(f"Invalid chain ID: {chain_id!r}")

        self._atom = <_Atom*> malloc(sizeof(_Atom))
        if self._atom is NULL:
            raise MemoryError("Failed to allocate atom")

        self._atom.serial = serial
        self._atom.altLoc = ord(altloc)
        self._atom.chainID1 = ord(chain_id[0]) if len(chain_id) > 0 else 0
        self._atom.chainID2 = ord(chain_id[1]) if len(chain_id) > 1 else ord(' ')
        self._atom.resSeq = residue_number
        self._atom.iCode = ord(insertion_code)
        self._atom.x[0] = x
        self._atom.x[1] = y
        self._atom.x[2] = z
        self._atom.occupancy = occupancy
        self._atom.tempFactor = temperature_factor
        self._atom.charge = charge
        copy_token(self._atom.resName, residue_name.encode('ascii').ljust(3, b'\0'), 3)
        copy_token(self._atom.segID, segment.encode('ascii').ljust(3, b'\0'), 3)
        copy_token(self._atom.element, element.encode('ascii').ljust(2, b'\0'), 2)

        _name = bytearray(name, 'ascii')
        if len(_name) < 4:
            _name.insert(0, ord('_'))
        copy_token(self._atom.name, _name.ljust(4, b'\0'), 4)

    def __copy__(self):
        return self.copy()

    cdef dict _state(self):
        return {
            "serial": self.serial,
            "name": self.name,
            "altloc": self.altloc,
            "residue_name": self.residue_name,
            "chain_id": self.chain_id,
            "residue_number": self.residue_number,
            "insertion_code": self.insertion_code,
            "x": self.x,
            "y": self.y,
            "z": self.z,
            "temperature_factor": self.temperature_factor,
            "occupancy": self.occupancy,
            "segment": self.segment,
            "element": self.element,
            "charge": self.charge,
        }

    def __reduce__(self):
        cdef dict state = self._state()
        return functools.partial(type(self), **state), ()

    def __repr__(self):
        cdef str  ty   = type(self).__name__
        cdef list args = []
        for k,v in self._state().items():
            if v is not None:
                args.append(f"{k}={v!r}")
        return f"{ty}({', '.join(args)})"

    def __sizeof__(self):
        cdef size_t size = sizeof(self)
        if not self.owned:
            size += sizeof(_Atom)
        return size

    def __eq__(self, object other):
        cdef Atom other_
        if not isinstance(other, Atom):
            return NotImplemented
        other_ = other
        # FIXME: it should be possible to do a memcmp here.
        return self._state() == other_._state()

    def __hash__(self):
        return hash(tuple(self._state().values()))

    @property
    def serial(self):
        """`int`: The atom serial number.
        """
        assert self._atom is not NULL
        return self._atom.serial

    @property
    def altloc(self):
        """`str`: The alternate location indicator for the atom.
        """
        assert self._atom is not NULL
        return chr(self._atom.altLoc)

    @property
    def name(self):
        """`str`: The atom name.
        """
        assert self._atom is not NULL
        return self._atom.name[:4].decode('ascii').strip("_")

    @property
    def residue_name(self):
        """`str`: The residue name.
        """
        assert self._atom is not NULL
        return self._atom.resName[:3].decode('ascii').strip("_")

    @property
    def residue_number(self):
        """`int`: The residue sequence number.
        """
        assert self._atom is not NULL
        return self._atom.resSeq

    @property
    def segment(self):
        """`str`: The segment identifier.
        """
        assert self._atom is not NULL
        return self._atom.segID[:3].decode('ascii').strip('_')

    @property
    def element(self):
        """`str`: The element symbol.
        """
        assert self._atom is not NULL
        return self._atom.element[:2].decode('ascii').strip('_')

    @property
    def insertion_code(self):
        """`str`: The code for insertion of residues.
        """
        assert self._atom is not NULL
        return chr(self._atom.iCode)

    @property
    def chain_id(self):
        """`str`: The identifier of the chain the atom belongs to.
        """
        assert self._atom is not NULL
        return "{}{}".format(chr(self._atom.chainID1), chr(self._atom.chainID2)).strip()

    @property
    def occupancy(self):
        """`float`: The atom occupancy.
        """
        assert self._atom is not NULL
        return self._atom.occupancy

    @property
    def temperature_factor(self):
        """`float`: The atom temperature factor.
        """
        assert self._atom is not NULL
        return self._atom.tempFactor

    @property
    def charge(self):
        """`int`: The atom charge.
        """
        assert self._atom is not NULL
        return self._atom.charge

    @property
    def x(self):
        assert self._atom is not NULL
        return self._atom.x[0]

    @property
    def y(self):
        assert self._atom is not NULL
        return self._atom.x[1]

    @property
    def z(self):
        assert self._atom is not NULL
        return self._atom.x[2]

    cpdef Atom copy(self):
        """Create a copy of this atom.

        Returns:
            `~pyjess.Atom`: A newly allocated atom with identical attributes.

        .. versionadded:: 0.4.0

        """
        cdef Atom copy = Atom.__new__(Atom)
        copy._atom = <_Atom*> malloc(sizeof(_Atom))
        if copy._atom is NULL:
            raise MemoryError("Failed to allocate atom")
        memcpy(copy._atom, self._atom, sizeof(_Atom))
        return copy


cdef class TemplateAtom:
    """A single template atom.

    .. versionadded:: 0.4.0
       Equality, hashing and pickle protocol support.

    """
    cdef object     owner
    cdef bint       owned
    cdef _TessAtom* _atom

    @classmethod
    def load(cls, file):
        """Load a template atom from the given file.

        Arguments:
            file (str, os.PathLike or file-like object): A file-like object
                opened in text or binary mode to read the template atom from.

        """
        try:
            handle = open(file)
        except TypeError:
            handle = nullcontext(file)
        with handle as f:
            return cls.loads(f.read())

    @classmethod
    def loads(cls, text):
        """Load a template atom from the given string.

        Arguments:
            text (`str`, `bytes` or `bytearray`): The atom line to read the
                atom metadata from.

        """
        cdef bytearray    b
        cdef TemplateAtom atom

        if isinstance(text, str):
            b = bytearray(text, 'utf-8')
        else:
            b = bytearray(text)
        if not b.endswith(b'\n'):
            b.append(b'\n')
        b.append(b'\0')

        atom = TemplateAtom.__new__(TemplateAtom)
        atom._atom = jess.tess_atom.TessAtom_create(<const char*> b)
        if atom._atom == NULL:
            raise ValueError(f"Failed to parse template atom: {text!r}")

        # validate match mode *now* to avoid Jess exiting when it does so later
        if atom.match_mode not in range(-1, 9) and atom.match_mode not in range(100, 108):
            raise ValueError(f"Invalid match mode: {atom.match_mode!r}")

        return atom

    def __cinit__(self):
        self.owner = None
        self.owned = False
        self._atom = NULL

    def __dealloc__(self):
        if not self.owned:
            jess.tess_atom.TessAtom_free(self._atom)

    def __init__(
        self,
        *,
        str chain_id,
        int residue_number,
        double x,
        double y,
        double z,
        object residue_names,
        object atom_names,
        double distance_weight = 0.0,
        int match_mode = 0,
    ):
        """__init__(self, *, chain_id, residue_number, x, y, z, residue_names, atom_names, distance_weight=0.0, match_mode=0)\n--\n

        Create a new template atom.

        Raises:
            `MemoryError`: When the system allocator fails to allocate
                enough memory for the template atom storage.

        """
        cdef size_t m
        cdef char*  p
        cdef size_t ac
        cdef size_t rc
        cdef size_t alloc_size

        # validate match mode to avoid a potential hard exit later
        if match_mode not in range(-1, 9) and match_mode not in range(100, 108):
            raise ValueError(f"Invalid match mode: {match_mode!r}")
        if len(chain_id) > 2:
            raise ValueError(f"Invalid chain ID: {chain_id!r}")

        # compute total allocation
        ac = len(atom_names)
        rc = len(residue_names)
        alloc_size = sizeof(_TessAtom) + sizeof(char*) * (ac + rc) + sizeof(char) * (5*ac + 4*rc)

        # allocate base memory
        self._atom = <_TessAtom*> malloc(alloc_size)
        if self._atom is NULL:
            raise MemoryError("Failed to allocate template atom")

        # copy base data
        self._atom.code = match_mode
        self._atom.resSeq = residue_number
        self._atom.pos[0] = x
        self._atom.pos[1] = y
        self._atom.pos[2] = z
        self._atom.chainID1, self._atom.chainID2 = map(ord, chain_id.ljust(2))
        self._atom.nameCount = ac
        self._atom.resNameCount = rc
        self._atom.distWeight = distance_weight

        # setup string pointers
        p = <char*> &self._atom[1]
        self._atom.name = <char**> p
        p += sizeof(char*)*ac
        for m in range(ac):
            self._atom.name[m] = <char*> p
            p += 5
        self._atom.resName = <char**> p
        p += sizeof(char*)*rc
        for m in range(rc):
            self._atom.resName[m] = <char*> p
            p += 4

        # copy atom names
        for m, name in enumerate(atom_names):
            if isinstance(name, str):
                _name = bytearray(name, 'ascii')
            else:
                _name = bytearray(name)
            if len(_name) > 4:
                raise ValueError(f"Invalid atom name: {name!r}")
            elif len(_name) < 3:
                _name.insert(0, ord('_'))
            copy_token(self._atom.name[m], _name.ljust(4, b'\0'), 4)

        # copy residue names
        for m, name in enumerate(residue_names):
            _name = name.encode('ascii') if isinstance(name, str) else name
            if len(_name) > 3:
                raise ValueError(f"Invalid residue name: {name!r}")
            copy_token(self._atom.resName[m], _name.ljust(3, b'\0'), 3)

    cdef dict _state(self):
        return {
            "chain_id": self.chain_id,
            "residue_number": self.residue_number,
            "x": self.x,
            "y": self.y,
            "z": self.z,
            "residue_names": self.residue_names,
            "atom_names": self.atom_names,
            "distance_weight": self.distance_weight,
            "match_mode": self.match_mode,
        }

    def __repr__(self):
        cdef str  ty   = type(self).__name__
        cdef list args = []
        for k, v in self._state().items():
            args.append(f"{k}={v!r}")
        return f"{ty}({', '.join(args)})"

    def __copy__(self):
        return self.copy()

    def __eq__(self, object other):
        cdef TemplateAtom other_
        if not isinstance(other, TemplateAtom):
            return NotImplemented
        other_ = other
        return self._state() == other_._state()

    def __hash__(self):
        return hash(tuple(self._state().values()))

    def __reduce__(self):
        return functools.partial(type(self), **self._state()), ()

    def __sizeof__(self):
        assert self._atom is not NULL

        cdef size_t ac   = self._atom.nameCount
        cdef size_t rc   = self._atom.resNameCount
        cdef size_t size = sizeof(self)

        if not self.owned:
            size += (
                sizeof(_TessAtom)
                + sizeof(char*) * (ac + rc)
                + sizeof(char) * (5*ac + 4*rc)
            )
        return size

    @property
    def match_mode(self):
        """`int`: The match mode for this particular atom.
        """
        assert self._atom is not NULL
        return self._atom.code

    @property
    def residue_number(self):
        """`int`: The residue sequence number.
        """
        assert self._atom is not NULL
        return self._atom.resSeq

    @property
    def chain_id(self):
        """`str`: The identifier of the chain the atom belongs to.
        """
        assert self._atom is not NULL
        cdef char c1 = jess.tess_atom.TessAtom_chainID1(self._atom)
        cdef char c2 = jess.tess_atom.TessAtom_chainID2(self._atom)
        return "{}{}".format(chr(c1), chr(c2)).strip()

    @property
    def x(self):
        """`float`: The x coordinate of the atom.
        """
        assert self._atom is not NULL
        return self._atom.pos[0]

    @property
    def y(self):
        """`float`: The y coordinate of the atom.
        """
        assert self._atom is not NULL
        return self._atom.pos[1]

    @property
    def z(self):
        """`float`: The z coordinate of the atom.
        """
        assert self._atom is not NULL
        return self._atom.pos[2]

    @property
    def atom_names(self):
        """`tuple` of `str`: The different atom names for this atom.

        .. versionchanged:: 0.4.1
           Property now returns a `tuple` rather than a `list`.

        """
        assert self._atom is not NULL

        cdef int  i
        cdef list l = []

        for i in range(self._atom.nameCount):
            l.append(self._atom.name[i].replace(b'_', b'').decode())
        return tuple(l)

    @property
    def residue_names(self):
        """`tuple` of `str`: The different residue names for this atom.

        .. versionchanged:: 0.4.1
           Property now returns a `tuple` rather than a `list`.

        """
        assert self._atom is not NULL

        cdef int  i
        cdef list l = []

        for i in range(self._atom.resNameCount):
            l.append(self._atom.resName[i].replace(b'_', b'').decode())
        return tuple(l)

    @property
    def distance_weight(self):
        """`float`: The distance weight for this atom.
        """
        assert self._atom is not NULL
        return self._atom.distWeight

    cpdef TemplateAtom copy(self):
        """Create a copy of this template atom.

        Returns:
            `~pyjess.TemplateAtom`: A new template atom object with
            identical attributes.

        .. versionadded:: 0.4.0

        """
        return type(self)(**self._state())


cdef class Template:
    """A template, as a sequence of `TemplateAtom` objects.

    .. versionadded:: 0.4.0
       Equality, hashing and pickle protocol support.

    """
    cdef object         owner
    cdef bint           owned
    cdef _Template*     _tpl
    cdef _TessTemplate* _tess

    @classmethod
    def loads(cls, text, str id = None):
        """Load a template from a string.

        Arguments:
            file (`str`, `os.PathLike`, or file-like object): Either the path
                to a file, or a file-like object opened in **text mode**
                containing the template.
            id (`str`, optional): The identifier of the template. By default,
                the parser will take the one from the ``PDB_ID`` remark if
                found in the header.

        Returns:
            `~pyjess.Template`: The template parsed from the given string.

        See Also:
            `Template.load` to load a template from a file-like object or
            from a path.

        """
        return cls.load(io.StringIO(text), id=id)

    @classmethod
    def load(cls, file, str id = None):
        """Load a template from the given file.

        Arguments:
            file (`str`, `os.PathLike` or file-like object): Either the
                path to a file, or a file-like object opened in **text mode**
                to read the template from.
            id (`str`, optional): The identifier of the template. By default,
                the parser will take the one from the ``PDB_ID`` remark if
                found in the header.

        Returns:
            `~pyjess.Template`: The template parsed from the given file.

        """
        try:
            handle = open(file)
        except TypeError:
            handle = nullcontext(file)
        with handle as f:
            atoms = []
            for line in f:
                if line.startswith("ATOM"):
                    atoms.append(TemplateAtom.loads(line))
                elif id is None and line.startswith("REMARK PDB_ID"):
                    id = line.split(" ", maxsplit=2)[2].strip()
        return cls(atoms, id=id)

    def __cinit__(self):
        self._tpl  = NULL
        self._tess = NULL
        self.owner = None
        self.owned = False

    def __dealloc__(self):
        if not self.owned:
            jess.tess_template.TessTemplate_free(self._tpl)

    def __init__(self, object atoms = (), str id = None):
        """__init__(self, atoms=(), id=None)\n--\n

        Create a new template.

        Arguments:
            atoms (sequence of `~pyjess.TemplateAtom`): The atoms of the
                templates.
            id (`str`, optional): The identifier of the template.

        Raises:
            `MemoryError`: When the system allocator fails to allocate
                enough memory for the template storage.

        """
        cdef int          i
        cdef int          j
        cdef double       dist
        cdef TemplateAtom atom
        cdef size_t       alloc_size
        cdef int          count      = len(atoms)

        alloc_size = (
            sizeof(_Template) + sizeof(_TessTemplate)
            + count * sizeof(_TessAtom*)
            + count * sizeof(double*)
            + count * count * sizeof(double)
        )

        self._tpl = <_Template*> calloc(1, alloc_size)
        if self._tpl is NULL:
            raise MemoryError("Failed to allocate template")

        # setup memory for atoms
        self._tess = <_TessTemplate*> &self._tpl[1]
        self._tess.atom = <_TessAtom**> &self._tess[1]
        for i in range(count):
            self._tess.atom[i] = NULL

        # setup memory and pointers for distances
        self._tess.distance = <double**> &self._tess.atom[count]
        if count > 0:
            self._tess.distance[0] = <double*> &self._tess.distance[count]
        for i in range(1, count):
            self._tess.distance[i] = <double*> &self._tess.distance[i-1][count]

        # setup template function pointers
        self._tpl.free = jess.tess_template.TessTemplate_free
        self._tpl.match = jess.tess_template.TessTemplate_match
        self._tpl.position = jess.tess_template.TessTemplate_position
        self._tpl.count = jess.tess_template.TessTemplate_count
        self._tpl.range = jess.tess_template.TessTemplate_range
        self._tpl.check = jess.tess_template.TessTemplate_check
        self._tpl.candidates = jess.tess_template.TessTemplate_candidates
        self._tpl.name = jess.tess_template.TessTemplate_name
        self._tpl.logE = jess.tess_template.TessTemplate_logE
        self._tpl.distWeight = jess.tess_template.TessTemplate_distWeight
        self._tpl.copy = jess.tess_template.TessTemplate_copy

        # copy name and atom count
        self._tess.count = count
        self._tess.symbol = NULL if id is None else strdup(id.encode())

        # copy atom data
        for i, atom in enumerate(atoms):
            assert i < count
            self._tess.atom[i] = jess.tess_atom.TessAtom_copy(atom._atom)
            if self._tess.atom[i] is NULL:
                raise MemoryError("Failed to allocate template atom")

        # compute distances
        for i in range(count):
            self._tess.distance[i][i] = 0.0
            for j in range(i+1, count):
                dist = jess.tess_atom.TessAtom_distance(self._tess.atom[i], self._tess.atom[j])
                self._tess.distance[i][j] = dist
                self._tess.distance[j][i] = dist

        # compute dimension
        residues = {
            (
                self._tess.atom[i].resSeq ,
                self._tess.atom[i].chainID1,
                self._tess.atom[i].chainID2,
            )
            for i in range(count)
        }
        self._tess.dim = len(residues)

    def __copy__(self):
        return self.copy()

    def __len__(self):
        assert self._tpl is not NULL
        return self._tess.count

    def __getitem__(self, object index):
        assert self._tess is not NULL

        cdef TemplateAtom atom
        cdef ssize_t      length = self._tess.count
        cdef ssize_t      index_

        if isinstance(index, slice):
            indices = range(*index.indices(length))
            return type(self)(atoms=[self[i] for i in indices], id=self.id)
        else:
            index_ = index
            if index_ < 0:
                index_ += length
            if index_ < 0 or index_ >= length:
                raise IndexError(index)
            atom = TemplateAtom.__new__(TemplateAtom)
            atom.owner = self
            atom.owned = True
            atom._atom = self._tess.atom[index_]
            return atom

    def __eq__(self, object other):
        cdef Template other_
        if not isinstance(other, Template):
            return NotImplemented
        other_ = other
        if self.id != other_.id:
            return False
        if self.dimension != other_.dimension:
            return False
        if len(self) != len(other_):
            return False
        return all(x == y for x,y in zip(self, other_))

    def __hash__(self):
        return hash((
            self.id,
            *(hash(x) for x in self)
        ))

    def __reduce__(self):
        return type(self), (list(self), self.id)

    def __sizeof__(self):
        assert self._tess is not NULL

        cdef size_t     i
        cdef size_t     ac
        cdef size_t     rc
        cdef _TessAtom* atom
        cdef size_t     size = sizeof(self)

        size = (
            sizeof(_Template)
            + sizeof(_TessTemplate)
            + self._tess.count * sizeof(_TessAtom*)
            + self._tess.count * sizeof(double*)
            + self._tess.count * self._tess.count * sizeof(double)
        )
        for i in range(self._tess.count):
            atom = self._tess.atom[i]
            ac = atom.nameCount
            rc = atom.resNameCount
            size += (
                sizeof(_TessAtom)
                + sizeof(char*) * (ac + rc)
                + sizeof(char) * (5*ac + 4*rc)
            )
        return size

    @property
    def id(self):
        assert self._tpl is not NULL

        cdef const char* name = self._tpl.name(self._tpl)
        if name is NULL:
            return None
        return name.decode()

    @property
    def dimension(self):
        """`int`: The dimension of the template (i.e. number of residues).
        """
        assert self._tess is not NULL
        return self._tess.dim

    cpdef Template copy(self):
        return Template(
            self,
            self.id
        )


cdef class Query:
    """A query over templates with a given molecule.

    Jess iterates over the templates and attempt matches the query
    molecule, so the hits can actually be generated iteratively. This
    class allows accessing the hits as a Python iterator.

    Attributes:
        jess (`~pyjess.Jess`): The templates this object is currently
            scanning.
        molecule (`~pyjess.Molecule`): The query molecule to align to
            the templates.
        rmsd_threshold (`float`): The RMSD threshold for reporting
            results.
        best_match (`bool`): Whether the query will return only the
            best match to each template.

    """
    cdef _JessQuery* _jq
    cdef bint        _partial
    cdef int         _candidates
    cdef uintptr_t   _prev_tpl
    cdef int         _max_candidates
    cdef _IgnoreType _ignore_chain

    cdef readonly Jess     jess
    cdef readonly Molecule molecule
    cdef readonly bint     best_match
    cdef readonly double   rmsd_threshold

    def __cinit__(self):
        self._jq = NULL
        self._candidates = 0
        self._partial = False
        self._prev_tpl = 0

    def __dealloc__(self):
        jess.jess.JessQuery_free(self._jq)

    def __iter__(self):
        return self

    @property
    def ignore_chain(self):
        """`str` or `None`: The way atom chains are considered or discarded.
        """
        if self._ignore_chain == _IgnoreType.ignoreNone:
            return None
        elif self._ignore_chain == _IgnoreType.ignoreResidues:
            return "residues"
        elif self._ignore_chain == _IgnoreType.ignoreAtoms:
            return "atoms"

    @ignore_chain.setter
    def ignore_chain(self, ignore_chain):
        if ignore_chain is None:
            self._ignore_chain = _IgnoreType.ignoreNone
        elif ignore_chain == "residues":
            self._ignore_chain = _IgnoreType.ignoreResidues
        elif ignore_chain == "atoms":
            self._ignore_chain = _IgnoreType.ignoreAtoms
        else:
            raise ValueError(f"invalid value for `ignore_chain`: {ignore_chain!r}")

    @property
    def max_candidates(self):
        """`int`: The maximum number of candidate hits to report *by template*.
        """
        return None if self._max_candidates == -1 else self._max_candidates

    @max_candidates.setter
    def max_candidates(self, max_candidates):
        if max_candidates is None:
            self._max_candidates = -1
        elif max_candidates >= 0:
            self._max_candidates = max_candidates
        else:
            raise ValueError(f"invalid value for `max_candidates` argument: {max_candidates!r}")

    cdef bint _advance(self) noexcept nogil:
        if self._partial:
            self._partial = False
            return True
        return jess.jess.JessQuery_next(self._jq, self._ignore_chain)

    cdef bint _rewind(self) noexcept nogil:
        self._partial = True

    cdef int _copy_atoms(self, const _Template* tpl, Hit hit) except -1 nogil:
        cdef _Atom** atoms = jess.jess.JessQuery_atoms(self._jq)
        cdef int     count = tpl.count(tpl)

        hit._atoms = <_Atom*> realloc(hit._atoms, count * sizeof(_Atom))
        if hit._atoms is NULL:
            raise MemoryError("Failed to allocate hit atoms")
        for i in range(count):
            memcpy(&hit._atoms[i], atoms[i], sizeof(_Atom))
        return count

    cdef int _copy_superposition(self, _Superposition* sup, Hit hit) noexcept nogil:
        cdef const double* M = jess.super.Superposition_rotation(sup)
        cdef const double* c = jess.super.Superposition_centroid(sup, 0)
        cdef const double* v = jess.super.Superposition_centroid(sup, 1)
        memcpy(hit._rotation, M, 9*sizeof(double))
        memcpy(hit._centre[0], c, 3*sizeof(double))
        memcpy(hit._centre[1], v, 3*sizeof(double))
        return 0

    def __next__(self):
        assert self._jq is not NULL

        cdef double          rmsd
        cdef const double*   rot
        cdef _Template*      tpl     = NULL
        cdef _Template*      hit_tpl = NULL
        cdef _Superposition* sup     = NULL
        cdef Hit             hit     = Hit.__new__(Hit)

        # prepare the hit to be returned
        hit.rmsd = INFINITY
        hit._atoms = NULL
        hit._molecule = self.molecule
        hit_tpl = NULL
        hit_found = False

        # search the next hit without the GIL to allow parallel queries.
        with nogil:
            while self._advance():
                # load current iteration template, and check that the hit
                # was obtained with the current template and not with the
                # previous one
                self._prev_tpl = <uintptr_t> tpl
                tpl = jess.jess.JessQuery_template(self._jq)
                if hit_found and hit_tpl != tpl:
                    self._rewind()
                    break

                # load superposition and compute RMSD for the current iteration
                sup = jess.jess.JessQuery_superposition(self._jq)
                rmsd = jess.super.Superposition_rmsd(sup)

                # NB(@althonos): we don't need to compute the E-value to get the
                #                best match by molecule/template pair since the
                #                logE-value for a fixed pair varies by the RMSD
                #                term only (see `TessTemplate_logE`)

                # check that the candidate passes threshold, and return it
                # if not in best match, otherwise record it until the next
                # template is reached (or the iterator finished)
                if rmsd <= self.rmsd_threshold and rmsd < hit.rmsd:
                    # check if the rotation matrix contains NaN values
                    rot = jess.super.Superposition_rotation(sup)
                    nan = False
                    for i in range(9):
                        nan |= isnan(rot[i])

                    if nan:
                        with gil:
                            warnings.warn(
                                "Jess returned a superposition matrix with NaN values",
                                UserWarning,
                                stacklevel=2,
                            )
                    else:
                        self._copy_atoms(tpl, hit)
                        self._copy_superposition(sup, hit)
                        hit.rmsd = rmsd
                        hit_tpl = tpl
                        hit_found = True

                # check if we already made it to the next template,
                # or if we need to short-circuit the iteration and
                # force the query to move to the next template as
                # we found too many candidates already.
                if <uintptr_t> tpl != self._prev_tpl:
                    self._candidates = 0
                else:
                    self._candidates += 1
                if self._max_candidates != -1 and self._candidates > self._max_candidates:
                    self._candidates = 0
                    jess.jess.JessQuery_nextTemplate(self._jq)

                # free superposition items (as relevant data was copied in
                # the Hit if needed) and return hits immediately if we are
                # not in best match mode
                jess.super.Superposition_free(sup)
                if hit_found and not self.best_match:
                    break

        if not hit_found:
            raise StopIteration

        # get the template object for the hit
        hit.template = self.jess._templates[self.jess._indices[<size_t> hit_tpl]]
        return hit


cdef class Hit:
    """A hit identified between a query molecule and a target template.

    Attributes:
        rmsd (`float`): The RMSD between the aligned structures.
        template (`~pyjess.Template`): The template that matched the
            query molecule.
        molecule (`~pyjess.Molecule`): The query molecule.

    """
    cdef double[9]       _rotation
    cdef double[2][3]    _centre
    cdef _Atom*          _atoms

    cdef readonly double   rmsd
    cdef readonly Template template
    cdef          Molecule _molecule

    def __dealloc__(self):
        free(self._atoms)

    def __getstate__(self):
        return {
            "rotation": list(self._rotation),
            "centre": list(self._centre),
            "atoms": self.atoms(transform=False),
            "rmsd": self.rmsd,
            "template": self.template,
            "molecule": self.molecule(transform=False),
        }

    def __setstate__(self, state):
        cdef size_t i
        cdef size_t count
        cdef Atom   atom

        self.rmsd = state["rmsd"]
        self.template = state["template"]
        self._molecule = state["molecule"]
        self._rotation = state["rotation"]
        self._centre = state["centre"]

        # check number of atoms is consistent
        count = len(self.template)
        if len(state["atoms"]) != count:
            raise ValueError(f"unexpected number of atoms: {len(state['atoms'])!r} (expected {count!r})")
        # allocate or reallocate memory for atoms
        self._atoms = <_Atom*> realloc(self._atoms, count * sizeof(_Atom))
        if self._atoms is NULL:
            raise MemoryError("Failed to allocate hit atoms")
        # copy atom data
        for i, atom in enumerate(state["atoms"]):
            memcpy(&self._atoms[i], atom._atom, sizeof(_Atom))

    @property
    def determinant(self):
        """`float`: The determinant of the rotation matrix.
        """
        cdef const double* p   = self._rotation
        cdef double        det = 0.0

        with nogil:
            det += p[0] * (p[4] * p[8] - p[5] * p[7])
            det -= p[1] * (p[3] * p[8] - p[5] * p[6])
            det += p[2] * (p[3] * p[7] - p[4] * p[6])
        return det

    @property
    def log_evalue(self):
        """`float`: The logarithm of the E-value estimated for the hit.
        """
        assert self.template._tpl is not NULL

        cdef int    n
        cdef double e

        with nogil:
            n = jess.molecule.Molecule_count(self._molecule._mol)
            e = self.template._tpl.logE(self.template._tpl, self.rmsd, n)
        return e

    @property
    def evalue(self):
        """`float`: The E-value estimated for the hit.
        """
        cdef int    n
        cdef double e

        with nogil:
            n = jess.molecule.Molecule_count(self._molecule._mol)
            e = exp(self.template._tpl.logE(self.template._tpl, self.rmsd, n))
        return e

    cpdef list atoms(self, bint transform=True):
        """Get the list of query atoms matching the template.

        Arguments:
            transform (`bool`): Whether or not to transform coordinates
                of hits into template frame.

        Returns:
            `list` of `~pyjess.Atom`: The list of matching atoms.

        """
        assert self.template._tpl is not NULL

        cdef Atom atom
        cdef int  i
        cdef int  j
        cdef int  k
        cdef int  count = self.template._tpl.count(self.template._tpl)
        cdef list atoms = []

        cdef const double* M = self._rotation
        cdef const double* c = self._centre[0]
        cdef const double* v = self._centre[1]

        for k in range(count):
            atom = Atom.__new__(Atom)
            if transform:
                atom._atom = <_Atom*> malloc(sizeof(_Atom))
                memcpy(atom._atom, &self._atoms[k], sizeof(_Atom))
                for i in range(3):
                    atom._atom.x[i] = v[i]
                    for j in range(3):
                        atom._atom.x[i] += M[3*i + j] * (self._atoms[k].x[j] - c[j])
            else:
                atom.owned = True
                atom.owner = self
                atom._atom = &self._atoms[k]

            atoms.append(atom)

        return atoms

    cpdef Molecule molecule(self, bint transform=False):
        """Get the molecule matching the template.

        Arguments:
            transform (`bool`): Whether or not to transform coordinates
                of the molecule atoms into template frame.

        Returns:
            `~pyjess.Molecule`: The matching molecule, optionally
            rotated to match the template coordinate.

        .. versionadded:: 0.5.0

        """
        assert self.template._tpl is not NULL

        cdef _Atom*        atom
        cdef Molecule      mol
        cdef size_t        i
        cdef size_t        j
        cdef size_t        k
        cdef const double* M = self._rotation
        cdef const double* c = self._centre[0]
        cdef const double* v = self._centre[1]

        if not transform:
            return self._molecule

        mol = self._molecule.copy()
        for k in range(mol._mol.count):
            atom = mol._mol.atom[k]
            for i in range(3):
                atom.x[i] = v[i]
                for j in range(3):
                    atom.x[i] += M[3*i + j] * (self._molecule._mol.atom[k].x[j] - c[j])

        return mol


cdef class Jess:
    """A handle to run Jess over a list of templates.

    .. versionadded:: 0.4.0
       Equality, hashing and pickle protocol support.

    """
    cdef _Jess* _jess
    cdef dict   _indices
    cdef tuple  _templates
    cdef size_t length

    def __cinit__(self):
        self._jess = NULL
        self.length = 0

    def __dealloc__(self):
        jess.jess.Jess_free(self._jess)

    def __init__(self, object templates = ()):
        """__init__(self, templates=())\n--\n

        Create a new Jess database containing the given templates.

        Arguments:
            templates (sequence of `~pyjess.Template`): The templates to
                index in the database for further querying.

        Caution:
            The `~pyjess.Template` objects given in argument will be copied
            because the internal C data structure requires ownership of the
            data. Modification to the original `~pyjess.Template` objects will
            not have an effect on the newly created `~pyjess.Jess` templates.

        """
        cdef Template   template
        cdef _Template* tpl
        cdef list       _templates = []

        self._jess = jess.jess.Jess_create()
        self._indices = {}

        for template in templates:
            # NOTE: the Jess storage owns the data, so we make a copy of the
            #       template given as argument to avoid a double-free.
            tpl = template._tpl.copy(template._tpl)
            jess.jess.Jess_addTemplate(self._jess, tpl)
            self._indices[<size_t> tpl] = self.length
            _templates.append(template)
            self.length += 1

        self._templates = tuple(_templates)

    def __copy__(self):
        return self.copy()

    def __reduce__(self):
        return type(self), (self._templates,)

    def __eq__(self, object other):
        cdef Jess other_
        if not isinstance(other, Jess):
            return NotImplemented
        other_ = other
        return self._templates == other_._templates

    def __hash__(self):
        return hash((Jess, self._templates))

    def __len__(self):
        return self.length

    def __getitem__(self, object index):
        cdef ssize_t index_

        if isinstance(index, slice):
            indices = range(*index.indices(self.length))
            return type(self)(map(self.__getitem__, indices))
        else:
            index_ = index
            if index_ < 0:
                index_ += self.length
            if index_ < 0 or index_ >= self.length:
                raise IndexError(index)
            return self._templates[index_]

    cpdef Jess copy(self):
        """Create a copy of the `Jess` object.

        Returns:
            `~pyjess.Jess`: A `Jess` object containing the same templates.

        .. versionadded:: 0.4.0

        """
        return type(self)(self._templates)

    def query(
        self,
        Molecule molecule,
        double rmsd_threshold,
        double distance_cutoff,
        double max_dynamic_distance,
        *,
        object max_candidates = None,
        object ignore_chain = None,
        bint best_match = False,
        bint reorder = True,
    ):
        """Scan for templates matching the given molecule.

        Arguments:
            molecule (`~pyjess.Molecule`): The protein to match the
                templates to.
            rmsd_threshold (`float`): The RMSD threshold for reporting
                results.
            distance_cutoff (`float`): The global distance cutoff
                used to guide the search.
            max_dynamic_distance (`float`): The maximum template/query
                dynamic distance after adding the global distance cutoff
                and the individual atom distance cutoff defined for each
                atom of the template.
            max_candidates (`int` or `None`): The maximum number of candidate
                hits to report by template. If a non-`None` value is given,
                it may speed up querying for unspecific templates, but also
                produce results potentially inconsistent with Jess.
            ignore_chain (`str` or `None`): Whether to check or ignore the
                chain of the atoms to match. The different supported modes
                are:

                - `None`: Force the atoms in the molecule to belong
                  to different (resp. same) chains if so is the case
                  in the template.
                - ``residues``: Allow atoms to belong to different
                  (resp. same) chains even if it is not the case in
                  the template, but force all atoms of a residue to
                  belong to the same chain.
                - ``atoms``: Allow atoms to belong to any chain,
                  independently to the template or the residue they
                  belong to.

            best_match (`bool`): Pass `True` to return only the best match
                to each template, based on RMSD. In case of ties, the
                first match is returned. Note that a match must still
                be passing the RMSD threshold given in ``rmsd_threshold``
                to be returned.
            reorder (`bool`): Whether to enable template atom reordering
                to accelerate matching in the scanner algorithm. Pass
                `False` to reverse to the original, slower algorithm
                which matches atoms in the same order as they appear in
                the template, at the cost of longer run times.

        Returns:
            `~pyjess.Query`: An iterator over the query hits.

        Caution:
            Since ``v0.6.0``, this function uses an optimized variant of
            the Jess scanning algorithm which minimized the number of steps
            needed to generate matches, by re-ordering the order the
            template atoms are iterated upon. Because of this change,
            the query may return *exactly* the same matches but in an order
            that *differs* from the original Jess version. If you really
            need results in the original order, set ``reorder`` to `False`.

        .. versionadded:: 0.6.0
            The ``reorder`` argument, defaulting to `True`.

        .. versionchanged:: 0.7.0
            Default value of ``max_candidates`` argument to `None`.

        .. versionchanged:: 0.7.0
            ``ignore_chain`` now expects string variants rather than `bool`.

        """

        if ignore_chain is True:
            warnings.warn(
                "`ignore_chain` parameter expects string parameters "
                "to specificy the mode since PyJess v0.7.0. "
                "Use `ignore_chain='atoms'` instead of `ignore_chain=True`",
                DeprecationWarning,
            )
            ignore_chain="atoms"
        elif ignore_chain is False:
            warnings.warn(
                "`ignore_chain` parameter expects string parameters "
                "to specificy the mode since PyJess v0.7.0. "
                "Use `ignore_chain=None` instead of `ignore_chain=False`",
                DeprecationWarning,
            )
            ignore_chain=None

        cdef Query query = Query.__new__(Query)
        query.max_candidates = max_candidates
        query.ignore_chain = ignore_chain
        query.rmsd_threshold = rmsd_threshold
        query.best_match = best_match
        query.molecule = molecule
        query.jess = self
        query._jq = jess.jess.Jess_query(
            self._jess,
            molecule._mol,
            distance_cutoff,
            max_dynamic_distance,
            reorder,
        )
        return query
