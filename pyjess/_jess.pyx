# coding: utf-8
# cython: language_level=3, linetrace=True
"""Bindings to Jess, a 3D template matching software.

References:
    - Barker, J. A., & Thornton, J. M. (2003). &An algorithm for
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

from libc.math cimport exp
from libc.stdio cimport FILE, fclose, fdopen, printf
from libc.stdlib cimport calloc, free, malloc
from libc.string cimport memcpy, memset, strncpy, strdup

cimport jess.atom
cimport jess.jess
cimport jess.molecule
cimport jess.super
cimport jess.tess_template
cimport jess.tess_atom
from jess.atom cimport Atom as _Atom
from jess.jess cimport Jess as _Jess
from jess.jess cimport JessQuery as _JessQuery
from jess.molecule cimport Molecule as _Molecule
from jess.super cimport Superposition as _Superposition
from jess.template cimport Template as _Template
from jess.tess_template cimport TessTemplate as _TessTemplate
from jess.tess_atom cimport TessAtom as _TessAtom

# --- Python imports ---------------------------------------------------------

import contextlib
import io
import os

# --- Utils ------------------------------------------------------------------

cdef inline void copy_token(char* dst, const char* src, size_t n) noexcept nogil:
    cdef size_t i
    strncpy(dst, src, n)
    for i in range(n):
        if dst[i] == ord(' ') or dst[i] == 0:
            dst[i] = ord('_')
    dst[n] = 0

@contextlib.contextmanager
def nullcontext(return_value=None):
    yield return_value

# --- Classes ----------------------------------------------------------------

cdef class Molecule:
    """A molecule structure, as a sequence of `Atom` objects.
    """
    cdef _Molecule* _mol

    @classmethod
    def loads(cls, text):
        return cls.load(io.StringIO(text))

    @classmethod
    def load(cls, file):
        try:
            handle = open(file)
        except TypeError:
            handle = nullcontext(file)
        with handle as f:
            atoms = []
            id = None
            for line in f:
                if line.startswith("HEADER"):
                    id = line[62:66].strip()
                    if not id:
                        id = None
                elif line.startswith(("ATOM", "HETATM")):
                    atoms.append(Atom.loads(line))
        return cls(atoms, id=id)

    def __cinit__(self):
        self._mol = NULL

    def __dealloc__(self):
        jess.molecule.Molecule_free(self._mol)

    def __init__(self, object atoms = (), str id = None):
        cdef Atom atom
        cdef int i
        cdef int count = len(atoms)

        if id is not None and len(id) != 4:
            raise ValueError(f"Invalid PDB ID: {id!r}")

        self._mol = <_Molecule*> malloc(sizeof(_Molecule) + count * sizeof(_Atom*))
        if self._mol is NULL:
            raise MemoryError("Failed to allocate molecule")

        self._mol.count = count
        for i in range(count):
            self._mol.atom[i] = NULL
        if id is None:
            memset(self._mol.id, 0, 5)
        else:
            strncpy(self._mol.id, id.encode('ascii'), 5)

        for i, atom in enumerate(atoms):
            self._mol.atom[i] = <_Atom*> malloc(sizeof(_Atom))
            if self._mol.atom[i] is NULL:
                raise MemoryError("Failed to allocate atom")
            memcpy(self._mol.atom[i], atom._atom, sizeof(_Atom))

    def __len__(self):
        assert self._mol is not NULL
        return self._mol.count

    def __getitem__(self, ssize_t index):
        assert self._mol is not NULL

        cdef Atom    atom
        cdef ssize_t length = self._mol.count
        cdef ssize_t index_ = index

        if index_ < 0:
            index_ += length
        if index_ < 0 or index_ >= length:
            raise IndexError(index)

        atom = Atom.__new__(Atom)
        atom.owner = self
        atom.owned = True
        atom._atom = <_Atom*> jess.molecule.Molecule_atom(self._mol, index_)
        return atom

    def __sizeof__(self):
        assert self._mol is not NULL
        return (
            sizeof(self)
            + sizeof(_Molecule)
            + self._mol.count*(sizeof(_Atom*) + sizeof(_Atom))
        )

    @property
    def id(self):
        assert self._mol is not NULL

        cdef const char* pdb_id

        pdb_id = jess.molecule.Molecule_id(self._mol)
        if pdb_id is NULL:
            return None
        return pdb_id.decode()

    cpdef Molecule conserved(self, double cutoff = 0.0):
        assert self._mol is not NULL
        cdef Atom atom
        return Molecule(
            id=self.id,
            atoms=[
                atom
                for atom in self
                if cutoff <= 0.0
                or atom._atom.tempFactor >= cutoff
            ]
        )


cdef class Atom:
    """A single atom in a molecule.
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
        cdef bytearray b
        cdef Atom      atom

        if isinstance(text, str):
            b = bytearray(text, 'utf-8')
        else:
            b = bytearray(text)
        if not b.endswith(b'\n'):
            b.append(b'\n')
        b.append(b'\0')

        atom = cls.__new__(cls)
        atom._atom = <_Atom*> malloc(sizeof(_Atom))
        if atom._atom == NULL:
            raise MemoryError("Failed to allocate atom")

        if not jess.atom.Atom_parse(atom._atom, b):
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
        float x,
        float y,
        float z,
        float occupancy = 0.0,
        float temperature_factor = 0.0,
        str segment = '',
        str element = '',
        int charge = 0,
    ):
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
        self._atom.chainID2 = ord(chain_id[1]) if len(chain_id) > 1 else ord('0')
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

    def __repr__(self):
        cdef str  ty   = type(self).__name__
        cdef list args = [
            f"serial={self.serial!r}",
            f"name={self.name!r}",
            f"altloc={self.altloc!r}",
            f"residue_name={self.residue_name!r}",
            f"chain_id={self.chain_id!r}",
            f"residue_number={self.residue_number!r}",
            f"x={self.x!r}",
            f"y={self.y!r}",
            f"z={self.z!r}",
            f"segment={self.segment!r}",
            f"insertion_code={self.insertion_code!r}",
        ]
        if self.occupancy:
            args.append(f"occupancy={self.occupancy!r}")
        if self.temperature_factor:
            args.append(f"temperature_factor={self.temperature_factor!r}")
        if self.element:
            args.append(f"element={self.element!r}")
        if self.charge:
            args.append(f"charge={self.charge!r}")
        return f"{ty}({', '.join(args)})"

    def __sizeof__(self):
        cdef size_t size = sizeof(self)
        if not self.owned:
            size += sizeof(_Atom)
        return size

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


cdef class TemplateAtom:
    """A single template atom.
    """
    cdef object     owner
    cdef bint       owned
    cdef _TessAtom* _atom

    @classmethod
    def load(cls, file):
        """Load a template atom from the given file.

        Arguments:
            file (file-like object): A file-like object opened in
                text mode to read the template atom from.

        """
        return cls.loads(file.read())

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

    def __repr__(self):
        cdef str  ty   = type(self).__name__
        cdef list args = [
            f"chain_id={self.chain_id!r}",
            f"residue_number={self.residue_number!r}",
            f"x={self.x!r}",
            f"y={self.y!r}",
            f"z={self.z!r}",
            f"residue_names={self.residue_names!r}",
            f"atom_names={self.atom_names!r}",
        ]
        if self.distance_weight:
            args.append(f"distance_weight={self.distance_weight!r}")
        if self.match_mode:
            args.append(f"match_mode={self.match_mode!r}") 
        return f"{ty}({', '.join(args)})"

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
        """`list` of `str`: The different atom names for this atom.
        """
        assert self._atom is not NULL

        cdef int  i
        cdef list l = []

        for i in range(self._atom.nameCount):
            l.append(self._atom.name[i].replace(b'_', b'').decode())
        return l

    @property
    def residue_names(self):
        """`list` of `str`: The different residue names for this atom.
        """
        assert self._atom is not NULL

        cdef int  i
        cdef list l = []

        for i in range(self._atom.resNameCount):
            l.append(self._atom.resName[i].replace(b'_', b'').decode())
        return l

    @property
    def distance_weight(self):
        """`float`: The distance weight for this atom.
        """
        assert self._atom is not NULL
        return self._atom.distWeight


cdef class Template:
    """A template, as a sequence of `TemplateAtom` objects.
    """
    cdef object         owner
    cdef bint           owned
    cdef _Template*     _tpl
    cdef _TessTemplate* _tess

    @classmethod
    def loads(cls, text, str id = None):
        return cls.load(io.StringIO(text), id=id)

    @classmethod
    def load(cls, file, str id = None):
        try:
            handle = open(file)
        except TypeError:
            handle = nullcontext(file)
        with handle as f:
            atoms = []
            for line in f:
                if line.startswith("ATOM"):
                    atoms.append(TemplateAtom.loads(line))
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
        residues = { self._tess.atom[i].resSeq for i in range(count) }
        self._tess.dim = len(residues)

    def __len__(self):
        assert self._tpl is not NULL
        return self._tess.count

    def __getitem__(self, ssize_t index):
        assert self._tess is not NULL

        cdef TemplateAtom atom
        cdef ssize_t      length = self._tess.count
        cdef ssize_t      index_ = index

        if index_ < 0:
            index_ += length
        if index_ < 0 or index_ >= length:
            raise IndexError(index)

        atom = TemplateAtom.__new__(TemplateAtom)
        atom.owner = self
        atom.owned = True
        atom._atom = self._tess.atom[index_]
        return atom

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
        max_candidates (`int`): The maximum number of candidate hits
            to report.

    """
    cdef _JessQuery* _jq
    cdef int         _candidates

    cdef readonly Jess     jess
    cdef readonly Molecule molecule
    cdef readonly bint     ignore_chain
    cdef readonly double   rmsd_threshold
    cdef readonly int      max_candidates

    def __cinit__(self):
        self._jq = NULL
        self._candidates = 0

    def __dealloc__(self):
        jess.jess.JessQuery_free(self._jq)

    def __iter__(self):
        return self

    def __next__(self):
        assert self._jq is not NULL

        cdef _Template*      tpl
        cdef _Superposition* sup
        cdef _Atom**         atoms
        cdef double          rmsd
        cdef Hit             hit

        while jess.jess.JessQuery_next(self._jq, self.ignore_chain) and self._candidates < self.max_candidates:
            tpl = jess.jess.JessQuery_template(self._jq)
            count = tpl.count(tpl)
            sup = jess.jess.JessQuery_superposition(self._jq)
            atoms = jess.jess.JessQuery_atoms(self._jq)
            rmsd = jess.super.Superposition_rmsd(sup)
            # check that the candidate passes threshold
            if rmsd <= self.rmsd_threshold:
                # create a new hit
                hit = Hit.__new__(Hit)
                hit.rmsd = rmsd
                # record superposition object
                # (TODO: expose as a Superposition object?)
                hit._sup = sup
                # copy atoms (the pointer will be invalidated on the next
                # call of JessQuery_next)
                hit._atoms = <_Atom*> calloc(count, sizeof(_Atom))
                for i in range(count):
                    memcpy(&hit._atoms[i], atoms[i], sizeof(_Atom))
                # record molecule from which we got the hit (the query molecule)
                hit.molecule = self.molecule
                # create a new object to wrap the template that got a hit
                # (no need to copy, we keep a reference to the template
                # list from the original Jess object).
                hit.template = Template.__new__(Template)
                hit.template._tpl = tpl
                hit.template.owner = self.jess
                hit.template.owned = True
                return hit
            # free superposition items that are not used for hits
            jess.super.Superposition_free(sup)
            self._candidates += 1

        raise StopIteration


cdef class Hit:
    """A hit identified between a query molecule and a target template.

    Attributes:
        rmsd (`float`): The RMSD between the aligned structures.
        template (`~pyjess.Template`): The template that matched the
            query molecule.
        molecule (`~pyjess.Molecule`): The query molecule.

    """
    cdef _Superposition* _sup
    cdef _Atom*          _atoms

    cdef readonly double   rmsd
    cdef readonly Template template
    cdef readonly Molecule molecule

    def __dealloc__(self):
        jess.super.Superposition_free(self._sup)
        free(self._atoms)

    @property
    def determinant(self):
        assert self._sup is not NULL
        cdef const double* p = jess.super.Superposition_rotation(self._sup)
        cdef double det = 0.0
        det += p[0] * (p[4] * p[8] - p[5] * p[7])
        det -= p[1] * (p[3] * p[8] - p[5] * p[6])
        det += p[2] * (p[3] * p[7] - p[4] * p[6])
        return det

    @property
    def log_evalue(self):
        assert self.template._tpl is not NULL
        cdef int n = jess.molecule.Molecule_count(self.molecule._mol)
        return self.template._tpl.logE(self.template._tpl, self.rmsd, n)

    @property
    def evalue(self):
        assert self.template._tpl is not NULL
        cdef int n = jess.molecule.Molecule_count(self.molecule._mol)
        return exp(self.template._tpl.logE(self.template._tpl, self.rmsd, n))

    cpdef list atoms(self, bint transform=True):
        """Get the list of query atoms matching the template.

        Arguments:
            transform (`bool`): Whether or not to transform coordinates
                of hits into template frame.

        Returns:
            `list` of `~pyjess.Atom`: The list of matching atoms.

        """
        assert self.template._tpl is not NULL
        assert self._sup is not NULL

        cdef Atom atom
        cdef int  i
        cdef int  j
        cdef int  k
        cdef int  count = self.template._tpl.count(self.template._tpl)
        cdef list atoms = []

        cdef const double* M = jess.super.Superposition_rotation(self._sup)
        cdef const double* c = jess.super.Superposition_centroid(self._sup, 0)
        cdef const double* v = jess.super.Superposition_centroid(self._sup, 1)

        for k in range(count):

            atom = Atom.__new__(Atom)
            atom._atom = <_Atom*> malloc(sizeof(_Atom))
            memcpy(atom._atom, &self._atoms[k], sizeof(_Atom))

            if transform:
                for i in range(3):
                    atom._atom.x[i] = v[i]
                    for j in range(3):
                        atom._atom.x[i] += M[3*i + j] * (self._atoms[k].x[j] - c[j])

            atoms.append(atom)

        return atoms


cdef class Jess:
    """A handle to run Jess over a list of templates.
    """
    cdef _Jess* _jess

    def __cinit__(self):
        self._jess = NULL

    def __dealloc__(self):
        jess.jess.Jess_free(self._jess)

    def __init__(self, object templates):
        cdef Template template
        self._jess = jess.jess.Jess_create()
        for template in templates:
            # NOTE: the Jess storage owns the data, so we make a copy of the
            #       template given as argument to avoid a double-free.
            jess.jess.Jess_addTemplate(self._jess, template._tpl.copy(template._tpl))

    cpdef Query query(
        self,
        Molecule molecule,
        double rmsd_threshold,
        double distance_cutoff,
        double max_dynamic_distance,
        int max_candidates = 1000,
        bint ignore_chain = False,
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

        Returns:
            `~pyjess.Query`: An iterator over the query hits.

        """
        cdef Query query = Query.__new__(Query)
        query.ignore_chain = ignore_chain
        query.max_candidates = max_candidates
        query.rmsd_threshold = rmsd_threshold
        query.molecule = molecule
        query.jess = self
        query._jq = jess.jess.Jess_query(
            self._jess,
            molecule._mol,
            distance_cutoff,
            max_dynamic_distance,
        )
        return query
