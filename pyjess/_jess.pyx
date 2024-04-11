cimport cython
from cpython.unicode cimport PyUnicode_FromStringAndSize

from libc.math cimport exp
from libc.stdio cimport FILE, fclose, fdopen, printf
from libc.stdlib cimport calloc, free, malloc
from libc.string cimport memcpy, strncpy

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

import os


cdef class Molecule:
    cdef _Molecule* _mol

    def __cinit__(self):
        self._mol = NULL

    def __dealloc__(self):
        jess.molecule.Molecule_free(self._mol)

    def __init__(
        self,
        object file,
        bint ignore_endmdl = False,
        double conservation_cutoff = 0.0,
    ):
        """Create a new molecule by loading the given file.
        """
        cdef int   fd = os.open(file, os.O_RDONLY)
        cdef FILE* f  = fdopen(fd, "r")

        try:
            self._mol = jess.molecule.Molecule_create(f, ignore_endmdl, conservation_cutoff)
            if self._mol is NULL:
                raise ValueError(f"Failed to parse molecule file from {file!r}")
        finally:
            fclose(f)

    def __len__(self):
        assert self._mol is not NULL
        return jess.molecule.Molecule_count(self._mol)

    def __getitem__(self, ssize_t index):
        assert self._mol is not NULL

        cdef Atom    atom
        cdef ssize_t length = jess.molecule.Molecule_count(self._mol)
        cdef ssize_t index_ = index

        if index_ < 0:
            index_ += length
        if index_ < 0 or index_ >= length:
            raise IndexError(index)

        atom = Atom.__new__(Atom)
        atom.owner = self
        atom._atom = jess.molecule.Molecule_atom(self._mol, index_)
        return atom

    @property
    def id(self):
        assert self._mol is not NULL

        cdef const char* pdb_id = NULL

        pdb_id = jess.molecule.Molecule_id(self._mol)
        if pdb_id is NULL:
            return None
        return pdb_id.decode()


@cython.no_gc_clear
cdef class Atom:
    cdef object owner
    cdef _Atom* _atom

    def __cinit__(self):
        self._atom = NULL
        self.owner = None

    def __dealloc__(self):
        if self.owner is None:
            free(self._atom)

    @property
    def serial(self):
        """`int`: The atom serial number.
        """
        assert self._atom is not NULL
        return self._atom.serial

    @property
    def altloc(self):
        """`str`: The alternate locaiton indicator for the atom.
        """
        assert self._atom is not NULL
        return chr(self._atom.altLoc)

    @property
    def name(self):
        """`str`: The atom name.
        """
        assert self._atom is not NULL
        return PyUnicode_FromStringAndSize(self._atom.name, 4).strip("_")

    @property
    def residue_name(self):
        """`str`: The residue name.
        """
        assert self._atom is not NULL
        return PyUnicode_FromStringAndSize(self._atom.resName, 3).strip("_")

    @property
    def residue_number(self):
        """`int`: The residue sequence number.
        """
        assert self._atom is not NULL
        return self._atom.resSeq

    @property
    def segment(self):
        """`str`: The segment identifier
        """
        assert self._atom is not NULL
        return PyUnicode_FromStringAndSize(self._atom.segID, 4).strip()

    @property
    def element(self):
        """`str`: The element symbol.
        """
        assert self._atom is not NULL
        return PyUnicode_FromStringAndSize(self._atom.element, 3)

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


@cython.no_gc_clear
cdef class TessAtom:
    cdef object owner
    cdef _TessAtom* _atom

    @classmethod
    def load(cls, file):
        return cls.loads(file.read())

    @classmethod
    def loads(cls, text):
        cdef bytearray b
        cdef TessAtom  atom
        
        if isinstance(text, str):
            b = bytearray(text, 'utf-8')
        else:
            b = bytearray(text)
        if not b.endswith(b'\n'):
            b.append(b'\n')
        
        atom = TessAtom.__new__(TessAtom)
        atom._atom = jess.tess_atom.TessAtom_create(<const char*> b)
        if atom._atom == NULL:
            return ValueError("Invalid atom")

        # validate match mode *now* to avoid Jess exiting when it does so later
        if atom.match_mode not in range(-1, 9) and atom.match_mode not in range(100, 108):
            raise ValueError(f"Invalid match mode: {atom.match_mode!r}")

        return atom

    def __cinit__(self):
        self.owner = None
        self._atom = NULL

    def __dealloc__(self):
        if self.owner is None:
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
        cdef void*  p
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
        self.owner = None
        self._atom = <_TessAtom*> malloc(alloc_size)
        if self._atom is NULL:
            raise MemoryError("failed to allocate atom")
        
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
        p = &self._atom[1]
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
        for i, name in enumerate(atom_names):
            _name = name.encode('ascii') if isinstance(name, str) else name
            if len(_name) > 4:
                raise ValueError(f"Invalid atom name: {name!r}")
            strncpy(self._atom.name[i], b'___\0', 5)
            for j, c in enumerate(_name):
                self._atom.name[i][j] = c

        # copy residue names
        for i, name in enumerate(residue_names):
            _name = name.encode('ascii') if isinstance(name, str) else name
            if len(_name) > 3:
                raise ValueError(f"Invalid residue name: {name!r}")
            strncpy(self._atom.resName[i], b'___\0', 4)
            for j, c in enumerate(_name):
                self._atom.resName[i][j] = c

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"{ty}(chain_id={self.chain_id!r}, residue_number={self.residue_number!r}, x={self.x!r}, y={self.y!r}, z={self.z!r}, residue_names={self.residue_names!r}, atom_names={self.atom_names!r}, distance_weight={self.distance_weight!r}, match_mode={self.match_mode!r})"

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


@cython.no_gc_clear
cdef class Template:
    cdef object         owner
    cdef _Template*     _tpl
    cdef _TessTemplate* _tess

    def __cinit__(self):
        self._tpl  = NULL
        self._tess = NULL
        self.owner = None

    def __dealloc__(self):
        if self.owner is None:
            jess.tess_template.TessTemplate_free(self._tpl)

    def __init__(self, str name, object file):
        """Create a new template by loading the given file.

        Arguments:
            name (`str`): The name of the template.
            file (`str`, `bytes` or `os.PathLike`): The filename of the file
                containing the template atoms.

        """
        cdef bytes name_  = name.encode()
        cdef int   fd     = os.open(file, os.O_RDONLY)
        cdef FILE* f      = fdopen(fd, "r")

        try:
            self._tpl = jess.tess_template.TessTemplate_create(f, name_)
            if self._tpl is NULL:
                raise ValueError(f"Failed to parse template file from {file!r}")
            self._tess = <_TessTemplate*> &self._tpl[1]
        finally:
            fclose(f)

    def __len__(self):
        assert self._tpl is not NULL
        return self._tess.count

    def __getitem__(self, ssize_t index):
        assert self._tess is not NULL

        cdef TessAtom atom
        cdef ssize_t  length = self._tess.count
        cdef ssize_t  index_ = index

        if index_ < 0:
            index_ += length
        if index_ < 0 or index_ >= length:
            raise IndexError(index)

        atom = TessAtom.__new__(TessAtom)
        atom.owner = self
        atom._atom = self._tess.atom[index_]
        return atom

    @property
    def name(self):
        assert self._tpl is not NULL
        return self._tpl.name(self._tpl).decode()

    @property
    def dimension(self):
        """`int`: The dimension of the template (i.e. number of residues).
        """
        assert self._tess is not NULL
        return self._tess.dim


cdef class JessQuery:
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
        self._jq = NULL

    def __iter__(self):
        return self

    def __next__(self):
        assert self._jq is not NULL

        cdef _Template*      tpl   = NULL
        cdef _Superposition* sup   = NULL
        cdef _Atom**         atoms = NULL
        cdef double*         p     = NULL
        cdef const double*   c     = NULL
        cdef const double*   v     = NULL
        cdef double          det
        cdef double          logE
        cdef double          rmsd
        cdef int             kill_switch = 0
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
                # (TODO: expose as a Superposition object)
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
                return hit
            # free superposition items that are not used for hits
            jess.super.Superposition_free(sup)
            self._candidates += 1

        raise StopIteration


cdef class Hit:
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
        """Get the list of atoms.
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

    cpdef JessQuery query(
        self,
        Molecule molecule,
        double rmsd_threshold,
        double distance_cutoff,
        double max_dynamic_distance,
        int max_candidates = 1000,
        bint ignore_chain = False,
    ):
        cdef JessQuery query = JessQuery.__new__(JessQuery)
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
