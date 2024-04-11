cimport cython
from cpython.unicode cimport PyUnicode_FromStringAndSize

from libc.math cimport exp
from libc.stdio cimport FILE, fclose, fdopen, printf
from libc.stdlib cimport calloc, free, malloc
from libc.string cimport memcpy

cimport jess.jess
cimport jess.molecule
cimport jess.super
cimport jess.tess_template
from jess.atom cimport Atom as _Atom
from jess.jess cimport Jess as _Jess
from jess.jess cimport JessQuery as _JessQuery
from jess.molecule cimport Molecule as _Molecule
from jess.super cimport Superposition as _Superposition
from jess.template cimport Template as _Template

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
cdef class Template:
    cdef object     owner
    cdef _Template* _tpl

    def __cinit__(self):
        self._tpl  = NULL
        self.owner = None

    def __dealloc__(self):
        if self.owner is None:
            self._tpl.free(self._tpl)

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
        finally:
            fclose(f)

    def __len__(self):
        assert self._tpl is not NULL
        return self._tpl.count(self._tpl)

    @property
    def name(self):
        assert self._tpl is not NULL
        return self._tpl.name(self._tpl).decode()


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
