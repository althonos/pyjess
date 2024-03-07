cimport cython

from libc.stdio cimport FILE, fdopen, fclose, printf
from libc.stdlib cimport free 
from libc.math cimport exp

cimport jess.jess
cimport jess.molecule
cimport jess.tess_template
cimport jess.super
from jess.atom cimport Atom as _Atom
from jess.jess cimport Jess as _Jess, JessQuery as _JessQuery
from jess.molecule cimport Molecule as _Molecule
from jess.template cimport Template as _Template
from jess.super cimport Superposition as _Superposition

import os

@cython.no_gc_clear
cdef class Template:
    cdef object     owner
    cdef _Template* _tpl

    def __cinit__(self):
        self._tpl  = NULL
        self.owner = None

    def __dealloc__(self):
        # if self.owner is None:
            # self._tpl.free(self._tpl)
        pass

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

    @property
    def name(self):
        assert self._tpl is not NULL
        return self._tpl.name(self._tpl).decode()


cdef class Molecule:
    cdef _Molecule* _mol

    def __cinit__(self):
        self._mol = NULL

    def __dealloc__(self):
        jess.molecule.Molecule_free(self._mol)

    def __init__(self, object file, bint ignore_endmdl = False):
        """Create a new molecule by loading the given file.
        """
        cdef int   fd = os.open(file, os.O_RDONLY)
        cdef FILE* f  = fdopen(fd, "r")

        try:
            self._mol = jess.molecule.Molecule_create(f, ignore_endmdl)
            if self._mol is NULL:
                raise ValueError(f"Failed to parse molecule file from {file!r}")
        finally:
            fclose(f)

    @property
    def id(self):
        assert self._mol is not NULL
        return jess.molecule.Molecule_id(self._mol).decode()


cdef class JessQuery:
    cdef _JessQuery* _jq
    cdef int         _candidates

    cdef readonly Jess     jess

    cdef readonly Molecule molecule    
    cdef readonly bint     ignore_chain
    cdef readonly double   rmsd_threshold

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

        while jess.jess.JessQuery_next(self._jq, self.ignore_chain) and self._candidates < 200:

            self._candidates += 1

            tpl = jess.jess.JessQuery_template(self._jq)
            count = tpl.count(tpl)
            sup = jess.jess.JessQuery_superposition(self._jq)
            atoms = jess.jess.JessQuery_atoms(self._jq)

            rmsd = jess.super.Superposition_rmsd(sup)
            if rmsd <= self.rmsd_threshold:

                hit = Hit.__new__(Hit)
                hit.rmsd = rmsd
                hit._sup = sup

                hit.molecule = self.molecule
                hit.template = Template.__new__(Template)
                hit.template._tpl = tpl
                hit.template.owner = self.jess

                return hit

            # 
            jess.super.Superposition_free(sup)

        raise StopIteration


cdef class Hit:
    cdef _Superposition* _sup

    cdef readonly double   rmsd
    cdef readonly Template template
    cdef readonly Molecule molecule

    def __dealloc__(self):
        jess.super.Superposition_free(self._sup)

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
            # FIXME: copy templates here so that the Jess storage owns the data
            jess.jess.Jess_addTemplate(self._jess, template._tpl.copy(template._tpl))

    cpdef JessQuery query(
        self, 
        Molecule molecule, 
        double rmsd_threshold, 
        double distance_cutoff,
        double max_total_threshold,
        bint transform = True,
        bint ignore_chain = False,
    ):
        cdef JessQuery query = JessQuery.__new__(JessQuery)
        query.ignore_chain = ignore_chain
        query.rmsd_threshold = rmsd_threshold
        query.molecule = molecule
        query.jess = self
        query._jq = jess.jess.Jess_query(
            self._jess, 
            molecule._mol, 
            distance_cutoff, 
            max_total_threshold
        )
        return query
