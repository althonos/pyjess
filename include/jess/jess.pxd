from .atom cimport Atom
from .molecule cimport Molecule
from .super cimport Superposition
from .template cimport Template, IgnoreType


cdef extern from "Jess.h" nogil:

    ctypedef struct _Jess:
        pass

    ctypedef struct _JessQuery:
        pass

    ctypedef _Jess Jess
    ctypedef _JessQuery JessQuery

    Jess* Jess_create()
    void Jess_free(Jess*)
    void Jess_addTemplate(Jess*, Template*)
    JessQuery* Jess_query(Jess*, Molecule*, double, double, bint)

    void JessQuery_free(JessQuery*)
    int JessQuery_next(JessQuery*, IgnoreType)
    int JessQuery_nextTemplate(JessQuery*)
    Template* JessQuery_template(JessQuery*)
    const Molecule* JessQuery_molecule(JessQuery*)
    Atom** JessQuery_atoms(JessQuery*)
    Superposition* JessQuery_superposition(JessQuery*)


