from libc.stdio cimport FILE

from .atom cimport Atom
from .tess_atom cimport TessAtom
from .template cimport Template


cdef extern from "TessTemplate.h" nogil:

    cdef struct _TessTemplate:
        int count
        TessAtom** atom
        double** distance
        char* symbol
        int dim

    ctypedef _TessTemplate TessTemplate

    int TessTemplate_count(const Template *T)
    int TessTemplate_match(const Template *T,int k,const Atom *A)
    int TessTemplate_range(const Template *T,int i,int j,double *a,double *b)
    const double *TessTemplate_position(const Template *T, int k)
    double TessTemplate_distWeight(const Template *T, int k)
    int TessTemplate_check(const Template *T, Atom **A, int k, int ignore_chain)
    const char *TessTemplate_name(const Template *T)
    double TessTemplate_logE(const Template *T,double rmsd, int n)
    void TessTemplate_free(Template *T)

    Template* TessTemplate_create(FILE*, const char*)
    Template* TessTemplate_copy(Template*)

