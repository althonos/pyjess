from .atom cimport Atom

cdef extern from "Template.h" nogil:

    ctypedef _Template Template
    cdef struct _Template:
        void (*free)(Template*)
        int (*count)(const Template*)
        int (*match)(const Template*, int, const Atom*)
        int (*range)(const Template*, int, int, double*, double*)
        int (*check)(const Template*, Atom**, int, int)
        const double* (*position)(const Template*, int)
        const char* (*name)(const Template*)
        double (*logE)(const Template*, double, int)
        double (*distWeight)(const Template*, int)
