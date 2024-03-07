from .atom cimport Atom


cdef extern from "Template.h" nogil:

    ctypedef _Template Template
    cdef struct _Template:
        void (*free)(Template*) nogil
        int (*count)(const Template*) nogil
        int (*match)(const Template*, int, const Atom*) nogil
        int (*range)(const Template*, int, int, double*, double*) nogil
        int (*check)(const Template*, Atom**, int, int) nogil
        const double* (*position)(const Template*, int) nogil
        const char* (*name)(const Template*) nogil
        double (*logE)(const Template*, double, int) nogil
        double (*distWeight)(const Template*, int) nogil
        Template* (*copy)(const Template*) nogil
