from .atom cimport Atom


cdef extern from "ResIndex.h" nogil:

    struct _ResIndex:
        int count
        Atom **atom
        double **coord

    ctypedef _ResIndex ResIndex

    extern ResIndex* ResIndex_create(Atom**, int)
    extern Atom** ResIndex_get(ResIndex*, const char[4])
    extern void ResIndex_free(ResIndex*)