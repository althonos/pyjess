from .region cimport Region


cdef extern from "Annulus.h" nogil:
    Region* Annulus_create(double*, double, double, int)
    void Annulus_free(Region*)