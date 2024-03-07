cdef extern from "Super.h" nogil:

    cdef struct _Superposition:
        pass
    ctypedef _Superposition Superposition

    Superposition* Superposition_create()
    void Superposition_free(Superposition*)
    void Superposition_align(Superposition*, const double*, const double*)
    int Superposition_count(const Superposition*)
    double Superposition_rmsd(Superposition*)
    double Superposition_rmsd100(Superposition*)
    const double* Superposition_centroid(Superposition*, int)
    const double* Superposition_rotation(Superposition*)