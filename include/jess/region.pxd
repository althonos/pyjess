cdef extern from "Region.h" nogil:

    struct _Region:
        int (*intersectionQ)(Region*, double*, double*, int)
        int (*inclusionQ)(Region*, double*, int)
        void (*free)(Region*)
    ctypedef _Region Region

    double Region_volume(Region*, double, double*, double*, int)