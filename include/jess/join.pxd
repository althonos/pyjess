from .region cimport Region


cdef extern from "Join.h" nogil:

    ctypedef enum JoinType:
        innerJoin
        outerJoin

    Region* Join_create(Region**, int, JoinType)
    void Join_free(Region*)

