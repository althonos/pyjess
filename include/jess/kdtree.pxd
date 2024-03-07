from .region cimport Region


cdef extern from "KdTree.h" nogil:

    cdef struct _KdTree:
        pass
    cdef struct _KdTreeQuery:
        pass

    ctypedef _KdTree KdTree
    ctypedef _KdTreeQuery KdTreeQuery

    KdTree* KdTree_create(double**, int, int)
    void KdTree_free(KdTree**)
    KdTreeQuery* KdTree_query(KdTree*, Region*)

    void KdTreeQuery_free(KdTreeQuery*)
    int KdTreeQuery_next(KdTreeQuery*)