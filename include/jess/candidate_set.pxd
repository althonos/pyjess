from .atom cimport Atom
from .molecule cimport Molecule


cdef extern from "CandidateSet.h" nogil:

    struct _CandidateSet:
        int count
        Atom **atom
        double **coord

    ctypedef _CandidateSet CandidateSet

    CandidateSet *CandidateSet_create(const Molecule *M)
    void CandidateSet_free(CandidateSet*)