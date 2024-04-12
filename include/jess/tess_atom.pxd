from .atom cimport Atom

cdef extern from "TessAtom.h" nogil:

    struct _TessAtom:
        int code
        int resSeq
        int nameCount
        int resNameCount
        char chainID1
        char chainID2
        char **name
        char **resName
        double pos[3]
        double distWeight

    ctypedef _TessAtom TessAtom

    TessAtom* TessAtom_create(const char*)
    void TessAtom_free(TessAtom*)
    const double* TessAtom_position(const TessAtom*)
    int TessAtom_match(const TessAtom*, const Atom*)
    int TessAtom_resSeq(const TessAtom*)
    char TessAtom_chainID1(const TessAtom*)
    char TessAtom_chainID2(const TessAtom*)
    double TessAtom_distWeight(const TessAtom*)
    
    TessAtom* TessAtom_copy(const TessAtom*)
    double TessAtom_distance(const TessAtom*, const TessAtom*) 