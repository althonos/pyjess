cdef extern from "Atom.h" nogil:

    struct _Atom:
        int serial
        char[5] name
        char altLoc
        char[4] resName
        char chainID1
        char chainID2
        int resSeq
        char iCode
        double[3] x
        double occupancy
        double tempFactor
        char[4] segID
        char[3] element
        int charge

    ctypedef _Atom Atom

    int Atom_parse(Atom*, const char*)