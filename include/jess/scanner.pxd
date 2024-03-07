from .molecule import Molecule
from .template import Template


cdef extern from "Scanner.h" nogil:

    cdef struct _Scanner:
        pass
    ctypedef _Scanner Scanner

    Scanner* Scanner_create(Molecule*, Template*, double, double)
    void Scanner_free(Scanner*)
    Atom** Scanner_next(Scanner*, int)
    double Scanner_rmsd(Scanner*)

