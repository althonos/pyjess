from libc.stdio cimport FILE

from .atom cimport Atom


cdef extern from "Molecule.h" nogil:

    cdef struct _Molecule:
        int count
        char[5] id
        Atom** atom

    ctypedef _Molecule Molecule

    Molecule* Molecule_create(FILE*, int, double)
    void Molecule_free(Molecule*) 
    int Molecule_count(const Molecule*)
    Atom* Molecule_atom(const Molecule*, int)
    const char* Molecule_id(const Molecule*)