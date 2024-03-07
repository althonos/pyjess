from libc.stdio cimport FILE

from .template cimport Template

cdef extern from "TessTemplate.h" nogil:

    Template* TessTemplate_create(FILE*, const char*)
    Template* TessTemplate_copy(Template*)