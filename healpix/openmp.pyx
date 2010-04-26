"""
Simple interface to some OpenMP control routines.
"""
cdef extern from "omp.h":
    void omp_set_num_threads(int)
    int omp_get_num_threads()
    int omp_get_max_threads()
    double omp_get_wtime()
    double omp_get_wtick()

def set_num_threads(n):
    omp_set_num_threads(n)

def get_num_threads():
    return omp_get_num_threads()

def get_max_threads():
    return omp_get_max_threads()

def get_wtime():
    return omp_get_wtime()

def get_wtick():
    return omp_get_wtick()

def force_segfault():
    cdef int* i = NULL
    return i[0]
