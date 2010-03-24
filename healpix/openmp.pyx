"""
Simple interface to some OpenMP control routines.
"""
cdef extern from "omp.h":
    void omp_set_num_threads(int)
    double omp_get_wtime()
    double omp_get_wtick()

def set_num_threads(n):
    omp_set_num_threads(n)

def get_wtime():
    return omp_get_wtime()

def get_wtick():
    return omp_get_wtick()
