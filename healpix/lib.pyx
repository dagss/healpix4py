from __future__ import division
##############################################################################
#    Copyright (C) 2010 Dag Sverre Seljebotn <dagss@student.matnat.uio.no>
#  Distributed under the terms of the GNU General Public License (GPL),
#  either version 2 of the License, or (at your option) any later version.
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################

cimport numpy as np
import numpy as np

np.import_array()

nocopy = False

def log(*args):
    import sys
    sys.stderr.write("%s\n"%repr(args))

def must_copy(arr):
    while len(arr.shape) > 0 and arr.shape[0] == 1:
        arr = arr[0,:]
    while len(arr.shape) > 0 and arr.shape[-1] == 1:
        arr = arr[:,0]
    result = not arr.flags.f_contiguous
    if result and nocopy:
        raise ValueError("Copying disabled but have to copy!")
    return result

real_dtype = np.double
complex_dtype = np.complex

cdef extern from *:
    ctypedef int size_t

cdef struct complex_double:
    double real
    double imag

#ctypedef np.int32_t i4b

ctypedef struct planck_rng:
    i4b x, y, z, w
    double gset
    i4b empty

cdef extern:
    np.int32_t nside2npix_ "pix_tools_mp_nside2npix_"(np.int32_t*)
    np.int32_t npix2nside_ "pix_tools_mp_npix2nside_"(np.int32_t*)

    void cywrap_vec2ang_(double* vector, double* theta, double* phi)

    void cywrap_pix2vec_nest_(i4b* nside, i4b* ipix, double* vector)
    void cywrap_pix2vec_ring_(i4b* nside, i4b* ipix, double* vector)
    void cywrap_vec2pix_ring_(i4b* nside, double* vector, i4b* ipix)
    
    void alm2map_sc_d_ "alm_tools_mp_alm2map_sc_d_"(
        i4b* nsmax,
        i4b* nlmax,
        i4b* nmmax,
        complex_double* alm,
        double* map) nogil
    
    void map2alm_sc_d_ "alm_tools_mp_map2alm_sc_d_"(
        i4b* nsmax,
        i4b* nlmax,
        i4b* nmmax,
        double* map,
        complex_double* alm,
        double* zbound,
        double* w8ring) nogil

    void cywrap_output_map_d_(
        double* map,
        char* outfile,
        int* shapes)

    void cywrap_convert_ring2nest_d_(
        i4b* nside,
        double* map,
        i4b* nd)

    void cywrap_convert_nest2ring_d_(
        i4b* nside,
        double* map,
        i4b* nd)

    void cywrap_alms2fits_(
        char* filename,
        i4b* nalms,
        double* alms,
        i4b* ncl,
        i4b* next,
        i4b* shapes)

    void cywrap_fits2alms_(
        char* filename,
        i4b* nalms,
        double* alms,
        i4b* ncl,
        i4b* next,
        i4b* shapes)

    void cywrap_sub_udgrade_nest_d_(double* map_in, i4b* nside_in, double* map_out, i4b* nside_out)


    void cywrap_read_unformatted_2d_complex_d_(double* data, i4b* n, i4b* m,
                                               char* filename, i4b* filename_len)
                                               
    void cywrap_read_unformatted_1d_real_d_(double* data, i4b* n, char* filename, i4b* filename_len)


    void cywrap_rotate_alm_d_(i4b* lmax, double* alm, double* psi, double* theta, double* phi,
                              i4b* alm_s0, i4b* alm_s1, i4b* alm_s2) nogil

    void cywrap_remove_dipole_double_(i4b* nside, double* map, i4b* ordering,
                                      i4b* degree, double* multipoles,
                                      double* mask)
    
cpdef np.int32_t nside2npix(np.int32_t nside):
    return nside2npix_(&nside)

cpdef np.int32_t npix2nside(np.int32_t npix):
    return npix2nside_(&npix)

cpdef int alm2map_sc_d(np.int32_t nsmax,
                       np.int32_t nlmax,
                       np.int32_t nmmax,
                       alm,
                       map) except -1:
    cdef np.ndarray[complex_double, ndim=3] alm_work = alm
    cdef np.ndarray[double, ndim=1] map_work = map
    assert map.shape[0] == 12*nsmax*nsmax
    assert (alm.shape[0] == 1 and
            alm.shape[1] == nlmax+1 and
            alm.shape[2] == nmmax+1)

    if must_copy(alm):
        alm_work = alm.copy('F')
    if must_copy(map):
        map_work = map.copy('F')
    with nogil:
        alm2map_sc_d_(&nsmax,
                      &nlmax,
                      &nmmax,
                      <complex_double*>alm_work.data,
                      <double*>map_work.data)
    if map_work is not map:
        map[...] = map_work
    return 0

cpdef int map2alm_sc_d(np.int32_t nsmax,
                       np.int32_t nlmax,
                       np.int32_t nmmax,
                       map,
                       alm,
                       weight_ring) except -1:
    cdef np.ndarray[double, ndim=1] map_work = map
    cdef np.ndarray[complex_double, ndim=3] alm_work = alm
    cdef np.ndarray[double, ndim=1] weight_ring_work = weight_ring
    assert map.shape[0] == 12*nsmax*nsmax
    assert alm.shape[0] == 1
    assert alm.shape[1] == nlmax+1
    assert alm.shape[2] == nmmax+1
    assert weight_ring.shape[0] == 2 * nsmax

    if must_copy(alm):
        alm_work = alm.copy('F')
    if must_copy(map):
        map_work = map.copy('F')
    if must_copy(weight_ring):
        weight_ring_work = weight_ring.copy('F')
    cdef np.ndarray zbounds = np.array([-1, 1], real_dtype)
    with nogil:
        map2alm_sc_d_(&nsmax,
                      &nlmax,
                      &nmmax,
                      <double*>map_work.data,
                      <complex_double*>alm_work.data,
                      #                  NULL, NULL)
                      <double*>zbounds.data,
                      <double*>weight_ring_work.data)
    if alm_work is not alm:
        alm[...] = alm_work
    return 0

cpdef int output_map_d(object map_, bytes filename) except -1:
    cdef np.ndarray[double, ndim=2] map_work
    if not must_copy(map_):
        map_work = map_
    else:
        map_work = map_.copy('F')

    print map_work.shape[0], map_work.shape[1]
    cdef int* mapshape = [map_work.shape[0],
                          map_work.shape[1],
                          len(filename)]
    cywrap_output_map_d_(<double*>map_work.data,
                      <char*>filename, mapshape)
    return 0

cpdef int convert_ring2nest_d(i4b nside, map_, i4b nd=-1) except -1:
    cdef np.ndarray[double, ndim=2] map_work
    assert map_.shape[0] == 12*nside*nside
    if not must_copy(map_):
        map_work = map_
    else:
        map_work = map_.copy('F')
    if nd == -1:
        nd = map_work.shape[1]

    cywrap_convert_ring2nest_d_(&nside, <double*>map_work.data,
                                &nd)
    if map_work is not map_:
        map_[...] = map_work
    
    return 0

cpdef int alms2fits(bytes filename, alms) except -1:
    cdef i4b nalms, ncl, next
    cdef i4b* shapes = [len(filename)]

    nalms, ncl, next = alms.shape
    ncl -= 1
    cdef np.ndarray alms_work
    if must_copy(alms):
        alms_work = alms.copy('F')
        
    print ncl
    cywrap_alms2fits_(<char*>filename,
                      &nalms,
                      <double*>alms_work.data,
                      &ncl,
                      &next,
                      shapes)                      
    return 0

cpdef int fits2alms(bytes filename, alms) except -1:
    cdef i4b nalms, ncl, next
    cdef i4b* shapes = [len(filename)]
    cdef np.ndarray alms_work

    nalms, ncl, next = alms.shape
    ncl -= 1

    if must_copy(alms):
        alms_work = alms.copy('F')
    else:
        alms_work = alms

    cywrap_fits2alms_(<char*>filename,
                      &nalms,
                      <double*>alms_work.data,
                      &ncl,
                      &next,
                      shapes)

    if alms_work is not alms:
        alms[...] = alms_work
    
    return 0

def pix2vec_nest(nside, ipix):
    """
    >>> pix2vec_nest(16, 32)
    array([ 0.81322541,  0.54337985,  0.20833333])
    >>> pix2vec_nest(16, 40000)
    Traceback (most recent call last):
        ...
    ValueError: ipix too large
    >>> x = pix2vec_nest([32,64], [[40],[400]]); x
    array([[[ 0.79462058,  0.58933079,  0.14583333],
            [ 0.75519318,  0.65143412,  0.07291667]],
    <BLANKLINE>
           [[ 0.37563625,  0.702766  ,  0.60416667],
            [ 0.56786914,  0.7656829 ,  0.30208333]]])
    >>> x.shape
    (2, 2, 3)
    """
    nside = np.asarray(nside, np.int32)
    ipix = np.asarray(ipix, np.int32)
    cdef np.broadcast input_it = np.broadcast(nside, ipix)
    out_shape = (<object>input_it).shape + (3,)
    out = np.empty(out_shape, dtype=np.double)
    cdef int vecaxis = len(out_shape) - 1
    cdef i4b nside_val, ipix_val
    cdef double* vector
    cdef np.flatiter out_it = np.PyArray_IterAllButAxis(out, &vecaxis)
    import sys
    assert out.strides[vecaxis] == sizeof(double)
    while np.PyArray_ITER_NOTDONE(out_it):
        nside_val = (<i4b*>np.PyArray_MultiIter_DATA(input_it, 0))[0]
        ipix_val = (<i4b*>np.PyArray_MultiIter_DATA(input_it, 1))[0]
        if ipix_val >= nside2npix(nside_val):
            raise ValueError("ipix too large")
        vector = <double*>np.PyArray_ITER_DATA(out_it)
        
        cywrap_pix2vec_nest_(&nside_val, &ipix_val, vector)

        np.PyArray_ITER_NEXT(out_it)
        np.PyArray_MultiIter_NEXT(input_it)
    assert not np.PyArray_MultiIter_NOTDONE(input_it)
    return out
    
def pix2vec_ring(nside, ipix):
    """
    >>> pix2vec_ring(16, 32)
    array([-0.19915651, -0.03961469,  0.97916667])
    >>> pix2vec_ring(16, 40000)
    Traceback (most recent call last):
        ...
    ValueError: ipix too large
    >>> x = pix2vec_ring([32,64], [[40],[400]]); x
    array([[[ 0.12575028,  0.01991689,  0.99186198],
            [ 0.0629714 ,  0.00997369,  0.99796549]],
    <BLANKLINE>
           [[-0.20338749, -0.28664785,  0.93619792],
            [-0.10294272, -0.14508418,  0.98404948]]])
    >>> x.shape
    (2, 2, 3)
    """
    nside = np.asarray(nside, np.int32)
    ipix = np.asarray(ipix, np.int32)
    cdef np.broadcast input_it = np.broadcast(nside, ipix)
    out_shape = (<object>input_it).shape + (3,)
    out = np.empty(out_shape, dtype=np.double)
    cdef int vecaxis = len(out_shape) - 1
    cdef i4b nside_val, ipix_val
    cdef double* vector
    cdef np.flatiter out_it = np.PyArray_IterAllButAxis(out, &vecaxis)
    import sys
    assert out.strides[vecaxis] == sizeof(double)
    while np.PyArray_ITER_NOTDONE(out_it):
        nside_val = (<i4b*>np.PyArray_MultiIter_DATA(input_it, 0))[0]
        ipix_val = (<i4b*>np.PyArray_MultiIter_DATA(input_it, 1))[0]
        if ipix_val >= nside2npix(nside_val):
            raise ValueError("ipix too large")
        vector = <double*>np.PyArray_ITER_DATA(out_it)
        
        cywrap_pix2vec_ring_(&nside_val, &ipix_val, vector)

        np.PyArray_ITER_NEXT(out_it)
        np.PyArray_MultiIter_NEXT(input_it)
    assert not np.PyArray_MultiIter_NOTDONE(input_it)
    return out

def vec2pix_ring(i4b nside, np.ndarray[double, ndim=2] vecs):
    """
    >>> vec2pix_ring(16, (-0.19915651, -0.03961469,  0.97916667))
    32
    """
    cdef Py_ssize_t i
    cdef i4b ipix
    cdef np.ndarray[i4b] out = np.zeros(vecs.shape[0], np.int32)
    cdef double* vec = [0, 0, 0]

    for i in range(vecs.shape[0]):
        vec[0] = vecs[i, 0]
        vec[1] = vecs[i, 1]
        vec[2] = vecs[i, 2]
        cywrap_vec2pix_ring_(&nside, vec, &ipix)
        out[i] = ipix
    return out

## def vec2pix_ring(i4b nside, double x, double y, double z):
##     """
##     >>> vec2pix_ring(4, 0, 0, 1)
    
##     """
##     cdef double* vec = [x, y, z]
##     cdef i4b ipix
##     cywrap_vec2pix_ring_(&nside, vec, &ipix)
##     return ipix


def convert_ring2nest(i4b nside, np.ndarray[double, ndim=2, mode='fortran'] map):
    cdef i4b nd = map.shape[1]
    if map.shape[0] != 12 * nside**2:
        raise ValueError("Wrong number of pixels for given nside")
    cywrap_convert_ring2nest_d_(&nside, <double*>map.data, &nd)

def convert_nest2ring(i4b nside, np.ndarray[double, ndim=2, mode='fortran'] map):
    cdef i4b nd = map.shape[1]
    if map.shape[0] != 12 * nside**2:
        raise ValueError("Wrong number of pixels for given nside")
    cywrap_convert_nest2ring_d_(&nside, <double*>map.data, &nd)

#
# comp_normalised_Plm
#
cdef extern:
    void comp_normalised_plm_(i4b* nlmax, i4b* m, double* theta, double* lambda_)

cpdef comp_normalised_Plm(i4b nlmax, i4b m, double theta, out=None):
    """
    >>> out = np.zeros(5, dtype=np.double)
    >>> lambda_ = comp_normalised_Plm(4, 0, np.pi/2, out)
    >>> Plm = lambda_ / np.sqrt((2*np.r_[0:5] + 1) / 4.0 / np.pi)
    >>> ' '.join(['%.2f' % x for x in Plm])
    '1.00 0.00 -0.50 -0.00 0.38'
    """
    cdef np.ndarray out_work = out#[double, ndim=1] out_work
#    if out is None:
#        out_work = np.empty(nlmax+1, dtype=np.double)
#    else:
#        out_work = np.asfortranarray(out)
    comp_normalised_plm_(&nlmax, &m, &theta, <double*>out_work.data)
#    if out is not out_work:
#        out[...] = out_work[...]
#    return out

cdef extern from "stdlib.h":
    void* malloc(size_t)

def vec2ang(double x, double y, double z):
    """
    >>> vec2ang(0, 1, np.pi/2)
    >>> vec2ang(1, 0, 0)
    >>> vec2ang(0, 1, 0)
    """
    
    cdef double* vector = [x, y, z]
    cdef double theta, phi
    cywrap_vec2ang_(vector, &theta, &phi)
    return (theta, phi)

def sub_udgrade_nest(np.ndarray[double, mode='fortran'] map_in,
                     np.ndarray[double, mode='fortran'] map_out):
    cdef i4b nside_in = npix2nside(map_in.shape[0]), nside_out = npix2nside(map_out.shape[0])
    cywrap_sub_udgrade_nest_d_(<double*>map_in.data, &nside_in,
                               <double*>map_out.data, &nside_out)

def read_unformatted_2d_complex(np.ndarray[np.complex128_t, mode='fortran', ndim=2] data, filename):
    cdef bytes filename_bytes = filename.encode('ASCII')
    cdef char* filename_buf = filename_bytes
    cdef i4b n = data.shape[0], m=data.shape[1], flen = len(filename_bytes)
    cywrap_read_unformatted_2d_complex_d_(<double*>data.data, &n, &m, filename_buf, &flen)

def read_unformatted_1d_real(np.ndarray[np.float64_t, mode='fortran', ndim=1] data, filename):
    cdef bytes filename_bytes = filename.encode('ASCII')
    cdef char* filename_buf = filename_bytes
    cdef i4b n = data.shape[0], flen = len(filename_bytes)
    cywrap_read_unformatted_1d_real_d_(<double*>data.data, &n, filename_buf, &flen)

def rotate_alm_d(i4b lmax, np.ndarray[np.complex128_t, ndim=3, mode='fortran'] alm,
                 double psi, double theta, double phi):
    cdef i4b i, j, k
    i, j, k = alm.shape[0], alm.shape[1], alm.shape[2]
    if single is None:
        single = False
    with nogil:
        cywrap_rotate_alm_d_(&lmax, <double*>alm.data, &psi, &theta, &phi, &i, &j, &k)

def remove_dipole_double(i4b nside,
                         np.ndarray[np.double_t, ndim=1, mode='fortran'] map,
                         ordering,
                         i4b degree,
                         np.ndarray[np.double_t, ndim=1, mode='fortran'] mask):
    cdef double* multipoles = [0, 0, 0, 0]
    cdef i4b npix = nside2npix(nside)
    cdef i4b orderingint
    if map.shape[0] != npix or mask.shape[0] != npix:
        raise ValueError("Arrays do not match nside")
    if ordering == 'ring':
        orderingint = 1
    elif ordering == 'nested':
        orderingint = 2
    else:
        raise ValueError("ordering must be ring or nested")
    if degree < 0 or degree > 2:
        raise ValueError("illegal degree")
    cywrap_remove_dipole_double_(&nside, <double*>map.data, &orderingint,
                                 &degree, multipoles, <double*>mask.data)
    if degree == 0:
        return ()
    elif degree == 1:
        return (multipoles[0],)
    elif degree == 2:
        return (multipoles[0], multipoles[1], multipoles[2], multipoles[3])

#cdef extern from *:
#    void doit_()

#def doit_yeah():
#    
