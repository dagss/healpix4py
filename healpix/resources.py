"""
Load resources from the Healpix library (pixel window and quadrature
ring weights).

The resources are cached in memory.
"""
import pyfits
import os
import numpy as np

class HealpixData(object):
    def __init__(self, data_path):
        self.data_path = data_path
        self.weight_ring_cache = {}
        self.pixel_window_cache = {}

    def weight_ring(self, Nside):
        """
        Returns the ring weights to use. NOTE that 1 is automatically
        added; these weights can be passed directly to HEALPix.

        >>> res = get_default()
        >>> x = res.weight_ring(Nside=8)
        >>> x.shape
        (3, 16)
        >>> x[0,0]
        1.16567683090317
        >>> y = res.weight_ring(Nside=8)
        >>> y is x
        True

        FITS files don't always have the same format, so try for different
        Nside::

            >>> res.weight_ring(Nside=512).shape
            (3, 1024)

        """
        result = self.weight_ring_cache.get(Nside)
        if result is None:
            hdulist = pyfits.open(os.path.join(self.data_path,
                                               'weight_ring_n%05d.fits' % Nside))
            try:
                # Get data to plain 2D float array
                data = hdulist[1].data.view(np.ndarray)
                temp = data['TEMPERATURE WEIGHTS'].ravel()
                qpol = data['Q-POLARISATION WEIGHTS'].ravel()
                upol = data['U-POLARISATION WEIGHTS'].ravel()
                # Convert to native endian...
                data = np.asarray([temp, qpol, upol], np.double, order='C')
                # Add 1
                data += 1
                result = data
                self.weight_ring_cache[Nside] = result
            finally:
                hdulist.close()
        return result

    def pixel_window(self, Nside, dtype=np.double):
        """
        >>> res = get_default()
        >>> t, p = res.pixel_window(Nside=8)
        >>> t.shape, p.shape
        ((33,), (33,))
        >>> t.dtype
        dtype('float64')
        >>> np.round(t[0], 4)
        1.0
        >>> t2, p2 = res.pixel_window(Nside=8)
        >>> t is t2, p is p2
        (True, True)

        >>> res.pixel_window(Nside=512)[0].shape
        (2049,)

        """
        result = self.pixel_window_cache.get(Nside)
        if result is None:
            hdulist = pyfits.open(os.path.join(self.data_path,
                                               'pixel_window_n%04d.fits' % Nside))
            try:
                # Important to use astype, in order to convert to native endia
                data = hdulist[1].data
                result = (data.field('TEMPERATURE').astype(dtype),
                          data.field('POLARIZATION').astype(dtype))
                self.pixel_window_cache[Nside] = result
            finally:
                hdulist.close()
        return result
        


_default = None
def get_default():
    global _default
    if _default is None:
        import os
        try:
            path = os.environ['HEALPIX']
        except KeyError:
            raise RuntimeError('"HEALPIX" environment variable not set; cannot get default '
                               'HEALPix data')
        path = os.path.join(path, 'data')
        if not os.path.isdir(path):
            raise RuntimeError('"%s" is not a directory' % path)
        _default = HealpixData(path)
    return _default
