from healpix.lib import npix2nside

def nestedmap_to_fits(map, filename, map_units='K', output_units='mK'):
    import pyfits
    import os

    if map.ndim != 1:
        raise ValueError('Array must be 1D')
    Npix = map.shape[0]
    Nside = 12*Npix**2

    if map_units != output_units:
        # TODO: Do a proper output/input unit conversion, for now
        # we very stupidly only support one
        assert output_units == 'mK'
        if map_units == 'K':
            map = map * 1e3
        elif map_units == 'mK':
            pass
        elif map_units == 'uK':
            map = map * 1e6
        elif map_units == 'raw':
            pass
        else:
            raise ValueError('Illegal map unit')

    pri = pyfits.PrimaryHDU()
    col = pyfits.Column(name='TEMPERATURE', format='D',
                        unit='%s,thermodynamic' % output_units, array=map)
    sec = pyfits.new_table(pyfits.ColDefs([col]))
    sec.header.update('PIXTYPE', 'HEALPIX')
    sec.header.update('ORDERING', 'NESTED')
    sec.header.update('NSIDE', Nside)
    sec.header.update('FIRSTPIX', 0)
    sec.header.update('LASTPIX', Npix)
    if os.path.isfile(filename):
        os.unlink(filename)
    pyfits.HDUList([pri, sec]).writeto(filename)

def map2gif(map, outfile, bar=True, title=None, col=5,
                max=None, min=None):
    import os
    import tempfile
    fd, infile = tempfile.mkstemp(suffix='.fits')
    os.close(fd)
    nestedmap_to_fits(map, infile)
    if os.path.isfile(outfile):
        os.remove(outfile)
    flags = []
    if title is not None: flags.append('-ttl "%s"' % title)
    if max is not None: flags.append('-max %f' % max)
    if min is not None: flags.append('-min %f' % min)
    if bar: flags.append('-bar .true.')
    
    cmd = ('map2gif %s -col %d -inp "%s" -out "%s"' %
                (' '.join(flags), col, infile, outfile))
    os.system(cmd)
    os.remove(infile)



