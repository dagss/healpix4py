import os
import numpy
import distutils.sysconfig

if not os.environ.has_key('HEALPIX'):
    raise ValueError("HEALPIX environment variable must be set")

settings = Environment(HEALPIX=os.environ['HEALPIX'],
                       CFITSIO='/mn/corcaroli/d1/dagss/cfitsio')
#settings = Environment(HEALPIX='/home/dagss/local/Healpix_2.11c',
#                       CFITSIO='/home/dagss/local/cfitsio')

env = settings.Environment(
    FORTRAN="ifort",
    F90="ifort",
    FORTRANFLAGS=["-g", "-vec_report0"],
    PYEXT_USE_DISTUTILS=True)


env.Tool("pyext")
env.Tool("cython")

env.Append(LIBPATH=[settings.subst('$HEALPIX/lib'), settings.subst('$CFITSIO')],
           F90PATH=[settings.subst('$HEALPIX/include')],
           PYEXTINCPATH=[numpy.get_include()])
env.Replace(PYEXTCFLAGS=['-fno-strict-aliasing', '-DNDEBUG', '-Wall',
                         '-fwrapv', '-g', '-Wstrict-prototypes'],#, '-DCYTHON_REFNANNY'],
            CYTHON="python /uio/arkimedes/s07/dagss/cython/devel/cython.py",
            CYTHONFLAGS=['-a'])
env['ENV']['PATH'] = os.environ['PATH']

Export('env')

SConscript(['healpix/SConscript'])

