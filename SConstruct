import os
import numpy
import distutils.sysconfig
from pprint import pprint

v = Variables('config')
v.Add('HEALPIX_INCLUDE')
v.Add('HEALPIX_LIB')
v.Add('F90')
v.Add('F90FLAGS')
v.Add('CFITSIO_LIB')
v.Add('PYEXTCFLAGS')

env = Environment(variables=v, PYEXT_USE_DISTUTILS=True)

env.Tool("pyext")
env.Tool("cython")
for x in ['PATH', 'INTEL_LICENSE_FILE', 'LIBRARY_PATH', 'LD_LIBRARY_PATH']:
    if x in os.environ:
        env['ENV'][x] = os.environ[x]

env.Append(PYEXTINCPATH=[numpy.get_include()],
           PYEXTCFLAGS=['-O0'],
           CYTHONFLAGS=['-a'],
           )

env['ENV']['PATH'] = os.environ['PATH']
if 'LD_LIBRARY_PATY' in os.environ:
    env['ENV']['LIBRARY_PATH'] = env['ENV']['LD_LIBRARY_PATH'] = os.environ['LD_LIBRARY_PATH']

Export('env')

SConscript(['healpix/SConscript'])

