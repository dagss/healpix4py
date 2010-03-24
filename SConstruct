import os
import numpy
import distutils.sysconfig

if not os.environ.has_key('HEALPIX'):
    raise ValueError("HEALPIX environment variable must be set")

healpix_template_builder_dd = Builder(
    action="sed -e 's/KLOAD/d/g' -e 's/KMAP/DP/g' "
    "-e 's/KALMC/DPC/g'  -e 's/KALM/DP/g' $SOURCE > $TARGET")
healpix_template_builder_ss = Builder(
    action="sed -e 's/KLOAD/s/g' -e 's/KMAP/SP/g' "
    "-e 's/KALMC/SPC/g'  -e 's/KALM/SP/g' $SOURCE > $TARGET")

ip = '%s/include' % os.environ['HEALPIX']
if not os.path.isdir(ip):
    ip = '%s/include_ifort' % os.environ['HEALPIX']
    if not os.path.isdir(ip):
        raise Exception()

env = Environment(
    FORTRAN="ifort",
    F90="ifort",
    F90PATH=[ip],
    FORTRANFLAGS=["-g", "-vec_report0"],
    F90FLAGS=["-O3", '-openmp', '-parallel', '-cm', '-w', '-sox'],
    PYEXT_USE_DISTUTILS=True,
    BUILDERS=dict(HealpixTemplateDD=healpix_template_builder_dd,
                  HealpixTemplateSS=healpix_template_builder_ss))

env.Tool("pyext")
env.Tool("cython")

env.Append(PYEXTINCPATH=[numpy.get_include()])
env.Replace(PYEXTCFLAGS=['-fno-strict-aliasing', '-DNDEBUG', '-Wall',
                         '-fwrapv', '-g', '-Wstrict-prototypes'],#, '-DCYTHON_REFNANNY'],
            #CYTHON="python /uio/arkimedes/s07/dagss/cython/stable/cython.py",
            CYTHONFLAGS=['-a'])

env['ENV']['PATH'] = os.environ['PATH']
env['ENV']['LIBRARY_PATH'] = env['ENV']['LD_LIBRARY_PATH'] = os.environ['LD_LIBRARY_PATH']

Export('env')

SConscript(['healpix/SConscript'])

