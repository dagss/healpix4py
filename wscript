#------------------------------------------------------------------------------
# Copyright (c) 2010, Kurt W. Smith
# All rights reserved. See LICENSE.txt.
#------------------------------------------------------------------------------

top = '.'
out = 'build'

def options(ctx):
    for lib in ['healpix', 'cfitsio']:
        ctx.add_option('--%s' % lib, action='store', default=None)
        ctx.add_option('--%s-libpath' % lib, action='store')
        ctx.add_option('--%s-lib' % lib, action='store', default=lib)
        ctx.add_option('--%s-include' % lib, action='store')
    ctx.load('compiler_c')
    ctx.load('compiler_fc')
    ctx.load('python')

def configure(ctx):
    ctx.load('compiler_c')
    ctx.load('compiler_fc')
    ctx.check_fortran()
    ctx.check_fortran_verbose_flag()
    ctx.check_fortran_clib()

    ctx.load('python')
    ctx.check_python_version((2,5))
    ctx.check_python_headers()
    ctx.check_python_module('numpy')
    ctx.check_numpy_version(minver=(1,3))
    ctx.get_numpy_includes()

    ctx.find_program('cython', var='CYTHON')
    ctx.check_cython_version(minver=(0,11,1))

    for lib in ['healpix', 'cfitsio']:
        LIB = lib.upper()
        if getattr(ctx.options, lib) is not None:
            libpath = os.path.join(getattr(ctx.options, lib), 'lib')
            includepath = os.path.join(getattr(ctx.options, lib), 'include')
        else:
            libpath = getattr(ctx.options, '%s_lib' % lib)
            includepath = getattr(ctx.options, '%s_include' % lib)

        setattr(ctx.env, 'LIBPATH_%s' % LIB, libpath)
        setattr(ctx.env, 'INCLUDES_%s' % LIB, includepath)
        setattr(ctx.env, 'LIB_%s' % LIB, getattr(ctx.options, '%s_lib' % lib))

    ctx.add_os_flags('INCLUDES')
    ctx.add_os_flags('LIB')
    ctx.add_os_flags('LIBPATH')
    ctx.add_os_flags('STLIB')
    ctx.add_os_flags('STLIBPATH')

    ctx.check_healpix()

def build(ctx):
    ctx.recurse('healpix')



from waflib.Configure import conf
@conf
def check_numpy_version(conf, minver, maxver=None):
    conf.start_msg("Checking numpy version")
    minver = tuple(minver)
    if maxver: maxver = tuple(maxver)
    (np_ver_str,) = conf.get_python_variables(#conf.env['PYTHON'],
            ['numpy.version.short_version'], ['import numpy'])
    np_ver = tuple([int(x) for x in np_ver_str.split('.')])
    if np_ver < minver or (maxver and np_ver > maxver):
        conf.end_msg(False)
        conf.fatal("numpy version %s is not in the "
                "range of supported versions: minimum=%s, maximum=%s" % (np_ver_str, minver, maxver))
    conf.end_msg(str(np_ver))

@conf
def get_numpy_includes(conf):
    conf.start_msg("Checking numpy includes")
    (np_includes,) = conf.get_python_variables(#conf.env['PYTHON'],
            ['numpy.get_include()'], ['import numpy'])
    conf.env.INCLUDES_NUMPY = np_includes
    conf.end_msg('ok (%s)' % np_includes)

@conf
def get_healpix_includes(conf):
    conf.start_msg("Checking HEALPix includes")
    conf.end_msg("(TODO)")

@conf
def check_cython_version(conf, minver):
    conf.start_msg("Checking cython version")
    minver = tuple(minver)
    import re
    version_re = re.compile(r'cython\s*version\s*(?P<major>\d*)\.(?P<minor>\d*)(?:\.(?P<micro>\d*))?', re.I).search
    cmd = conf.cmd_to_list(conf.env['CYTHON'])
    cmd = cmd + ['--version']
    from waflib.Tools import fc_config
    stdout, stderr = fc_config.getoutput(conf, cmd)
    if stdout:
        match = version_re(stdout)
    else:
        match = version_re(stderr)
    if not match:
        conf.fatal("cannot determine the Cython version")
    cy_ver = [match.group('major'), match.group('minor')]
    if match.group('micro'):
        cy_ver.append(match.group('micro'))
    else:
        cy_ver.append('0')
    cy_ver = tuple([int(x) for x in cy_ver])
    if cy_ver < minver:
        conf.end_msg(False)
        conf.fatal("cython version %s < %s" % (cy_ver, minver))
    conf.end_msg(str(cy_ver))

@conf
def check_healpix(conf):
    from waflib.Errors import BuildError
    
    FRAG = '''\
program test
  use pix_tools
  use healpix_types
  real(dp) :: x
  if (nside2npix(32) .ne. 12288) then
    write (*,*) "ERROR"
    stop
  endif
end program
'''
    conf.start_msg('Checking for HEALPix')
    for flags in ['', '-fopenmp', '-liomp5']:
        conf.env.LINKFLAGS_HEALPIX = flags
        try:
            kw = dict(uselib='HEALPIX',
                      fragment=FRAG,
                      compile_filename = 'test.f90',
                      features         = 'fc fcprogram')
            conf.validate_c(kw)
            conf.run_c_code(**kw)
        except conf.errors.ConfigurationError:
            pass
        else:
            conf.end_msg('ok')
            return
    conf.fatal('not found, wrong FC, or wrong link flags')

    

import os
from waflib import Logs, Build, Utils

from waflib import TaskGen, Task

TaskGen.declare_chain(
        name = "cython",
        rule = "${CYTHON} ${CYTHONFLAGS} ${CPPPATH_ST:INCPATHS} ${SRC} -o ${TGT}",
        ext_in = ['.pyx'],
        ext_out = ['.c'],
        reentrant = True,
        )

# vim:ft=python
