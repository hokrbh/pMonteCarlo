from distutils.core import setup, Extension

mcml = Extension('pMonteCarlo', sources = ['mcmlModule.c','mcmlSingle.c','allocate.c','hybridTaus.c','vector.c'], extra_compile_args = ['-std=c11'])

setup (name = 'pMonteCarlo',
        version = '1.0',
        description = 'A Monte Carlo library for python',
        author = 'Brett H. Hokr',
        author_email = 'brett.hokr@gmail.com',
        ext_modules = [mcml])
