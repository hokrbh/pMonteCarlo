from setuptools import setup, find_packages, Extension

mcml = Extension('mcml',
        sources = ['pmontecarlo/mcmlModule.c',\
                'pmontecarlo/mcmlSingle.c',\
                'pmontecarlo/allocate.c',\
                'pmontecarlo/hybridTaus.c',\
                'pmontecarlo/vector.c'],
        extra_compile_args = ['-std=c11'])

setup (name = 'pmontecarlo',
        version = '1.0',
        description = 'A Monte Carlo radiation transport library for python',
        author = 'Brett H. Hokr',
        author_email = 'brett.hokr@gmail.com',
        license = 'GNU GPL v2',
        #py_modules=['pmontecarlo/functions'],
        packages = find_packages(),
        #packages = ['pmontecarlo'],
        #ext_package = 'pmontecarlo/montecarlo',
        #ext_modules = [mcml]
        )
