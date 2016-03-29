from setuptools import setup, find_packages, Extension

mcml = Extension('mcml.mcml_run',
        sources = ['pmontecarlo/mcml/mcmlModule.c',\
                'pmontecarlo/mcml/mcmlSingle.c',\
                'pmontecarlo/allocate.c',\
                'pmontecarlo/hybridTaus.c',\
                'pmontecarlo/vector.c'],
        extra_compile_args = ['-std=c11'])
        
fmcml = Extension('fmcml.fmcml_run',
        sources = ['pmontecarlo/fmcml/fmcmlModule.c',\
                'pmontecarlo/fmcml/fmcmlSingle.c',\
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
        packages = find_packages(),
        ext_package = 'pmontecarlo',
        ext_modules = [mcml,fmcml]
        )
