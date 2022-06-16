from setuptools import setup, find_packages
from distutils.extension import Extension
from pathlib import Path

import numpy

# get the information about numpy header files
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

# handle cython modules: pyqcprot module for the rotation matrix
#try:
#    from Cython.Distutils import build_ext
#    use_cython = True
#    cmdclass = {'build_ext': build_ext}
#except ImportError:
use_cython = False
cmdclass = {}
#finally:
#    print (f'use_cython: {use_cython}')

ext_modules = [Extension("ties/pyqcprot", [f"ties/pyqcprot/pyqcprot.{'pyx' if use_cython else 'c'}"],
                         include_dirs=[numpy_include],
                         extra_compile_args=["-O3","-ffast-math"])]

setup(
    name='ties',
    version='0.0.1.dev1',
    description='TIES: Thermodynamic Integratin wtih Enhanced Sampling',
    long_description='Copy from README file',
    url='http://ccs.chem.ucl.ac.uk',
    author='Mateusz K. Bieniek',
    author_email='bieniekmat@gmail.com',
    packages=find_packages(),
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'ties = ties:cmdties.command_line_script'
        ]
    },
    cmdclass = cmdclass,
    ext_modules = ext_modules,
)
