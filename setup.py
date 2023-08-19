from setuptools import setup, find_packages
from distutils.extension import Extension

# pyqcprot module for the rotation matrix
import numpy
ext_modules = [Extension("ties.pyqcprotext.pyqcprot", ["ties/pyqcprotext/pyqcprot.pyx"],
                         include_dirs=[numpy.get_include()],
                         extra_compile_args=["-O3","-ffast-math"])]

setup(
    name='ties',
    version=open("ties/version.txt").read().strip(),
    description='TIES: Thermodynamic Integration with Enhanced Sampling',
    long_description='Copy from README file',
    url='http://ccs.chem.ucl.ac.uk',
    author='Mateusz K. Bieniek',
    author_email='bieniekmat@gmail.com',
    packages=find_packages(),
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'ties = ties:cli.command_line_script'
        ]
    },
    ext_modules = ext_modules,
)
