from setuptools import setup, find_packages
from setuptools.command.develop import develop
from setuptools.command.install import install

class PostInstallCommand(install):
    def run(self):
        install.run(self)
        # importing ties prints a message about how to use it
        import ties

setup(
    name='ties',
    version='0.0.1.dev1',
    description='TIES: Thermodynamic Integratin wtih Enhanced Sampling',
    long_description='Copy from README file',
    url='http://ccs.chem.ucl.ac.uk',
    author='Mateusz K. Bieniek',
    author_email='bieniekmat@gmail.com',
    install_requires=['numpy', 'mdanalysis', 'cython', 'setuptools', 'matplotlib', 'networkx'],
    packages=find_packages(),
    include_package_data=True,
    cmdclass={
        'install': PostInstallCommand,
    },
)