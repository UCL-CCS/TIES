from setuptools import setup, find_packages

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
)
