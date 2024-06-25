from setuptools import setup, find_packages

setup(
    name='tethercad',
    version='0.1.0',
    author='Austen Goddu',
    py_modules=['tethercad'],
    author_email='austen.j.goddu@jpl.nasa.gov',
    packages=find_packages(include=['calculation_libraries', 'tether_analysis', 'databases']),
    install_requires=[
        'matplotlib>=3.7.1',
        'pandas>=1.5.3',
        'numpy>=1.24.2',
        'pytest>=7.2.2',
        'scipy>=1.10.0',
        'mpmath>=1.2.1'
    ],
    package_data={'databases':['*.csv']}
)