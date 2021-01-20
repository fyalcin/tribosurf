
import os

from setuptools import setup


if __name__ == "__main__":
    setup(
        name='triboflow',
        version='0.0.1',
        description='triboflow contains complex workflows for interface '
                    'tribology based on atomate and FireWorks.',
        long_description=open('README.md').read(),
        url='https://gitlab.com/triboteam/fireflow',
        author='Michael Wolloch, Gabriele Losi',
        author_email='michael.wolloch@univie.ac.at, gabriele.losi@unimore.it',
        license='Not public yet',
        install_requires=['atomate@git+https://github.com/hackingmaterials/atomate.git',
                          'mpinterfaces@git+https://github.com/henniggroup/MPInterfaces.git'],
        classifiers=["Programming Language :: Python :: 3",
                     "Programming Language :: Python :: 3.6",
                     "Programming Language :: Python :: 3.7",
                     'Development Status :: 3 - Alpha',
                     'Intended Audience :: Science/Research',
                     'Intended Audience :: System Administrators',
                     'Operating System :: OS Independent',
                     'Topic :: Other/Nonlisted Topic',
                     'Topic :: Scientific/Engineering'],
    )
