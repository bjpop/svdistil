#!/usr/bin/env python

from distutils.core import setup

LONG_DESCRIPTION = \
'''Convert DNA structural variants in VCF files into BED format'''


setup(
    name='svdistil',
    version='0.1.0.0',
    author='Bernie Pope',
    author_email='bjpope@unimelb.edu.au',
    packages=['svdistil'],
    package_dir={'svdistil': 'svdistil'},
    entry_points={
        'console_scripts': ['svdistil = svdistil.svdistil:main']
    },
    url='https://github.com/bjpop/svdistil',
    license='LICENSE',
    description=('Convert DNA structural variants in VCF files into BED format'),
    long_description=(LONG_DESCRIPTION),
    install_requires=["cyvcf2"],
)
