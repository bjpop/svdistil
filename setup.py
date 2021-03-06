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
        'console_scripts': ['svdistil = svdistil.svdistil:main',
            'svqualfilter = svdistil.svqualfilter:main',
            'svmerge = svdistil.svmerge:main',
            'cnvmerge = svdistil.cnvmerge:main',
            'svannotate = svdistil.svannotate:main',
            'cnvannotate = svdistil.cnvannotate:main',
            'snvdistil = svdistil.snvdistil:main',
            'snvannotate = svdistil.snvannotate:main',
            'snvinsight = svdistil.snvinsight:main',
            'snvhsf = svdistil.snvhsf:main',
            'snvfilter = svdistil.snvfilter:main',
            'snvgene = svdistil.snvgene:main',
            ]
    },
    url='https://github.com/bjpop/svdistil',
    license='LICENSE',
    description=('Convert DNA structural variants in VCF files into BED format'),
    long_description=(LONG_DESCRIPTION),
    install_requires=["cyvcf2", "networkx", "quicksect", "numpy", "scipy", "plotly"],
)
