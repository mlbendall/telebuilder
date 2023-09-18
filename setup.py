# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
import versioneer

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2023 Matthew L. Bendall"

setup(
    name='telebuilder',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description='Identify transposable element loci from repeat annotations',
    url='http://github.com/mlbendall/telebuilder',
    author='Matthew L. Bendall',
    author_email='bendall@gwu.edu',
    license='MIT',
    packages=find_packages(),
    install_requires=[
        'intervaltree',
        'biopython',
    ],
    entry_points={
        'console_scripts': [
            'buildERV=telebuilder.buildERV:console',
            'analyzeERV=telebuilder.analyzeERV:console',
            'buildL1=telebuilder.buildL1:console',
            'gtftools=telebuilder.gtftools:console',
        ],
    },
    zip_safe=False,
)

