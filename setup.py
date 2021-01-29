#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from setuptools import setup, find_packages

SOURCE_DIR = "src/mitohifi"
version = '0.1'

PROJECT_NAME = "MitoHiFi"
PROJECT_OWNER = PROJECT_USERAME = "marcelauliano"
PROJECT_URL = "https://github.com/marcelauliano/MitoHiFi"
PROJECT_AUTHOR = "Marcela Uliano-Silva"
PROJECT_EMAIL = "mu2@sanger.ac.uk"
RAW_CONTENT_URL = "https://raw.github.com/%s/%s/master/" % (
    PROJECT_USERAME, PROJECT_NAME
)

TEST_DIR = 'exampleFiles'
PROJECT_DESCRIPTION = '''MitoHiFi circularises, cuts and annotates the mitogenome from contigs assembled with
                         PacBio HiFi reads and softwares such as HiCanu or Hifiasm.'''
ENTRY_POINTS = '''
        [console_scripts]
        circularization-check=mitohifi.circularizationCheck_modified:main
        cut-coords=mitohifi.cut_coords:main
        filter-fasta=mitohifi.filterfasta:main
        get-larger=mitohifi.get_Larger:main
        parse-blast=mitohifi.parse_blast:main
        parse-blast-subject=mitohifi.parse_blast_subject:main
        rotate=mitohifi.rotate:main
        parse-blast-all-query=mitohifi.parse_blast_allQueryPerc:main
'''

PACKAGE_DATA = {
    # Be sure to update MANIFEST.in for source dist.
}
PACKAGE_DIR = {
    SOURCE_DIR: SOURCE_DIR,
}

readme = open('README.md').read()

if os.path.exists("requirements.txt"):
    requirements = open("requirements.txt").read().split("\n")
else:
    # In tox, it will cover them anyway.
    requirements = []


test_requirements = [
    # TODO: put package test requirements here
]


setup(
    name=PROJECT_NAME,
    version=version,
    description=PROJECT_DESCRIPTION,
    long_description=readme + '\n',
    author=PROJECT_AUTHOR,
    author_email=PROJECT_EMAIL,
    url=PROJECT_URL,
    packages=find_packages('src'),
    entry_points=ENTRY_POINTS,
    package_data=PACKAGE_DATA,
    package_dir={'': 'src'},
    include_package_data=True,
    install_requires=requirements,
    license="GPL-3",
    zip_safe=False,
    keywords='hifi',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'Environment :: Console',
        'License :: OSI Approved :: Academic Free License (AFL)',
        'Operating System :: POSIX',
        'Topic :: Software Development',
        'Topic :: Software Development :: Code Generators',
        'Topic :: Software Development :: Testing',
        'Natural Language :: English',
        "Programming Language :: Python :: 3",
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    test_suite=TEST_DIR,
    tests_require=test_requirements
)
