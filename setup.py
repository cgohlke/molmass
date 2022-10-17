# molmass/setup.py

"""Molmass package Setuptools script."""

import sys
import re

from setuptools import setup


def search(pattern, code, flags=0):
    # return first match for pattern in code
    match = re.search(pattern, code, flags)
    if match is None:
        raise ValueError(f'{pattern!r} not found')
    return match.groups()[0]


with open('molmass/molmass.py') as fh:
    code = fh.read()

version = search(r"__version__ = '(.*?)'", code)

description = search(r'"""(.*)\.(?:\r\n|\r|\n)', code)

readme = search(
    r'(?:\r\n|\r|\n){2}"""(.*)"""(?:\r\n|\r|\n){2}from __future__',
    code,
    re.MULTILINE | re.DOTALL,
)

readme = '\n'.join(
    [description, '=' * len(description)] + readme.splitlines()[1:]
)

license = search(
    r'(# Copyright.*?(?:\r\n|\r|\n))(?:\r\n|\r|\n)+""',
    code,
    re.MULTILINE | re.DOTALL,
)

license = license.replace('# ', '').replace('#', '')

if 'sdist' in sys.argv:
    with open('LICENSE', 'w') as fh:
        fh.write('BSD 3-Clause License\n\n')
        fh.write(license)
    with open('README.rst', 'w') as fh:
        fh.write(readme)

setup(
    name='molmass',
    version=version,
    license='BSD',
    description=description,
    long_description=readme,
    author='Christoph Gohlke',
    author_email='cgohlke@cgohlke.com',
    url='https://www.cgohlke.com',
    project_urls={
        'Bug Tracker': 'https://github.com/cgohlke/molmass/issues',
        'Source Code': 'https://github.com/cgohlke/molmass',
        # 'Documentation': 'https://',
    },
    packages=['molmass'],
    package_data={'molmass': ['icon.png']},
    entry_points={
        'console_scripts': [
            'molmass = molmass.molmass:main',
            'molmass_web = molmass.web:main',
        ],
        'gui_scripts': ['elements_gui = molmass.elements_gui:main'],
    },
    extras_require={
        'all': ['Flask', 'pandas'],
        'gui': ['wxPython>=4.0'],
    },
    platforms=['any'],
    python_requires='>=3.8',
    classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: BSD License',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
)
