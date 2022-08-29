# This file is part of 'python-fmask' - a cloud masking module
# Copyright (C) 2015  Neil Flood
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

import fmask

from setuptools import setup, Extension

try:
    from numpy import get_include as numpy_get_include
    withExtensions = True
except ImportError:
    withExtensions = False

# use the latest numpy API
NUMPY_MACROS = ('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')

if withExtensions:
    # This is for a normal build
    fillminimaC = Extension(name='_fillminima',
        define_macros=[NUMPY_MACROS],
        sources=['src/fillminima.c'],
        include_dirs=[numpy_get_include()])
    valueIndexesC = Extension(name='_valueindexes',
        define_macros=[NUMPY_MACROS],
        sources=['src/valueindexes.c'],
        include_dirs=[numpy_get_include()])
    extensionsList = [fillminimaC, valueIndexesC]
else:
    # This would be for a ReadTheDocs build. 
    extensionsList = []

# do the setup
setup(name='python-fmask',
    version=fmask.__version__,
    description='Module to implement the fmask cloud masking algorithm (Zhu, Wang & Woodcock 2015)',
    author='Neil Flood',
    author_email='neil.flood@des.qld.gov.au',
    entry_points={
        'console_scripts': [
            'fmask_sentinel2makeAnglesImage.py = fmask.cmdline.sentinel2makeAnglesImage:mainRoutine',
            'fmask_sentinel2Stacked.py = fmask.cmdline.sentinel2Stacked:mainRoutine',
            'fmask_usgsLandsatMakeAnglesImage.py = fmask.cmdline.usgsLandsatMakeAnglesImage:mainRoutine',
            'fmask_usgsLandsatSaturationMask.py = fmask.cmdline.usgsLandsatSaturationMask:mainRoutine',
            'fmask_usgsLandsatStacked.py = fmask.cmdline.usgsLandsatStacked:mainRoutine',
            'fmask_usgsLandsatTOA.py = fmask.cmdline.usgsLandsatTOA:mainRoutine'
        ]
    },
    packages=['fmask', 'fmask/cmdline'],
    ext_package='fmask',
    ext_modules=extensionsList,
    license='LICENSE.txt',
    url='https://www.pythonfmask.org/',
    classifiers=['Intended Audience :: Developers',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10'])
          
