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

import os
import glob
import fmask

from setuptools import setup, Extension

try:
    from numpy import get_include as numpy_get_include
    withExtensions = True
except ImportError:
    withExtensions = False

# Are we installing the command line scripts?
# this is an experimental option for users who are
# using the Python entry point feature of setuptools and Conda instead
NO_INSTALL_CMDLINE = int(os.getenv('FMASK_NOCMDLINE', '0')) > 0

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

if NO_INSTALL_CMDLINE:
    scriptList = None
else:
    scriptList = glob.glob("bin/*.py")
    
# do the setup
setup(name='python-fmask',
    version=fmask.__version__,
    description='Module to implement the fmask cloud masking algorithm (Zhu, Wang & Woodcock 2015)',
    author='Neil Flood',
    author_email='neil.flood@des.qld.gov.au',
    scripts=scriptList,
    packages=['fmask', 'fmask/cmdline'],
    ext_package='fmask',
    ext_modules=extensionsList,
    license='LICENSE.txt',
    data_files=[('', ['LICENSE.txt'])],  # add this to tarball
    url='https://www.pythonfmask.org/',
    classifiers=['Intended Audience :: Developers',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6'])
          
