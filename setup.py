# This file is part of 'python-fmask' - a cloud masking module
# Copyright (C) 2015  Neil Flood
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
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
import sys
import glob
import fmask

# If we fail to import the numpy version of setup(), still try to proceed, as it is possibly
# because we are being run by ReadTheDocs, and so we just need to be able to generate documentation. 
try:
    from numpy.distutils.core import setup, Extension
    withExtensions = True
except ImportError:
    from distutils.core import setup
    withExtensions = False

# When building the sdist on Linux we want the extra .bat
# files that are need for the Windows install. 
INCLUDE_WINDOWS_BAT = int(os.getenv('FMASK_INCLUDEBAT', '0')) > 0

# use the latest numpy API
NUMPY_MACROS = ('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')

if withExtensions:
    # This is for a normal build
    fillminimaC = Extension(name='_fillminima', 
                define_macros=[NUMPY_MACROS],
                sources=['src/fillminima.c'])
    valueIndexesC = Extension(name='_valueindexes',
                define_macros=[NUMPY_MACROS],
                sources=['src/valueindexes.c'])
    extensionsList = [fillminimaC, valueIndexesC]
else:
    # This would be for a ReadTheDocs build. 
    from distutils.core import setup
    extensionsList = []

scriptList = glob.glob("bin/*.py")
if sys.platform == 'win32' or INCLUDE_WINDOWS_BAT:
    # include any .bat file helpers also (just one at this stage)
    batList = glob.glob("bin/*.bat")
    scriptList.extend(batList)
    
# do the setup
setup( name = 'python-fmask',
        version = fmask.__version__,
        description = 'Module to implement the fmask cloud masking algorithm (Zhu, Wang & Woodcock 2015)',
        author = 'Neil Flood',
        author_email = 'neil.flood@dsiti.qld.gov.au',
        scripts = scriptList,
        packages = ['fmask'],
        ext_package = 'fmask',
        ext_modules = extensionsList,
        license='LICENSE.txt',
        url='https://bitbucket.org/chchrsc/python-fmask',
        classifiers=['Intended Audience :: Developers',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 2',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.4',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6'])
          
