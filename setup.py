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

"""
This setup.py now only covers how to compile the C extension modules. All
other information has been moved into pyprojects.toml.
"""

import sys

from setuptools import setup, Extension


# So I can import fmask itself, which will allow access to fmask.__version__
# to support the dynamic version attribute specified in pyproject.toml
sys.path.append(".")

try:
    from numpy import get_include as numpy_get_include
    withExtensions = True
except ImportError:
    # I suspect that this can no longer happen. We had to add numpy as a
    # build-time dependency inside the pyproject.toml file, which I think
    # means that it will always be present during any build process, even
    # on ReadTheDocs or similar.
    withExtensions = False

# Trigger the most up-to-date numpy compile-time deprecation warnings.
# See https://numpy.org/doc/stable/reference/c-api/deprecations.html
# for a not-very-clear explanation.
# The messages are not visible during the install process (i.e. using pip),
# but ARE visible when building a wheel file (e.g. when using "python -m build")
NUMPY_DEPR_WARN = ('NPY_NO_DEPRECATED_API', 'NPY_2_0_API_VERSION')

if withExtensions:
    # This is for a normal build
    fillminimaC = Extension(name='fmask._fillminima',
        define_macros=[NUMPY_DEPR_WARN],
        sources=['c_src/fillminima.c'],
        include_dirs=[numpy_get_include()])
    valueIndexesC = Extension(name='fmask._valueindexes',
        define_macros=[NUMPY_DEPR_WARN],
        sources=['c_src/valueindexes.c'],
        include_dirs=[numpy_get_include()])
    extensionsList = [fillminimaC, valueIndexesC]
else:
    # This would be for a ReadTheDocs build. As noted above, this probably never
    # happens anymore, since the addition of pyproject.toml
    extensionsList = []

# do the setup
setup(ext_modules=extensionsList)
          
