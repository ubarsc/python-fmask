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

from numpy.distutils.core import setup, Extension

fillminimaC = Extension(name='_fillminima', 
                sources=['src/fillminima.c'])
valueIndexesC = Extension(name='_valueindexes',
                sources=['src/valueindexes.c'])

# do the setup
setup( name = 'python-fmask',
        version = '0.1',
        description = 'Module to implement the fmask cloud masking algorithm (Zhu, Wang & Woodcock 2015)',
        author = 'Neil Flood',
        author_email = 'neil.flood@science.dsiti.qld.gov.au',
        packages = ['fmask'],
        ext_package = 'fmask',
        ext_modules = [fillminimaC, valueIndexesC],
        license='LICENSE.txt',
        url='https://bitbucket.org/chchrsc/python-fmask',
        classifiers=['Intended Audience :: Developers',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 2',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.2',
          'Programming Language :: Python :: 3.3',
          'Programming Language :: Python :: 3.4',
          'Programming Language :: Python :: 3.5'])
          