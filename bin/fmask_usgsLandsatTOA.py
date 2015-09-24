#!/usr/bin/env python

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
from __future__ import print_function, division

import sys
import optparse
from fmask import landsatTOA

class CmdArgs(object):
    """
    Class for processing command line arguments
    """
    def __init__(self):
        self.parser = optparse.OptionParser()
        self.parser.add_option('-i', '--infile', dest='infile',
            help='Input raw DN reflective image')
        self.parser.add_option('-m', '--mtl', dest='mtl', 
            help='.MTL  file')
        self.parser.add_option('-o', '--output', dest='output',
            help='Output TOA reflectance file')

        (options, self.args) = self.parser.parse_args()
        self.__dict__.update(options.__dict__)

        if (self.infile is None or self.mtl is None or  
                self.output is None):
            self.parser.print_help()
            sys.exit()

def mainRoutine():
    cmdargs = CmdArgs()
    
    landsatTOA.makeTOAReflectance(cmdargs.infile, cmdargs.mtl, cmdargs.output)
    
if __name__ == '__main__':
    mainRoutine()
    
    