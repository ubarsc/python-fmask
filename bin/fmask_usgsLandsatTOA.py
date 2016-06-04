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
import argparse
from fmask import landsatTOA

def getCmdargs():
    """
    Get command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', help='Input raw DN radiance image')
    parser.add_argument('-m', '--mtl', help='.MTL  file')
    parser.add_argument("-z", "--anglesfile", 
        help="Image of sun and satellite angles (see fmask_usgsLandsatMakeAnglesImage.py)")
    parser.add_argument('-o', '--output', help='Output TOA reflectance file')

    cmdargs = parser.parse_args()

    if (cmdargs.infile is None or cmdargs.mtl is None or  
            cmdargs.output is None or cmdargs.anglesfile is None):
        parser.print_help()
        sys.exit()
    return cmdargs


def mainRoutine():
    cmdargs = getCmdargs()
    
    landsatTOA.makeTOAReflectance(cmdargs.infile, cmdargs.mtl, cmdargs.anglesfile, cmdargs.output)
    
if __name__ == '__main__':
    mainRoutine()
    
    
