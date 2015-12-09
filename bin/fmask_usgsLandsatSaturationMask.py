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
from fmask import saturationcheck
from fmask import config

class CmdArgs(object):
    """
    Class for processing command line arguments
    """
    def __init__(self):
        self.parser = optparse.OptionParser()
        self.parser.add_option('-i', '--infile', dest='infile',
            help='Input raw DN radiance image')
        self.parser.add_option('-m', '--mtl', dest='mtl', 
            help='.MTL  file')
        self.parser.add_option('-o', '--output', dest='output',
            help='Output saturation mask file')

        (options, self.args) = self.parser.parse_args()
        self.__dict__.update(options.__dict__)

        if (self.infile is None or self.mtl is None or  
                self.output is None):
            self.parser.print_help()
            sys.exit()

def mainRoutine():
    cmdargs = CmdArgs()
    
    mtlInfo = config.readMTLFile(cmdargs.mtl)
    landsat = mtlInfo['SPACECRAFT_ID'][-1]
    
    if landsat == '4':
        sensor = config.FMASK_LANDSAT47
    elif landsat == '5':
        sensor = config.FMASK_LANDSAT47
    elif landsat == '7':
        sensor = config.FMASK_LANDSAT47
    elif landsat == '8':
        sensor = config.FMASK_LANDSAT8
    else:
        raise SystemExit('Unsupported Landsat sensor')

    # needed so the saturation function knows which
    # bands are visible etc.
    fmaskConfig = config.FmaskConfig(sensor)
    
    saturationcheck.makeSaturationMask(fmaskConfig, cmdargs.infile, 
            cmdargs.output)
    
if __name__ == '__main__':
    mainRoutine()
    

