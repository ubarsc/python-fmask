#!/usr/bin/env python

"""
Script that takes USGS landsat stacked separately for 
reflective and thermal and runs the fmask on it.
"""
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
from fmask import fmask

class CmdArgs(object):
    """
    Class for processing command line arguments
    """
    def __init__(self):
        self.parser = optparse.OptionParser()
        self.parser.add_option('-v', '--visible', dest='visible',
            help='Input stack of visible bands')
        self.parser.add_option('-t', '--thermal', dest='thermal',
            help='Input stack of thermal bands')
        self.parser.add_option('-a', '--toa', dest='toa',
            help='Input stack of TOA reflectance (see fmask_usgsLandsatTOA.py)')
        self.parser.add_option('-m', '--mtl', dest='mtl',
            help='Input .MTL file')
        self.parser.add_option('-o', '--output', dest='output',
            help='output cloud mask')
        self.parser.add_option('-V', '--verbose', dest='verbose', default=False,
            action='store_true', help='verbose output')
        self.parser.add_option('-k', '--keepintermediates', dest='keepintermediates', 
            default=False, action='store_true', help='verbose output')
            
        (options, self.args) = self.parser.parse_args()
        self.__dict__.update(options.__dict__)
        
        if (self.visible is None or self.thermal is None or 
                self.mtl is None is None or self.output is None
                or self.toa is None):
            self.parser.print_help()
            sys.exit(1)
            
def mainRoutine():
    """
    Main routine that calls fmask
    """
    cmdargs = CmdArgs()
    
    # 1040nm thermal band should always be the first (or only) band in a
    # stack of Landsat thermal bands
    thermalInfo = fmask.readThermalInfoFromLandsatMTL(cmdargs.mtl, 
                        cmdargs.thermal, 0)
                        
    anglesInfo = fmask.readAnglesFromLandsatMTL(cmdargs.mtl)
    
    mtlInfo = fmask.readMTLFile(cmdargs.mtl)
    landsat = mtlInfo['SPACECRAFT_ID'][-1]
    
    # first 3 are equivalent anyway, but for completeness..
    if landsat == '4':
        bandInfo = fmask.LANDSAT4_TM_BANDS
    elif landsat == '5':
        bandInfo = fmask.LANDSAT5_TM_BANDS
    elif landsat == '7':
        bandInfo = fmask.LANDSAT7_ETM_BANDS
    elif landsat == '8':
        bandInfo = fmask.LANDSAT8_OLI_BANDS
    
    fmask.doFmask(cmdargs.visible, bandInfo, cmdargs.toa, anglesInfo, 
                cmdargs.output, cmdargs.thermal, thermalInfo=thermalInfo,
                verbose=cmdargs.verbose, keepintermediates=cmdargs.keepintermediates)

if __name__ == '__main__':
    mainRoutine()

