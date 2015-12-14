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
from fmask import config

class CmdArgs(object):
    """
    Class for processing command line arguments
    """
    def __init__(self):
        self.parser = optparse.OptionParser()
        self.parser.add_option('-t', '--thermal', dest='thermal',
            help='Input stack of thermal bands')
        self.parser.add_option('-a', '--toa', dest='toa',
            help='Input stack of TOA reflectance (see fmask_usgsLandsatTOA.py)')
        self.parser.add_option('-m', '--mtl', dest='mtl',
            help='Input .MTL file')
        self.parser.add_option('-s', '--saturation', dest='saturation',
            help='Input saturation mask (see fmask_usgsLandsatSaturationMask.py)')
        self.parser.add_option('-o', '--output', dest='output',
            help='output cloud mask')
        self.parser.add_option('-v', '--verbose', dest='verbose', default=False,
            action='store_true', help='verbose output')
        self.parser.add_option('-k', '--keepintermediates', dest='keepintermediates', 
            default=False, action='store_true', help='verbose output')
        self.parser.add_option('-e', '--tempdir', dest='tempdir',
            default='.', help="Temp directory to use (default=%default)")
            
        (options, self.args) = self.parser.parse_args()
        self.__dict__.update(options.__dict__)
        
        if (self.thermal is None or 
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
    thermalInfo = config.readThermalInfoFromLandsatMTL(cmdargs.mtl)
                        
    anglesInfo = config.readAnglesFromLandsatMTL(cmdargs.mtl)
    
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
        
    fmaskFilenames = config.FmaskFilenames()
    fmaskFilenames.setTOAReflectanceFile(cmdargs.toa)
    fmaskFilenames.setThermalFile(cmdargs.thermal)
    fmaskFilenames.setOutputCloudMaskFile(cmdargs.output)
    if cmdargs.saturation is not None:
        fmaskFilenames.setSaturationMask(cmdargs.saturation)
    else:
        print('saturation mask not supplied - see fmask_usgsLandsatSaturationMask.py')
    
    fmaskConfig = config.FmaskConfig(sensor)
    fmaskConfig.setThermalInfo(thermalInfo)
    fmaskConfig.setAnglesInfo(anglesInfo)
    fmaskConfig.setKeepIntermediates(cmdargs.keepintermediates)
    fmaskConfig.setVerbose(cmdargs.verbose)
    fmaskConfig.setTempDir(cmdargs.tempdir)
    
    fmask.doFmask(fmaskFilenames, fmaskConfig)
    
if __name__ == '__main__':
    mainRoutine()

