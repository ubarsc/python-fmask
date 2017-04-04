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
import argparse
from fmask import fmask
from fmask import config

from rios import fileinfo


def getCmdargs():
    """
    Get command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--thermal', help='Input stack of thermal bands')
    parser.add_argument('-a', '--toa', 
        help='Input stack of TOA reflectance (see fmask_usgsLandsatTOA.py)')
    parser.add_argument('-m', '--mtl', help='Input .MTL file')
    parser.add_argument('-s', '--saturation', 
        help='Input saturation mask (see fmask_usgsLandsatSaturationMask.py)')
    parser.add_argument("-z", "--anglesfile", 
        help="Image of sun and satellite angles (see fmask_usgsLandsatMakeAnglesImage.py)")
    parser.add_argument('-o', '--output', dest='output',
        help='output cloud mask')
    parser.add_argument('-v', '--verbose', default=False,
        action='store_true', help='verbose output')
    parser.add_argument('-k', '--keepintermediates', dest='keepintermediates', 
        default=False, action='store_true', help='verbose output')
    parser.add_argument('-e', '--tempdir', 
        default='.', help="Temp directory to use (default=%(default)s)")

    params = parser.add_argument_group(title="Configurable parameters", description="""
        Changing these parameters will affect the way the algorithm works, and thus the 
        quality of the final output masks. 
        """)
    params.add_argument("--mincloudsize", type=int, default=0, 
        help="Mininum cloud size (in pixels) to retain, before any buffering. Default=%(default)s)")
    params.add_argument("--cloudbufferdistance", type=float, default=150,
        help="Distance (in metres) to buffer final cloud objects (default=%(default)s)")
    params.add_argument("--shadowbufferdistance", type=float, default=300,
        help="Distance (in metres) to buffer final cloud shadow objects (default=%(default)s)")
    defaultCloudProbThresh = 100 * config.FmaskConfig.Eqn17CloudProbThresh
    params.add_argument("--cloudprobthreshold", type=float, default=defaultCloudProbThresh,
        help=("Cloud probability threshold (percentage) (default=%(default)s). This is "+
            "the constant term at the end of equation 17, given in the paper as 0.2 (i.e. 20%%). "+
            "To reduce commission errors, increase this value, but this will also increase "+
            "omission errors. "))
    dfltNirSnowThresh = config.FmaskConfig.Eqn20NirSnowThresh
    params.add_argument("--nirsnowthreshold", default=dfltNirSnowThresh, type=float,
        help=("Threshold for NIR reflectance (range [0-1]) for snow detection "+
            "(default=%(default)s). Increase this to reduce snow commission errors"))
    dfltGreenSnowThresh = config.FmaskConfig.Eqn20GreenSnowThresh
    params.add_argument("--greensnowthreshold", default=dfltGreenSnowThresh, type=float,
        help=("Threshold for Green reflectance (range [0-1]) for snow detection "+
            "(default=%(default)s). Increase this to reduce snow commission errors"))

    cmdargs = parser.parse_args()

    if (cmdargs.thermal is None or cmdargs.anglesfile is None or 
            cmdargs.mtl is None is None or cmdargs.output is None
            or cmdargs.toa is None):
        parser.print_help()
        sys.exit(1)
    
    return cmdargs

            
def mainRoutine():
    """
    Main routine that calls fmask
    """
    cmdargs = getCmdargs()
    
    # 1040nm thermal band should always be the first (or only) band in a
    # stack of Landsat thermal bands
    thermalInfo = config.readThermalInfoFromLandsatMTL(cmdargs.mtl)
                        
    anglesfile = cmdargs.anglesfile
    anglesInfo = config.AnglesFileInfo(anglesfile, 3, anglesfile, 2, anglesfile, 1, anglesfile, 0)
    
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
    fmaskConfig.setMinCloudSize(cmdargs.mincloudsize)
    fmaskConfig.setEqn17CloudProbThresh(cmdargs.cloudprobthreshold / 100)    # Note conversion from percentage
    fmaskConfig.setEqn20NirSnowThresh(cmdargs.nirsnowthreshold)
    fmaskConfig.setEqn20GreenSnowThresh(cmdargs.greensnowthreshold)

    # Work out a suitable buffer size, in pixels, dependent on the resolution of the input TOA image
    toaImgInfo = fileinfo.ImageInfo(cmdargs.toa)
    fmaskConfig.setCloudBufferSize(int(cmdargs.cloudbufferdistance / toaImgInfo.xRes))
    fmaskConfig.setShadowBufferSize(int(cmdargs.shadowbufferdistance / toaImgInfo.xRes))
    
    fmask.doFmask(fmaskFilenames, fmaskConfig)
    
if __name__ == '__main__':
    mainRoutine()

