#!/usr/bin/env python
"""
Script that takes a stacked Sentinel 2 Level 1C image and runs
fmask on it.
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
import os
import argparse
import numpy
import tempfile

from rios import fileinfo

from fmask import config
from fmask import fmask

def getCmdargs():
    """
    Get command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--toa', 
        help='Input stack of TOA reflectance (as supplied by ESA)')
    parser.add_argument('-z', '--anglesfile', 
        help=("Input angles file containing satellite and sun azimuth and zenith. " +
            "See fmask_sentinel2makeAnglesImage.py for assistance in creating this"))
    parser.add_argument('-o', '--output', help='Output cloud mask')
    parser.add_argument('-v', '--verbose', dest='verbose', default=False,
        action='store_true', help='verbose output')
    parser.add_argument('-k', '--keepintermediates', 
        default=False, action='store_true', help='Keep intermediate temporary files (normally deleted)')
    parser.add_argument('-e', '--tempdir', 
        default='.', help="Temp directory to use (default='%(default)s')")
    
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

    if cmdargs.output is None or cmdargs.toa is None or cmdargs.anglesfile is None:
        parser.print_help()
        sys.exit(1)
    
    return cmdargs


def checkAnglesFile(inputAnglesFile, toafile):
    """
    Check that the resolution of the input angles file matches that of the input
    TOA reflectance file. If not, make a VRT file which will resample it 
    on-the-fly. Only checks the resolution, assumes that if these match, then everything
    else will match too. 
    
    Return the name of the angles file to use. 
    
    """
    toaImgInfo = fileinfo.ImageInfo(toafile)
    anglesImgInfo = fileinfo.ImageInfo(inputAnglesFile)
    
    outputAnglesFile = inputAnglesFile
    if (toaImgInfo.xRes != anglesImgInfo.xRes) or (toaImgInfo.yRes != anglesImgInfo.yRes):
        (fd, vrtName) = tempfile.mkstemp(prefix='angles', suffix='.vrt')
        os.close(fd)
        cmdFmt = ("gdalwarp -q -of VRT -tr {xres} {yres} -te {xmin} {ymin} {xmax} {ymax} "+
            "-r near {infile}  {outfile} ")
        cmd = cmdFmt.format(xres=toaImgInfo.xRes, yres=toaImgInfo.yRes, xmin=toaImgInfo.xMin,
            ymin=toaImgInfo.yMin, xmax=toaImgInfo.xMax, ymax=toaImgInfo.yMax, 
            outfile=vrtName, infile=inputAnglesFile)
        os.system(cmd)
        outputAnglesFile = vrtName
    
    return outputAnglesFile
            

def mainRoutine():
    """
    Main routine that calls fmask
    """
    cmdargs = getCmdargs()
    
    anglesfile = checkAnglesFile(cmdargs.anglesfile, cmdargs.toa)
    anglesInfo = config.AnglesFileInfo(anglesfile, 3, anglesfile, 2, anglesfile, 1, anglesfile, 0)
    
    fmaskFilenames = config.FmaskFilenames()
    fmaskFilenames.setTOAReflectanceFile(cmdargs.toa)
    fmaskFilenames.setOutputCloudMaskFile(cmdargs.output)
    
    fmaskConfig = config.FmaskConfig(config.FMASK_SENTINEL2)
    fmaskConfig.setAnglesInfo(anglesInfo)
    fmaskConfig.setKeepIntermediates(cmdargs.keepintermediates)
    fmaskConfig.setVerbose(cmdargs.verbose)
    fmaskConfig.setTempDir(cmdargs.tempdir)
    fmaskConfig.setTOARefScaling(10000.0)
    fmaskConfig.setMinCloudSize(cmdargs.mincloudsize)
    fmaskConfig.setEqn17CloudProbThresh(cmdargs.cloudprobthreshold / 100)    # Note conversion from percentage
    fmaskConfig.setEqn20NirSnowThresh(cmdargs.nirsnowthreshold)
    fmaskConfig.setEqn20GreenSnowThresh(cmdargs.greensnowthreshold)
    
    # Work out a suitable buffer size, in pixels, dependent on the resolution of the input TOA image
    toaImgInfo = fileinfo.ImageInfo(cmdargs.toa)
    fmaskConfig.setCloudBufferSize(int(cmdargs.cloudbufferdistance / toaImgInfo.xRes))
    fmaskConfig.setShadowBufferSize(int(cmdargs.shadowbufferdistance / toaImgInfo.xRes))
    
    fmask.doFmask(fmaskFilenames, fmaskConfig)
    
    if anglesfile != cmdargs.anglesfile:
        # Must have been a temporary vrt, so remove it
        os.remove(anglesfile)
    

if __name__ == '__main__':
    mainRoutine()

