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
import glob

from rios import fileinfo
from rios.parallel.jobmanager import find_executable

from fmask import config
from fmask import fmaskerrors
from fmask import fmask

def getCmdargs():
    """
    Get command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--safedir", help=("Name of .SAFE directory, as unzipped from " +
        "a standard ESA L1C zip file. Using this option will automatically create intermediate " +
        "stacks of the input bands, and so does NOT require --toa or --anglesfile. "))
    parser.add_argument('-a', '--toa', 
        help=('Input stack of TOA reflectance (as supplied by ESA). This is only required if NOT' +
            'using the --safedir option. '))
    parser.add_argument('-z', '--anglesfile', 
        help=("Input angles file containing satellite and sun azimuth and zenith. " +
            "See fmask_sentinel2makeAnglesImage.py for assistance in creating this. " +
            "This option is only required if NOT using the --safedir option. "))
    parser.add_argument('-o', '--output', help='Output cloud mask')
    parser.add_argument('-v', '--verbose', dest='verbose', default=False,
        action='store_true', help='verbose output')
    parser.add_argument("--pixsize", default=20, type=int, 
        help="Output pixel size in metres (default=%(default)s)")
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
    params.add_argument("--parallaxtest", default=False, action="store_true",
        help="Turn on the parallax displacement test from Frantz (2018) (default will not use this test)")

    cmdargs = parser.parse_args()

    safeDirGiven = (cmdargs.safedir is not None)
    stackAnglesGiven = (cmdargs.toa is not None and cmdargs.anglesfile is not None)
    multipleInputGiven = safeDirGiven and stackAnglesGiven
    inputGiven = safeDirGiven or stackAnglesGiven
    if cmdargs.output is None or multipleInputGiven or not inputGiven:
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


def makeStackAndAngles(cmdargs):
    """
    Make an intermediate stack of all the TOA reflectance bands. Also make an image
    of the angles. Fill in the names of these in the cmdargs object. 
        
    """
    # Find the commands we need, even under Windoze
    anglesScript = find_executable("fmask_sentinel2makeAnglesImage.py")
    gdalTranslateCmd = find_executable("gdal_translate")
    gdalWarpCmd = find_executable("gdalwarp")
    buildvrtCmd = find_executable("gdalbuildvrt")

    # Make the angles file
    (fd, anglesfile) = tempfile.mkstemp(dir=cmdargs.tempdir, prefix="angles_tmp_", 
        suffix=".tif")
    os.close(fd)
    xmlfile = findGranuleXml(cmdargs.safedir)
    cmd = "{} -i {} -o {}".format(anglesScript, xmlfile, anglesfile)
    if cmdargs.verbose:
        print("Making angles image")
    os.system(cmd)
    cmdargs.anglesfile = anglesfile
    
    # Make a stack of the reflectance bands. Not that we do an explicit resample to the
    # output pixel size, to avoid picking up the overview layers with the ESA jpg files. 
    # According to @vincentschut, these are shifted slightly, and should be avoided.
    bandList = ['B01', 'B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A',
        'B09', 'B10', 'B11', 'B12']
    imgDir = "{}/GRANULE/L1C_*/IMG_DATA".format(cmdargs.safedir)
    resampledBands = []
    (fd, tmpBand1) = tempfile.mkstemp(dir=cmdargs.tempdir, prefix="tmp1_",
        suffix=".tif")
    if cmdargs.verbose:
        print("Resampling explicitly, to avoid shifted overviews")
    for band in bandList:
        (fd, tmpBand2) = tempfile.mkstemp(dir=cmdargs.tempdir, prefix="tmp2_{}_".format(band),
            suffix=".tif")
        inBandImgList = glob.glob("{}/*_{}.jp2".format(imgDir, band))
        if len(inBandImgList) != 1:
            raise fmaskerrors.FmaskFileError("Cannot find input band {}".format(band))
        inBandImg = inBandImgList[0]
        # First make a plain copy with no overviews
        cmd = ("{gdaltranslate} -co TILED=YES -of GTiff {inimg} {outimg}").format(
            gdaltranslate=gdalTranslateCmd, inimg=inBandImg, outimg=tmpBand1)
        os.system(cmd)
        # Now make a resampled copy to the desired pixel size
        cmd = ("{gdalwarp} -tr {pixsize} {pixsize} -co TILED=YES -of GTiff "+
            "-r average {inimg} {outimg}").format(gdalwarp=gdalWarpCmd, 
            pixsize=cmdargs.pixsize, inimg=tmpBand1, outimg=tmpBand2)
        os.system(cmd)
        os.remove(tmpBand1)
        
        resampledBands.append(tmpBand2)
    
    # Now make a vrt stack of these
    cmdargs.toa = "allbands.vrt"
    cmd = "{buildvrt} -separate {outvrt} {inimgs}".format(buildvrt=buildvrtCmd,
        outvrt=cmdargs.toa, inimgs=' '.join(resampledBands))
    if cmdargs.verbose:
        print(cmd)
    os.system(cmd)
    
    return resampledBands


def findGranuleXml(safedir):
    """
    Search the given .SAFE directory, and find the main XML file at the GRANULE level.
    
    Note that this currently only works for the new-format zip files, with one 
    tile per zipfile. The old ones are being removed from service, so we won't 
    cope with them. 
    
    """
    granuleDirPattern = "{}/GRANULE/L1C_*".format(safedir)
    granuleDirList = glob.glob(granuleDirPattern)
    if len(granuleDirList) == 0:
        raise fmaskerrors.FmaskFileError("Unable to find GRANULE sub-directory {}".format(granuleDirPattern))
    elif len(granuleDirList) > 1:
        dirstring = ','.join(granuleDirList)
        msg = "Found multiple GRANULE sub-directories: {}".format(dirstring)
        raise fmaskerrors.FmaskFileError(msg)
    
    granuleDir = granuleDirList[0]
    xmlfile = "{}/MTD_TL.xml".format(granuleDir)
    if not os.path.exists(xmlfile):
        raise fmaskerrors.FmaskFileError("Unable to find XML file {}".format(xmlfile))
    return xmlfile


def mainRoutine():
    """
    Main routine that calls fmask
    """
    cmdargs = getCmdargs()
    if cmdargs.safedir is not None:
        resampledBands = makeStackAndAngles(cmdargs)
    
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
    fmaskConfig.setSen2displacementTest(cmdargs.parallaxtest)
    
    # Work out a suitable buffer size, in pixels, dependent on the resolution of the input TOA image
    toaImgInfo = fileinfo.ImageInfo(cmdargs.toa)
    fmaskConfig.setCloudBufferSize(int(cmdargs.cloudbufferdistance / toaImgInfo.xRes))
    fmaskConfig.setShadowBufferSize(int(cmdargs.shadowbufferdistance / toaImgInfo.xRes))
    
    fmask.doFmask(fmaskFilenames, fmaskConfig)
    
    if anglesfile != cmdargs.anglesfile:
        # Must have been a temporary vrt, so remove it
        os.remove(anglesfile)
    
    if cmdargs.safedir is not None:
        # These must be temporary files created by this script
        for fn in resampledBands:
            if os.path.exists(fn):
                os.remove(fn)
        if os.path.exists(cmdargs.toa):
            os.remove(cmdargs.toa)
    

if __name__ == '__main__':
    mainRoutine()

