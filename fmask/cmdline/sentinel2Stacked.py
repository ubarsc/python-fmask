"""
Script that takes a stacked Sentinel 2 Level 1C image and runs
fmask on it.
"""
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
from __future__ import print_function, division

import sys
import os
import argparse
import tempfile
import glob

from osgeo import gdal
from osgeo_utils import gdal_merge

from rios import fileinfo
from rios.imagewriter import DEFAULTDRIVERNAME, dfltDriverOptions

from fmask import config
from fmask import fmaskerrors
from fmask.cmdline import sentinel2makeAnglesImage
from fmask import fmask
from fmask import sen2meta

# for GDAL.
CMDLINECREATIONOPTIONS = []
if DEFAULTDRIVERNAME in dfltDriverOptions:
    for opt in dfltDriverOptions[DEFAULTDRIVERNAME]:
        CMDLINECREATIONOPTIONS.append('-co')
        CMDLINECREATIONOPTIONS.append(opt)


def getCmdargs(argv=None):
    """
    Get command line arguments
    
    If argv is given, it should be a list of pairs of parameter and arguments like in command line.
    See getCmdargs(['-h']) for details on available parameters.
    Example: getCmdargs(['--safedir', '<.SAFE directory>', '-o', '<output file>'])
    
    If argv is None or not given, command line sys.args are used, see argparse.parse_args.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--safedir", help=("Name of .SAFE directory, as unzipped from " +
        "a standard ESA L1C zip file. Using this option will automatically create intermediate " +
        "stacks of the input bands, and so does NOT require --toa or --anglesfile. "))
    parser.add_argument("--granuledir", help=("Name of granule sub-directory within the " +
        ".SAFE directory, as unzipped from a standard ESA L1C zip file. This option is an " +
        "alternative to --safedir, for use with ESA's old format zipfiles which had multiple " +
        "granules in each zipfile. Specify the subdirectory of the single tile, under the " +
        "<safedir>/GRANULE/ directory. " +
        "Using this option will automatically create intermediate " +
        "stacks of the input bands, and so does NOT require --toa or --anglesfile. "))
    parser.add_argument('-a', '--toa', 
        help=('Input stack of TOA reflectance (as supplied by ESA). This is obsolete, and is ' +
            'only required if NOT using the --safedir or --granuledir option. '))
    parser.add_argument('-z', '--anglesfile', 
        help=("Input angles file containing satellite and sun azimuth and zenith. " +
            "See fmask_sentinel2makeAnglesImage.py for assistance in creating this. " +
            "This option is obsolete, and is only required if NOT using the --safedir " +
            "or --granuledir option. "))
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

    cmdargs = parser.parse_args(argv)

    # Do some sanity checks on what was given
    safeDirGiven = (cmdargs.safedir is not None)
    granuleDirGiven = (cmdargs.granuledir is not None)
    if granuleDirGiven and safeDirGiven:
        print("Only give one of --safedir or --granuledir. The --granuledir is only ")
        print("required for multi-tile zipfiles in the old ESA format")
        sys.exit(1)
    stackAnglesGiven = (cmdargs.toa is not None and cmdargs.anglesfile is not None)
    multipleInputGiven = (safeDirGiven or granuleDirGiven) and stackAnglesGiven
    inputGiven = safeDirGiven or granuleDirGiven or stackAnglesGiven
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
        
        options = gdal.WarpOptions(format='VRT', xRes=toaImgInfo.xRes, 
            yRes=toaImgInfo.yRes, outputBounds=[toaImgInfo.xMin, toaImgInfo.yMin,
            toaImgInfo.xMax, toaImgInfo.yMax], resampleAlg='near')
        gdal.Warp(vrtName, inputAnglesFile, options=options)

        outputAnglesFile = vrtName
    
    return outputAnglesFile


def makeStackAndAngles(cmdargs):
    """
    Make an intermediate stack of all the TOA reflectance bands. Also make an image
    of the angles. Fill in the names of these in the cmdargs object. 
        
    """
    if cmdargs.granuledir is None and cmdargs.safedir is not None:
        cmdargs.granuledir = findGranuleDir(cmdargs.safedir)

    # Make the angles file
    (fd, anglesfile) = tempfile.mkstemp(dir=cmdargs.tempdir, prefix="angles_tmp_", 
        suffix=".img")
    os.close(fd)
    xmlfile = findGranuleXml(cmdargs.granuledir)
    if cmdargs.verbose:
        print("Making angles image")
    sentinel2makeAnglesImage.makeAngles(xmlfile, anglesfile)
    cmdargs.anglesfile = anglesfile
    
    # Make a stack of the reflectance bands. Not that we do an explicit resample to the
    # output pixel size, to avoid picking up the overview layers with the ESA jpg files. 
    # According to @vincentschut, these are shifted slightly, and should be avoided.
    bandList = ['B01', 'B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A',
        'B09', 'B10', 'B11', 'B12']
    imgDir = "{}/IMG_DATA".format(cmdargs.granuledir)
    resampledBands = []
    for band in bandList:
        (fd, tmpBand) = tempfile.mkstemp(dir=cmdargs.tempdir, prefix="tmp_{}_".format(band),
            suffix=".vrt")
        os.close(fd)
        inBandImgList = glob.glob("{}/*_{}.jp2".format(imgDir, band))
        if len(inBandImgList) != 1:
            raise fmaskerrors.FmaskFileError("Cannot find input band {}".format(band))
        inBandImg = inBandImgList[0]

        # Now make a resampled copy to the desired pixel size, using the right resample method
        resampleMethod = chooseResampleMethod(cmdargs.pixsize, inBandImg)
        
        options = gdal.WarpOptions(format='VRT', resampleAlg=resampleMethod,
            xRes=cmdargs.pixsize, yRes=cmdargs.pixsize)
        gdal.Warp(tmpBand, inBandImg, options=options)
        
        resampledBands.append(tmpBand)
    
    # Now make a stack of these
    if cmdargs.verbose:
        print("Making stack of all bands, at {}m pixel size".format(cmdargs.pixsize))
    (fd, tmpStack) = tempfile.mkstemp(dir=cmdargs.tempdir, prefix="tmp_allbands_",
        suffix=".img")
    os.close(fd)
    cmdargs.toa = tmpStack

    # We need to turn off exceptions while using gdal_merge, as it doesn't cope
    usingExceptions = gdal.GetUseExceptions()
    gdal.DontUseExceptions()
    gdal_merge.main(['-q', '-of', DEFAULTDRIVERNAME] + CMDLINECREATIONOPTIONS + 
        ['-separate', '-o', cmdargs.toa] + resampledBands)
    if usingExceptions:
        gdal.UseExceptions()
    
    for fn in resampledBands:
        fmask.deleteRaster(fn)

    return resampledBands


def chooseResampleMethod(outpixsize, inBandImg):
    """
    Choose the right resample method, given the image and the desired output pixel size
    """
    imginfo = fileinfo.ImageInfo(inBandImg)
    inPixsize = imginfo.xRes
    
    if outpixsize == inPixsize:
        resample = "near"
    elif outpixsize > inPixsize:
        resample = "average"
    else:
        resample = "cubic"
    
    return resample


def findGranuleDir(safedir):
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
    return granuleDir


def findGranuleXml(granuleDir):
    """
    Find the granule-level XML file, given the granule dir
    """
    xmlfile = "{}/MTD_TL.xml".format(granuleDir)
    if not os.path.exists(xmlfile):
        # Might be old-format zipfile, so search for *.xml
        xmlfilePattern = "{}/*.xml".format(granuleDir)
        xmlfileList = glob.glob(xmlfilePattern)
        if len(xmlfileList) == 1:
            xmlfile = xmlfileList[0]
        else:
            raise fmaskerrors.FmaskFileError("Unable to find XML file {}".format(xmlfile))
    return xmlfile


def readTopLevelMeta(cmdargs):
    """
    Since ESA introduced the radiometric offsets in version 04.00
    of their processing, we now need to read the XML metadata from 
    the top level of the SAFE directory. 

    Return a sen2meta.Sen2ZipfileMeta object. 

    """
    safeDir = cmdargs.safedir
    if safeDir is None:
        safeDir = os.path.join(cmdargs.granuledir, os.pardir, os.pardir)
        if not os.path.exists(os.path.join(safeDir, 'GRANULE')):
            msg = "Cannot identify the SAFE-level directory, which is now required (since ESA version 04.00)"
            raise fmaskerrors.FmaskFileError(msg)

    # If we can find the top-level XML file, then we can check it for
    # offset values. First look for the stndard named file
    topLevelXML = os.path.join(safeDir, 'MTD_MSIL1C.xml')
    if not os.path.exists(topLevelXML):
        # We may have an old-format zip file, in which the XML is
        # named for the date of acquisition. It should be the only 
        # .xml file in that directory. 
        topLevelXML = None
        xmlList = [f for f in glob.glob(os.path.join(safeDir, "*.xml")) if "INSPIRE.xml" not in f]
        if len(xmlList) == 1:
            topLevelXML = xmlList[0]
    if topLevelXML is None:
        msg = "Unable to find top-level XML file, which is now required"
        raise fmaskerrors.FmaskFileError(msg)

    topMeta = sen2meta.Sen2ZipfileMeta(xmlfilename=topLevelXML)
    return topMeta


def makeRefOffsetDict(topMeta):
    """
    Take the given sen2meta.Sen2ZipfileMeta object and convert it
    into a dictionary suitable to give to FmaskConfig.setTOARefOffsetDict.

    """
    # This dictionary established a correspondance between the string
    # given in sen2meta.nameFromBandId and the internal index values used 
    # for bands within python-fmask. 
    # Note that this should include every band used within the Fmask code, 
    # although not necessarily every Sentinel-2 band. 
    bandIndexNameDict = {config.BAND_BLUE: "B02", config.BAND_GREEN: "B03", 
        config.BAND_RED: "B04", config.BAND_NIR: "B08", 
        config.BAND_SWIR1: "B11", config.BAND_SWIR2: "B12", 
        config.BAND_CIRRUS: "B10", config.BAND_S2CDI_NIR8A: "B08A", 
        config.BAND_S2CDI_NIR7: "B07", config.BAND_WATERVAPOUR: "B09"}

    offsetDict = {}
    for bandNdx in bandIndexNameDict:
        bandNameStr = bandIndexNameDict[bandNdx]
        offsetVal = topMeta.offsetValDict[bandNameStr]
        offsetDict[bandNdx] = offsetVal
    return offsetDict


def mainRoutine(argv=None):
    """
    Main routine that calls fmask
    
    If argv is given, it should be a list of pairs of parameter and arguments like in command line.
    See mainRoutine(['-h']) for details on available parameters.
    Example: mainRoutine(['--safedir', '<.SAFE directory>', '-o', '<output file>'])
    
    If argv is None or not given, command line sys.args are used, see argparse.parse_args.
    """
    cmdargs = getCmdargs(argv)
    tempStack = False
    if cmdargs.safedir is not None or cmdargs.granuledir is not None:
        tempStack = True
        makeStackAndAngles(cmdargs)
    topMeta = readTopLevelMeta(cmdargs)
    
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
    fmaskConfig.setTOARefScaling(topMeta.scaleVal)
    offsetDict = makeRefOffsetDict(topMeta)
    fmaskConfig.setTOARefOffsetDict(offsetDict)
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
    
    if (anglesfile != cmdargs.anglesfile):
        # Must have been a temporary, so remove it
        fmask.deleteRaster(anglesfile)
    
    if tempStack and not cmdargs.keepintermediates:
        for fn in [cmdargs.toa, cmdargs.anglesfile]:
            if os.path.exists(fn):
                fmask.deleteRaster(fn)
    
