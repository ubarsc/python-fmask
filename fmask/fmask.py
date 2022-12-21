
"""
Implement the cloud and shadow algorithms known collectively as Fmask, 
as published in 

Zhu, Z. and Woodcock, C.E. (2012). 
Object-based cloud and cloud shadow detection in Landsat imagery
Remote Sensing of Environment 118 (2012) 83-94. 
    
and
    
Zhu, Z., Wang, S. and Woodcock, C.E. (2015).
Improvement and expansion of the Fmask algorithm: cloud, cloud
shadow, and snow detection for Landsats 4-7, 8, and Sentinel 2 images
Remote Sensing of Environment 159 (2015) 269-277.
    
Taken from Neil Flood's implementation by permission.

The notation and variable names are largely taken from the paper. Equation
numbers are also from the paper. 

Input is a top of atmosphere (TOA) reflectance file. 

The output file is a single thematic raster layer with codes representing
null, clear, cloud, shadow, snow and water. These are the values 0-5
respectively, but there are constants defined for the different codes, 
as fmask.fmask.OUTCODE_*

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

import os
import tempfile

import numpy
from osgeo import gdal
from scipy.ndimage import uniform_filter, maximum_filter, label
import scipy.stats

# We use RIOS intensively here
from rios import applier
from rios import pixelgrid
from rios import rat
from rios import imageio
from rios import fileinfo

# our wrappers for bits of C that are installed with this package
from . import fillminima
from . import valueindexes
# configuration classes
from . import config
# exceptions
from . import fmaskerrors
# so we can check if thermal all zeroes
from . import zerocheck

numpy.seterr(all='raise')

# Bands in the saturation mask, if supplied
SATURATION_BLUE = 0
SATURATION_GREEN = 1
SATURATION_RED = 2

# The values used in the final output raster
#: Output pixel value for null
OUTCODE_NULL = 0
#: Output pixel value for clear land
OUTCODE_CLEAR = 1
#: Output pixel value for cloud
OUTCODE_CLOUD = 2
#: Output pixel value for cloud shadow
OUTCODE_SHADOW = 3
#: Output pixel value for snow
OUTCODE_SNOW = 4
#: Output pixel value for water
OUTCODE_WATER = 5


def doFmask(fmaskFilenames, fmaskConfig):
    """
    Main routine for whole Fmask algorithm. Calls all other routines in sequence. 
    Parameters:
    
    * **fmaskFilenames** an instance of :class:`fmask.config.FmaskFilenames` that contains the files to use
    * **fmaskConfig** an instance of :class:`fmask.config.FmaskConfig` that contains the parameters to use
    
    If :func:`fmask.config.FmaskConfig.setKeepIntermediates` has been called with True, then
    a dictionary of intermediate files will be returned. Otherwise None is returned.
    
    """
    
    # check config.thermal and filenames.thermal both set or unset
    if (fmaskFilenames.thermal is None) != (fmaskConfig.thermalInfo is None):
        msg = 'Either both thermal filename and thermal info should be set, or neither'
        raise fmaskerrors.FmaskParameterError(msg)
        
    # do we have thermal?
    missingThermal = fmaskFilenames.thermal is None
    if not missingThermal:
        # check that it is not all zeros
        if zerocheck.isBandAllZeroes(fmaskFilenames.thermal, 
                fmaskConfig.thermalInfo.thermalBand1040um):
            if fmaskConfig.verbose:
                print('Ignoring thermal data since file is empty')
            missingThermal = True

    # do some basic checking of inputs
    if fmaskFilenames.toaRef is None:
        msg = 'Must provide input TOA reflectance file via fmaskFilenames parameter'
        raise fmaskerrors.FmaskParameterError(msg)
        
    if fmaskConfig.anglesInfo is None:
        msg = 'Must provide Angles information via fmaskConfig.setAnglesInfo'
        raise fmaskerrors.FmaskParameterError(msg)

    if fmaskFilenames.outputMask is None:
        msg = 'Output filename must be provided via fmaskFilenames parameter'
        raise fmaskerrors.FmaskParameterError(msg)

    if ((fmaskConfig.sensor == config.FMASK_SENTINEL2) and 
            (fmaskConfig.TOARefDNoffsetDict is None)):
        msg = """
            When using Fmask with Sentinel-2, it is now a requirement that
            reflectance offsets be explicitly set. This is due to ESA's 
            breaking changes in their processing version 04.00 (Nov 2021),
            which added offsets to the imagery. 
            See fmask.config.setTOARefOffsetDict() for more details. 
            Also, sen2meta.Sen2ZipfileMeta can read the necessary XML file,
            and fmask.cmdline.sentinel2Stacked.makeRefOffsetDict for further
            sample code. 
        """
        raise fmaskerrors.Sen2MetaError(msg)
        
    if fmaskConfig.strictFmask:
        # change these values back to match the paper
        fmaskConfig.setCloudBufferSize(0)
        fmaskConfig.setShadowBufferSize(3)
    
    if fmaskConfig.verbose:
        print("Cloud layer, pass 1")
    (pass1file, Twater, Tlow, Thigh, NIR_17, nonNullCount) = doPotentialCloudFirstPass(
        fmaskFilenames, fmaskConfig, missingThermal)
    if fmaskConfig.verbose:
        print("  Twater=", Twater, "Tlow=", Tlow, "Thigh=", Thigh, "NIR_17=", 
            NIR_17, "nonNullCount=", nonNullCount)
    
    if fmaskConfig.verbose:
        print("Cloud layer, pass 2")
    (pass2file, landThreshold) = doPotentialCloudSecondPass(fmaskFilenames, 
        fmaskConfig, pass1file, Twater, Tlow, Thigh, missingThermal, nonNullCount)
    if fmaskConfig.verbose:
        print("  landThreshold=", landThreshold)

    if fmaskConfig.verbose:
        print("Cloud layer, pass 3")
    interimCloudmask = doCloudLayerFinalPass(fmaskFilenames, fmaskConfig, 
        pass1file, pass2file, landThreshold, Tlow, missingThermal)
        
    if fmaskConfig.verbose:
        print("Potential shadows")
    potentialShadowsFile = doPotentialShadows(fmaskFilenames, fmaskConfig, NIR_17)
    
    if fmaskConfig.verbose:
        print("Clumping clouds")
    (clumps, numClumps) = clumpClouds(interimCloudmask)
    
    if fmaskConfig.verbose:
        print("Making 3d clouds")
    (cloudShape, cloudBaseTemp, cloudClumpNdx) = make3Dclouds(fmaskFilenames, 
        fmaskConfig, clumps, numClumps, missingThermal)
    
    if fmaskConfig.verbose:
        print("Making cloud shadow shapes")
    shadowShapesDict = makeCloudShadowShapes(fmaskFilenames, fmaskConfig,
        cloudShape, cloudClumpNdx)
    
    if fmaskConfig.verbose:
        print("Matching shadows")
    interimShadowmask = matchShadows(fmaskConfig, interimCloudmask, 
        potentialShadowsFile, shadowShapesDict, cloudBaseTemp, Tlow, Thigh, 
        pass1file)
    
    if fmaskConfig.verbose:
        print("Doing final tidy up")
    finalizeAll(fmaskFilenames, fmaskConfig, interimCloudmask, interimShadowmask, 
        pass1file)
    
    # Remove temporary files
    retVal = None
    if not fmaskConfig.keepIntermediates:
        for filename in [pass1file, pass2file, interimCloudmask, potentialShadowsFile,
                interimShadowmask]:
            deleteRaster(filename)
    else:
        # create a dictionary with the intermediate filenames so we can return them.
        retVal = {'pass1': pass1file, 'pass2': pass2file, 
            'interimCloud': interimCloudmask, 
            'potentialShadows': potentialShadowsFile, 
            'interimShadow': interimShadowmask}

    if fmaskConfig.verbose:
        print('finished fmask')
    
    return retVal


#: An offset so we can scale brightness temperature (BT, in deg C) to the range 0-255, for use in histograms.
BT_OFFSET = 176    

BT_HISTSIZE = 256
BYTE_MIN = 0
BYTE_MAX = 255
#: Gain to scale b4 reflectances to 0-255 for histograms
B4_SCALE = 500.0

#: Global RIOS window size
RIOS_WINDOW_SIZE = 512


def doPotentialCloudFirstPass(fmaskFilenames, fmaskConfig, missingThermal):
    """
    Run the first pass of the potential cloud layer. Also
    finds the temperature thresholds which will be needed 
    in the second pass, because it has the relevant data handy. 
    
    """
    infiles = applier.FilenameAssociations()
    outfiles = applier.FilenameAssociations()
    otherargs = applier.OtherInputs()
    controls = applier.ApplierControls()
    
    infiles.toaref = fmaskFilenames.toaRef
    if not missingThermal:
        infiles.thermal = fmaskFilenames.thermal
    if fmaskFilenames.saturationMask is not None:
        infiles.saturationMask = fmaskFilenames.saturationMask
    elif fmaskConfig.verbose:
        print('Saturation mask not supplied - saturated areas may not be detected')
    
    (fd, outfiles.pass1) = tempfile.mkstemp(prefix='pass1', dir=fmaskConfig.tempDir, 
                                suffix=fmaskConfig.defaultExtension)
    os.close(fd)
    if (fmaskConfig.sensor == config.FMASK_SENTINEL2) and fmaskConfig.sen2displacementTest:
        # needs overlap because of focalVariance
        overlap = int((fmaskConfig.sen2cdiWindow - 1) / 2)
        controls.setOverlap(overlap)
    controls.setWindowXsize(RIOS_WINDOW_SIZE)
    controls.setWindowYsize(RIOS_WINDOW_SIZE)
    controls.setReferenceImage(infiles.toaref)
    controls.setCalcStats(False)

    otherargs.refBands = fmaskConfig.bands  
    otherargs.thermalInfo = fmaskConfig.thermalInfo
    otherargs.waterBT_hist = numpy.zeros(BT_HISTSIZE, dtype=numpy.uint32)
    otherargs.clearLandBT_hist = numpy.zeros(BT_HISTSIZE, dtype=numpy.uint32)
    otherargs.clearLandB4_hist = numpy.zeros(BT_HISTSIZE, dtype=numpy.uint32)
    otherargs.fmaskConfig = fmaskConfig
    refImgInfo = fileinfo.ImageInfo(fmaskFilenames.toaRef)
    otherargs.refNull = refImgInfo.nodataval[0]
    if otherargs.refNull is None:
        # The null value used by USGS is 0, but is not recorded in the TIF files
        otherargs.refNull = 0
    if not missingThermal:
        thermalImgInfo = fileinfo.ImageInfo(fmaskFilenames.thermal)
        otherargs.thermalNull = thermalImgInfo.nodataval[0]
        if otherargs.thermalNull is None:
            otherargs.thermalNull = 0
    otherargs.nonNullCount = 0
    
    # Which reflective bands do we use to make a null mask. The numbers being set here 
    # are zero-based index numbers for use as array indexes. It should be just all bands, 
    # but because of the Sentinel-2 madness, we have to make special cases. 
    if fmaskConfig.sensor == config.FMASK_LANDSAT47:
        nullBandNdx = [config.BAND_BLUE, config.BAND_GREEN, config.BAND_RED, config.BAND_NIR, 
            config.BAND_SWIR1, config.BAND_SWIR2]

    elif fmaskConfig.sensor in (config.FMASK_LANDSAT8, config.FMASK_LANDSATOLI):
        nullBandNdx = [config.BAND_BLUE, config.BAND_GREEN, config.BAND_RED, config.BAND_NIR, 
            config.BAND_SWIR1, config.BAND_SWIR2, config.BAND_CIRRUS]

    elif fmaskConfig.sensor == config.FMASK_SENTINEL2:
        # For Sentinel-2, only use the visible bands to define the null mask. This is because ESA
        # are leaving a lot of spurious nulls in their imagery, most particularly in the IR bands
        # and the cirrus band. 
        nullBandNdx = [config.BAND_BLUE, config.BAND_GREEN, config.BAND_RED]

    else:
        msg = 'Unknown sensor'
        raise fmaskerrors.FmaskParameterError(msg)

    otherargs.bandsForRefNull = numpy.array([fmaskConfig.bands[i] for i in nullBandNdx])

    applier.apply(potentialCloudFirstPass, infiles, outfiles, otherargs, controls=controls)
    
    (Twater, Tlow, Thigh) = calcBTthresholds(otherargs)
    
    # 17.5 percentile of band 4, for clear land pixels. Used later in shadow masking. 
    b4_17 = scoreatpcnt(otherargs.clearLandB4_hist, 17.5)
    if b4_17 is not None:
        b4_17 = b4_17 / B4_SCALE
    else:
        # Not enough land to work this out, so guess a low value. 
        b4_17 = 0.01
    
    return (outfiles.pass1, Twater, Tlow, Thigh, b4_17, otherargs.nonNullCount)


def potentialCloudFirstPass(info, inputs, outputs, otherargs):
    """
    Called from RIOS. 
    
    Calculate the first pass potential cloud layer (equation 6)
        
    """
    fmaskConfig = otherargs.fmaskConfig

    ref = refDNtoUnits(inputs.toaref, fmaskConfig)
    # Clamp off any reflectance <= 0
    ref[ref<=0] = 0.00001

    # Extract the bands we need
    blue = otherargs.refBands[config.BAND_BLUE]
    green = otherargs.refBands[config.BAND_GREEN]
    red = otherargs.refBands[config.BAND_RED]
    nir = otherargs.refBands[config.BAND_NIR]
    swir1 = otherargs.refBands[config.BAND_SWIR1]
    swir2 = otherargs.refBands[config.BAND_SWIR2]
    if hasattr(inputs, 'thermal'):
        THERM = otherargs.thermalInfo.thermalBand1040um
    
    # Special mask needed only for resets in final pass
    refNullmask = (inputs.toaref[otherargs.bandsForRefNull] == otherargs.refNull).any(axis=0)
    if hasattr(inputs, 'thermal'):
        thermNullmask = (inputs.thermal[THERM] == otherargs.thermalNull)
        nullmask = (refNullmask | thermNullmask)
        # Brightness temperature in degrees C
        bt = otherargs.thermalInfo.scaleThermalDNtoC(inputs.thermal)
    else:
        thermNullmask = numpy.zeros_like(ref[0], dtype=bool)
        nullmask = refNullmask
    
    # Equation 1
    ndsi = (ref[green] - ref[swir1]) / (ref[green] + ref[swir1])
    ndvi = (ref[nir] - ref[red]) / (ref[nir] + ref[red])
    # In two parts, in case we have no thermal.
    basicTest = (ref[swir2] > fmaskConfig.Eqn1Swir2Thresh) & (ndsi < 0.8) & (ndvi < 0.8)
    if hasattr(inputs, 'thermal'):
        basicTest = (basicTest & (bt < fmaskConfig.Eqn1ThermThresh))
    
    # Equation 2
    meanVis = (ref[blue] + ref[green] + ref[red]) / 3.0
    whiteness = numpy.zeros(ref[0].shape)
    for n in [blue, green, red]:
        whiteness = whiteness + numpy.absolute((ref[n] - meanVis) / meanVis)

    whitenessTest = (whiteness < fmaskConfig.Eqn2WhitenessThresh)
    
    # Haze test, equation 3
    hazeTest = ((ref[blue] - 0.5 * ref[red] - 0.08) > 0)
    
    # Equation 4
    b45test = ((ref[nir] / ref[swir1]) > 0.75)
    
    # Equation 5
    waterTest = numpy.logical_or(
        numpy.logical_and(ndvi < 0.01, ref[nir] < 0.11),
        numpy.logical_and(ndvi < 0.1, ref[nir] < 0.05)
    )

    waterTest[nullmask] = False
    
    if config.BAND_CIRRUS in otherargs.refBands:
        # Zhu et al 2015, section 2.2.1. 
        cirrus = otherargs.refBands[config.BAND_CIRRUS]
        cirrusBandTest = (ref[cirrus] > fmaskConfig.cirrusBandTestThresh)
    
    # Equation 6. Potential cloud pixels (first pass)
    pcp = basicTest & whitenessTest & hazeTest & b45test
    
    # If Sentinel-2, we can use the Frantz 2018 displacement test
    if (fmaskConfig.sensor == config.FMASK_SENTINEL2) and fmaskConfig.sen2displacementTest:
        (ratio8a8, ratio8a7, v8a8, v8a7, cdi) = calcCDI(ref, fmaskConfig, otherargs.refBands)
        selection = pcp & (cdi < -0.5)
        # erode selection with 1 px
        selection = scipy.ndimage.binary_erosion(selection)
        # region grow within (cdi < -0.25)
        rg_mask = pcp & (cdi < -0.25)
        selection = scipy.ndimage.binary_dilation(selection, mask=rg_mask, iterations=0)
        pcp[~selection] = False

    # Include cirrusBandTest, from 2015 paper. Zhu et al. are not clear whether it is
    # supposed to be combined with previous tests using AND or OR, so I tried both
    # and picked what seemed best. 
    if config.BAND_CIRRUS in otherargs.refBands:
        pcp = (pcp | cirrusBandTest)
    
    # This is an extra saturation test added by DERM, and is not part of the Fmask algorithm. 
    # However, some cloud centres are saturated, and thus fail the whiteness and haze tests
    if hasattr(inputs, 'saturationMask'):
        saturatedVis = (inputs.saturationMask != 0).any(axis=0)
        veryBright = (meanVis > 0.45)
        saturatedAndBright = saturatedVis & veryBright
        pcp[saturatedAndBright] = True
        whiteness[saturatedAndBright] = 0
    
    pcp[nullmask] = False
    
    # Equation 7
    clearSkyWater = numpy.logical_and(waterTest, ref[swir2] < fmaskConfig.Eqn7Swir2Thresh)
    clearSkyWater[nullmask] = False
    
    # Equation 12
    clearLand = numpy.logical_and(numpy.logical_not(pcp), numpy.logical_not(waterTest))
    clearLand[nullmask] = False
    
    # Equation 15
    # Need to modify ndvi/ndsi by saturation......
    if hasattr(inputs, 'saturationMask'):
        modNdvi = numpy.where((inputs.saturationMask[SATURATION_GREEN] != 0), 0, ndvi)
        modNdsi = numpy.where((inputs.saturationMask[SATURATION_RED] != 0), 0, ndsi)
    else:
        modNdvi = ndvi
        modNdsi = ndsi
    # Maximum of three indices
    maxNdx = numpy.absolute(modNdvi)
    maxNdx = numpy.maximum(maxNdx, numpy.absolute(modNdsi))
    maxNdx = numpy.maximum(maxNdx, whiteness)
    variabilityProb = 1 - maxNdx
    variabilityProb[nullmask] = 0
    variabilityProbPcnt = numpy.round(variabilityProb * PROB_SCALE)
    variabilityProbPcnt = variabilityProbPcnt.clip(BYTE_MIN, BYTE_MAX).astype(numpy.uint8)
    
    # Equation 20
    # In two parts, in case we are missing thermal
    snowmask = ((ndsi > 0.15) & (ref[nir] > fmaskConfig.Eqn20NirSnowThresh) & 
        (ref[green] > fmaskConfig.Eqn20GreenSnowThresh))
    if hasattr(inputs, 'thermal'):
        snowmask = snowmask & (bt < fmaskConfig.Eqn20ThermThresh)
    snowmask[nullmask] = False
    
    # Output the pcp and water test layers. 
    outputs.pass1 = numpy.array([pcp, waterTest, clearLand, variabilityProbPcnt, 
        nullmask, snowmask, refNullmask, thermNullmask])
    
    # Accumulate histograms of temperature for land and water separately
    if hasattr(inputs, 'thermal'):
        scaledBT = (bt + BT_OFFSET).clip(0, BT_HISTSIZE)
        otherargs.waterBT_hist = accumHist(otherargs.waterBT_hist, scaledBT[clearSkyWater])
        otherargs.clearLandBT_hist = accumHist(otherargs.clearLandBT_hist, scaledBT[clearLand])
    scaledB4 = (ref[nir] * B4_SCALE).astype(numpy.uint8)
    otherargs.clearLandB4_hist = accumHist(otherargs.clearLandB4_hist, scaledB4[clearLand])
    otherargs.nonNullCount += numpy.count_nonzero(~nullmask)


def accumHist(counts, vals):
    """
    Accumulate the given values into the given (partial) counts
    """
    (valsHist, edges) = numpy.histogram(vals, bins=BT_HISTSIZE, range=(0, BT_HISTSIZE))
    # some versions of numpy seem to give an error if dtypes don't match here
    counts += valsHist.astype(counts.dtype)
    return counts


def scoreatpcnt(counts, pcnt):
    """
    Given histogram counts (binned on the range 0-255), find the value
    which corresponds to the given percentile value (0-100). 
    
    """
    n = None
    
    total = counts.sum()
    if total > 0:
        # Cumulative counts as percentages
        cumHist = numpy.cumsum(counts) * 100.0 / total
        (gtNdx, ) = numpy.where(cumHist >= pcnt)
        if len(gtNdx) > 0:
            n = gtNdx[0]
        else:
            n = 255
    return n


def refDNtoUnits(refDN, fmaskConfig):
    """
    Convert the given reflectance pixel value array to physical units,
    using parameters given in fmaskConfig. 

    Scaling is ref = (dn+offset)/scaleVal

    """
    scaleVal = float(fmaskConfig.TOARefScaling)

    bandNdxLookup = {}
    for idNum in fmaskConfig.bands:
        ndx = fmaskConfig.bands[idNum]
        bandNdxLookup[ndx] = idNum

    refUnits = numpy.zeros(refDN.shape, dtype=numpy.float32)
    numBands = refDN.shape[0]
    for bandNdx in range(numBands):
        offset = 0
        # Convert the band index number into the ID number
        # we use as key for the offset dictionary
        bandIdNum = None
        if bandNdx in bandNdxLookup:
            bandIdNum = bandNdxLookup[bandNdx]

        if bandIdNum is not None and fmaskConfig.TOARefDNoffsetDict is not None:
            offset = fmaskConfig.TOARefDNoffsetDict[bandIdNum]

        refUnits[bandNdx] = singleRefDNtoUnits(refDN[bandNdx], scaleVal, offset)
    return refUnits


def singleRefDNtoUnits(refDN, scaleVal, offset):
    """
    Apply the given scale and offset to transform a single band
    of reflectance from digital number (DN) to reflectance units.

    Calculation is
        ref = (refDN + offset) / scaleVal
    """
    ref = (numpy.float32(refDN) + offset) / scaleVal
    return ref


def calcBTthresholds(otherargs):
    """
    Calculate some global thresholds based on the results of the first pass
    """
    # Equation 8
    Twater = scoreatpcnt(otherargs.waterBT_hist, 82.5)
    if Twater is not None:
        Twater = Twater - BT_OFFSET
    # Equation 13
    Tlow = scoreatpcnt(otherargs.clearLandBT_hist, 17.5)
    if Tlow is not None:
        Tlow = Tlow - BT_OFFSET
    Thigh = scoreatpcnt(otherargs.clearLandBT_hist, 82.5)
    if Thigh is not None:
        Thigh = Thigh - BT_OFFSET
    return (Twater, Tlow, Thigh)


#: For scaling probability values so I can store them in 8 bits
PROB_SCALE = 100.0


def doPotentialCloudSecondPass(fmaskFilenames, fmaskConfig, pass1file, 
                Twater, Tlow, Thigh, missingThermal, nonNullCount):
    """
    Second pass for potential cloud layer
    """
    infiles = applier.FilenameAssociations()
    outfiles = applier.FilenameAssociations()
    otherargs = applier.OtherInputs()
    controls = applier.ApplierControls()
    
    infiles.pass1 = pass1file
    infiles.toaref = fmaskFilenames.toaRef
    if not missingThermal:
        infiles.thermal = fmaskFilenames.thermal
    (fd, outfiles.pass2) = tempfile.mkstemp(prefix='pass2', dir=fmaskConfig.tempDir, 
                                    suffix=fmaskConfig.defaultExtension)
    os.close(fd)
    otherargs.refBands = fmaskConfig.bands
    otherargs.thermalInfo = fmaskConfig.thermalInfo
    
    otherargs.Twater = Twater
    otherargs.Tlow = Tlow
    otherargs.Thigh = Thigh
    otherargs.lCloudProb_hist = numpy.zeros(BT_HISTSIZE, dtype=numpy.uint32)
    otherargs.fmaskConfig = fmaskConfig
    
    controls.setWindowXsize(RIOS_WINDOW_SIZE)
    controls.setWindowYsize(RIOS_WINDOW_SIZE)
    controls.setReferenceImage(fmaskFilenames.toaRef)
    controls.setCalcStats(False)
    
    applier.apply(potentialCloudSecondPass, infiles, outfiles, otherargs, controls=controls)
    
    # Equation 17
    # Need at least 3% of nonnull pixels as clear land for this to be reliable. 
    minPixelsReqd = 0.03 * nonNullCount
    if otherargs.lCloudProb_hist.sum() < minPixelsReqd:
        # Almost no clear land pixels
        landThreshold = None
    else:
        landThreshold = scoreatpcnt(otherargs.lCloudProb_hist, 82.5)
    if landThreshold is not None:
        landThreshold = landThreshold / PROB_SCALE + fmaskConfig.Eqn17CloudProbThresh
    else:
        landThreshold = fmaskConfig.Eqn17CloudProbThresh
    return (outfiles.pass2, landThreshold)


def potentialCloudSecondPass(info, inputs, outputs, otherargs):
    """
    Called from RIOS
    
    Second pass of potential cloud layer
    """
    fmaskConfig = otherargs.fmaskConfig
    
    ref = refDNtoUnits(inputs.toaref, fmaskConfig)
    # Clamp off any reflectance <= 0
    ref[ref<=0] = 0.00001
    
    if hasattr(inputs, 'thermal'):
        # Brightness temperature in degrees C
        bt = otherargs.thermalInfo.scaleThermalDNtoC(inputs.thermal)
        
    Twater = otherargs.Twater
    (Tlow, Thigh) = (otherargs.Tlow, otherargs.Thigh)
    # Values from first pass
    clearLand = inputs.pass1[2].astype(bool)
    variabilityProbPcnt = inputs.pass1[3]
    variability_prob = variabilityProbPcnt / PROB_SCALE
    
    # Cirrus band. From Zhu et al 2015, equation 1
    if config.BAND_CIRRUS in otherargs.refBands:
        cirrus = otherargs.refBands[config.BAND_CIRRUS]
        cirrusProb = ref[cirrus] / fmaskConfig.cirrusProbRatio

    # Equation 9
    if Twater is not None:
        wTemperature_prob = (Twater - bt) / 4.0
    else:
        # There is no water, so who cares. 
        wTemperature_prob = 1
    
    swir1 = otherargs.refBands[config.BAND_SWIR1]
    # Equation 10
    brightness_prob = numpy.minimum(ref[swir1], 0.11) / 0.11
    
    # Equation 11
    wCloud_prob = wTemperature_prob * brightness_prob
    # Zhu et al 2015, equation 2
    if config.BAND_CIRRUS in otherargs.refBands:
        wCloud_prob += cirrusProb
    
    # Equation 14
    if Thigh is not None and Tlow is not None:
        lTemperature_prob = (Thigh + 4 - bt) / (Thigh + 4 - (Tlow - 4))
    else:
        # there was no land available for temperature thresholds, so it is probably cloud. 
        lTemperature_prob = 1
    
    # Equation 16
    lCloud_prob = lTemperature_prob * variability_prob
    if config.BAND_CIRRUS in otherargs.refBands:
        lCloud_prob += cirrusProb
    
    outstack = numpy.array([
        (wCloud_prob * PROB_SCALE).clip(BYTE_MIN, BYTE_MAX), 
        (lCloud_prob * PROB_SCALE).clip(BYTE_MIN, BYTE_MAX)], dtype=numpy.uint8)
    outputs.pass2 = outstack
    
    # Accumulate histogram of lCloud_prob
    scaledProb = (lCloud_prob * PROB_SCALE).clip(BYTE_MIN, BYTE_MAX).astype(numpy.uint8)
    otherargs.lCloudProb_hist = accumHist(otherargs.lCloudProb_hist, scaledProb[clearLand])


def doCloudLayerFinalPass(fmaskFilenames, fmaskConfig, pass1file, pass2file, 
                    landThreshold, Tlow, missingThermal):
    """
    Final pass
    """
    infiles = applier.FilenameAssociations()
    outfiles = applier.FilenameAssociations()
    otherargs = applier.OtherInputs()
    controls = applier.ApplierControls()
    
    infiles.pass1 = pass1file
    infiles.pass2 = pass2file
    if not missingThermal:
        infiles.thermal = fmaskFilenames.thermal
    otherargs.landThreshold = landThreshold
    otherargs.Tlow = Tlow
    otherargs.thermalInfo = fmaskConfig.thermalInfo
    otherargs.minCloudSize = fmaskConfig.minCloudSize_pixels
    otherargs.sensor = fmaskConfig.sensor

    (fd, outfiles.cloudmask) = tempfile.mkstemp(prefix='interimcloud', 
        dir=fmaskConfig.tempDir, suffix=fmaskConfig.defaultExtension)
    os.close(fd)
    # Need overlap so we can do Fmask's 3x3 fill-in
    overlap = 1
    # Also need overlap for cloud size filter
    overlap = max(overlap, fmaskConfig.minCloudSize_pixels)
        
    controls.setOverlap(overlap)
    controls.setWindowXsize(RIOS_WINDOW_SIZE)
    controls.setWindowYsize(RIOS_WINDOW_SIZE)
    controls.setReferenceImage(pass1file)
    controls.setCalcStats(False)
    
    applier.apply(cloudFinalPass, infiles, outfiles, otherargs, controls=controls)
    
    return outfiles.cloudmask


def cloudFinalPass(info, inputs, outputs, otherargs):
    """
    Called from RIOS
    
    Final pass of cloud mask layer
    """
    nullmask = inputs.pass1[4].astype(bool)
    pcp = inputs.pass1[0].astype(bool)
    waterTest = inputs.pass1[1].astype(bool)
    notWater = numpy.logical_not(waterTest)
    notWater[nullmask] = False
    wCloud_prob = inputs.pass2[0] / PROB_SCALE
    lCloud_prob = inputs.pass2[1] / PROB_SCALE
    if hasattr(inputs, 'thermal'):
        # Brightness temperature in degrees C
        bt = otherargs.thermalInfo.scaleThermalDNtoC(inputs.thermal)
        
    landThreshold = otherargs.landThreshold
    Tlow = otherargs.Tlow
    
    cloudmask1 = pcp & waterTest & (wCloud_prob>0.5)
    cloudmask2 = pcp & notWater & (lCloud_prob>landThreshold)
    # according to [Zhu 2015] the lCloudprob > 0.99 test should be removed.
    # For now I only disabled it for S2, because it gives a lot of false
    # positives due to missing a thermal band.
    if (otherargs.sensor == config.FMASK_SENTINEL2):
        cloudmask3 = numpy.zeros(cloudmask1.shape, dtype=bool)
    else:
        cloudmask3 = (lCloud_prob > 0.99) & notWater
    if Tlow is not None:
        cloudmask4 = (bt < (Tlow - 35))
    else:
        # Not enough land for final test. Also come here when missing thermal.
        cloudmask4 = numpy.zeros(cloudmask1.shape, dtype=bool)
        
    # Equation 18
    cloudmask = cloudmask1 | cloudmask2 | cloudmask3 | cloudmask4
    cloudmask[nullmask] = 0
    
    # If required, filter out small clouds. 
    if otherargs.minCloudSize > 1:
        (clumps, numClumps) = label(cloudmask)
        clumpSizes = numpy.bincount(clumps.flatten())
        clumpSizes[0] = 0       # Knock out the size of the null area
        sizeImg = clumpSizes[clumps]
        cloudmask[sizeImg < otherargs.minCloudSize] = 0

    # Apply the prescribed 3x3 buffer. According to Zhu&Woodcock (page 87, end of section 3.1.2) 
    # they set a pixel to cloud if 5 or more of its 3x3 neighbours is cloud. 
    # This little incantation will do exactly the same. 
    bufferedCloudmask = (uniform_filter(cloudmask * 2.0, size=3) >= 1.0)

    bufferedCloudmask[nullmask] = 0
    
    outputs.cloudmask = numpy.array([bufferedCloudmask])


def doPotentialShadows(fmaskFilenames, fmaskConfig, NIR_17):
    """
    Make potential shadow layer, as per section 3.1.3 of Zhu&Woodcock. 
    """
    (fd, potentialShadowsFile) = tempfile.mkstemp(prefix='shadows', dir=fmaskConfig.tempDir, 
                                        suffix=fmaskConfig.defaultExtension)
    os.close(fd)

    # convert from numpy (0 based) to GDAL (1 based) indexing
    NIR_lyr = fmaskConfig.bands[config.BAND_NIR] + 1
    
    # Read in whole of band 4
    ds = gdal.Open(fmaskFilenames.toaRef)
    band = ds.GetRasterBand(NIR_lyr)
    nullval = band.GetNoDataValue()
    if nullval is None:
        nullval = 0
    # Sentinel2 is uint16 which causes problems...
    scaledNIR = band.ReadAsArray().astype(numpy.int16)
    # Check for ESA's stoopid offset, for NIR band only
    NIRoffset = 0
    if (fmaskConfig.TOARefDNoffsetDict is not None and 
            config.BAND_NIR in fmaskConfig.TOARefDNoffsetDict):
        NIRoffset = fmaskConfig.TOARefDNoffsetDict[config.BAND_NIR]
    scaleVal = fmaskConfig.TOARefScaling
    NIR_17_dn = NIR_17 * scaleVal - NIRoffset
    
    scaledNIR_filled = fillminima.fillMinima(scaledNIR, nullval, NIR_17_dn)

    NIR = singleRefDNtoUnits(scaledNIR, scaleVal, NIRoffset)
    NIR_filled = singleRefDNtoUnits(scaledNIR_filled, scaleVal, NIRoffset)
    del scaledNIR, scaledNIR_filled
    
    # Equation 19
    potentialShadows = ((NIR_filled - NIR) > fmaskConfig.Eqn19NIRFillThresh)
    
    driver = gdal.GetDriverByName(applier.DEFAULTDRIVERNAME)
    creationOptions = applier.dfltDriverOptions[applier.DEFAULTDRIVERNAME]
    outds = driver.Create(potentialShadowsFile, ds.RasterXSize, ds.RasterYSize, 
                    1, gdal.GDT_Byte, creationOptions)
    proj = ds.GetProjection()
    outds.SetProjection(proj)
    transform = ds.GetGeoTransform()
    outds.SetGeoTransform(transform)
    outband = outds.GetRasterBand(1)
    outband.WriteArray(potentialShadows)
    outband.SetNoDataValue(0)
    del outds

    return potentialShadowsFile


def clumpClouds(cloudmaskfile):
    """
    Clump cloud pixels to make a layer of cloud objects. Currently assumes
    that the cloud mask contains only zeros and ones. 
    """
    ds = gdal.Open(cloudmaskfile)
    band = ds.GetRasterBand(1)
    cloudmask = band.ReadAsArray()
    
    (clumps, numClumps) = label(cloudmask, structure=numpy.ones((3, 3)))

    return (clumps, numClumps)


CLOUD_HEIGHT_SCALE = 10


def make3Dclouds(fmaskFilenames, fmaskConfig, clumps, numClumps, missingThermal):
    """
    Create 3-dimensional cloud objects from the cloud mask, and the thermal 
    information. Assumes a constant lapse rate to convert temperature into height.
    Resulting cloud heights are relative to cloud base. 
    
    Returns an image of relative cloud height (relative to cloud base for 
    each cloud object), and a dictionary of cloud base temperature, for each 
    cloud object, and valueindexes.ValueIndexes object for use in extracting the location of 
    every pixel for a given cloud object. 
    
    """
    # Find out the pixel grid of the toareffile, so we can use that for RIOS.
    # this is necessary because the thermal might be on a different grid,
    # and we can use RIOS to resample that. 
    referencePixgrid = pixelgrid.pixelGridFromFile(fmaskFilenames.toaRef)
    
    infiles = applier.FilenameAssociations()
    outfiles = applier.FilenameAssociations()
    otherargs = applier.OtherInputs()
    controls = applier.ApplierControls()
    
    # if we have thermal, run against that 
    # otherwise we are just 
    if not missingThermal:
        infiles.thermal = fmaskFilenames.thermal
    else:
        infiles.toaRef = fmaskFilenames.toaRef
        
    otherargs.clumps = clumps
    otherargs.cloudClumpNdx = valueindexes.ValueIndexes(clumps, nullVals=[0])
    otherargs.numClumps = numClumps
    otherargs.thermalInfo = fmaskConfig.thermalInfo
    
    # Run RIOS on whole image as one block
    (nRows, nCols) = referencePixgrid.getDimensions()
    controls.setWindowXsize(nCols)
    controls.setWindowYsize(nRows)
    controls.setReferencePixgrid(referencePixgrid)
    controls.setCalcStats(False)
    
    applier.apply(cloudShapeFunc, infiles, outfiles, otherargs, controls=controls)
    
    return (otherargs.cloudShape, otherargs.cloudBaseTemp, otherargs.cloudClumpNdx)


def cloudShapeFunc(info, inputs, outputs, otherargs):
    """
    Called from RIOS.
    
    Calculate the 3d cloud shape image. Requires that RIOS be run on the whole
    image at once, as substantial spatial structure is required. Needed to
    use RIOS because of the possibility of the thermal input being on 
    a different resolution to the cloud input
    
    Returns the result as whole arrays in otherargs, rather than writing them to 
    output files, as they would just be read in again as whole arrays immediately. 
    
    """
    
    cloudBaseTemp = {}
    
    # If we are missing the thermal, then the clouds are flat 2-d shapes.
    if hasattr(inputs, 'thermal'):
        bt = otherargs.thermalInfo.scaleThermalDNtoC(inputs.thermal)
        cloudShape = numpy.zeros(bt.shape, dtype=numpy.uint8)
        
        cloudIDlist = otherargs.cloudClumpNdx.values
        for cloudID in cloudIDlist:
            cloudNdx = otherargs.cloudClumpNdx.getIndexes(cloudID)
            btCloud = bt[cloudNdx]
        
            numPixInCloud = len(cloudNdx[0])
        
            # Equation 22, in several pieces
            R = numpy.sqrt(numPixInCloud / (2 * numpy.pi))
            if R >= 8:
                percentile = 100.0 * (R - 8.0)**2 / (R**2)
                Tcloudbase = scipy.stats.scoreatpercentile(btCloud, percentile)
            else:
                Tcloudbase = btCloud.min()
        
            # Equation 23
            btCloud[btCloud>Tcloudbase] = Tcloudbase
        
            # Equation 24 (relative to cloud base). 
            # N.B. Equation given in paper appears to be wrong, it multiplies by lapse
            # rate instead of dividing by it. 
            LAPSE_RATE_WET = 6.5        # degrees/km
            Htop_relative = (Tcloudbase - btCloud) / LAPSE_RATE_WET
        
            # Put this back into the cloudShape array at the right place
            cloudShape[cloudNdx] = numpy.round(Htop_relative * CLOUD_HEIGHT_SCALE).astype(numpy.uint8)
        
            # Save the Tcloudbase for this cloudID
            cloudBaseTemp[cloudID] = Tcloudbase
    else:
        # fake it
        cloudShape = numpy.zeros(inputs.toaRef[0].shape, dtype=numpy.uint8)
    
    otherargs.cloudShape = cloudShape
    otherargs.cloudBaseTemp = cloudBaseTemp


METRES_PER_KM = 1000.0
BYTES_PER_VOXEL = 4
SOLIDCLOUD_MAXMEM = float(1024 * 1024 * 1024)


def makeCloudShadowShapes(fmaskFilenames, fmaskConfig,
        cloudShape, cloudClumpNdx):
    """
    Project the 3d cloud shapes onto horizontal surface, along the sun vector, to
    make the 2d shape of the shadow. 
    """
    # Read in the two solar angles. Assumes that the angles file is on the same 
    # pixel grid as the cloud, which should always be the case. 
    ds = gdal.Open(fmaskFilenames.toaRef)
    geotrans = ds.GetGeoTransform()
    (xRes, yRes) = (float(geotrans[1]), float(geotrans[5]))
    (nrows, ncols) = (ds.RasterYSize, ds.RasterXSize)
    del ds

    # tell anglesInfo it may need to read data into memory
    fmaskConfig.anglesInfo.prepareForQuerying()
    
    shadowShapesDict = {}
    
    cloudIDlist = cloudClumpNdx.values
    for cloudID in cloudIDlist:
        cloudNdx = cloudClumpNdx.getIndexes(cloudID)
        
        sunAz = fmaskConfig.anglesInfo.getSolarAzimuthAngle(cloudNdx)
        sunZen = fmaskConfig.anglesInfo.getSolarZenithAngle(cloudNdx)
        satAz = fmaskConfig.anglesInfo.getViewAzimuthAngle(cloudNdx)
        satZen = fmaskConfig.anglesInfo.getViewZenithAngle(cloudNdx)
        
        # Cloudtop height of each pixel in cloud, in metres
        cloudHgt = METRES_PER_KM * cloudShape[cloudNdx] / CLOUD_HEIGHT_SCALE
        
        # Relative (x, y) positions of each pixel in the cloud, in metres. Note 
        # that the negative yRes flips the Y axis (which is what we want)
        x = (cloudNdx[1] * xRes)
        y = (cloudNdx[0] * yRes)
        
        maxHgt = numpy.ceil(cloudHgt.max() / xRes) * xRes
        # Ensure we have at least one layer of voxels
        maxHgt = max(maxHgt, xRes)
 
# The following commented-out lines are the rigorous approach to constructing
# the 3-dimensional shape of the cloud. The original paper is not clear about 
# exactly how they code this, so I initially went for the most rigorous
# version, but it turns out to be a memory hog for large individual clouds.
# It constructs a solid 3-d representation of the cloud object, and then projects every 
# voxel along the sun vector to its 2-d position at cloudbase height. This will
# capture the whole shadow shape, but takes up a lot of memory to do it. 
#        # A cubical voxel has the potential to take up large amounts of memory, if the
#        # single cloud is very large, so reduce the voxel height sufficiently to keep 
#        # the memory requirements down
#        maxLayers = int(numpy.ceil(SOLIDCLOUD_MAXMEM / (numPix * BYTES_PER_VOXEL)))
#        voxelHeight = max(maxHgt / maxLayers, xRes)
#
#        # Make a "solid" cloud, with a stack of (cubical) voxels
#        solidCloud = numpy.mgrid[:maxHgt:voxelHeight].astype(numpy.float32)[:, None].repeat(numPix, axis=1)
#        # The z coordinate of every voxel which is inside the cloud, and zero for above cloud
#        z = solidCloud * (solidCloud <= cloudHgt)
#        del solidCloud, cloudHgt
#        
#        d = z * numpy.tan(sunZen, dtype=numpy.float32)
#        del z

        # This is the much less rigorous approach to calculating the projected position
        # of each part of the cloud. It only uses the top of the cloud on each pixel, 
        # and assumes that this is sufficient to capture the whole cloud. For a very 
        # tall thin cloud this might not be true, but I have yet to see an example of 
        # it failing. It uses substantially less memory, so I am going with this for now. 
        d = cloudHgt * numpy.tan(sunZen).astype(numpy.float32)

        # (x', y') are coordinates of each voxel projected onto the plane of the cloud base,
        # for every voxel in the solid cloud
        xDash = x - d * float(numpy.sin(sunAz))
        yDash = y - d * float(numpy.cos(sunAz))
        del d, x, y
        
        # Turn these back into row/col coordinates
        rows = (yDash / yRes).astype(numpy.uint32).clip(0, nrows - 1)
        cols = (xDash / xRes).astype(numpy.uint32).clip(0, ncols - 1)
        
        # Make the row/col arrays have the right shape, and store as a single tuple
        shadowNdx = (rows.flatten(), cols.flatten())
        
        # The row/cols can contain duplicates, as many 3-d points will project
        # into the same 2-d location at cloudbase height. 
        # So, we must remove the duplicates and give them the right shape. It seems
        # that the most efficient way of doing this is to use them as indexes to 
        # set pixels, then get back the indexes of the pixels which were set. Seems
        # a bit roundabout, but I can't think of a better way. 
        # Sadly, I have had to comment this bit out, as it makes the whole
        # thing take several times as long. Obviously I need a better method of removing
        # the duplicates. Sigh.....
        # blankImg[shadowNdx] = True
        # shadowNdx = numpy.where(blankImg)
        # blankImg[shadowNdx] = False
        
        # Stash these shapes in a dictionary, along with the corresponding sun and satellite angles
        shadowShapesDict[cloudID] = (shadowNdx, satAz, satZen, sunAz, sunZen)
    
    # no more querying needed
    fmaskConfig.anglesInfo.releaseMemory()
    
    return shadowShapesDict


def getIntersectionCoords(filelist):
    """
    Use the RIOS utilities to get the correct area of intersection
    for a set of files, although we are not going to read the files using
    RIOS itself, but with GDAL directly. 
    """
    pixgridList = [pixelgrid.pixelGridFromFile(filename) for filename in filelist]
    # Just assume that they are all on the same grid, so the first one is 
    # fine to use as reference grid
    intersectionPixgrid = pixelgrid.findCommonRegion(pixgridList, pixgridList[0])
    # Get top-left in pixel coords for each file
    
    (tlX, tlY) = (intersectionPixgrid.xMin, intersectionPixgrid.yMax)
    tlList = [imageio.wld2pix(pixgrid.makeGeoTransform(), tlX, tlY) 
        for pixgrid in pixgridList]
    tlDict = {}
    for i in range(len(filelist)):
        filename = filelist[i]
        tl = tlList[i]
        tlDict[filename] = (int(tl.x), int(tl.y))
    return (tlDict, intersectionPixgrid)


def makeBufferKernel(buffsize):
    """
    Make a 2-d array for buffering. It represents a circle of 
    radius buffsize pixels, with 1 inside the circle, and zero outside.
    """
    bufferkernel = None
    if buffsize > 0:
        n = 2 * buffsize + 1
        (r, c) = numpy.mgrid[:n, :n]
        radius = numpy.sqrt((r - buffsize)**2 + (c - buffsize)**2)
        bufferkernel = (radius <= buffsize).astype(numpy.uint8)
    return bufferkernel


def matchShadows(fmaskConfig, interimCloudmask, potentialShadowsFile, 
        shadowShapesDict, cloudBaseTemp, Tlow, Thigh, pass1file):
    """
    Match the cloud shadow shapes to the potential cloud shadows. 
    Write an output file of the resulting shadow layer. 
    Includes a 3-pixel buffer on the final shadows. 
    """
    # Do a bunch of fancy footwork to read the same region from the whole
    # raster, of each of three separate rasters. Really RIOS should be able to do this
    # better, but for now I think I need to do it this way. This is only necessary
    # because the potentialShadow file can, on rare occasions, come out a different
    # shape to the other two, due to the thermal having a slightly different number of 
    # rows and columns. Don't know why this happens - less than 0.5% of cases, but still. 
    (topLeftDict, intersectionPixgrid) = getIntersectionCoords([potentialShadowsFile, 
        interimCloudmask, pass1file])
    (nrows, ncols) = intersectionPixgrid.getDimensions()
    
    # Read in whole rasters from the three relevant files. 
    ds = gdal.Open(potentialShadowsFile)
    band = ds.GetRasterBand(1)
    (xoff, yoff) = topLeftDict[potentialShadowsFile]
    potentialShadow = band.ReadAsArray(xoff, yoff, ncols, nrows).astype(bool)
    del ds
    ds = gdal.Open(interimCloudmask)
    band = ds.GetRasterBand(1)
    (xoff, yoff) = topLeftDict[interimCloudmask]
    cloudmask = band.ReadAsArray(xoff, yoff, ncols, nrows).astype(bool)
    geotrans = ds.GetGeoTransform()
    (xRes, yRes) = (geotrans[1], geotrans[5])
    (xsize, ysize) = (ds.RasterXSize, ds.RasterYSize)
    proj = ds.GetProjection()
    del ds
    ds = gdal.Open(pass1file)
    band = ds.GetRasterBand(5)
    (xoff, yoff) = topLeftDict[pass1file]
    nullmask = band.ReadAsArray(xoff, yoff, ncols, nrows).astype(bool)
    del ds

    (fd, interimShadowmask) = tempfile.mkstemp(prefix='matchedshadows', dir=fmaskConfig.tempDir, 
                                        suffix=fmaskConfig.defaultExtension)
    os.close(fd)
    
    shadowmask = numpy.zeros(potentialShadow.shape, dtype=bool)
    
    unmatchedCount = 0
    cloudIDlist = shadowShapesDict.keys()
    for cloudID in cloudIDlist:
        shadowEntry = shadowShapesDict[cloudID]
        if cloudID in cloudBaseTemp:
            Tcloudbase = cloudBaseTemp[cloudID]
        else:
            Tcloudbase = 0

        matchedShadowNdx = matchOneShadow(cloudmask, shadowEntry, potentialShadow, Tcloudbase, 
            Tlow, Thigh, xRes, yRes, cloudID, nullmask)
        
        if matchedShadowNdx is not None:
            shadowmask[matchedShadowNdx] = True
        else:
            unmatchedCount += 1

    if fmaskConfig.verbose:
        print("No shadow found for %s of %s clouds " % (unmatchedCount, len(cloudIDlist)))

    del potentialShadow, cloudmask, nullmask
    
    # Now apply a 3-pixel buffer, as per section 3.2 (2nd-last paragraph)
    # I have the buffer size settable from the commandline, with our default
    # being larger than the original. 
    if fmaskConfig.shadowBufferSize > 0:
        kernel = makeBufferKernel(fmaskConfig.shadowBufferSize)
        shadowmaskBuffered = maximum_filter(shadowmask, footprint=kernel)
    else:
        shadowmaskBuffered = shadowmask

    driver = gdal.GetDriverByName(applier.DEFAULTDRIVERNAME)
    creationOptions = applier.dfltDriverOptions[applier.DEFAULTDRIVERNAME]
    ds = driver.Create(interimShadowmask, xsize, ysize, 1, gdal.GDT_Byte,
                creationOptions)
    ds.SetProjection(proj)
    ds.SetGeoTransform(geotrans)
    band = ds.GetRasterBand(1)
    band.WriteArray(shadowmaskBuffered)
    del ds
    
    return interimShadowmask


def matchOneShadow(cloudmask, shadowEntry, potentialShadow, Tcloudbase, Tlow, Thigh, 
        xRes, yRes, cloudID, nullmask):
    """
    Given the temperatures and sun angles for a single cloud object, and a shadow
    shape, search along the sun vector for a matching shadow object. 
    
    """
    (imgNrows, imgNcols) = cloudmask.shape

    # Not enough clear land to work out temperature thresholds, so guess. 
    if Tlow is None:
        Tlow = 0.0
    if Thigh is None:
        Thigh = 10.0
    
    # Equation 21. Convert these to metres instead of kilometres
    Hcloudbase_min = max(0.2, (Tlow - 4 - Tcloudbase) / 9.8) * METRES_PER_KM
    Hcloudbase_max = min(12, (Thigh + 4 - Tcloudbase)) * METRES_PER_KM
    
    # Entry for this cloud shadow object
    (shapeNdx, satAz, satZen, sunAz, sunZen) = shadowEntry
    
    tanSunZen = numpy.tan(sunZen)
    sinSunAz = numpy.sin(sunAz)
    cosSunAz = numpy.cos(sunAz)
    tanSatZen = numpy.tan(satZen)
    sinSatAz = numpy.sin(satAz)
    cosSatAz = numpy.cos(satAz)
    
    # We want to shift the cloud up, from Hcloudbase_min to Hcloudbase_max.
    # Given the sun angles, this corresponds to shifting the shadow along
    # the ground from Dmin to Dmax. 
    Dmin = Hcloudbase_min * tanSunZen
    Dmax = Hcloudbase_max * tanSunZen
    
    # This corresponds to the following offsets in X and Y
    Xoff_min = Dmin * sinSunAz
    Xoff_max = Dmax * sinSunAz
    Yoff_min = Dmin * cosSunAz
    Yoff_max = Dmax * cosSunAz
    
    # We want the step to be xRes in at least one direction. 
    longestShift = max(abs(Xoff_max - Xoff_min), abs(Yoff_max - Yoff_min))
    numSteps = max(1, int(numpy.ceil(longestShift / xRes)))      # Assumes square pixels
    Xstep = (Xoff_max - Xoff_min) / numSteps
    Ystep = (Yoff_max - Yoff_min) / numSteps
    
    # shadowTemplate is a rectangle containing just the shadow shape to be shifted
    row0 = shapeNdx[0].min()
    rowN = shapeNdx[0].max()
    col0 = shapeNdx[1].min()
    colN = shapeNdx[1].max()
    (nrows, ncols) = ((rowN - row0 + 1), (colN - col0 + 1))
    shadowTemplate = numpy.zeros((nrows, ncols), dtype=bool)
    shadowTemplate[shapeNdx[0] - row0, shapeNdx[1] - col0] = True
    
    # Step this template across the potential shadows until we match. 
    i = 0
    bestSimilarity = 0
    bestRC = (0, 0)
    bestOverlapRegion = None
    for i in range(numSteps):
        # Cloudbase height for this step
        H = (Xoff_min + i * Xstep) / (tanSunZen * sinSunAz)
        # Calculate the shift in the cloud position due to the view angle and the cloud elevation
        D_viewoffset = H * tanSatZen
        X_viewoffset = D_viewoffset * sinSatAz
        Y_viewoffset = D_viewoffset * cosSatAz
        
        # Shadow shift in metres
        Xoff = Xoff_min + i * Xstep - X_viewoffset
        Yoff = Yoff_min + i * Ystep - Y_viewoffset
        
        # Shift in pixels. Note that negative yRes inverts the row axis
        rowOff = int(Yoff / yRes)
        colOff = int(Xoff / xRes)
        
        # Extract the potential shadow, and also the cloud, from the shifted region of
        # the full images
        r = row0 - rowOff
        c = col0 - colOff
        if r >= 0 and r + nrows <= imgNrows and c >= 0 and c + ncols <= imgNcols:
            cloud = cloudmask[r:r + nrows, c:c + ncols]
            potShadow = potentialShadow[r:r + nrows, c:c + ncols]
            null = nullmask[r:r + nrows, c:c + ncols]
            # mask the potential shadow layer, so we don't include anything we think is cloud
            potShadow[cloud] = 0
            # Similarly with the areas which are null in the imagery
            potShadow[null] = 0

            # Mask the shadow template with the cloud from this area
            shadowTemplateMasked = shadowTemplate.copy()
            shadowTemplateMasked[cloud] = False
            shadowTemplateMasked[null] = False

            similarity = 0
            overlap = numpy.logical_and(potShadow, shadowTemplateMasked)
            # Calculate overlap area (by counting pixels)
            overlapArea = overlap.sum()
            # Remaining area of shadow shape
            shadowArea = shadowTemplateMasked.sum()
            if shadowArea > 0:
                similarity = float(overlapArea) / shadowArea

            # We don't use the Zhu & Woodcock termination condition, as this
            # very often results in stopping search too soon. We just check the whole
            # transect, and save the best position. 
            # TODO: strict version should use new threshold
            if similarity > bestSimilarity:
                bestRC = (r, c)
                bestSimilarity = similarity
                bestOverlapRegion = overlap
    
    if bestSimilarity > 0.3:
        # We accept the match, now save the index for the pixels in the overlap region
        overlapNdx = numpy.where(bestOverlapRegion)
        matchedShadowNdx = (bestRC[0] + overlapNdx[0], bestRC[1] + overlapNdx[1])
    else:
        matchedShadowNdx = None
    
    return matchedShadowNdx


def finalizeAll(fmaskFilenames, fmaskConfig, interimCloudmask, interimShadowmask, 
        pass1file):
    """
    Use the cloud and shadow masks to mask the snow layer (as per Zhu & Woodcock). 
    Apply the optional extra buffer to the cloud mask, and write to final file. 
    """
    infiles = applier.FilenameAssociations()
    outfiles = applier.FilenameAssociations()
    otherargs = applier.OtherInputs()
    controls = applier.ApplierControls()
    
    infiles.cloud = interimCloudmask
    infiles.shadow = interimShadowmask
    infiles.pass1 = pass1file
    outfiles.out = fmaskFilenames.outputMask
    controls.setOverlap(fmaskConfig.cloudBufferSize)
    controls.setThematic(True)
    controls.setStatsIgnore(OUTCODE_NULL)
    controls.setWindowXsize(RIOS_WINDOW_SIZE)
    controls.setWindowYsize(RIOS_WINDOW_SIZE)
    controls.setOutputDriverName(fmaskConfig.gdalDriverName)
    
    if fmaskConfig.cloudBufferSize > 0:
        otherargs.bufferkernel = makeBufferKernel(fmaskConfig.cloudBufferSize)

    applier.apply(maskAndBuffer, infiles, outfiles, otherargs, controls=controls)
    
    rat.setColorTable(outfiles.out, numpy.array([[2, 255, 0, 255, 255],
                                                 [3, 255, 255, 0, 255],
                                                 [4, 85, 255, 255, 255],
                                                 [5, 0, 0, 255, 255]]))

    usingExceptions = gdal.GetUseExceptions()
    gdal.UseExceptions()
    try:
        rat.writeColumn(outfiles.out, "Classification", [b"Null", b"Valid", b"Cloud", 
                                                    b"Cloud Shadow", b"Snow", b"Water"])
    except Exception:
        # Failed to write the RAT, probably because the selected format does not support it. 
        # Just ignore it silently
        pass
    finally:
        if not usingExceptions:
            gdal.DontUseExceptions()


def maskAndBuffer(info, inputs, outputs, otherargs):
    """
    Called from RIOS
    
    Apply cloud and shadow masks to snow layer, and buffer cloud layer
    
    The main aims of all this re-masking are:
        1) A pixel should be either cloud, shadow, snow or not, but never 
           more than one
        2) Areas which are null in the input imagery should be null in the 
           mask, even after buffering, etc. 
    
    """
    snow = inputs.pass1[5].astype(bool)
    nullmask = inputs.pass1[4].astype(bool)
    resetNullmask = nullmask

    cloud = inputs.cloud[0].astype(bool)
    shadow = inputs.shadow[0].astype(bool)
    water = inputs.pass1[1].astype(bool)
    
    # Buffer the cloud
    if hasattr(otherargs, 'bufferkernel'):
        cloud = maximum_filter(cloud, footprint=otherargs.bufferkernel)
    
    # now convert these masks to
    # 0 - null
    # 1 - not null and not mask
    # 2 - cloud
    # 3 - cloud shadow
    # 4 - snow
    # 5 - water
    out = numpy.full(cloud.shape, fill_value=OUTCODE_CLEAR, dtype=numpy.uint8)
    out[water] = OUTCODE_WATER
    out[snow] = OUTCODE_SNOW
    out[shadow] = OUTCODE_SHADOW
    out[cloud] = OUTCODE_CLOUD
    out[resetNullmask] = OUTCODE_NULL
    
    outputs.out = numpy.array([out])


def focalVariance(img, winSize):
    """
    Calculate the focal variance of the given 2-d image, over a moving window of
    size winSize pixels.

    """
    img32 = img.astype(numpy.float32)
    focalMean = uniform_filter(img32, size=winSize)
    meanSq = uniform_filter(img32**2, size=winSize)
    variance = meanSq - focalMean**2
    return variance


def calcCDI(ref, fmaskConfig, refBands):
    """
    Calculate the Cloud Displacement Index, as per Frantz et al (2018).
    """
    # Equations 5 & 6
    ratio8a8 = ref[refBands[config.BAND_NIR]] / ref[refBands[config.BAND_S2CDI_NIR8A]]
    ratio8a7 = ref[refBands[config.BAND_S2CDI_NIR7]] / ref[refBands[config.BAND_S2CDI_NIR8A]]

    # Equation 7
    v8a8 = focalVariance(ratio8a8, fmaskConfig.sen2cdiWindow)
    v8a7 = focalVariance(ratio8a7, fmaskConfig.sen2cdiWindow)
    
    # Mask out where we would divide by zero
    cdi = numpy.zeros(v8a7.shape, dtype=numpy.float32)
    divOK = ((v8a7 + v8a8) != 0)
    cdi[divOK] = (v8a7[divOK] - v8a8[divOK]) / (v8a7[divOK] + v8a8[divOK])

    return (ratio8a8, ratio8a7, v8a8, v8a7, cdi)


def deleteRaster(filename):
    """
    Use GDAL's Driver.Delete() method to fully delete the given raster file
    """
    ds = gdal.Open(filename)
    if ds is not None:
        drvr = ds.GetDriver()
        del ds
        drvr.Delete(filename)
