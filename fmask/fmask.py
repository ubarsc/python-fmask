#!/usr/bin/env python
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
    
    
The notation and variable names are largely taken from the paper. Equation
numbers are also from the paper. 

Input is a radiance file (the da1 file, in our terms), and outputs
are cloud, cloud shadow and snow mask files. 

Taken from Neil Flood's implementation by permission.
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
import abc
import subprocess
import tempfile

import numpy
numpy.seterr(all='raise')
from osgeo import gdal
from scipy.ndimage import uniform_filter, maximum_filter, label
import scipy.stats
import scipy.constants

from rios import applier
from rios import pixelgrid
from rios import rat
from rios import imageio

#from . import fillminima

"""
Some constants for the various visible bands used in fmask.
The wavelength numbers are approximate - they move about
for different sensors.
"""
VIS_BAND_045um = 0
VIS_BAND_052um = 1
VIS_BAND_063um = 2
VIS_BAND_076um = 3
VIS_BAND_136um = 4 # Sentinel2 only
VIS_BAND_155um = 5
VIS_BAND_208um = 6

"""
Some standard file configurations for different sensors.
Assumed that panchromatic + thermal bands stored in separate files.
"""
LANDSAT5_TM_BANDS = {VIS_BAND_045um:0, VIS_BAND_052um:1, VIS_BAND_063um:2,
                    VIS_BAND_076um:3, VIS_BAND_155um:4, VIS_BAND_208um:5}
LANDSAT7_ETM_BANDS = LANDSAT5_TM_BANDS
LANDSAT8_OLI_BANDS = {VIS_BAND_045um:1, VIS_BAND_052um:2, VIS_BAND_063um:3,
                    VIS_BAND_076um:4, VIS_BAND_155um:5, VIS_BAND_208um:6}
SENTINEL2_BANDS = {VIS_BAND_045um:1, VIS_BAND_052um:2, VIS_BAND_063um:3,
                    VIS_BAND_076um:8, VIS_BAND_136um:10, VIS_BAND_155um:11,
                    VIS_BAND_208um:12}

class FmaskException(Exception):
    "An exception rasied by Fmask"
    
class FmaskParameterError(FmaskException):
    "Something is wrong with a parameter"

class InfileError(FmaskException): pass

def setDefaultDriver():
    """
    Sets some default values into global variables, defining
    what defaults we should use for GDAL driver. On any given
    output file these can be over-ridden, and can be over-ridden globally
    using the environment variables 
    * $FMASK_DFLT_DRIVER
    * $FMASK_DFLT_DRIVEROPTIONS
    
    If FMASK_DFLT_DRIVER is set, then it should be a gdal short driver name
    If FMASK_DFLT_DRIVEROPTIONS is set, it should be a space-separated list
    of driver creation options, e.g. "COMPRESS=LZW TILED=YES", and should
    be appropriate for the selected GDAL driver. This can also be 'None'
    in which case an empty list of creation options is passed to the driver.
    
    If not otherwise supplied, the default is to use the HFA driver, with compression. 
    """                                                            
    global DEFAULTDRIVERNAME, DEFAULTCREATIONOPTIONS, DEFAULTEXTENSION
    DEFAULTDRIVERNAME = os.getenv('FMASK_DFLT_DRIVER', default='HFA')
    DEFAULTCREATIONOPTIONS = ['COMPRESSED=TRUE','IGNOREUTM=TRUE']
    creationOptionsStr = os.getenv('FMASK_DFLT_DRIVEROPTIONS')
    if creationOptionsStr is not None:
        DEFAULTCREATIONOPTIONS = creationOptionsStr.split()
        
    driver = gdal.GetDriverByName(DEFAULTDRIVERNAME)
    if driver is None:
        msg = 'Cannot find GDAL driver %s' % DEFAULTDRIVERNAME
        raise FmaskParameterError(msg)
        
    DEFAULTEXTENSION = '.tmp'
    drivermeta = driver.GetMetadata()
    if gdal.DMD_EXTENSION in drivermeta:
        DEFAULTEXTENSION = '.' + drivermeta[gdal.DMD_EXTENSION]
        
setDefaultDriver()

class ThermalFileInfo(object):
    """
    Contains parameters for interpreting thermal file.
    
    """
    # TODO: way of reading in/defaulting for different sensors
    thermalFile = None
    thermalBand1040um = None
    thermalGain1040um = None
    thermalOffset1040um = None
    thermalK1_1040um = None
    thermalK2_1040um = None

    def scaleThermalDNtoC(self, scaledBT):
        """
        Use the given params to unscale the thermal, and then 
        convert it from K to C. Return a single 2-d array of the 
        temperature in deg C. 
        """
        KELVIN_ZERO_DEGC = scipy.constants.zero_Celsius
        rad = (scaledBT[self.thermalBand1040um].astype(float) * 
                    self.thermalGain1040um + self.thermalOffset1040um)
        # see http://www.yale.edu/ceo/Documentation/Landsat_DN_to_Kelvin.pdf
        # and https://landsat.usgs.gov/Landsat8_Using_Product.php
        rad[rad <= 0] = 0.00001 # to stop errors below
        temp = self.thermalK2_1040um / numpy.log(self.thermalK1_1040um / rad + 1.0)
        bt = temp - KELVIN_ZERO_DEGC
        return bt
    
class AnglesInfo(object):
    """
    Abstract base class that Contains view and solar angle 
    information for file (in radians).
    
    """
    __metaclass__ = abc.ABCMeta
    
    def prepareForQuerying(self):
        """
        Called when fmask is about to query this object for angles.
        Derived class should do any reading of files into memory required here.
        """
        
    def releaseMemory(self):
        """
        Called when fmask has finished querying this object.
        Can release any allocated memory.
        """
    
    @abc.abstractmethod
    def getSolarZenithAngle(self, indices):
        """
        Return the average solar zenith angle for the given indices
        """

    @abc.abstractmethod
    def getSolarAzimuthAngle(self, indices):
        """
        Return the average solar azimuth angle for the given indices
        """
    
    @abc.abstractmethod
    def getViewZenithAngle(self, indices):
        """
        Return the average view zenith angle for the given indices
        """

    @abc.abstractmethod
    def getViewAzimuthAngle(self, indices):
        """
        Return the average view azimuth angle for the given indices
        """
    
class AnglesFileInfo(AnglesInfo):
    """
    An implementation of AnglesInfo that reads the information from
    GDAL supported files.
    """
    def __init__(self, solarZenithFilename, solarZenithBand, solarAzimuthFilename,
            solarAzimuthBand, viewZenithFilename, viewZenithBand, 
            viewAzimuthFilename, viewAzimuthBand):
        """
        Initialises the object with the names and band numbers of the angles.
        band numbers should be 0 based - ie first band is 0.
        """
        self.solarZenithFilename = solarZenithFilename
        self.solarZenithBand = solarZenithBand
        self.solarAzimuthFilename = solarAzimuthFilename
        self.solarAzimuthBand = solarAzimuthBand
        self.viewZenithFilename = viewZenithFilename
        self.viewZenithBand = viewZenithBand
        self.viewAzimuthFilename = viewAzimuthFilename
        self.viewAzimuthBand = viewAzimuthBand
        # these will contain the actual image data once read
        # by prepareForQuerying()
        self.solarZenithData = None
        self.solarAzimuthData = None
        self.viewZenithData = None
        self.viewAzimuthData = None
    
    @staticmethod
    def readData(filename, bandNum):
        ds = gdal.Open(filename)
        band = ds.GetRasterBand(bandNum + 1)
        data = band.ReadAsArray()
        del ds
        return data
    
    def prepareForQuerying(self):
        """
        Called when fmask is about to query this object for angles.
        """
        self.solarZenithData = self.readData(self.solarZenithFilename, 
                                self.solarZenithFilename)
        self.solarAzimuthData = self.readData(self.solarAzimuthFilename, 
                                self.solarAzimuthFilename)
        self.viewZenithData = self.readData(self.viewZenithFilename, 
                                self.viewZenithFilename)
        self.viewAzimuthData = self.readData(self.viewAzimuthFilename, 
                                self.viewAzimuthFilename)
        
    def releaseMemory(self):
        """
        Called when fmask has finished querying this object.
        """
        del self.solarZenithData
        del self.solarAzimuthData
        del self.viewZenithData
        del self.viewAzimuthData
    
    def getSolarZenithAngle(self, indices):
        """
        Return the average solar zenith angle for the given indices
        """
        return self.solarZenithData[indices].mean()

    def getSolarAzimuthAngle(self, indices):
        """
        Return the average solar azimuth angle for the given indices
        """
        return self.solarAzimuthData[indices].mean()
    
    def getViewZenithAngle(self, indices):
        """
        Return the average view zenith angle for the given indices
        """
        return self.viewZenithData[indices].mean()

    def getViewAzimuthAngle(self, indices):
        """
        Return the average view azimuth angle for the given indices
        """
        return self.viewAzimuthData[indices].mean()

class AngleConstantInfo(AnglesInfo):
    """
    An implementation of AnglesInfo that uses constant
    angles accross the scene. 
    """
    def __init__(self, solarZenithAngle, solarAzimuthAngle, viewZenithAngle,
                    viewAzimuthAngle):
        self.solarZenithAngle = solarZenithAngle
        self.solarAzimuthAngle = solarAzimuthAngle
        self.viewZenithAngle = viewZenithAngle
        self.viewAzimuthAngle = viewAzimuthAngle

    def getSolarZenithAngle(self, indices):
        """
        Return the solar zenith angle
        """
        return self.solarZenithAngle

    def getSolarAzimuthAngle(self, indices):
        """
        Return the solar azimuth angle
        """
        return self.solarAzimuthAngle
    
    def getViewZenithAngle(self, indices):
        """
        Return the view zenith angle
        """
        return self.viewZenithAngle

    def getViewAzimuthAngle(self, indices):
        """
        Return the view azimuth angle
        """
        return self.viewAzimuthAngle

# TODO: way of reading in for different sensors.
def readAnglesFromLandsatMTL(mtlfile):
    raise NotImplementedError()
    
    
def doFmask(radianceFile, radianceBands, toaRefFile, anglesInfo, outMask, 
                outDriver=DEFAULTDRIVERNAME,
                outCreationOpions=DEFAULTCREATIONOPTIONS, thermalFile=None, 
                thermalInfo=None,
                keepintermediates=False, cloudbuffersize=5, shadowbuffersize=10,
                nosaturationtest=False, verbose=False, strictfmask=False, 
                tempdir='.'):
    """
    Main routine for whole Fmask algorithm. Calls all other routines in sequence. 
    Parameters:
    
    * **radianceFile** must be path to a GDAL readable file with all the visible bands.
    * **radianceBands** should be a dictionary of wavelengths and bands. See VIS_BAND_* and LANDSAT5_TM_BANDS etc.
    * **toaRefFile** a file containing Top of Atmosphere reflectance * 1000 with the same bands as radianceFile. See fmask.toaref
    * **anglesInfo** an instance of AnglesInfo.
    * **outMask** name of output cloud mask
    * **outDriver** name of GDAL driver to use for output - see setDefaultDriver() for discussion of defaults.
    * **outCreationOpions** list of creation options for output - see setDefaultDriver() for discussion of defaults.
    * **thermalFile** input thermal file. Set to None if not available.
    * **thermalInfo** an instance of ThermalFileInfo.
    * **keepintermediates** set to True to prevent intermediate files from being deleted. A dictionary will be returned.
    * **cloudbuffersize** Extra buffer of this many pixels on cloud layer
    * **shadowbuffersize** Buffer of this many pixels on shadow layer
    * **nosaturationtest** Omit extra test for saturated cloud pixels (i.e. use only strict Fmask algorithm)
    * **verbose** Print informative messages
    * **strictfmask** Set whatever options are necessary to run strictly as per Fmask paper (Zhu & Woodcock)
    * **tempdir** Temp directory to use
    
    """
    if radianceFile is None or radianceBands is None:
        msg = 'Must provide input radianceFile and radianceBands'
        raise FmaskParameterError(msg)
        
    if toaRefFile is None:
        msg = 'Must provide input toaRefFile'
        raise FmaskParameterError(msg)
        
    if anglesInfo is None:
        msg = 'Must provide input anglesInfo'
        raise FmaskParameterError(msg)
        
    if outMask is None:
        msg = 'Output filename must be provided'
        raise FmaskParameterError(msg)
        
    if strictfmask:
        cloudbuffersize = 0
        shadowbuffersize = 3
        nosaturationtest = True
    
    if verbose: print("Cloud layer, pass 1")
    (pass1file, Twater, Tlow, Thigh, b4_17) = doPotentialCloudFirstPass(radiancefile, 
        radianceBands, toaRefFile, thermalFile, thermalInfo, nosaturationtest, tempdir)
 
    if verbose: print("Cloud layer, pass 2")
    (pass2file, landThreshold) = doPotentialCloudSecondPass(toaRefFile, radianceBands, 
        thermalFile, thermalInfo, pass1file, Twater, Tlow, Thigh, tempdir)

    if verbose: print("Cloud layer, pass 3")
    interimCloudmask = doCloudLayerFinalPass(thermalFile, thermalInfo, pass1file, 
        pass2file, landThreshold, Tlow, tempdir)
        
    if verbose: print("Potential shadows")
    potentialShadowsFile = doPotentialShadows(toaRefFile, radianceBands, b4_17, tempdir)
    
    if verbose: print("Clumping clouds")
    (clumps, numClumps) = clumpClouds(interimCloudmask)
    
    if verbose: print("Making 3d clouds")
    (cloudShape, cloudBaseTemp, cloudClumpNdx) = make3Dclouds(clumps, numClumps, 
        thermalFile, thermalInfo, toaRefFile, radianceBands)
    
    if verbose: print("Making cloud shadow shapes")
    shadowShapesDict = makeCloudShadowShapes(toaRefFile, radianceBands, cloudShape, 
        cloudClumpNdx, anglesInfo)
    
    if verbose: print("Matching shadows")
    interimShadowmask = matchShadows(interimCloudmask, potentialShadowsFile, shadowShapesDict, 
        cloudBaseTemp, Tlow, Thigh, cmdargs, pass1file, shadowbuffersize)
    
    finalizeAll(interimCloudmask, interimShadowmask, pass1file, outfile, cloudbuffersize)
    
    # Remove temporary files
    if not cmdargs.keepintermediates:
        for filename in [pass1file, pass2file, interimCloudmask, potentialShadowsFile,
                interimShadowmask, toareffile]:
            os.remove(filename)

    if verbose: print('finished fmask')

def getFiles(cmdargs):
    """
    Work out the other files we need
    """
    sensor = lcrfs.lcrfs('sensor', cmdargs.radiancefile)
    thsensor = sensor[:-2] + 'th'
    thermalfile = lcrfs.change(cmdargs.radiancefile, 'sensor', thsensor)
    print('thermalfile', thermalfile)
    
    anglesfile = lcrfs.change(cmdargs.radiancefile, 'stage', 'rect')
    anglesfile = lcrfs.change(anglesfile, 'ext', 'ang')
    resensor = sensor[:-2] + 're'
    anglesfile = lcrfs.change(anglesfile, 'sensor', resensor)
    for zone in (58, 59, 60, 61):
        anglesfile = lcrfs.change(anglesfile, 'proj', 'utm%d' % zone)
        if os.path.exists(anglesfile):
            break
    calfile = lcrfs.change(anglesfile, 'ext', 'cal')
    mtlfile = lcrfs.change(anglesfile, 'ext', 'mtl')
    
    return (thermalfile, anglesfile, calfile, mtlfile)

RADIANCE_MULT = 'RADIANCE_MULT_BAND_%s'
RADIANCE_ADD = 'RADIANCE_ADD_BAND_%s'
K1_CONST = 'K1_CONSTANT_BAND_%s'
K2_CONST = 'K2_CONSTANT_BAND_%s'

# band numbers in mtl file for gain and offset for thermal
TH_BAND_NUM_DICT = {'landsat4_tmth' : '6', 
        'landsat5_tmth' : '6',
        'landsat7_etmth' : '6_VCID_1',
        'landsat8_olith' : '10'}

# for some reason L4, 5, and 7 don't
# have these numbers in the mtl file, but L8 does
# from http://www.yale.edu/ceo/Documentation/Landsat_DN_to_Kelvin.pdf
K1_DICT = {'tmth' : 607.76, 'etmth' : 666.09}
K2_DICT = {'tmth' : 1260.56, 'etmth' : 1282.71}

def readMTLFile(mtl):
    """
    Very simple reader that just creates a dictionary
    of key and values
    """
    dict = {}
    for line in open(mtl):
        arr = line.split('=')
        if len(arr) == 2:
            (key, value) = arr
            dict[key.strip()] = value.replace('"', '').strip()

    return dict

def getThermalGainAndOffset(mtl, thermalFile):
    mtlData = readMTLFile(mtl)
    satellite = lcrfs.lcrfs('satellite', thermalFile)
    sensor = lcrfs.lcrfs('sensor', thermalFile)
    satellite_sensor = satellite + '_' + sensor

    band = TH_BAND_NUM_DICT[satellite_sensor]

    s = RADIANCE_MULT % band
    gain = float(mtlData[s])

    s = RADIANCE_ADD % band
    offset = float(mtlData[s])

    s = K1_CONST % band
    if s in mtlData:
        k1 = float(mtlData[s])
    else:
        k1 = K1_DICT[sensor]

    s = K2_CONST % band
    if s in mtlData:
        k2 = float(mtlData[s])
    else:
        k2 = K2_DICT[sensor]

    return gain, offset, k1, k2

# An offset so we can scale BT(deg C) to the range 0-255, for use in histograms. 
BT_OFFSET = 176
BT_HISTSIZE = 256
BYTE_MIN = 0
BYTE_MAX = 255
# Gain to scale b4 reflectances to 0-255 for histograms
# TODO: check this ok for Sentinel
B4_SCALE = 500.0

# Global RIOS window size
RIOS_WINDOW_SIZE = 512

def doPotentialCloudFirstPass(radiancefile, radianceBands, toaRefFile, 
                thermalFile, thermalInfo, nosaturationtest, tempdir):
    """
    Run the first pass of the potential cloud layer. Also
    finds the temperature thresholds which will be needed 
    in the second pass, because it has the relevant data handy. 
    
    """
    infiles = applier.FilenameAssociations()
    outfiles = applier.FilenameAssociations()
    otherargs = applier.OtherInputs()
    controls = applier.ApplierControls()
    
    # TODO: handle no thermal
    infiles.radiance = radiancefile
    infiles.toaref = toareffile
    infiles.thermal = thermalfile
    (fd, outfiles.pass1) = tempfile.mkstemp(prefix='pass1', dir=tempdir, 
                                suffix=DEFAULTEXTENSION)
    os.close(fd)
    controls.setWindowXsize(RIOS_WINDOW_SIZE)
    controls.setWindowYsize(RIOS_WINDOW_SIZE)
    controls.setReferenceImage(toareffile)
    controls.setCalcStats(False)
    controls.setOutputDriverName(DEFAULTDRIVERNAME)
    controls.setCreationOptions(DEFAULTCREATIONOPTIONS)

    otherargs.radianceBands = radianceBands    
    otherargs.thermalInfo = thermalInfo
    otherargs.waterBT_hist = numpy.zeros(BT_HISTSIZE, dtype=numpy.uint32)
    otherargs.clearLandBT_hist = numpy.zeros(BT_HISTSIZE, dtype=numpy.uint32)
    otherargs.clearLandB4_hist = numpy.zeros(BT_HISTSIZE, dtype=numpy.uint32)
    otherargs.saturationtest = not nosaturationtest

    applier.apply(potentialCloudFirstPass, infiles, outfiles, otherargs, controls=controls)
    
    (Twater, Tlow, Thigh) = calcBTthresholds(otherargs)
    
    # 17.5 percentile of band 4, for clear land pixels. Used later in shadow masking. 
    b4_17 = scoreatpcnt(otherargs.clearLandB4_hist, 17.5)
    if b4_17 is not None:
        b4_17 = b4_17 / B4_SCALE
    else:
        # Not enough land to work this out, so guess a low value. 
        b4_17 = 0.01
    
    return (outfiles.pass1, Twater, Tlow, Thigh, b4_17)


def potentialCloudFirstPass(info, inputs, outputs, otherargs):
    """
    Called from RIOS. 
    
    Calculate the first pass potential cloud layer (equation 6)
        
    """
    ref = inputs.toaref.astype(numpy.float) / 1000.0
    # Clamp off any reflectance <= 0
    ref[ref<=0] = 0.00001
    # Brightness temperature in degrees C
    bt = otherargs.thermalInfo.scaleThermalDNtoC(inputs.thermal)

    # Extract the bands we need - expressed as Landsat Bands since
    # this is how the original paper was.
    TM1 = otherinputs.radianceInfo[VIS_BAND_045um]
    TM2 = otherinputs.radianceInfo[VIS_BAND_052um]
    TM3 = otherinputs.radianceInfo[VIS_BAND_063um]
    TM4 = otherinputs.radianceInfo[VIS_BAND_076um]
    TM5 = otherinputs.radianceInfo[VIS_BAND_155um]
    TM7 = otherinputs.radianceInfo[VIS_BAND_208um]
    THERM = otherargs.thermalInfo.thermalBand1040um
    
    radNull = info.getNoDataValueFor(inputs.radiance)
    thermalNull = info.getNoDataValueFor(inputs.thermal)
    nullmask = ((inputs.radiance[TM1] == radNull) | (inputs.thermal[THERM] == thermalNull))
    # Special mask needed only for resets in final pass
    refNullmask = (inputs.radiance[TM1] == radNull)
    thermNullmask = (inputs.thermal[THERM] == thermalNull)
    
    imgshape = bt.shape
    
    # Equation 1
    ndsi = (ref[TM2] - ref[TM5]) / (ref[TM2] + ref[TM5])
    ndvi = (ref[TM4] - ref[TM3]) / (ref[TM4] + ref[TM3])
    basicTest = (ref[TM7] > 0.03) & (bt < 27) & (ndsi < 0.8) & (ndvi < 0.8)
    
    # Equation 2
    meanVis = (ref[TM1] + ref[TM2] + ref[TM3]) / 3.0
    whiteness = numpy.zeros(imgshape)
    for n in [TM1, TM2, TM3]:
        whiteness = whiteness + numpy.absolute((ref[n] - meanVis) / meanVis)
    whitenessTest = (whiteness < 0.7)
    
    # Haze test, equation 3
    hazeTest = ((ref[TM1] - 0.5 * ref[TM3] - 0.08) > 0)
    
    # Equation 4
    b45test = ref[TM4] / ref[TM5] > 0.75
    
    # Equation 5
    waterTest = numpy.logical_or(
        numpy.logical_and(ndvi < 0.01, ref[TM4] < 0.11),
        numpy.logical_and(ndvi < 0.1, ref[TM4] < 0.05)
    )
    waterTest[nullmask] = False
    
    # Equation 6. Potential cloud pixels (first pass)
    pcp = basicTest & whitenessTest & hazeTest & b45test
    
    # This is an extra saturation test added by DERM, and is not part of the Fmask algorithm. 
    # However, some cloud centres are saturated,and thus fail the whiteness and haze tests
    if otherargs.saturationtest:
        saturatedVis = reduce(numpy.logical_or, [(inputs.radiance[n] == 255) for n in [TM1, TM2, TM3]])
        veryBright = (meanVis > 0.45)
        saturatedAndBright = saturatedVis & veryBright
        pcp[saturatedAndBright] = True
        whiteness[saturatedAndBright] = 0
    
    pcp[nullmask] = False
    
    # Equation 7
    clearSkyWater = numpy.logical_and(waterTest, ref[TM7] < 0.03)
    clearSkyWater[nullmask] = False
    
    # Equation 12
    clearLand = numpy.logical_and(numpy.logical_not(pcp), numpy.logical_not(waterTest))
    clearLand[nullmask] = False
    
    # Equation 15
    # Need to modify ndvi/ndsi by saturation......
    saturatedBand2 = (inputs.radiance[TM2] == 255)
    saturatedBand3 = (inputs.radiance[TM3] == 255)
    modNdvi = numpy.where(saturatedBand3, 0, ndvi)
    modNdsi = numpy.where(saturatedBand2, 0, ndsi)
    # Maximum of three indices
    maxNdx = numpy.absolute(modNdvi)
    maxNdx = numpy.maximum(maxNdx, numpy.absolute(modNdsi))
    maxNdx = numpy.maximum(maxNdx, whiteness)
    variabilityProb = 1 - maxNdx
    variabilityProb[nullmask] = 0
    variabilityProbPcnt = numpy.round(variabilityProb * PROB_SCALE)
    variabilityProbPcnt = variabilityProbPcnt.clip(BYTE_MIN, BYTE_MAX).astype(numpy.uint8)
    
    # Equation 20
    snowmask = (ndsi > 0.15) & (bt < 3.8) & (ref[TM4] > 0.11) & (ref[TM2] > 0.1)
    snowmask[nullmask] = False
    
    # Output the pcp and water test layers. 
    outputs.pass1 = numpy.array([pcp, waterTest, clearLand, variabilityProbPcnt, 
        nullmask, snowmask, refNullmask, thermNullmask])
    
    # Accumulate histograms of temperature for land and water separately
    scaledBT = (bt + BT_OFFSET).clip(0, BT_HISTSIZE)
    otherargs.waterBT_hist = accumHist(otherargs.waterBT_hist, scaledBT[clearSkyWater])
    otherargs.clearLandBT_hist = accumHist(otherargs.clearLandBT_hist, scaledBT[clearLand])
    scaledB4 = (ref[TM4] * B4_SCALE).astype(numpy.uint8)
    otherargs.clearLandB4_hist = accumHist(otherargs.clearLandB4_hist, scaledB4[clearLand])


def accumHist(counts, vals):
    """
    Accumulate the given values into the given (partial) counts
    """
    (valsHist, edges) = numpy.histogram(vals, bins=BT_HISTSIZE, range=(0, BT_HISTSIZE))
    counts += valsHist
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

# For scaling probability values so I can store them in 8 bits
PROB_SCALE = 100.0

def doPotentialCloudSecondPass(toareffile, radianceBands, thermalfile, 
        thermalInfo, pass1file, Twater, Tlow, Thigh, tempdir):
    """
    Second pass for potential cloud layer
    """
    infiles = applier.FilenameAssociations()
    outfiles = applier.FilenameAssociations()
    otherargs = applier.OtherInputs()
    controls = applier.ApplierControls()
    
    infiles.pass1 = pass1file
    infiles.toaref = toareffile
    infiles.thermal = thermalfile
    (fd, outfiles.pass2) = tempfile.mkstemp(prefix='pass2', dir=tempdir, 
                                suffix=DEFAULTEXTENSION)
    os.close(fd)
    otherargs.radianceBands = radianceBands
    otherargs.thermalInfo = thermalInfo
    
    otherargs.Twater = Twater
    otherargs.Tlow = Tlow
    otherargs.Thigh = Thigh
    otherargs.lCloudProb_hist = numpy.zeros(BT_HISTSIZE, dtype=numpy.uint32)
    
    controls.setWindowXsize(RIOS_WINDOW_SIZE)
    controls.setWindowYsize(RIOS_WINDOW_SIZE)
    controls.setReferenceImage(toareffile)
    controls.setCalcStats(False)
    controls.setOutputDriverName(DEFAULTDRIVERNAME)
    controls.setCreationOptions(DEFAULTCREATIONOPTIONS)
    
    applier.apply(potentialCloudSecondPass, infiles, outfiles, otherargs, controls=controls)
    
    # Equation 17
    landThreshold = scoreatpcnt(otherargs.lCloudProb_hist, 82.5)
    if landThreshold is not None:
        landThreshold = landThreshold / PROB_SCALE + 0.2
    else:
        landThreshold = 0.2
    return (outfiles.pass2, landThreshold)


def potentialCloudSecondPass(info, inputs, outputs, otherargs):
    """
    Called from RIOS
    
    Second pass of potential cloud layer
    """
    ref = inputs.toaref.astype(numpy.float) / 1000.0
    # Clamp off any reflectance <= 0
    ref[ref<=0] = 0.00001
    # Brightness temperature in degrees C
    bt = otherargs.thermalInfo.scaleThermalDNtoC(inputs.thermal)
    Twater = otherargs.Twater
    (Tlow, Thigh) = (otherargs.Tlow, otherargs.Thigh)
    # Values from first pass
    clearLand = inputs.pass1[2].astype(numpy.bool)
    variabilityProbPcnt = inputs.pass1[3]
    variability_prob = variabilityProbPcnt / PROB_SCALE

    # Equation 9
    if Twater is not None:
        wTemperature_prob = (Twater - bt) / 4.0
    else:
        # There is no water, so who cares. 
        wTemperature_prob = 1
    
    TM5 = otherinputs.radianceInfo[VIS_BAND_155um]
    # Equation 10
    brightness_prob = numpy.minimum(ref[TM5], 0.11) / 0.11
    
    # Equation 11
    wCloud_prob = wTemperature_prob * brightness_prob
    
    # Equation 14
    if Thigh is not None and Tlow is not None:
        lTemperature_prob = (Thigh + 4 - bt) / (Thigh + 4 - (Tlow - 4))
    else:
        # there was no land available for temperature thresholds, so it is probably cloud. 
        lTemperature_prob = 1
    
    # Equation 16
    lCloud_prob = lTemperature_prob * variability_prob
    
    outstack = numpy.array([
        (wCloud_prob * PROB_SCALE).clip(BYTE_MIN, BYTE_MAX), 
        (lCloud_prob * PROB_SCALE).clip(BYTE_MIN, BYTE_MAX)], dtype=numpy.uint8)
    outputs.pass2 = outstack
    
    # Accumulate histogram of lCloud_prob
    scaledProb = (lCloud_prob * PROB_SCALE).clip(BYTE_MIN, BYTE_MAX).astype(numpy.uint8)
    otherargs.lCloudProb_hist = accumHist(otherargs.lCloudProb_hist, scaledProb[clearLand])


def doCloudLayerFinalPass(thermalfile, thermalInfo, pass1file, pass2file, 
                    landThreshold, Tlow, tempdir):
    """
    Final pass
    """
    infiles = applier.FilenameAssociations()
    outfiles = applier.FilenameAssociations()
    otherargs = applier.OtherInputs()
    controls = applier.ApplierControls()
    
    infiles.pass1 = pass1file
    infiles.pass2 = pass2file
    infiles.thermal = thermalfile
    otherargs.landThreshold = landThreshold
    otherargs.Tlow = Tlow
    otherargs.thermalInfo = thermalInfo

    (fd, outfiles.cloudmask) = tempfile.mkstemp(prefix='interimcloud', dir=tempdir, 
                                suffix=DEFAULTEXTENSION)
    os.close(fd)
    # Need overlap so we can do Fmask's 3x3 fill-in
    overlap = 1
        
    controls.setOverlap(overlap)
    controls.setWindowXsize(RIOS_WINDOW_SIZE)
    controls.setWindowYsize(RIOS_WINDOW_SIZE)
    controls.setReferenceImage(pass1file)
    controls.setCalcStats(False)
    controls.setOutputDriverName(DEFAULTDRIVERNAME)
    controls.setCreationOptions(DEFAULTCREATIONOPTIONS)
    
    applier.apply(cloudFinalPass, infiles, outfiles, otherargs, controls=controls)
    
    return outfiles.cloudmask


def cloudFinalPass(info, inputs, outputs, otherargs):
    """
    Called from RIOS
    
    Final pass of cloud mask layer
    """
    nullmask = inputs.pass1[4].astype(numpy.bool)
    pcp = inputs.pass1[0]
    waterTest = inputs.pass1[1]
    notWater = numpy.logical_not(waterTest)
    notWater[nullmask] = False
    wCloud_prob = inputs.pass2[0] / PROB_SCALE
    lCloud_prob = inputs.pass2[1] / PROB_SCALE
    # Brightness temperature in degrees C
    bt = otherargs.thermalInfo.scaleThermalDNtoC(inputs.thermal)
    landThreshold = otherargs.landThreshold
    Tlow = otherargs.Tlow
    
    cloudmask1 = pcp & waterTest & (wCloud_prob>0.5)
    cloudmask2 = pcp & notWater & (lCloud_prob>landThreshold)
    cloudmask3 = (lCloud_prob > 0.99) & notWater
    if Tlow is not None:
        cloudmask4 = (bt < (Tlow-35))
    else:
        # Not enough land for final test
        cloudmask4 = True
        
    # Equation 18
    cloudmask = cloudmask1 | cloudmask2 | cloudmask3 | cloudmask4
    cloudmask[nullmask] = 0

    # Apply the prescribed 3x3 buffer. According to Zhu&Woodcock (page 87, end of section 3.1.2) 
    # they set a pixel to cloud if 5 or more of its 3x3 neighbours is cloud. 
    # This little incantation will do exactly the same. 
    bufferedCloudmask = (uniform_filter(cloudmask*2.0, size=3) >= 1.0)

    bufferedCloudmask[nullmask] = 0
    
    outputs.cloudmask = numpy.array([bufferedCloudmask])

def doPotentialShadows(toareffile, radianceBands, b4_17, tempdir):
    """
    Make potential shadow layer, as per section 3.1.3 of Zhu&Woodcock. 
    """
    (fd, potentialShadowsFile) = tempfile.mkstemp(prefix='shadows', dir=tempdir, 
                                suffix=DEFAULTEXTENSION)
    os.close(fd)

    # convert from numpy (0 based) to GDAL (1 based) indexing
    TM4_lyr = otherinputs.radianceInfo[VIS_BAND_076um] + 1
    
    # Read in whole of band 4
    ds = gdal.Open(toareffile)
    band = ds.GetRasterBand(TM4_lyr)
    nullval = band.GetNoDataValue()
    scaledb4 = band.ReadAsArray()
    b4_17_dn = b4_17 * 1000
    
    scaledb4_filled = fillminima.fillMinima(scaledb4, nullval, b4_17_dn)

    b4 = scaledb4.astype(numpy.float) / 1000.0
    b4_filled = scaledb4_filled.astype(numpy.float) / 1000.0
    del scaledb4, scaledb4_filled
    
    # Equation 19
    potentialShadows = ((b4_filled - b4) > 0.02)
    
    driver = gdal.GetDriverByName(DEFAULTDRIVERNAME)
    outds = driver.Create(potentialShadowsFile, ds.RasterXSize, ds.RasterYSize, 
                    1, gdal.GDT_Byte, DEFAULTCREATIONOPTIONS)
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
    
    (clumps, numClumps) = label(cloudmask, structure=numpy.ones((3,3)))
    
    return (clumps, numClumps)


CLOUD_HEIGHT_SCALE = 10
def make3Dclouds(clumps, numClumps, thermalFile, thermalInfo, toaRefFile, 
                        radianceBands):
    """
    Create 3-dimensional cloud objects from the cloud mask, and the thermal 
    information. Assumes a constant lapse rate to convert temperature into height.
    Resulting cloud heights are relative to cloud base. 
    
    Returns an image of relative cloud height (relative to cloud base for 
    each cloud object), and a dictionary of cloud base temperature, for each 
    cloud object, and mdl.ValueIndexes object for use in extracting the location of 
    every pixel for a given cloud object. 
    
    """
    # Find out the pixel grid of the toareffile, so we can use that for RIOS.
    # this is necessary because the thermal might be on a different grid,
    # and we can use RIOS to resample that. 
    referencePixgrid = pixelgrid.pixelGridFromFile(toareffile)
    
    infiles = applier.FilenameAssociations()
    outfiles = applier.FilenameAssociations()
    otherargs = applier.OtherInputs()
    controls = applier.ApplierControls()
    
    infiles.thermal = thermalfile
    otherargs.clumps = clumps
    # TODO:
    otherargs.cloudClumpNdx = mdl.ValueIndexes(clumps, nullVals=[0])
    otherargs.numClumps = numClumps
    otherargs.thermalInfo = thermalInfo
    
    # Run RIOS on whole image as one block
    (nRows, nCols) = referencePixgrid.getDimensions()
    controls.setWindowXsize(nCols)
    controls.setWindowYsize(nRows)
    controls.setReferencePixgrid(referencePixgrid)
    controls.setCalcStats(False)
    controls.setOutputDriverName(DEFAULTDRIVERNAME)
    controls.setCreationOptions(DEFAULTCREATIONOPTIONS)
    
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
    bt = otherargs.thermalInfo.scaleThermalDNtoC(inputs.thermal)
    
    cloudShape = numpy.zeros(bt.shape, dtype=numpy.uint8)
    cloudBaseTemp = {}
    
    cloudIDlist = otherargs.cloudClumpNdx.values
    for cloudID in cloudIDlist:
        cloudNdx = otherargs.cloudClumpNdx.getIndexes(cloudID)
        btCloud = bt[cloudNdx]
        
        numPixInCloud = len(cloudNdx[0])
        
        # Equation 22, in several pieces
        R = numpy.sqrt(numPixInCloud/(2*numpy.pi))
        if R >= 8:
            percentile = 100.0 * (R-8.0)**2 / (R**2)
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
    
    otherargs.cloudShape = cloudShape
    otherargs.cloudBaseTemp = cloudBaseTemp

METRES_PER_KM = 1000.0
BYTES_PER_VOXEL = 4
SOLIDCLOUD_MAXMEM = float(1024*1024*1024)

def makeCloudShadowShapes(toareffile, radianceBands, cloudShape, cloudClumpNdx, anglesInfo):
    """
    Project the 3d cloud shapes onto horizontal surface, along the sun vector, to
    make the 2d shape of the shadow. 
    """
    # Read in the two solar angles. Assumes that the angles file is on the same 
    # pixel grid as the cloud, which should always be the case. 
    ds = gdal.Open(toareffile)
    geotrans = ds.GetGeoTransform()
    (xRes, yRes) = (float(geotrans[1]), float(geotrans[5]))
    (nrows, ncols) = (ds.RasterYSize, ds.RasterXSize)
    del ds

    # tell anglesInfo it may need to read data into memory
    anglesInfo.prepareForQuerying()
    
    shadowShapesDict = {}
    
    cloudIDlist = cloudClumpNdx.values
    for cloudID in cloudIDlist:
        cloudNdx = cloudClumpNdx.getIndexes(cloudID)
        numPix = len(cloudNdx[0])
        
        sunAz = anglesInfo.getSolarAzimuthAngle(cloudNdx)
        sunZen = anglesInfo.getSolarAzimuthAngle(cloudNdx)
        satAz = anglesInfo.getViewAzimuthAngle(cloudNdx)
        satZen = anglesInfo.getViewZenithAngle(cloudNdx)
        
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
        rows = (yDash / yRes).astype(numpy.uint32).clip(0, nrows-1)
        cols = (xDash / xRes).astype(numpy.uint32).clip(0, ncols-1)
        
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
        #blankImg[shadowNdx] = True
        #shadowNdx = numpy.where(blankImg)
        #blankImg[shadowNdx] = False
        
        # Stash these shapes in a dictionary, along with the corresponding sun and satellite angles
        shadowShapesDict[cloudID] = (shadowNdx, satAz, satZen, sunAz, sunZen)
    
    # no more querying needed
    anglesInfo.releaseMemory()
    
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


def matchShadows(interimCloudmask, potentialShadowsFile, shadowShapesDict, cloudBaseTemp, 
        Tlow, Thigh, cmdargs, pass1file, shadowbuffersize):
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
    potentialShadow = band.ReadAsArray(xoff, yoff, ncols, nrows).astype(numpy.bool)
    del ds
    ds = gdal.Open(interimCloudmask)
    band = ds.GetRasterBand(1)
    (xoff, yoff) = topLeftDict[interimCloudmask]
    cloudmask = band.ReadAsArray(xoff, yoff, ncols, nrows).astype(numpy.bool)
    geotrans = ds.GetGeoTransform()
    (xRes, yRes) = (geotrans[1], geotrans[5])
    (xsize, ysize) = (ds.RasterXSize, ds.RasterYSize)
    proj = ds.GetProjection()
    del ds
    ds = gdal.Open(pass1file)
    band = ds.GetRasterBand(5)
    (xoff, yoff) = topLeftDict[pass1file]
    nullmask = band.ReadAsArray(xoff, yoff, ncols, nrows).astype(numpy.bool)
    del ds

    (fd, interimShadowmask) = tempfile.mkstemp(prefix='matchedshadows', dir=tempdir, 
                                suffix=DEFAULTEXTENSION)
    os.close(fd)
    
    shadowmask = numpy.zeros(potentialShadow.shape, dtype=numpy.bool)
    
    unmatchedCount = 0
    cloudIDlist = shadowShapesDict.keys()
    for cloudID in cloudIDlist:
        shadowEntry = shadowShapesDict[cloudID]
        Tcloudbase = cloudBaseTemp[cloudID]

        matchedShadowNdx = matchOneShadow(cloudmask, shadowEntry, potentialShadow, Tcloudbase, 
            Tlow, Thigh, xRes, yRes, cloudID, nullmask)
        
        if matchedShadowNdx is not None:
            shadowmask[matchedShadowNdx] = True
        else:
            unmatchedCount += 1

    if cmdargs.verbose:
        print("No shadow found for %s of %s clouds " % (unmatchedCount, len(cloudIDlist)))

    del potentialShadow, cloudmask, nullmask
    
    # Now apply a 3-pixel buffer, as per section 3.2 (2nd-last paragraph)
    # I have the buffer size settable from the commandline, with our default
    # being larger than the original. 
    kernel = mdl.makeBufferKernel(shadowbuffersize)
    shadowmaskBuffered = maximum_filter(shadowmask, footprint=kernel)

    driver = gdal.GetDriverByName(DEFAULTDRIVERNAME)
    ds = driver.Create(interimShadowmask, xsize, ysize, 1, gdal.GDT_Byte,
                DEFAULTCREATIONOPTIONS)
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
    Hcloudbase_min = max(0.2, (Tlow - 4 - Tcloudbase)/9.8) * METRES_PER_KM
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
    (nrows, ncols) = ((rowN-row0+1), (colN-col0+1))
    shadowTemplate = numpy.zeros((nrows, ncols), dtype=numpy.bool)
    shadowTemplate[shapeNdx[0]-row0, shapeNdx[1]-col0] = True
    
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
        if r >= 0 and r+nrows <= imgNrows and c >= 0 and c+ncols <= imgNcols:
            cloud = cloudmask[r:r+nrows, c:c+ncols]
            potShadow = potentialShadow[r:r+nrows, c:c+ncols]
            null = nullmask[r:r+nrows, c:c+ncols]
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


def finalizeAll(interimCloudmask, interimShadowmask, pass1file, outputfile,
            cloudbuffersize):
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
    outfiles.out = outputfile
    controls.setOverlap(cloudbuffersize)
    controls.setThematic(True)
    controls.setStatsIgnore(0)
    controls.setWindowXsize(RIOS_WINDOW_SIZE)
    controls.setWindowYsize(RIOS_WINDOW_SIZE)
    controls.setOutputDriverName(DEFAULTDRIVERNAME)
    controls.setCreationOptions(DEFAULTCREATIONOPTIONS)
    
    if cloudbuffersize > 0:
        otherargs.bufferkernel = mdl.makeBufferKernel(cloudbuffersize)

    # determine if we are SLC-off
    otherargs.isSLCoff = False
    if lcrfs.lcrfs('satellite', outputfile) == 'landsat7':
        date = lcrfs.lcrfs('date', outputfile)
        year = int(date[:2])
        if year < 50:
            year = 2000 + year
        else:
            year = 1900 + year
        month = int(date[2:4])
        day = int(date[4:])
        if year > 2003:
            otherargs.isSLCoff = True
        elif year == 2003:
            if month > 5:
                otherargs.isSLCoff = True
            elif month == 5:
                otherargs.isSLCoff = day >= 31
        print('SLC off', otherargs.isSLCoff)

    applier.apply(maskAndBuffer, infiles, outfiles, otherargs, controls=controls)
    
    rat.setColorTable(outfiles.out, numpy.array([[2, 255, 0, 255, 255],
                                                 [3, 255, 255, 0, 255],
                                                 [4, 85, 255, 255, 255]]))
    rat.writeColumn(outfiles.out, "Classification", [b"Null", b"Valid", b"Cloud", 
                                                    b"Cloud Shadow", b"Snow"])

def maskAndBuffer(info, inputs, outputs, otherargs):
    """
    Called from RIOS
    
    Apply cloud and shadow masks to snow layer, and buffer cloud layer
    
    The main aims of all this re-masking are:
        1) A pixel should be either cloud, shadow, snow or not, but never 
           more than one
        2) Areas which are null in the input imagery should be null in the 
           mask, even after buffering, etc. 
    
    As a concession to people who might be manually editing the resulting outputs,
    the term "null in the input imagery" is here interpreted to mean "null in the
    reflective bands", which means that any pixels which are non-null in reflective
    bands are non-null in the output masks, even if they didn't have thermal. This 
    leaves a few pixels non-null in the masks for which we didn't really have input
    data, but it is a very small number of pixels, and would (I estimate) be less
    than the error rate anyway, so I don't believe it has a noticeable effect on
    overall accuracy. 
    
    """
    snow = inputs.pass1[5].astype(numpy.bool)
    nullmask = inputs.pass1[4].astype(numpy.bool)
    refNullmask = inputs.pass1[6].astype(numpy.bool)
    thermNullmask = inputs.pass1[7].astype(numpy.bool)
    resetNullmask = (refNullmask | thermNullmask)
    # In the case of SLC-off images, we just use the reflective null mask, 
    # because otherwise we end up with a dusting of pixels with 
    # nulls in the cloud/shadow/snow masks along the edges of SLC-off gaps. 
    # This would be more correct, but much more painful for someone manually 
    # editing. These are generally the only pixels affected. 
    # With the way we trim the edges of the reflective image (to da1 stage),
    # the edges are already smaller in the reflective than the thermal so no 
    # problem. The same is not true for Landsat-8, in particular.

    # (James Shepherd) I have found this decision detrimental to our mosaicking of Landsat 7 SLC-off
    # We prefer to only have valid data where both reflective and thermal is present
    #
    #    if otherargs.isSLCoff:
    #        resetNullmask = refNullmask

    cloud = inputs.cloud[0].astype(numpy.bool)
    shadow = inputs.shadow[0].astype(numpy.bool)
    
    # Buffer the cloud
    if hasattr(otherargs, 'bufferkernel'):
        cloud = maximum_filter(cloud, footprint=otherargs.bufferkernel)

    # Mask the shadow, against the buffered cloud, and the nullmask
    shadow[cloud] = 0
    
    # Mask the snow
    mask = cloud | shadow
    snow[mask] = 0
    
    # now convert these masks to 0 - null
    # 1 - not null and not mask
    # 2 - cloud
    # 3 - cloud shadow
    # 4 - snow
    outNullval = 0
    out = (cloud.astype(numpy.uint8) + 1)
    out[shadow] = 3
    out[snow] = 4        
    out[resetNullmask] = outNullval
    outputs.out = numpy.array([out])
    
