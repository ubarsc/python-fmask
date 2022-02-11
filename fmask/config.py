"""
Configuration classes that define the inputs and parameters
for the fmask function.
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

import abc
import numpy
import scipy.constants

from osgeo import gdal
from rios import applier
from . import fmaskerrors

gdal.UseExceptions()

FMASK_LANDSAT47 = 0
"Landsat 4 to 7"
FMASK_LANDSAT8 = 1
"Landsat 8"
FMASK_SENTINEL2 = 2
"Sentinel 2"
FMASK_LANDSATOLI = 3
"Landsat OLI"

"""
Some constants for the various reflective bands used in fmask.
"""
#: ~475nm
BAND_BLUE = 0
#: ~560nm
BAND_GREEN = 1
#: ~660nm
BAND_RED = 2
#: ~780nm
BAND_NIR = 3
#: ~1360nm
BAND_CIRRUS = 4    # Sentinel2 + Landsat8 only
#: ~1610nm
BAND_SWIR1 = 5
#: ~2200nm
BAND_SWIR2 = 6

# These bands are used only for the Sentinel-2 Cloud Displacement Index code. They
# are NIR bands, with slightly different look angles, and used as per Frantz et al 2017.
#: ~865nm
BAND_S2CDI_NIR8A = 7
#: ~783nm
BAND_S2CDI_NIR7 = 8

BAND_WATERVAPOUR = 9


class FmaskConfig(object):
    """
    Class that contains the configuration parameters of the fmask
    run.
    """
    # some parameters for fmask operation
    keepIntermediates = False
    cloudBufferSize = 5
    shadowBufferSize = 10
    verbose = False
    strictFmask = False
    tempDir = '.'
    TOARefScaling = 10000.0
    TOARefDNoffsetDict = None
    # Minimum number of pixels in a single cloud (before buffering). A non-zero value
    # would allow filtering of very small clouds. 
    minCloudSize_pixels = 0
        
    # constants from the paper that could probably be tweaked
    # equation numbers are from the original paper.
    Eqn1Swir2Thresh = 0.03
    Eqn1ThermThresh = 27
    Eqn2WhitenessThresh = 0.7
    cirrusBandTestThresh = 0.01
    Eqn7Swir2Thresh = 0.03
    Eqn20ThermThresh = 3.8
    Eqn20NirSnowThresh = 0.11
    Eqn20GreenSnowThresh = 0.1
    cirrusProbRatio = 0.04
    Eqn19NIRFillThresh = 0.02
    
    # Constant term at the end of Equation 17. Zhu's MATLAB code now has this as a configurable
    # value, which they recommend as 22.5% (i.e. 0.225)
    Eqn17CloudProbThresh = 0.2
    
    # GDAL driver for final output file
    gdalDriverName = applier.DEFAULTDRIVERNAME
    
    # Do we do the Sentinel-2 Cloud Displacement Test ?
    sen2displacementTest = False
    sen2cdiWindow = 7

    def __init__(self, sensor):
        """
        Pass in the sensor (one of: FMASK_LANDSAT47, FMASK_LANDSAT8 or
        FMASK_SENTINEL2) and default of parameters will be set. These 
        can be overridden using the functions on this object.
        """
        self.sensor = sensor
        
        # Some standard file configurations for different sensors.
        # Assumed that panchromatic + thermal bands stored in separate files.
        # zero based indexing
        if sensor == FMASK_LANDSAT47:
            self.bands = {BAND_BLUE: 0, BAND_GREEN: 1, BAND_RED: 2, BAND_NIR: 3,
                BAND_SWIR1: 4, BAND_SWIR2: 5}
        elif sensor in (FMASK_LANDSAT8, FMASK_LANDSATOLI):
            self.bands = {BAND_BLUE: 1, BAND_GREEN: 2, BAND_RED: 3, BAND_NIR: 4,
                BAND_SWIR1: 5, BAND_SWIR2: 6, BAND_CIRRUS: 7}
        elif sensor == FMASK_SENTINEL2:
            # Assumes the input stack has ALL bands, in their numeric order (with 8A after 8)
            self.bands = {BAND_BLUE: 1, BAND_GREEN: 2, BAND_RED: 3, BAND_NIR: 7,
                    BAND_SWIR1: 11, BAND_SWIR2: 12, BAND_WATERVAPOUR: 9, BAND_CIRRUS: 10,
                    BAND_S2CDI_NIR7: 6, BAND_S2CDI_NIR8A: 8}
        else:
            msg = 'unrecognised sensor'
            raise fmaskerrors.FmaskParameterError(msg)

        # we can't do anything with the thermal yet since
        # we need a .mtl file or equivalent to get the gains etc
        self.thermalInfo = None
        
        # same with angles
        self.anglesInfo = None
        
        # obtain the usual extension for the GDAL driver used by RIOS
        # so we can create temporary files with this extension.
        driver = gdal.GetDriverByName(applier.DEFAULTDRIVERNAME)
        if driver is None:
            msg = 'Cannot find GDAL driver %s used by RIOS' 
            msg = msg % applier.DEFAULTDRIVERNAME
            raise fmaskerrors.FmaskParameterError(msg)

        ext = driver.GetMetadataItem('DMD_EXTENSION')
        if ext is None:
            self.defaultExtension = '.tmp'
        else:
            self.defaultExtension = '.' + ext
        
    def setReflectiveBand(self, band, index):
        """
        Tell fmask which band is in which index in the reflectance
        data stack file. band should be one of the BAND_* constants.
        index is zero based (ie 0 is first band in the file).
        
        These are set to default values for each sensor which are
        normally correct, but this function can be used to update.
        
        """
        self.bands[band] = index
        
    def setThermalInfo(self, info):
        """
        Set an instance of ThermalFileInfo. By default this is
        None and fmask assumes there is no thermal data available.
        
        The :func:`fmask.config.readThermalInfoFromLandsatMTL`
        function can be used to obtain this from a Landsat .mtl file.
        
        """
        self.thermalInfo = info
        
    def setAnglesInfo(self, info):
        """
        Set an instance of AnglesInfo. By default this is 
        None and will need to be set before fmask will run.
        
        The :func:`fmask.config.readAnglesFromLandsatMTL` 
        function can be used to obtain this from a Landsat .mtl file.

        """
        self.anglesInfo = info
        
    def setTOARefScaling(self, scaling):
        """
        Set the scaling used in the Top of Atmosphere reflectance
        image. The calculation is done as

        ref = (dn + dnOffset) / scaling

        and so is used in conjunction with the offset values 
        (see setTOARefOffsets). 

        The dnOffset was added in 2021 to cope with ESA's absurd 
        decision to suddenly introduce an offset in their Sentinel-2 
        TOA reflectance imagery. For Landsat, there is no need for 
        it ever to be non-zero. 

        """
        self.TOARefScaling = scaling

    def setTOARefOffsetDict(self, offsetDict):
        """
        Set the reflectance offsets to the given list. 
        This should contain an offset value for each band used with
        the Fmask code. The keys are the named constants in the 
        config module, BAND_*. 

        The offset is added to the corresponding band pixel values
        before dividing by the scaling value. 

        This facility is made available largely for use with Sentinel-2,
        after ESA unilaterally starting using non-zero offsets in their
        Level-1C imagery (Nov 2021). However, it can be used 
        with Landsat if required. 

        """
        if len(offsetDict) != len(self.bands):
            msg = "Must supply offsets for all bands being used"
            raise fmaskerrors.FmaskParameterError(msg)

        self.TOARefDNoffsetDict = offsetDict

    def setKeepIntermediates(self, keepIntermediates):
        """
        Set to True to keep the intermediate files created in
        processed. This is False by default.
        
        """
        self.keepIntermediates = keepIntermediates
    
    def setCloudBufferSize(self, bufferSize):
        """
        Extra buffer of this many pixels on cloud layer. Defaults to 5.
        
        """
        self.cloudBufferSize = bufferSize

    def setShadowBufferSize(self, bufferSize):
        """
        Extra buffer of this many pixels on cloud layer. Defaults to 10.
        
        """
        self.shadowBufferSize = bufferSize
    
    def setMinCloudSize(self, minCloudSize):
        """
        Set the minimum cloud size retained. This minimum is applied before any
        buffering of clouds. Size is specified as an area, in pixels. 
        """
        self.minCloudSize_pixels = minCloudSize
        
    def setVerbose(self, verbose):
        """
        Print informative messages. Defaults to False.
        
        """
        self.verbose = verbose
        
    def setStrictFmask(self, strictFmask):
        """
        Set whatever options are necessary to run strictly as per Fmask paper 
        (Zhu & Woodcock). Setting this will override the settings of other
        parameters on this object.
    
        """
        self.strictFmask = strictFmask
        
    def setTempDir(self, tempDir):
        """
        Temporary directory to use. Defaults to '.' (the current directory).
        
        """
        self.tempDir = tempDir
        
    def setDefaultExtension(self, extension):
        """
        Sets the default extension used by temporary files created by
        fmask. Defaults to the extension of the driver that RIOS
        is configured to use.
        
        Note that this should include the '.' - ie '.img'.
        
        """
        self.defaultExtension = extension
        
    def setEqn1Swir2Thresh(self, thresh):
        """
        Change the threshold used by Equation 1 for the SWIR2 band.
        This defaults to 0.03
        
        """
        self.Eqn1Swir2Thresh = thresh
        
    def setEqn1ThermThresh(self, thresh):
        """
        Change the threshold used by Equation one for BT.
        This defaults to 27.
        
        """
        self.Eqn1ThermThresh = thresh
        
    def setEqn2WhitenessThresh(self, thresh):
        """
        Change the threshold used by Equation 2 to determine
        whiteness from visible bands. This defaults to 0.7.
        
        """
        self.Eqn2WhitenessThresh = thresh
        
    def setCirrusBandTestThresh(self, thresh):
        """
        Change the threshold used by Zhu et al 2015, section 2.2.1
        for the cirrus band test. Defaults to 0.01.
        
        """
        self.cirrusBandTestThresh = thresh
        
    def setEqn7Swir2Thresh(self, thresh):
        """
        Change the threshold used by Equation 7 (water test)
        for the Swir2 band. This defaults to 0.03.
        
        """
        self.Eqn7Swir2Thresh = thresh
        
    def setEqn17CloudProbThresh(self, thresh):
        """
        Change the threshold used by Equation 17. The threshold
        given here is the constant term added to the end of the equation
        for the land probability threshold. Original paper had this as 0.2,
        although Zhu et al's MATLAB code now defaults it to 0.225 (i.e. 22.5%)
        
        """
        self.Eqn17CloudProbThresh = thresh
        
    def setEqn20ThermThresh(self, thresh):
        """
        Change the threshold used by Equation 20 (snow)
        for BT. This defaults to 3.8.
        
        """
        self.Eqn20ThermThresh = thresh
        
    def setEqn20NirSnowThresh(self, thresh):
        """
        Change the threshold used by Equation 20 (snow)
        for NIR reflectance. This defaults to 0.11
        
        """
        self.Eqn20NirSnowThresh = thresh
        
    def setEqn20GreenSnowThresh(self, thresh):
        """
        Change the threshold used by Equation 20 (snow)
        for green reflectance. This defaults to 0.1
        
        """
        self.Eqn20GreenSnowThresh = thresh
        
    def setCirrusProbRatio(self, ratio):
        """
        Change the ratio used by Zhu et al 2015 Equation 1
        to determine the cirrus cloud probability. Defaults
        to 0.04.
        
        """
        self.cirrusProbRatio = ratio
        
    def setEqn19NIRFillThresh(self, thresh):
        """
        Change the threshold used by Equation 19 to determine
        potential cloud shadow from the difference between NIR
        and flood filled NIR. Defaults to 0.02.
        
        """
        self.Eqn19NIRFillThresh = thresh
    
    def setSen2displacementTest(self, useDisplacementTest):
        """
        Set whether or not to use the Frantz (2018) parallax displacement test
        to remove false clouds. Pass True if the test is desired, False otherwise. 
        
        """
        self.sen2displacementTest = useDisplacementTest
    
    def setGdalDriverName(self, driverName):
        """
        Change the GDAL driver used for writing the final output file. Default
        value is taken from the default for the RIOS package, as per $RIOS_DFLT_DRIVER. 
        """
        self.gdalDriverName = driverName


class FmaskFilenames(object):
    """
    Class that contains the filenames used in the fmask run.
    """
    toaRef = None
    thermal = None
    saturationMask = None
    outputMask = None
    
    def __init__(self, toaRefFile=None, thermalFile=None, outputMask=None,
                saturationMask=None):
        self.toaRef = toaRefFile
        self.thermal = thermalFile
        self.saturationMask = saturationMask
        self.outputMask = outputMask
    
    def setThermalFile(self, thermalFile):
        """
        Set the path of the input thermal file. To make use
        of this, the :func:`fmask.config.FmaskConfig.setThermalInfo`
        function must also be called so that fmask knows how
        to use the file.
        
        This file should be in any GDAL readable format.
        
        """
        self.thermal = thermalFile
    
    def setTOAReflectanceFile(self, toaRefFile):
        """
        Set the path of the input top of atmosphere (TOA) file. It pays
        to check that the default set of bands match what fmask expects in 
        the :class:`fmask.config.FmaskConfig` class and update if necessary.
        
        This should have numbers which are reflectance * 1000
        
        Use the :func:`fmask.landsatTOA.makeTOAReflectance` function to create
        this file from raw Landsat radiance (or the fmask_usgsLandsatTOA.py
        command line program supplied with fmask).
        
        It is assumed that any values that are nulls in the original radiance
        image are set to the ignore values in the toaRefFile.

        This file should be in any GDAL readable format.
        
        """
        self.toaRef = toaRefFile
    
    def setSaturationMask(self, mask):
        """
        Set the mask to use for ignoring saturated pixels. By default
        no mask is used and all pixels are assumed to be unsaturated.
        This will cause problems for the whiteness test if some pixels
        are in fact saturated, but not masked out.
        
        Use the :func:`fmask.saturation.makeSaturationMask` function to
        create this from input radiance data.
        
        This mask should be 1 for pixels that are saturated, 0 otherwise.
        
        Note that this is not in the original paper so cannot be considered
        'strict', but if provided is used no matter the strict setting in 
        :class:`fmask.config.FmaskConfig`.

        This file should be in any GDAL readable format.
        
        """
        self.saturationMask = mask
    
    def setOutputCloudMaskFile(self, cloudMask):
        """
        Set the output cloud mask path. 
        
        Note that this file will be written in the format
        that RIOS is currently configured to use. See the 
        `RIOS documentation <http://rioshome.org/rios_imagewriter.html#rios.imagewriter.setDefaultDriver>`_
        for more details. Note that the default is HFA (.img) and can
        be overridden using environment variables.
        
        """
        self.outputMask = cloudMask


class ThermalFileInfo(object):
    """
    Contains parameters for interpreting thermal file.
    See :func:`fmask.config.readThermalInfoFromLandsatMTL`.
    
    """
    thermalBand1040um = None
    thermalGain1040um = None
    thermalOffset1040um = None
    thermalK1_1040um = None
    thermalK2_1040um = None
    
    def __init__(self, thermalBand1040um, thermalGain1040um,
            thermalOffset1040um, thermalK1_1040um, thermalK2_1040um):
        self.thermalBand1040um = thermalBand1040um
        self.thermalGain1040um = thermalGain1040um
        self.thermalOffset1040um = thermalOffset1040um
        self.thermalK1_1040um = thermalK1_1040um
        self.thermalK2_1040um = thermalK2_1040um

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
        rad[rad <= 0] = 0.00001  # to stop errors below
        temp = self.thermalK2_1040um / numpy.log(self.thermalK1_1040um / rad + 1.0)
        bt = temp - KELVIN_ZERO_DEGC
        return bt


# Keys within a .mtl file for each band
LANDSAT_RADIANCE_MULT = 'RADIANCE_MULT_BAND_%s'
LANDSAT_RADIANCE_ADD = 'RADIANCE_ADD_BAND_%s'
LANDSAT_K1_CONST = 'K1_CONSTANT_BAND_%s'
LANDSAT_K2_CONST = 'K2_CONSTANT_BAND_%s'

# Oldest format of MTL file has only min/max values
LANDSAT_LMAX_KEY = 'LMAX_BAND%s'
LANDSAT_LMIN_KEY = 'LMIN_BAND%s'
LANDSAT_QCALMAX_KEY = 'QCALMAX_BAND%s'
LANDSAT_QCALMIN_KEY = 'QCALMIN_BAND%s'

# band numbers in mtl file for gain and offset for thermal
LANDSAT_TH_BAND_NUM_DICT = {'LANDSAT_4': '6', 
        'LANDSAT_5': '6',
        'LANDSAT_7': '6_VCID_1',
        'LANDSAT_8': '10',
        'LANDSAT_9': '10'}
                        

# for some reason L4, 5, and 7 don't
# have these numbers in the mtl file, but L8 does
# from http://www.yale.edu/ceo/Documentation/Landsat_DN_to_Kelvin.pdf
LANDSAT_K1_DICT = {'TM': 607.76, 'ETM': 666.09, 'ETM+': 666.09}
LANDSAT_K2_DICT = {'TM': 1260.56, 'ETM': 1282.71, 'ETM+': 1282.71}


def readThermalInfoFromLandsatMTL(mtlfile, thermalBand1040um=0):
    """
    Returns an instance of ThermalFileInfo given a path to the mtl
    file and the index of the thermal band.
    
    """
    mtlData = readMTLFile(mtlfile)
    gain = None
    offset = None
    k1 = None
    k2 = None
    if 'SPACECRAFT_ID' in mtlData:
        # we can now grab the gain and offset
        spaceCraft = mtlData['SPACECRAFT_ID']
        band = LANDSAT_TH_BAND_NUM_DICT[spaceCraft]
        
        s = LANDSAT_RADIANCE_MULT % band

        oldestMtlFormat = (s not in mtlData)
        
        if not oldestMtlFormat:
            gain = float(mtlData[s])
            s = LANDSAT_RADIANCE_ADD % band
            offset = float(mtlData[s])
        else:
            # Oldest format MTL file
            if spaceCraft == "LANDSAT_7":
                band = "61"
            lMax = float(mtlData[LANDSAT_LMAX_KEY % band])
            lMin = float(mtlData[LANDSAT_LMIN_KEY % band])
            qcalMax = float(mtlData[LANDSAT_QCALMAX_KEY % band])
            qcalMin = float(mtlData[LANDSAT_QCALMIN_KEY % band])
            gain = (lMax - lMin) / (qcalMax - qcalMin)
            offset = lMin - qcalMin * gain
        
    if 'SENSOR_ID' in mtlData:
        # look for k1 and k2
        sensor = mtlData['SENSOR_ID']
        s = LANDSAT_K1_CONST % band
        if s in mtlData:
            k1 = float(mtlData[s])
        else:
            # drop back to our own values if not in file
            k1 = LANDSAT_K1_DICT[sensor]
                                    
        s = LANDSAT_K2_CONST % band
        if s in mtlData:
            k2 = float(mtlData[s])
        else:
            # drop back to our own values if not in file
            k2 = LANDSAT_K2_DICT[sensor]
            
    if gain is not None and offset is not None and k1 is not None and k2 is not None:
        thermalInfo = ThermalFileInfo(thermalBand1040um, gain, 
                        offset, k1, k2)
    else:
        msg = 'Cannot find SPACECRAFT_ID/SENSOR_ID in MTL file'
        raise fmaskerrors.FmaskFileError(msg)
        
    return thermalInfo


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
    
    @abc.abstractmethod
    def setScaleToRadians(self, scale):
        """
        Set scaling factor to get radians from angles image values. 
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
        
        # This default value matches the file produced by fmask_usgsLandsatMakeAnglesImage.py
        self.scaleToRadians = 0.01
    
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
                                self.solarZenithBand)
        self.solarAzimuthData = self.readData(self.solarAzimuthFilename, 
                                self.solarAzimuthBand)
        self.viewZenithData = self.readData(self.viewZenithFilename, 
                                self.viewZenithBand)
        self.viewAzimuthData = self.readData(self.viewAzimuthFilename, 
                                self.viewAzimuthBand)
        
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
        return self.solarZenithData[indices].mean() * self.scaleToRadians

    def getSolarAzimuthAngle(self, indices):
        """
        Return the average solar azimuth angle for the given indices
        """
        return self.solarAzimuthData[indices].mean() * self.scaleToRadians
    
    def getViewZenithAngle(self, indices):
        """
        Return the average view zenith angle for the given indices
        """
        return self.viewZenithData[indices].mean() * self.scaleToRadians

    def getViewAzimuthAngle(self, indices):
        """
        Return the average view azimuth angle for the given indices
        """
        return self.viewAzimuthData[indices].mean() * self.scaleToRadians

    def setScaleToRadians(self, scale):
        """
        Set scaling factor to get radians from angles image values. 
        """
        self.scaleToRadians = scale


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


def readMTLFile(mtl):
    """
    Very simple .mtl file reader that just creates a dictionary
    of key and values and returns it
    """
    dict = {}
    for line in open(mtl):
        arr = line.split('=')
        if len(arr) == 2:
            (key, value) = arr
            dict[key.strip()] = value.replace('"', '').strip()

    # For the older format of the MTL file, a few fields had different names. So, we fake the
    # new names, so that the rest of the code can just use those. 
    if 'ACQUISITION_DATE' in dict:
        dict['DATE_ACQUIRED'] = dict['ACQUISITION_DATE']
    if 'SCENE_CENTER_SCAN_TIME' in dict:
        dict['SCENE_CENTER_TIME'] = dict['SCENE_CENTER_SCAN_TIME']
    
    # Oldest format has spacecraft ID string formatted differently, so reformat it. 
    spaceCraft = dict['SPACECRAFT_ID']
    if spaceCraft.startswith('Landsat') and '_' not in spaceCraft:
        satNum = spaceCraft[-1]
        dict['SPACECRAFT_ID'] = "LANDSAT_" + satNum

    return dict
                                                                    

def readAnglesFromLandsatMTL(mtlfile):
    """
    Given the path to a Landsat USGS .MTL file, read the angles
    out and return an instance of AngleConstantInfo.
    
    This is no longer supported, and this routine now raises an exception. 
    
    """
    msg = ("The simplified option of assuming constant angles across the whole image is "+
        "no longer supported. You must use per-pixel angles. ")
    raise fmaskerrors.FmaskNotSupportedError(msg)
