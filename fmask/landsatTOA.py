#!/usr/bin/env python
"""
Module that handles convertion of scaled radiance (DN) 
values from USGS to Top of Atmosphere (TOA) reflectance (\*1000).
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

import numpy
from osgeo import gdal
from rios import applier, cuiprogress, fileinfo
from . import config

# Derived by Pete Bunting from 6S
LANDSAT8_ESUN = [1876.61, 1970.03, 1848.9, 1571.3, 967.66, 245.73, 82.03, 361.72]
# From Chander, G., Markham, B.L., Helder, D.L. (2008)
# Summary of current radiometric calibration coefficients for Landsat MSS, TM, ETM+, and EO-1 ALI sensors
# http://landsathandbook.gsfc.nasa.gov/pdfs/Landsat_Calibration_Summary_RSE.pdf
LANDSAT4_ESUN = [1983.0, 1795.0, 1539.0, 1028.0, 219.8, 83.49]
LANDSAT5_ESUN = [1983.0, 1796.0, 1536.0, 1031.0, 220.0, 83.44]
LANDSAT7_ESUN = [1997.0, 1812.0, 1533.0, 1039.0, 230.8, 84.90]
# Derived by RSC, using Landsat-9 RSR & Kurucz solar spectrum
LANDSAT9_ESUN = [1881.97, 1972.61, 1851.75, 1572.03, 969.37, 246.18, 82.11, 360.48]

ESUN_LOOKUP = {'LANDSAT_4': LANDSAT4_ESUN,
            'LANDSAT_5': LANDSAT5_ESUN,
            'LANDSAT_7': LANDSAT7_ESUN,
            'LANDSAT_8': LANDSAT8_ESUN,
            'LANDSAT_9': LANDSAT9_ESUN}

RADIANCE_MULT = 'RADIANCE_MULT_BAND_%d'
RADIANCE_ADD = 'RADIANCE_ADD_BAND_%d'

LMAX_KEY = 'LMAX_BAND%d'
LMIN_KEY = 'LMIN_BAND%d'
QCALMAX_KEY = 'QCALMAX_BAND%d'
QCALMIN_KEY = 'QCALMIN_BAND%d'


# band numbers in mtl file for gain and offset for reflective
BAND_NUM_DICT = {'LANDSAT_4': (1, 2, 3, 4, 5, 7), 
    'LANDSAT_5': (1, 2, 3, 4, 5, 7),
    'LANDSAT_7': (1, 2, 3, 4, 5, 7),
    'LANDSAT_8': (1, 2, 3, 4, 5, 6, 7, 9),
    'LANDSAT_9': (1, 2, 3, 4, 5, 6, 7, 9)}


def readGainsOffsets(mtlInfo):
    """
    Read the gains and offsets out of the .MTL file
    """
    spaceCraft = mtlInfo['SPACECRAFT_ID']
    nbands = len(BAND_NUM_DICT[spaceCraft])

    gains = numpy.zeros(nbands)
    offsets = numpy.zeros(nbands)
    
    if (RADIANCE_MULT % 1) in mtlInfo:
        # Newer format for MTL file
        for idx, band in enumerate(BAND_NUM_DICT[spaceCraft]):
            s = RADIANCE_MULT % band
            gain = float(mtlInfo[s])
            gains[idx] = gain

            s = RADIANCE_ADD % band
            offset = float(mtlInfo[s])
            offsets[idx] = offset
    else:
        # Old format, calculate gain and offset from min/max values
        for (idx, band) in enumerate(BAND_NUM_DICT[spaceCraft]):
            lMax = float(mtlInfo[LMAX_KEY % band])
            lMin = float(mtlInfo[LMIN_KEY % band])
            qcalMax = float(mtlInfo[QCALMAX_KEY % band])
            qcalMin = float(mtlInfo[QCALMIN_KEY % band])
            
            gain = (lMax - lMin) / (qcalMax - qcalMin)
            offset = lMin - qcalMin * gain
            
            gains[idx] = gain
            offsets[idx] = offset
                                        
    return gains, offsets


def earthSunDistance(date):
    """
    Given a date in YYYYMMDD will compute the earth sun distance in astronomical units
    """
    import datetime
    year = int(date[:4])
    month = int(date[4:6])
    day = int(date[6:])
    d1 = datetime.datetime(year, month, day)
    d2 = datetime.datetime(year, 1, 1)  # first day of year
    deltaT = d1 - d2
    jday = deltaT.days + 1  # julian day of year.
    ds = (1.0 - 0.01673 * numpy.cos(0.9856 * (jday - 4) * numpy.pi / 180.0))
    return ds


def riosTOA(info, inputs, outputs, otherinputs):
    """
    Called from RIOS
    """
    nbands = inputs.infile.shape[0]

    infile = inputs.infile.astype(numpy.float32)
    inIgnore = otherinputs.inNull
    if inIgnore is None:
        inIgnore = 0

    cosSunZen = numpy.cos(inputs.angles[3] * otherinputs.anglesToRadians)
    
    nullMask = (inputs.infile == inIgnore).any(axis=0)
    
    toaRefList = []
    for band in range(nbands):
        rtoa = infile[band] * otherinputs.gains[band] + otherinputs.offsets[band]

        p = numpy.pi * rtoa * otherinputs.earthSunDistanceSq / (otherinputs.esun[band] * cosSunZen)
        # clip to a sensible range
        numpy.clip(p, 0.0, 2.0, out=p)

        toaRefList.append(p)

    out = numpy.array(toaRefList) * 10000.0
    # convert to int16 
    outputs.outfile = out.astype(numpy.int16)
    # Mask out where input is null
    for i in range(len(outputs.outfile)):
        outputs.outfile[i][nullMask] = otherinputs.outNull


def makeTOAReflectance(infile, mtlFile, anglesfile, outfile):
    """
    Main routine - does the calculation

    The eqn for TOA reflectance, p, is
    p = pi * L * d^2 / E * cos(theta)
    
    d = earthSunDistance(date)
    L = image pixel (radiance)
    E = exoatmospheric irradiance for the band, and 
    theta = solar zenith angle.
    
    Assumes infile is radiance values in DN from USGS.
    mtlFile is the .mtl file.
    outfile will be created in the default format that RIOS
    is configured to use and will be top of atmosphere 
    reflectance values *10000. Also assumes that the 
    angles image file is scaled as radians*100, and has layers for
    satAzimuth, satZenith, sunAzimuth, sunZenith, in that order. 
    
    """
    mtlInfo = config.readMTLFile(mtlFile)
    spaceCraft = mtlInfo['SPACECRAFT_ID']
    date = mtlInfo['DATE_ACQUIRED']
    date = date.replace('-', '')
    
    inputs = applier.FilenameAssociations()
    inputs.infile = infile
    inputs.angles = anglesfile

    outputs = applier.FilenameAssociations()
    outputs.outfile = outfile

    otherinputs = applier.OtherInputs()
    otherinputs.earthSunDistance = earthSunDistance(date)
    otherinputs.earthSunDistanceSq = otherinputs.earthSunDistance * otherinputs.earthSunDistance
    otherinputs.esun = ESUN_LOOKUP[spaceCraft]
    gains, offsets = readGainsOffsets(mtlInfo)
    otherinputs.gains = gains
    otherinputs.offsets = offsets
    otherinputs.anglesToRadians = 0.01
    otherinputs.outNull = 32767
    imginfo = fileinfo.ImageInfo(infile)
    otherinputs.inNull = imginfo.nodataval[0]

    controls = applier.ApplierControls()
    controls.progress = cuiprogress.GDALProgressBar()
    controls.setStatsIgnore(otherinputs.outNull)
    controls.setCalcStats(False)
    
    applier.apply(riosTOA, inputs, outputs, otherinputs, controls=controls)
    
    # Explicitly set the null value in the output
    ds = gdal.Open(outfile, gdal.GA_Update)
    for i in range(ds.RasterCount):
        ds.GetRasterBand(i + 1).SetNoDataValue(otherinputs.outNull)
