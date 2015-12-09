#!/usr/bin/env python

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
import numpy
from rios import applier, cuiprogress
from . import fmask
from . import config

# Derived by Pete Bunting from 6S
LANDSAT8_ESUN = [1876.61, 1970.03, 1848.9, 1571.3, 967.66, 245.73, 82.03]
# From Chander, G., Markham, B.L., Helder, D.L. (2008)
# Summary of current radiometric calibration coefficients for Landsat MSS, TM, ETM+, and EO-1 ALI sensors
# http://landsathandbook.gsfc.nasa.gov/pdfs/Landsat_Calibration_Summary_RSE.pdf
LANDSAT4_ESUN = [1983.0, 1795.0, 1539.0, 1028.0, 219.8, 83.49]
LANDSAT5_ESUN = [1983.0, 1796.0, 1536.0, 1031.0, 220.0, 83.44]
LANDSAT7_ESUN = [1997.0, 1812.0, 1533.0, 1039.0, 230.8, 84.90]

ESUN_LOOKUP = {'LANDSAT_4' : LANDSAT4_ESUN,
            'LANDSAT_5' : LANDSAT5_ESUN,
            'LANDSAT_7' : LANDSAT7_ESUN,
            'LANDSAT_8' : LANDSAT8_ESUN}

RADIANCE_MULT = 'RADIANCE_MULT_BAND_%d'
RADIANCE_ADD = 'RADIANCE_ADD_BAND_%d'

# band numbers in mtl file for gain and offset for reflective
BAND_NUM_DICT = {'LANDSAT_4' : (1, 2, 3, 4, 5, 7), 
    'LANDSAT_5' : (1, 2, 3, 4, 5, 7),
    'LANDSAT_7' : (1, 2, 3, 4, 5, 7),
    'LANDSAT_8' : (1, 2, 3, 4, 5, 6, 7, 9)}
                                
def readGainsOffsets(mtlInfo):
    """
    Read the gains and offsets out of the .MTL file
    """
    spaceCraft = mtlInfo['SPACECRAFT_ID']
    nbands = len(BAND_NUM_DICT[spaceCraft])

    gains = numpy.zeros(nbands)
    offsets = numpy.zeros(nbands)
    
    for idx, band in enumerate(BAND_NUM_DICT[spaceCraft]):
        s = RADIANCE_MULT % band
        gain = float(mtlInfo[s])
        gains[idx] = gain
                        
        s = RADIANCE_ADD % band
        offset = float(mtlInfo[s])
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
    d2 = datetime.datetime(year, 1, 1) # first day of year
    deltaT = d1-d2
    jday = deltaT.days + 1 # julian day of year.
    ds = (1.0 - 0.01673*numpy.cos(0.9856*(jday-4)*numpy.pi/180.0))
    return ds

def riosTOA(info, inputs, outputs, otherinputs):
    """
    Called from RIOS
    """
    nbands = inputs.infile.shape[0]

    infile = inputs.infile.astype(numpy.float)

    toaRefList = []
    for band in range(nbands):
        rtoa = infile[band] * otherinputs.gains[band] + otherinputs.offsets[band]

        p = numpy.pi * rtoa * otherinputs.earthSunDistanceSq / (otherinputs.esun[band] * otherinputs.csza)
        # clip to a sensible range
        numpy.clip(p, -0.2, 2.0, out=p)
        # mask out where zero in input
        p[infile[band] == 0] = 0

        toaRefList.append(p)

    out = numpy.array(toaRefList) * 1000.0
    # convert to int16 
    outputs.outfile = out.astype(numpy.int16)


def makeTOAReflectance(infile, mtlFile, outfile):
    """
    Main routine - does the calculation

    The eqn for TOA reflectance, p, is
    p = pi * L * d^2 / E * cos(theta)
    
    d = geodesy.earthSunDistance(date)
    L = image pixel (radiance)
    E = exoatmospheric irradiance for the band, and 
    theta = solar zenith angle.
    """
    mtlInfo = config.readMTLFile(mtlFile)
    spaceCraft = mtlInfo['SPACECRAFT_ID']
    date = mtlInfo['DATE_ACQUIRED']
    date = date.replace('-', '')
    sza = 90.0 - float(mtlInfo['SUN_ELEVATION'])
    
    inputs = applier.FilenameAssociations()
    inputs.infile = infile

    outputs = applier.FilenameAssociations()
    outputs.outfile = outfile

    otherinputs = applier.OtherInputs()
    otherinputs.earthSunDistance = earthSunDistance(date)
    otherinputs.earthSunDistanceSq = otherinputs.earthSunDistance * otherinputs.earthSunDistance
    otherinputs.esun = ESUN_LOOKUP[spaceCraft]
    otherinputs.csza = numpy.cos(numpy.radians(sza))
    gains, offsets = readGainsOffsets(mtlInfo)
    otherinputs.gains = gains
    otherinputs.offsets = offsets

    controls = applier.ApplierControls()
    controls.progress = cuiprogress.GDALProgressBar()
    
    applier.apply(riosTOA, inputs, outputs, otherinputs, controls=controls)

if __name__ == '__main__':

    cmds = CmdArgs()

    makeTOAReflectance(cmds.infile, cmds.mtl, cmds.output)

