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

"""
Generate an output image file with estimates of per-pixel angles for
sun and satellite azimuth and zenith. These are rough estimates, using the generic
characteristics of the Landsat 5 platform, and are not particularly accurate,
but good enough for the current purposes. 
In the future, the USGS have plans to distribute a set of parameters which can 
be used to directly generate these angles, derived from the actual orbit ephemeris
of the pass. This program is intended for use when these parameters are not 
available. Quite possibly, when they do make these available, I will modify this 
program to use those parameters if they are present in the MTL file, but fall back 
to the old estimates if not. 

The general approach for satellite angles is to estimate the nadir line by running it
down the middle of the image data area. The satellite azimuth is assumed to be
at right angles to this nadir line, which is only roughly correct. For the whisk-broom 
sensors on Landsat-5 and Landsat-7, this angles is not 90 degrees, but is affected by 
earth rotation and is latitude dependent, while for Landsat-8, the scan line is at 
right angles, due to the compensation for earth rotation, but the push-broom is 
made up of sub-modules which point in slightly different directions, giving 
slightly different satellite azimuths along the scan line. None of these effects
are included in the current estimates. The satellite zenith is estimated based on the
nadir point, the scan-line, and the assumed satellite altitude, and includes the
appropriate allowance for earth curvature. 

Because this works by searching the imagery for the non-null area, and assumes that 
this represents a full-swath image, it would not work for a subsection. 

The sun angles are approximated using the Michalsky algorithm (Michalsky, 1988), 
which is accurate to 0.011 degree, for the period 1950-2050. 
(Actually, maybe I should just translate the 6S code from
Fortran.......)

"""
from __future__ import print_function, division

import argparse

import numpy
from osgeo import osr

from rios import applier
from rios import fileinfo

from fmask import config


def getCmdargs():
    """
    Get commandline arguments
    """
    p = argparse.ArgumentParser()
    p.add_argument("-m", "--mtl", help="MTL text file of USGS metadata")
    p.add_argument("-t", "--templateimg", 
        help="Image filename to use as template for output angles image")
    p.add_argument("-o", "--outfile", help="Output image file")
    cmdargs = p.parse_args()
    return cmdargs


def mainRoutine():
    """
    Main routine
    """
    cmdargs = getCmdargs()
    
    mtlInfo = config.readMTLFile(cmdargs.mtl)
    
    imgInfo = fileinfo.ImageInfo(cmdargs.templateimg)
    corners = findImgCorners(cmdargs.templateimg, imgInfo)
    nadirLine = findNadirLine(corners)
    
    cornerSunAngles = sunAnglesForPoints(corners, mtlInfo)
    satAzimuth = satAzLeftRight(nadirLine)
    
    makeAnglesImage(cmdargs, nadirLine, corners, cornerSunAngles, satAzimuth, imgInfo)
    

def findImgCorners(img, imgInfo):
    """
    Find the corners of the data within the given template image
    Return a numpy array of (x, y) coordinates. The array has 2 columns, for X and Y.
    Each row is a corner, in the order
        top-left, top-right, bottom-left, bottom-right.
        
    Uses RIOS to pass through the image searching for non-null data,
    and find the extremes. Assumes we are working with a full-swathe Landsat
    image. 
    
    Each list element is a numpy array of (x, y)
        
    """
    infiles = applier.FilenameAssociations()
    outfiles = applier.FilenameAssociations()
    otherargs = applier.OtherInputs()
    
    infiles.img = img
    otherargs.tl = None
    otherargs.tr = None
    otherargs.bl = None
    otherargs.br = None
    otherargs.nullVal = imgInfo.nodataval[0]
    if otherargs.nullVal is None:
        otherargs.nullVal = 0

    applier.apply(findCorners, infiles, outfiles, otherargs)
    
    corners = numpy.array([
        otherargs.tl, 
        otherargs.tr, 
        otherargs.bl, 
        otherargs.br, 
    ])
    return corners 


def findCorners(info, inputs, outputs, otherargs):
    """
    Called from RIOS
    
    Checks non-null area of image block. Finds extremes, records coords
    of extremes against those already in otherargs. 
    
    Note that the logic is very specific to the orientation of the usual Landsat
    descending pass imagery. The same logic should not be applied to swathes 
    oriented in other, e.g. for other satellites. 
    
    """
    (xblock, yblock) = info.getBlockCoordArrays()
    
    nonnull = (inputs.img != otherargs.nullVal).all(axis=0)
    
    xNonnull = xblock[nonnull]
    yNonnull = yblock[nonnull]
    
    if len(xNonnull) > 0:
        topNdx = numpy.argmax(yNonnull)
        topXY = (xNonnull[topNdx], yNonnull[topNdx])
        leftNdx = numpy.argmin(xNonnull)
        leftXY = (xNonnull[leftNdx], yNonnull[leftNdx])
        bottomNdx = numpy.argmin(yNonnull)
        bottomXY = (xNonnull[bottomNdx], yNonnull[bottomNdx])
        rightNdx = numpy.argmax(xNonnull)
        rightXY = (xNonnull[rightNdx], yNonnull[rightNdx])

        # If these are more extreme than those already in otherargs, replace them
        if otherargs.tl is None or topXY[1] > otherargs.tl[1]:
            otherargs.tl = topXY
        if otherargs.tr is None or rightXY[0] > otherargs.tr[0]:
            otherargs.tr = rightXY
        if otherargs.bl is None or leftXY[0] < otherargs.bl[0]:
            otherargs.bl = leftXY
        if otherargs.br is None or bottomXY[1] < otherargs.br[1]:
            otherargs.br = bottomXY
    

def findNadirLine(corners):
    """
    Return the equation of the nadir line, from the given corners of the swathe. 
    Returns a numpy array of [b, m], for the equation
        y = mx + b
    giving the y coordinate of the nadir as a function of the x coordinate. 

    """
    # Find the top and bottom mid-points. 
    topMid = (corners[0] + corners[1]) / 2.0
    bottomMid = (corners[2] + corners[3]) / 2.0
    
    slope = (topMid[1] - bottomMid[1]) / (topMid[0] - bottomMid[0])
    intercept = topMid[1] - slope * topMid[0]
    
    coeffs = numpy.array([intercept, slope])
    return coeffs
   

def satAzLeftRight(nadirLine):
    """
    Calculate the satellite azimuth for the left and right sides of the nadir line.
    Assume that the satellite azimuth vector is at right angles to the nadir line
    (which is not really true, but reasonably close), and that there are only
    two possibilities, as a pixel is either to the left or to the right of the nadir
    line. 
    
    Return a numpy array of [satAzLeft, satAzRight], with angles in radians, 
    in the range [-pi, pi]
    
    """
    slope = nadirLine[1]
    # Slope of a line perpendicular to the nadir
    slopePerp = -1 / slope
    
    # Azimuth for pixels to the left of the line
    azimuthLeft = numpy.pi / 2.0 - numpy.arctan(slopePerp)
    # Azimuth for pixels to the right is directly opposite
    azimuthRight = azimuthLeft - numpy.pi
    
    return numpy.array([azimuthLeft, azimuthRight])


def localRadius(latitude):
    """
    Calculate a local radius of curvature, for the given geodetic latitude. 
    This approximates the earth curvature at this latitude. The given 
    latitude is in degrees. This is probably overkill, given some of the other 
    approximations I am making....
    
    """
    latRadians = numpy.radians(latitude)
    
    # Earth axis lengths
    a = osr.SRS_WGS84_SEMIMAJOR
    invFlat = osr.SRS_WGS84_INVFLATTENING
    f = 1 / invFlat
    eSqd = 2*f - f**2

    # Radius of curvature
    R = a / numpy.sqrt(1 - eSqd * numpy.sin(latRadians)**2)
    return R


def sunAnglesForPoints(corners, mtlInfo):
    """
    Return array of sun azimuth and zenith for each of the (x, y) points given
    in corners. First column is azimuth, second column is zenith, both are in radians. 
    
    """
    # For now, just kludge using the scene centre sun angles. Come back to this
    # when I have the proper algorithm......
    sunAzimuth = numpy.radians(float(mtlInfo['SUN_AZIMUTH']))
    sunZenith = numpy.radians(90.0 - float(mtlInfo['SUN_ELEVATION']))
    
    sunAngles = numpy.array([(sunAzimuth, sunZenith) for i in range(len(corners))])
    return sunAngles

# def sunAnglesAtPoint(


def makeAnglesImage(cmdargs, nadirLine, corners, cornerSunAngles, satAzimuth, imgInfo):
    """
    Make a single output image file of the sun and satellite angles for every
    pixel in the template image. 
    
    """
    imgInfo  = fileinfo.ImageInfo(cmdargs.templateimg)
    
    infiles = applier.FilenameAssociations()
    outfiles = applier.FilenameAssociations()
    otherargs = applier.OtherInputs()
    controls = applier.ApplierControls()
    
    infiles.img = cmdargs.templateimg
    outfiles.angles = cmdargs.outfile
    
    (ctrLat, ctrLong) = getCtrLatLong(imgInfo)
    otherargs.R = localRadius(ctrLat)
    otherargs.nadirLine = nadirLine
    otherargs.corners = corners
    otherargs.cornerSunAngles = cornerSunAngles
    otherargs.satAltitude = 705000      # Landsat nominal altitude in metres
    otherargs.satAzimuth = satAzimuth
    otherargs.radianScale = 100        # Store pixel values as (radians * radianScale)
    controls.setStatsIgnore(500)
    
    applier.apply(makeAngles, infiles, outfiles, otherargs, controls=controls)


def makeAngles(info, inputs, outputs, otherargs):
    """
    Called from RIOS
    
    Make 4-layer sun and satellite angles for the image block
    
    """
    (xblock, yblock) = info.getBlockCoordArrays()
    # Nadir line coefficients of y=mx+b
    (b, m) = otherargs.nadirLine
    
    # Distance of each pixel from the nadir line
    dist = numpy.absolute((m * xblock - yblock + b) / numpy.sqrt(m**2 + 1))
    
    # Zenith angle assuming a flat earth
    satZenith = numpy.arctan(dist / otherargs.satAltitude)
    
    # Adjust satZenith for earth curvature. This is a very simple approximation, but 
    # the adjustment is less than one degree anyway, so this is accurate enough. 
    curveAngle = numpy.arctan(dist / otherargs.R)
    satZenith += curveAngle
    
    # Work out whether we are left or right of the nadir line
    isLeft = (yblock - (m * xblock + b)) > 0
    (satAzimuthLeft, satAzimuthRight) = otherargs.satAzimuth
    satAzimuth = numpy.where(isLeft, satAzimuthLeft, satAzimuthRight)
    
    # Sun angles are fixed for now, but we should be interpolating between the corner values....
    sunAzimuth = numpy.zeros(xblock.shape, dtype=numpy.float32) + otherargs.cornerSunAngles[0][0]
    sunZenith = numpy.zeros(xblock.shape, dtype=numpy.float32) + otherargs.cornerSunAngles[0][1]
    
    angleStack = numpy.array([satAzimuth, satZenith, sunAzimuth, sunZenith])
    angleStackDN = angleStack * otherargs.radianScale
    
    outputs.angles = numpy.round(angleStackDN).astype(numpy.int16)


def getCtrLatLong(imgInfo):
    """
    Return the lat/long of the centre of the image as
        (ctrLat, ctrLong)
        
    """
    cornerLatLong = imgInfo.getCorners(outEPSG=4326)
    (ul_long, ul_lat, ur_long, ur_lat, lr_long, lr_lat, ll_long, ll_lat) = cornerLatLong
    ctrLat = numpy.array([ul_lat, ur_lat, lr_lat, ll_lat]).mean()
    ctrLong = numpy.array([ul_long, ur_long, lr_long, ll_long]).mean()
    return (ctrLat, ctrLong)


if __name__ == "__main__":
    mainRoutine()
