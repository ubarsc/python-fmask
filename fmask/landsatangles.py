#!/usr/bin/env python
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
"""
Functions relating to estimating the per-pixel sun and satellite angles for 
a given Landsat image. These are rough estimates, using the generic
characteristics of the Landsat 5 platform, and are not particularly accurate,
but good enough for the current purposes. 

Historically, the USGS have not supplied satellite zenith/azimuth angles, and have only 
supplied scene-centre values for sun zenith/azimuth angles. Since the satellite
view geometry is important in correctly tracking a shadow when matching shadows
to their respective clouds, the Fmask algorithm requires good estimates of all these
angles. The routines contained here are used to derive per-pixel estimates of 
these angles. 

As of mid-2016, the USGS are planning to supply sufficient information to calculate
these angles directly from orbit ephemeris data. When that comes about, it seems likely
that the need for the routines here will diminish, but any data downloaded from USGS
prior to then will still require this approach, as the associated angle metadata will 
not be present. 

The core Fmask code in this package is adaptable enough to be configured for either 
approach. 

The general approach for satellite angles is to estimate the nadir line by running it
down the middle of the image data area. The satellite azimuth is assumed to be
at right angles to this nadir line, which is only roughly correct. For the whisk-broom 
sensors on Landsat-5 and Landsat-7, this angle is not 90 degrees, but is affected by 
earth rotation and is latitude dependent. For Landsat-8, the scan line is at 
right angles, due to the compensation for earth rotation, but the push-broom is 
made up of sub-modules which point in slightly different directions, giving 
slightly different satellite azimuths along the scan line. None of these effects
are included in the current estimates. The satellite zenith is estimated based on the
nadir point, the scan-line, and the assumed satellite altitude, and includes the
appropriate allowance for earth curvature. 

Because this works by searching the imagery for the non-null area, and assumes that 
this represents a full-swath image, it would not work for a subset of a full image. 

The sun angles are approximated using the algorithm found in the Fortran code with
6S (Second Simulation of the Satellite Signal in the Solar Spectrum). The subroutine
in question is the POSSOL() routine. I translated the Fortran code into Python for
inclusion here. 

"""
from __future__ import print_function, division

import datetime

import numpy
from osgeo import osr

from rios import applier
from rios import fileinfo


def findImgCorners(img, imgInfo):
    """
    Find the corners of the data within the given template image
    Return a numpy array of (x, y) coordinates. The array has 2 columns, for X and Y.
    Each row is a corner, in the order:

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
    eSqd = 2 * f - f**2

    # Radius of curvature
    R = a / numpy.sqrt(1 - eSqd * numpy.sin(latRadians)**2)
    return R


def sunAnglesForExtent(imgInfo, mtlInfo):
    """
    Return array of sun azimuth and zenith for each of the corners of the image
    extent. Note that this is the raster extent, not the corners of the swathe.

    The algorithm used here has been copied from the 6S possol() subroutine. The
    Fortran code I copied it from was .... up to the usual standard in 6S. So, the
    notation is not always clear.

    """
    cornerLatLong = imgInfo.getCorners(outEPSG=4326)
    (ul_long, ul_lat, ur_long, ur_lat, lr_long, lr_lat, ll_long, ll_lat) = cornerLatLong
    pts = numpy.array([
        [ul_long, ul_lat],
        [ur_long, ur_lat],
        [ll_long, ll_lat],
        [lr_long, lr_lat]
    ])
    longDeg = pts[:, 0]
    latDeg = pts[:, 1]

    # Date/time in UTC
    dateStr = mtlInfo['DATE_ACQUIRED']
    timeStr = mtlInfo['SCENE_CENTER_TIME'].replace('Z', '')
    ymd = [int(i) for i in dateStr.split('-')]
    dateObj = datetime.date(ymd[0], ymd[1], ymd[2])
    julianDay = (dateObj - datetime.date(ymd[0], 1, 1)).days + 1
    juldayYearEnd = (datetime.date(ymd[0], 12, 31) - datetime.date(ymd[0], 1, 1)).days + 1
    # Julian day as a proportion of the year
    jdp = julianDay / juldayYearEnd
    # Hour in UTC
    hms = [float(x) for x in timeStr.split(':')]
    hourGMT = hms[0] + hms[1] / 60.0 + hms[2] / 3600.0

    (sunAz, sunZen) = sunAnglesForPoints(latDeg, longDeg, hourGMT, jdp)

    sunAngles = numpy.vstack((sunAz, sunZen)).T
    return sunAngles


def sunAnglesForPoints(latDeg, longDeg, hourGMT, jdp):
    """
    Calculate sun azimuth and zenith for the given location(s), for the given
    time of year. jdp is the julian day as a proportion, ranging from 0 to 1, where
    Jan 1 is 1.0/365 and Dec 31 is 1.0.
    Location is given in latitude/longitude, in degrees, and can be arrays to
    calculate for multiple locations. hourGMT is a decimal hour number giving the time
    of day in GMT (i.e. UTC).

    Return a tuple of (sunAz, sunZen). If latDeg and longDeg are arrays, then returned
    values will be arrays of the same shape.

    """
    latRad = numpy.radians(latDeg)
    # Express jdp in radians
    jdpr = jdp * 2 * numpy.pi

    # Now work out the solar position. This is copied from the 6S code, but
    # is also documented in the 6S manual. The notation
    a = numpy.array([0.000075, 0.001868, 0.032077, 0.014615, 0.040849])
    meanSolarTime = hourGMT + longDeg / 15.0
    localSolarDiff = (a[0] + a[1] * numpy.cos(jdpr) - a[2] * numpy.sin(jdpr) -
        a[3] * numpy.cos(2 * jdpr) - a[4] * numpy.sin(2 * jdpr)) * 12 * 60 / numpy.pi
    trueSolarTime = meanSolarTime + localSolarDiff / 60 - 12.0
    # Hour as an angle
    ah = trueSolarTime * numpy.radians(15)

    b = numpy.array([0.006918, 0.399912, 0.070257, 0.006758, 0.000907, 0.002697, 0.001480])
    delta = (b[0] - b[1] * numpy.cos(jdpr) + b[2] * numpy.sin(jdpr) - 
        b[3] * numpy.cos(2. * jdpr) + b[4] * numpy.sin(2. * jdpr) - 
        b[5] * numpy.cos(3. * jdpr) + b[6] * numpy.sin(3. * jdpr))

    cosSunZen = (numpy.sin(latRad) * numpy.sin(delta) +
        numpy.cos(latRad) * numpy.cos(delta) * numpy.cos(ah))
    sunZen = numpy.arccos(cosSunZen)

    # sun azimuth from south, turning west (yeah, I know, weird, isn't it....)
    sinSunAzSW = numpy.cos(delta) * numpy.sin(ah) / numpy.sin(sunZen)
    sinSunAzSW = sinSunAzSW.clip(-1.0, 1.0)

    # This next bit seems to be to get the azimuth in the correct quadrant. I do
    # not really understand it.
    cosSunAzSW = (-numpy.cos(latRad) * numpy.sin(delta) +
        numpy.sin(latRad) * numpy.cos(delta) * numpy.cos(ah)) / numpy.sin(sunZen)
    sunAzSW = numpy.arcsin(sinSunAzSW)
    sunAzSW = numpy.where(cosSunAzSW <= 0, numpy.pi - sunAzSW, sunAzSW)
    sunAzSW = numpy.where((cosSunAzSW > 0) & (sinSunAzSW <= 0), 2 * numpy.pi + sunAzSW, sunAzSW)

    # Now convert to azimuth from north, turning east, as is usual convention
    sunAz = sunAzSW + numpy.pi
    # Keep within [0, 2pi] range
    sunAz = numpy.where(sunAz > 2 * numpy.pi, sunAz - 2 * numpy.pi, sunAz)

    return (sunAz, sunZen)


def makeAnglesImage(templateimg, outfile, nadirLine, extentSunAngles, satAzimuth, imgInfo):
    """
    Make a single output image file of the sun and satellite angles for every
    pixel in the template image.

    """
    imgInfo = fileinfo.ImageInfo(templateimg)

    infiles = applier.FilenameAssociations()
    outfiles = applier.FilenameAssociations()
    otherargs = applier.OtherInputs()
    controls = applier.ApplierControls()

    infiles.img = templateimg
    outfiles.angles = outfile

    (ctrLat, ctrLong) = getCtrLatLong(imgInfo)
    otherargs.R = localRadius(ctrLat)
    otherargs.nadirLine = nadirLine
    otherargs.xMin = imgInfo.xMin
    otherargs.xMax = imgInfo.xMax
    otherargs.yMin = imgInfo.yMin
    otherargs.yMax = imgInfo.yMax
    otherargs.extentSunAngles = extentSunAngles
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

    # Interpolate the sun angles from those calculated at the corners of the whole raster extent
    (xMin, xMax, yMin, yMax) = (otherargs.xMin, otherargs.xMax, otherargs.yMin, otherargs.yMax)
    sunAzimuth = bilinearInterp(xMin, xMax, yMin, yMax, otherargs.extentSunAngles[:, 0], xblock, yblock)
    sunZenith = bilinearInterp(xMin, xMax, yMin, yMax, otherargs.extentSunAngles[:, 1], xblock, yblock)

    angleStack = numpy.array([satAzimuth, satZenith, sunAzimuth, sunZenith])
    angleStackDN = angleStack * otherargs.radianScale

    outputs.angles = numpy.round(angleStackDN).astype(numpy.int16)


def bilinearInterp(xMin, xMax, yMin, yMax, cornerVals, x, y):
    """
    Evaluate the given value on a grid of (x, y) points. The exact value is given
    on a set of corner points (top-left, top-right, bottom-left, bottom-right).
    The corner locations are implied from xMin, xMax, yMin, yMax.

    """
    p = (y - yMin) / (yMax - yMin)
    q = (x - xMin) / (xMax - xMin)

    # Give the known corner values some simple names
    (tl, tr, bl, br) = cornerVals

    # Calculate the interpolated values
    vals = tr * p * q + tl * p * (1 - q) + br * (1 - p) * q + bl * (1 - p) * (1 - q)
    return vals


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
