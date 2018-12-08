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
this represents a full-swath image, it would not work for a subset of a full image. 

The sun angles are approximated using the algorithm found in the Fortran code with
6S (Second Simulation of the Satellite Signal in the Solar Spectrum). The subroutine
in question is the POSSOL() routine. I translated the Fortran code into Python for
inclusion here. 

"""
from __future__ import print_function, division

import sys
import argparse

from fmask import landsatangles
from fmask import config

from rios import fileinfo


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
    if (cmdargs.mtl is None or cmdargs.templateimg is None or 
            cmdargs.outfile is None):
        p.print_help()
        sys.exit(1)
    return cmdargs


def mainRoutine():
    """
    Main routine
    """
    cmdargs = getCmdargs()

    makeAngles(cmdargs.mtl, cmdargs.templateimg, cmdargs.outfile)


def makeAngles(mtlfile, templateimg, outfile):
    """
    Callable main routine
    """
    mtlInfo = config.readMTLFile(mtlfile)
    
    imgInfo = fileinfo.ImageInfo(templateimg)
    corners = landsatangles.findImgCorners(templateimg, imgInfo)
    nadirLine = landsatangles.findNadirLine(corners)
    
    extentSunAngles = landsatangles.sunAnglesForExtent(imgInfo, mtlInfo)
    satAzimuth = landsatangles.satAzLeftRight(nadirLine)
    
    landsatangles.makeAnglesImage(templateimg, outfile, 
        nadirLine, extentSunAngles, satAzimuth, imgInfo)
