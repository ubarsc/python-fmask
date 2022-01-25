"""
Module for doing checks of visible band and reporting
and saturation. Note that this works off the original 
radiance file, not the TOA reflectance.
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
from rios import applier, cuiprogress
from . import config


def makeSaturationMask(fmaskConfig, radiancefile, outMask):
    """
    Checks the radianceFile and creates a mask with
    1's where there is saturation in one of more visible bands.
    0 otherwise. 
    
    The fmaskConfig parameter should be an instance of
    :class:`fmask.config.FmaskConfig`. This is used to determine
    which bands are visible.
    
    This mask is advisible since the whiteness test Eqn 2. 
    and Equation 6 are affected by saturated pixels and
    may determine a pixel is not cloud when it is.
    
    The format of outMask will be the current 
    RIOS default format.
    
    It is assumed that the input radianceFile has values in the
    range 0-255 and saturated pixels are set to 255.
    
    """
    inputs = applier.FilenameAssociations()
    inputs.radiance = radiancefile
    
    outputs = applier.FilenameAssociations()
    outputs.mask = outMask
    
    otherargs = applier.OtherInputs()
    otherargs.radianceBands = fmaskConfig.bands
    
    controls = applier.ApplierControls()
    controls.progress = cuiprogress.GDALProgressBar()
    
    applier.apply(riosSaturationMask, inputs, outputs, 
                otherargs, controls=controls)


def riosSaturationMask(info, inputs, outputs, otherargs):
    """
    Called from RIOS. Does the actual saturation test. Currently assumes that
    only 8-bit radiance inputs can be saturated, but if this turns out
    not to be true, we can come back to this. 
    
    """
    if inputs.radiance.dtype == numpy.uint8:
        blue = otherargs.radianceBands[config.BAND_BLUE]
        green = otherargs.radianceBands[config.BAND_GREEN]
        red = otherargs.radianceBands[config.BAND_RED]

        satMaskList = []
        for band in [blue, green, red]:
            satMaskList.append(inputs.radiance[band] == 255)

        outputs.mask = numpy.array(satMaskList).astype(numpy.uint8)
    else:
        # Assume that anything larger than 8-bit is immune to saturation
        outShape = (3, ) + inputs.radiance[0].shape
        outputs.mask = numpy.zeros(outShape, dtype=numpy.uint8)
    
