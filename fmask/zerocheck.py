"""
Utility function to determine if a particular band in a 
file is all zeros. This is useful for Landsat8 since
although the file and metadata exist for thermal, it may
be all zeroes and hence unusable form cloud masking.
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

from rios import applier
from rios import fileinfo


def isBandAllZeroes(filename, band=0):
    """
    Checks the specified band within a file to see if it is all zeroes.
    band should be 0 based.
    
    Returns True if all zeroes, False otherwise.
    
    This function firstly checks the stats if present and uses this
    and assumes that they are correct.
    
    If they are not present, then the band is read and the values determined.
    
    """

    info = fileinfo.ImageFileStats(filename)
    maxVal = info[band].max
    if maxVal is not None:
        # we have valid stats
        return maxVal == 0
    
    # otherwise read the thing with rios
    infiles = applier.FilenameAssociations()
    infiles.input = filename
    
    outfiles = applier.FilenameAssociations()
    
    otherArgs = applier.OtherInputs()
    otherArgs.nonZeroFound = False
    otherArgs.band = band
    
    applier.apply(riosAllZeroes, infiles, outfiles, otherArgs)
    
    return not otherArgs.nonZeroFound
    
    
def riosAllZeroes(info, inputs, outputs, otherArgs):
    """
    Called from RIOS. Does the actual work.
    
    """
    if inputs.input[otherArgs.band].max() > 0:
        otherArgs.nonZeroFound = True
    
