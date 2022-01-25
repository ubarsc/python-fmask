"""
Module to implement filling of local minima in a raster surface. 

The algorithm is from 
    Soille, P., and Gratin, C. (1994). An efficient algorithm for drainage network
        extraction on DEMs. J. Visual Communication and Image Representation. 
        5(2). 181-189. 
        
The algorithm is intended for hydrological processing of a DEM, but is used by the 
Fmask cloud shadow algorithm as part of its process for finding local minima which 
represent potential shadow objects. 

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

import os
import numpy
from scipy.ndimage import grey_dilation

# Fail slightly less drastically when running from ReadTheDocs
if os.getenv('READTHEDOCS', default='False') != 'True':
    from . import _fillminima


def fillMinima(img, nullval, boundaryval):
    """
    Fill all local minima in the input img. The input
    array should be a numpy 2-d array. This function returns
    an array of the same shape and datatype, with the same contents, but
    with local minima filled using the reconstruction-by-erosion algorithm. 
    
    """
    (nrows, ncols) = img.shape
    dtype = img.dtype
    nullmask = (img == nullval)
    nonNullmask = numpy.logical_not(nullmask)
    (hMax, hMin) = (int(img[nonNullmask].max()), int(img[nonNullmask].min()))
    boundaryval = max(boundaryval, hMin)
    img2 = numpy.zeros((nrows, ncols), dtype=dtype)
    img2.fill(hMax)

    if nullmask.sum() > 0:
        nullmaskDilated = grey_dilation(nullmask, size=(3, 3))
        innerBoundary = nullmaskDilated ^ nullmask
        (boundaryRows, boundaryCols) = numpy.where(innerBoundary)
    else:
        img2[0, :] = img[0, :]
        img2[-1, :] = img[-1, :]
        img2[:, 0] = img[:, 0]
        img2[:, -1] = img[:, -1]
        (boundaryRows, boundaryCols) = numpy.where(img2!=hMax)

    # on some systems (32 bit only?) numpy.where returns int32
    # rather than int64. Convert so we don't have to handle both in C.
    boundaryRows = boundaryRows.astype(numpy.int64)
    boundaryCols = boundaryCols.astype(numpy.int64)

    _fillminima.fillMinima(img, img2, hMin, hMax, nullmask, boundaryval,
                        boundaryRows, boundaryCols)    
    
    img2[nullmask] = nullval
    
    return img2
