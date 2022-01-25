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

import os
import numpy

# Fail slightly less drastically when running from ReadTheDocs
if os.getenv('READTHEDOCS', default='False') != 'True':
    from . import _valueindexes


class ValueIndexesError(Exception):
    pass


class NonIntTypeError(ValueIndexesError):
    pass


class RangeError(ValueIndexesError):
    pass


class ValueIndexes(object):
    """
    An object which contains the indexes for every value in a given array.
    This class is intended to mimic the reverse_indices clause in IDL,
    only nicer. 
    
    Takes an array, works out what unique values exist in this array. Provides
    a method which will return all the indexes into the original array for 
    a given value. 
    
    The array must be of an integer-like type. Floating point arrays will
    not work. If one needs to look at ranges of values within a float array,
    it is possible to use numpy.digitize() to create an integer array 
    corresponding to a set of bins, and then use ValueIndexes with that. 
    
    Example usage, for a given array a::

        valIndexes = ValueIndexes(a)
        for val in valIndexes.values:
            ndx = valIndexes.getIndexes(val)
            # Do something with all the indexes
    
    
    This is a much faster and more efficient alternative to something like::

        values = numpy.unique(a)
        for val in values:
            mask = (a == val)
            # Do something with the mask

    The point is that when a is large, and/or the number of possible values 
    is large, this becomes very slow, and can take up lots of memory. Each 
    loop iteration involves searching through a again, looking for a different 
    value. This class provides a much more efficient way of doing the same 
    thing, requiring only one pass through a. When a is small, or when the 
    number of possible values is very small, it probably makes little difference. 
    
    If one or more null values are given to the constructor, these will not 
    be counted, and will not be available to the getIndexes() method. This 
    makes it more memory-efficient, so it doesn't store indexes of a whole 
    lot of nulls. 
    
    A ValueIndexes object has the following attributes:
    
    * **values**            Array of all values indexed
    * **counts**            Array of counts for each value
    * **nDims**             Number of dimensions of original array
    * **indexes**           Packed array of indexes
    * **start**             Starting points in indexes array for each value
    * **end**               End points in indexes for each value
    * **valLU**             Lookup table for each value, to find it in the values array without explicitly searching. 
    * **nullVals**          Array of the null values requested. 
    
    Limitations:
    The array index values are handled using unsigned 32bit int values, so 
    it won't work if the data array is larger than 4Gb. I don't think it would
    fail very gracefully, either. 
    
    """
    def __init__(self, a, nullVals=[]):
        """
        Creates a ValueIndexes object for the given array a. 
        A sequence of null values can be given, and these will not be included
        in the results, so that indexes for these cannot be determined. 
        
        """
        if not numpy.issubdtype(a.dtype, numpy.integer):
            raise NonIntTypeError("ValueIndexes only works on integer-like types. Array is %s"%a.dtype)
             
        if numpy.isscalar(nullVals):
            self.nullVals = [nullVals]
        else:
            self.nullVals = nullVals

        # Get counts of all values in a
        minval = a.min()
        maxval = a.max()
        (counts, binEdges) = numpy.histogram(a, range=(minval, maxval + 1), 
            bins=(maxval - minval + 1))
            
        # Mask counts for any requested null values. 
        maskedCounts = counts.copy()
        for val in self.nullVals:
            maskedCounts[binEdges[:-1]==val] = 0
        self.values = binEdges[:-1][maskedCounts>0].astype(a.dtype)
        self.counts = maskedCounts[maskedCounts>0]
        
        # Allocate space to store all indexes
        totalCounts = self.counts.sum()
        self.nDims = a.ndim
        self.indexes = numpy.zeros((totalCounts, a.ndim), dtype=numpy.uint32)
        self.end = self.counts.cumsum()
        self.start = self.end - self.counts
        
        if len(self.values) > 0:
            # A lookup table to make searching for a value very fast.
            valrange = numpy.array([self.values.min(), self.values.max()])
            numLookups = valrange[1] - valrange[0] + 1
            maxUint32 = 2**32 - 1
            if numLookups > maxUint32:
                raise RangeError("Range of different values is too great for uint32")
            self.valLU = numpy.zeros(numLookups, dtype=numpy.uint32)
            self.valLU.fill(maxUint32)     # A value to indicate "not found", must match _valndxFunc below
            self.valLU[self.values - self.values[0]] = range(len(self.values))

            # For use within C. For each value, the current index 
            # into the indexes array. A given element is incremented whenever it finds
            # a new element of that value. 
            currentIndex = self.start.copy().astype(numpy.uint32)

            _valueindexes.valndxFunc(a, self.indexes, valrange[0], valrange[1], 
                        self.valLU, currentIndex)

    def getIndexes(self, val):
        """
        Return a set of indexes into the original array, for which the
        value in the array is equal to val. 
        
        """
        # Find where this value is listed. 
        valNdx = (self.values == val).nonzero()[0]
        
        # If this value is not actually in those listed, then we 
        # must return empty indexes
        if len(valNdx) == 0:
            start = 0
            end = 0
        else:
            # The index into counts, etc. for this value. 
            valNdx = valNdx[0]
            start = self.start[valNdx]
            end = self.end[valNdx]
            
        # Create a tuple of index arrays, one for each index of the original array. 
        ndx = ()
        for i in range(self.nDims):
            ndx += (self.indexes[start:end, i], )
        return ndx
