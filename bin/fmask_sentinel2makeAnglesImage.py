#!/usr/bin/env python
"""
Make a 4-layer image of satellite and sun angles, from the given tile metadata
file. 

The sun angles are exactly as provided in the XML. The satellite angles are more 
complicated, because they vary slightly with each band. I have done some limited 
testing and concluded that the variation is generally quite small, and so the 
layers written here are the per-pixel averages over all bands. The satellite
zenith this appears to vary only a couple of degrees across bands, while the 
variation for satellite azimuth is somewhat larger, depending on which bands. I
think that the 60m bands vary a bit more, but not sure. This is also complicated
by the variation across the scan, due to the different view angles of each detector
module. The pushbroom appears to be made up of 12 such modules, each looking
in slightly different directions, and each band within each module looking slightly
differently. Complicated...... sigh......

A bit of a description of how the instrument is structured can be found at
    https://sentinels.copernicus.eu/web/sentinel/technical-guides/sentinel-2-msi/msi-instrument

"""
from __future__ import print_function, division

import sys
import argparse

import numpy
from osgeo import gdal
gdal.UseExceptions()
from osgeo import osr

from rios import applier
from rios import calcstats
from rios import cuiprogress

from fmask import sen2meta

# This scale value will convert between DN and radians in output image file, 
#    radians = dn * SCALE_TO_RADIANS
SCALE_TO_RADIANS = 0.01


def getCmdargs():
    """
    Get commandline arguments
    """
    p = argparse.ArgumentParser()
    p.add_argument("-i", "--infile", help="Input sentinel-2 tile metafile")
    p.add_argument("-o", "--outfile", help="Output angles image file")
    cmdargs = p.parse_args()
    if cmdargs.infile is None or cmdargs.outfile is None:
        p.print_help()
        sys.exit()
        
    return cmdargs


def mainRoutine():
    """
    Main routine
    """
    cmdargs = getCmdargs()
    
    info = sen2meta.Sen2TileMeta(filename=cmdargs.infile)
    
    ds = createOutfile(cmdargs.outfile, info)
    nullValDN = 1000
    
    # Get a sorted list of the Sentinel-2 band names. Note that sometimes this
    # is an incomplete list of band names, which appears to be due to a bug in 
    # earlier versions of ESA's processing software. I suspect it relates to 
    # Anomaly number 11 in the following page. 
    # https://sentinel.esa.int/web/sentinel/news/-/article/new-processing-baseline-for-sentinel-2-products
    bandNames = sorted(info.viewAzimuthDict.keys())
    
    # Mean over all bands
    satAzDeg = numpy.array([info.viewAzimuthDict[i] for i in bandNames])
    satAzDegMeanOverBands = satAzDeg.mean(axis=0)
    
    satZenDeg = numpy.array([info.viewZenithDict[i] for i in bandNames])
    satZenDegMeanOverBands = satZenDeg.mean(axis=0)

    sunAzDeg = info.sunAzimuthGrid
    
    sunZenDeg = info.sunZenithGrid
    
    stackDeg = numpy.array([satAzDegMeanOverBands, satZenDegMeanOverBands, sunAzDeg, sunZenDeg])
    stackRadians = numpy.radians(stackDeg)
    
    stackDN = numpy.round(stackRadians / SCALE_TO_RADIANS).astype(numpy.int16)
    nullmask = numpy.isnan(stackDeg)
    stackDN[nullmask] = nullValDN
    
    lnames = ['SatelliteAzimuth', 'SatelliteZenith', 'SunAzimuth', 'SunZenith']
    for i in range(ds.RasterCount):
        b = ds.GetRasterBand(i+1)
        b.WriteArray(stackDN[i])
        b.SetNoDataValue(nullValDN)
        b.SetDescription(lnames[i])
    calcstats.calcStats(ds, ignore=nullValDN, progress=cuiprogress.SilentProgress())
    del ds
        

def createOutfile(filename, info):
    """
    Create the empty output image file
    """
    drvr = gdal.GetDriverByName(applier.DEFAULTDRIVERNAME)
    (nrows, ncols) = info.anglesGridShape
    ds = drvr.Create(filename, ncols, nrows, 4, gdal.GDT_Int16, applier.DEFAULTCREATIONOPTIONS)
    gt = (info.anglesULXY[0], info.angleGridXres, 0, info.anglesULXY[1], 0.0, -info.angleGridYres)
    ds.SetGeoTransform(gt)
    
    sr = osr.SpatialReference()
    sr.ImportFromEPSG(int(info.epsg))
    ds.SetProjection(sr.ExportToWkt())
    return ds


if __name__ == "__main__":
    mainRoutine()

