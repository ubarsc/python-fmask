"""
Classes for handling the various metadata files which come with Sentinel-2.

Currently only has a class for the tile-based metadata file. 

"""
from __future__ import print_function, division

import datetime
from xml.etree import ElementTree

import numpy
from osgeo import osr

from . import fmaskerrors

class Sen2TileMeta(object):
    """
    Metadata for a single 100km tile
    """
    def __init__(self, filename=None):
        """
        Constructor takes a filename for the XML file of tile-based metadata. 
        
        """
        f = open(filename)
        
        root = ElementTree.fromstring(f.read())
        # Stoopid XML namespace prefix
        nsPrefix = root.tag[:root.tag.index('}')+1]
        nsDict = {'n1':nsPrefix[1:-1]}
        
        generalInfoNode = root.find('n1:General_Info', nsDict)
        # N.B. I am still not entirely convinced that this SENSING_TIME is really 
        # the acquisition time, but the documentation is rubbish. 
        sensingTimeNode = generalInfoNode.find('SENSING_TIME')
        sensingTimeStr = sensingTimeNode.text.strip()
        self.datetime = datetime.datetime.strptime(sensingTimeStr, "%Y-%m-%dT%H:%M:%S.%fZ")
        tileIdNode = generalInfoNode.find('TILE_ID')
        tileIdFullStr = tileIdNode.text.strip()
        self.tileId = tileIdFullStr.split('_')[-2]
        self.satId = tileIdFullStr[:3]
        self.procLevel = tileIdFullStr[13:16]    # Not sure whether to use absolute pos or split by '_'....
        
        geomInfoNode = root.find('n1:Geometric_Info', nsDict)
        geocodingNode = geomInfoNode.find('Tile_Geocoding')
        epsgNode = geocodingNode.find('HORIZONTAL_CS_CODE')
        self.epsg = epsgNode.text.split(':')[1]
        
        # Dimensions of images at different resolutions. 
        self.dimsByRes = {}
        sizeNodeList = geocodingNode.findall('Size')
        for sizeNode in sizeNodeList:
            res = sizeNode.attrib['resolution']
            nrows = int(sizeNode.find('NROWS').text)
            ncols = int(sizeNode.find('NCOLS').text)
            self.dimsByRes[res] = (nrows, ncols)

        # Upper-left corners of images at different resolutions. As far as I can
        # work out, these coords appear to be the upper left corner of the upper left
        # pixel, i.e. equivalent to GDAL's convention. This also means that they
        # are the same for the different resolutions, which is nice. 
        self.ulxyByRes = {}
        posNodeList = geocodingNode.findall('Geoposition')
        for posNode in posNodeList:
            res = posNode.attrib['resolution']
            ulx = float(posNode.find('ULX').text)
            uly = float(posNode.find('ULY').text)
            self.ulxyByRes[res] = (ulx, uly)
        
        # Sun and satellite angles. 
        tileAnglesNode = geomInfoNode.find('Tile_Angles')
        sunZenithNode = tileAnglesNode.find('Sun_Angles_Grid').find('Zenith')
        self.angleGridXres = float(sunZenithNode.find('COL_STEP').text)
        self.angleGridYres = float(sunZenithNode.find('ROW_STEP').text)
        self.sunZenithGrid = self.makeValueArray(sunZenithNode.find('Values_List'))
        sunAzimuthNode = tileAnglesNode.find('Sun_Angles_Grid').find('Azimuth')
        self.sunAzimuthGrid = self.makeValueArray(sunAzimuthNode.find('Values_List'))
        self.anglesGridShape = self.sunAzimuthGrid.shape
        
        # Now build up the viewing angle per grid cell, from the separate layers
        # given for each detector for each band. Initially I am going to keep
        # the bands separate, just to see how that looks. 
        # The names of things in the XML suggest that these are view angles,
        # but the numbers suggest that they are angles as seen from the pixel's 
        # frame of reference on the ground, i.e. they are in fact what we ultimately want. 
        viewingAngleNodeList = tileAnglesNode.findall('Viewing_Incidence_Angles_Grids')
        self.viewZenithDict = self.buildViewAngleArr(viewingAngleNodeList, 'Zenith')
        self.viewAzimuthDict = self.buildViewAngleArr(viewingAngleNodeList, 'Azimuth')
        
        # Make a guess at the coordinates of the angle grids. These are not given 
        # explicitly in the XML, and don't line up exactly with the other grids, so I am 
        # making a rough estimate. Because the angles don't change rapidly across these 
        # distances, it is not important if I am a bit wrong (although it would be nice
        # to be exactly correct!). 
        (ulx, uly) = self.ulxyByRes["10"]
        self.anglesULXY = (ulx - self.angleGridXres / 2.0, uly + self.angleGridYres / 2.0)
    
    @staticmethod
    def makeValueArray(valuesListNode):
        """
        Take a <Values_List> node from the XML, and return an array of the values contained
        within it. This will be a 2-d numpy array of float32 values (should I pass the dtype in??)
        
        """
        valuesList = valuesListNode.findall('VALUES')
        vals = []
        for valNode in valuesList:
            text = valNode.text
            vals.append([numpy.float32(x) for x in text.strip().split()])
        return numpy.array(vals)
    
    def buildViewAngleArr(self, viewingAngleNodeList, angleName):
        """
        Build up the named viewing angle array from the various detector strips given as
        separate arrays. I don't really understand this, and may need to re-write it once
        I have worked it out......
        
        The angleName is one of 'Zenith' or 'Azimuth'.
        Returns a dictionary of 2-d arrays, keyed by the bandId string. 
        """
        angleArrDict = {}
        for viewingAngleNode in viewingAngleNodeList:
            bandId = viewingAngleNode.attrib['bandId']
            angleNode = viewingAngleNode.find(angleName)
            angleArr = self.makeValueArray(angleNode.find('Values_List'))
            if bandId not in angleArrDict:
                angleArrDict[bandId] = angleArr
            else:
                mask = (~numpy.isnan(angleArr))
                angleArrDict[bandId][mask] = angleArr[mask]
        return angleArrDict

    def getUTMzone(self):
        """
        Return the UTM zone of the tile, as an integer
        """
        if not (self.epsg.startswith("327") or self.epsg.startswith("326")):
            raise fmaskerrors.Sen2MetaError("Cannot determine UTM zone from EPSG:{}".format(self.epsg))
        return int(self.epsg[3:])
    
    def getCtrXY(self):
        """
        Return the (X, Y) coordinates of the scene centre (in image projection, generally UTM)
        """
        (nrows, ncols) = self.dimsByRes['10']
        (ctrRow, ctrCol) = (nrows // 2, ncols // 2)
        (ulx, uly) = self.ulxyByRes['10']
        (ctrX, ctrY) = (ulx + ctrCol * 10, uly - ctrRow * 10)
        return (ctrX, ctrY)
    
    def getCtrLongLat(self):
        """
        Return the (longitude, latitude) of the scene centre
        """
        (ctrX, ctrY) = self.getCtrXY()
        srUTM = osr.SpatialReference()
        srUTM.ImportFromEPSG(int(self.epsg))
        srLL = osr.SpatialReference()
        srLL.ImportFromEPSG(4326)
        tr = osr.CoordinateTransformation(srUTM, srLL)
        (longitude, latitude, z) = tr.TransformPoint(ctrX, ctrY)
        return (longitude, latitude)

    


