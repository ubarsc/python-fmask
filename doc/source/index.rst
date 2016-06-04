.. _contents:

Python Fmask
========================================================

Introduction
------------

A set of command line utilities and Python modules that implement
the 'fmask' algorithm as published in:

Zhu, Z. and Woodcock, C.E. (2012). 
Object-based cloud and cloud shadow detection in Landsat imagery
Remote Sensing of Environment 118 (2012) 83-94. 
    
and
    
Zhu, Z., Wang, S. and Woodcock, C.E. (2015).
Improvement and expansion of the Fmask algorithm: cloud, cloud
shadow, and snow detection for Landsats 4-7, 8, and Sentinel 2 images
Remote Sensing of Environment 159 (2015) 269-277.

Installation requires `Python <http://python.org/>`_, `numpy <http://www.numpy.org/>`_, `scipy <http://www.scipy.org/>`_,
`GDAL <http://gdal.org/>`_ and `RIOS <http://rioshome.org/>`_ and the ability to compile C extensions for Python.
It is licensed under GPL 3.

Originally developed by Neil Flood at `DSITI <https://www.qld.gov.au/dsiti/>`_ 
and additional work funded by `Landcare Research <http://www.landcareresearch.co.nz>`_.


Philosophy
----------
This package implements the Fmask algorithm as a Python module. It is intended that this 
can be wrapped in a variety of main programs which can handle the local details of how
the image files are named and organised, and is intended to provide maximum flexibility. It
should not be tied to expecting the imagery to be layed out in a particular manner.

This modular design also simplifies the use of the same core algorithm on either Landsat and 
Sentinel imagery. The wrapper programs take care of differences in file organisation and 
metadata formats, while the core algorithm is the same for both. 

However, we have also supplied some example wrapper scripts, based around the image organisation
as supplied by the usual distributors of the imagery. In the case of Landsat, we have supplied
main programs which can cope with the data as it comes from USGS, while in the case of Sentinel-2
we have supplied wrappers to deal with the data as supplied by ESA. 

It is expected that some users will use these directly, while larger organisations will wish to
create their own wrappers specific to their own file naming and layout conventions. 

The output from the core algorithm module is a single thematic raster, with integer
codes representing null, clear, cloud, shadow, snow, water respectively.

The examples shown below use the given example wrappers. 

Command Line Examples
---------------------

All the commandline programs given use argparse to handle commandline arguments, and hence will 
respond sensibly to the \-\-help option. Some have options to modify their behaviour. 

Please note that the output format used is defined by `RIOS <http://rioshome.org/>`_. This defaults to HFA (.img). 
See `RIOS documentation <http://rioshome.org/rios_imagewriter.html#rios.imagewriter.setDefaultDriver>`_
for more information and how to change this using the environment variable $RIOS_DFLT_DRIVER.

**Note:** these examples are for use in a Unix/Linux shell so that the filename wildcards
get expanded properly. Windows users should prefix these commands with "python PATH\\TO\\fmask_expandWildcards.py" where 
PATH\\TO\\fmask_expandWildcards.py is the output of "where fmask_expandWildcards.py".

USGS Landsat
^^^^^^^^^^^^

The command line scripts supplied can process an untarred USGS Landsat scene. Firstly,
the reflective and thermal bands must be stacked separately. This needs to be done
in a different manner depending on the sensor.

Landsat 4&5::

    gdal_merge.py -separate -of HFA -co COMPRESSED=YES -o ref.img L*_B[1,2,3,4,5,7].TIF
    gdal_merge.py -separate -of HFA -co COMPRESSED=YES -o thermal.img L*_B6.TIF

Landsat 7::

    gdal_merge.py -separate -of HFA -co COMPRESSED=YES -o ref.img L*_B[1,2,3,4,5,7].TIF
    gdal_merge.py -separate -of HFA -co COMPRESSED=YES -o thermal.img L*_B6_VCID_?.TIF

Landsat 8::

    gdal_merge.py -separate -of HFA -co COMPRESSED=YES -o ref.img LC8*_B[1-7,9].TIF
    gdal_merge.py -separate -of HFA -co COMPRESSED=YES -o thermal.img LC8*_B1[0,1].TIF

The next step is to create an image of the relevant angles, per-pixel, for use in subsequent steps.
For Landsat, this can be done using::

    fmask_usgsLandsatMakeAnglesImage.py -m *_MTL.txt -t ref.img -o angles.img

Then mask and Top of Atmosphere reflectance must be calculated and finally the cloud mask itself::

    fmask_usgsLandsatSaturationMask.py -i ref.img -m *_MTL.txt -o saturationmask.img
    fmask_usgsLandsatTOA.py -i ref.img -m *_MTL.txt -z angles.img -o toa.img
    fmask_usgsLandsatStacked.py -t thermal.img -a toa.img -m *_MTL.txt -z angles.img -s saturationmask.img -o cloud.img 

If the thermal band is empty (for Landsat-8 with the SSM anomaly, after 2015-11-01) then it 
is ignored gracefully.

Sentinel2
^^^^^^^^^

The command line scripts supplied can process a Sentinel2 Level C granule from the image directory. 
Here is an example of how to do this. This example works at 20m resolution, but the 
recipe can be varied as required. 

This makes a stack of ALL the bands, at the 20m resolution (a compromise between speed and detail). 
Bands are in order of numeric band number::

    gdalbuildvrt -resolution user -tr 20 20 -separate allbands.vrt S2*_B0[1-8].jp2 S2*_B8A.jp2 S2*_B09.jp2 S2*_B1[0-2].jp2

Make a separate image of the per-pixel sun and satellite angles. ::

    fmask_sentinel2makeAnglesImage.py -i ../S2*.xml -o angles.img
    
Now create the cloud mask output image. 
Note that this assumes the bands are in a particular order (as created in the vrt, above)::

    fmask_sentinel2Stacked.py -a allbands.vrt -z angles.img -o cloud.img


Re-wrapping and Re-configuring
------------------------------
To build a different set of wrappers, and configure things differently, the default
wrappers are a good place to start. The configuration is mainly handled by the 
:class:`fmask.config.FmaskConfig` class. For example, one would 
call :func:`fmask.config.FmaskConfig.setReflectiveBand` to change which layer of the stack 
corresponds to which wavelength band. 

Downloads
---------
Get the source from `BitBucket <https://bitbucket.org/chchrsc/python-fmask/downloads>`_.

`Conda <http://conda.pydata.org/miniconda.html#miniconda>`_ packages are available under the 'rios' channel.
Once you have installed `Conda <http://conda.pydata.org/miniconda.html#miniconda>`_, run the following commands on the command line to install python-fmask: ::

    conda config --add channels conda-forge 
    conda config --add channels rios 
    conda create -n myenv python-fmask 
    source activate myenv # omit 'source' on Windows

Issues
------

Please log bugs encountered with the `Issue Tracker <https://bitbucket.org/chchrsc/python-fmask/issues?status=new&status=open>`_.


Python developer documentation
------------------------------

.. toctree::
    :maxdepth: 1

    Running the fmask algorithm  <fmask_fmask>
    Configuring the fmask run <fmask_config>
    Creating Top of Atmosphere rasters for Landsat <fmask_landsatTOA>
    fmask_saturationcheck
    fmask_zerocheck
    fmask_fillminima
    fmask_valueindexes
    fmask_fmaskerrors

* :ref:`modindex`
* :ref:`search`

.. codeauthor:: Neil Flood & Sam Gillingham
