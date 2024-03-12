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

Also includes optional extension for Sentinel-2 from
Frantz, D., Hass, E., Uhl, A., Stoffels, J., & Hill, J. (2018). 
Improvement of the Fmask algorithm for Sentinel-2 images: Separating clouds 
from bright surfaces based on parallax effects. 
Remote Sensing of Environment 215, 471-481.

Installation requires `Python <http://python.org/>`_, `numpy <http://www.numpy.org/>`_, `scipy <http://www.scipy.org/>`_,
`GDAL <http://gdal.org/>`_ and `RIOS <http://rioshome.org/>`_ and the ability to compile C extensions for Python.
It is licensed under GPL 3.

Originally developed by Neil Flood at `DSITI <https://www.qld.gov.au/dsiti/>`_
and additional work funded by `Landcare Research <http://www.landcareresearch.co.nz>`_.

Later Fmask4 Methods
--------------------
A later publication by Qiu et al, 2019, suggests a number of further changes to Fmask, and
gives this newer version the name of Fmask4. We have done some testing of these suggested 
methods, using Landsat and Sentinel-2 data over large areas of Australia. We found that 
overall, the losses outweighed the benefits, and so are unsure whether to 
implement these changes. Currently we have chosen not to. 

Qiu, S., Zhu, Z., & He, B. (2019). Fmask 4.0: Improved cloud and cloud shadow 
detection in Landsats 4-8 and Sentinel-2 imagery. Remote Sensing of Environment, 231, 111205.



Disclaimer and Acknowledgement
------------------------------
This Python implementation has followed the work of the authors cited above, and we
offer our thanks for their work. None of them were involved in the creation of
this Python implementation, and all errors made are our own responsibility. 


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
respond sensibly to the -h option by printing their own help.
Some have options to modify their behaviour.

Please note that the output format used is defined by `RIOS <http://rioshome.org/>`_. This defaults to HFA (.img).
See `RIOS documentation <http://rioshome.org/rios_imagewriter.html#rios.imagewriter.setDefaultDriver>`_
for more information and how to change this using the environment variable $RIOS_DFLT_DRIVER.

USGS Landsat
^^^^^^^^^^^^

**Update:** Since the USGS released their Collection-1 data set (completed globally in 2017),
they now distribute cloud, shadow
and snow masks included in their QA layer. These are calculated using CFMask, which should,
in principle, be equivalent to the results of this code. Therefore, when processing USGS
Collection-1 data, users may prefer the USGS-supplied masks.
See `USGS QA Layer <https://landsat.usgs.gov/collectionqualityband>`_.

The command line scripts supplied can process an untarred USGS Landsat scene. 
Here is an example of how to do this. This command will take a given scene directory, 
find the right images, and create an output file called cloud.img::

    fmask_usgsLandsatStacked -o cloud.img --scenedir LC08_L1TP_150033_20150413_20170410_01_T1

If the thermal band is empty (for Landsat-8 with the SSM anomaly, after 2015-11-01) then it
is ignored gracefully.

There are command line options to modify many aspects of the algorithm's behaviour. 

There are four options which are now obsolete, for manually specifying a pre-stacked file
of reflectance bands, thermal bands, a saturation mask and the image of angles. 
These should be considered obsolete, and are 
replaced with the --scenedir option, which takes care of all internally. 

Sentinel2
^^^^^^^^^

The command line scripts supplied can process a Sentinel2 Level C granule from the image directory.
Here is an example of how to do this. This example works at 20m resolution, but the
recipe can be varied as required. Be warned, processing at 10m resolution would be considerably
slower, and is unlikely to be any more accurate.

This command will take a given .SAFE directory, find the right images, and create an
output file called cloud.img::

    fmask_sentinel2Stacked -o cloud.img --safedir S2B_MSIL1C_20180918T235239_N0206_R130_T56JNQ_20180919T011001.SAFE

When working with the old ESA zipfile format, which packed multiple tiles into a single SAFE-format
zipfile, this approach will not work, as it won't know which tile to process. So, instead, use
the option to specify the granule directory, as follows::

    fmask_sentinel2Stacked -o cloud.img --granuledir S2A_OPER_PRD_MSIL1C_PDMC_20160111T072442_R030_V20160111T000425_20160111T000425.SAFE/GRANULE/S2A_OPER_MSI_L1C_TL_SGS__20160111T051031_A002887_T56JNQ_N02.01

This would also work on a new-format directory, but specifying the top .SAFE directory is easier. 

There are command line options to modify many aspects of the algorithm's behaviour. 

There are two options which are now obsolete, for manually specifying a pre-stacked file
of reflectance bands, and the image of angles. These should be considered obsolete, and are 
replaced with the --safedir or --granuledir option, which take care of all internally. 

Re-wrapping and Re-configuring
------------------------------
To build a different set of wrappers, and configure things differently, the default
wrappers are a good place to start. The configuration is mainly handled by the
:class:`fmask.config.FmaskConfig` class. For example, one would
call :func:`fmask.config.FmaskConfig.setReflectiveBand` to change which layer of the stack
corresponds to which wavelength band.

Downloads
---------
Get the source as a bundle from `GitHub <https://github.com/ubarsc/python-fmask/releases>`_.
Release notes for each version can be read in :doc:`releasenotes`. To install from source,
read the INSTALL.txt file included inside the source bundle.

Pre-built binary `Conda <https://github.com/conda-forge/miniforge>`_ packages are available
under the 'conda-forge' channel. Once you have installed
`Conda <https://github.com/conda-forge/miniforge>`_, run the following commands on the
command line to install python-fmask: ::

    conda create -n myenv python-fmask
    conda activate myenv

Applications that use python-fmask
----------------------------------

* `Cloud Masking <https://smbyc.bitbucket.io/qgisplugins/cloudmasking/>`_: It  is a Qgis plugin for cloud masking the Landsat (4, 5, 7 and 8) products using different process and filters such as Fmask, Blue Band, Cloud QA, Aerosol and Pixel QA.

Issues
------

Please log bugs encountered with the `Issue Tracker <https://github.com/ubarsc/python-fmask/issues>`_.


Python developer documentation
------------------------------

.. toctree::
    :maxdepth: 1

    Running the fmask algorithm  <fmask_fmask>
    Configuring the fmask run <fmask_config>
    Creating Top of Atmosphere rasters for Landsat <fmask_landsatTOA>
    fmask_saturationcheck
    fmask_landsatangles
    fmask_zerocheck
    fmask_fillminima
    fmask_valueindexes
    fmask_fmaskerrors

* :ref:`modindex`
* :ref:`search`

.. codeauthor:: Neil Flood & Sam Gillingham
