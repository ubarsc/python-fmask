Release Notes
=============

Version 0.5.9 (2024-03-13)
--------------------------

Bug Fixes:
  * Convert from bool to uint8 before writing a mask as latest RIOS no 
    longer supports writing bool arrays
    
Enhancements:
  * Add support for entrypoints using "pip install". These entrypoints 
    cannot have an extension so to reduce confusion the documentation
    has been updated to cover using these new entrypoints (no '.py').
    The conda package will be similarly updated. The old scripts remain 
    as they are, although they now emit a warning.

Version 0.5.8 (2022-12-22)
--------------------------

Bug Fixes
  * Cope with numpy-1.24 removal of deprecated type symbols like numpy.bool
  * Better handling of gdal exception use/non-use

Enhancements
  * Use gdal.Driver.Delete to remove temporary files in Sentinel-2 cmdline
    script
  * Use gdal internal function calls to avoid use of os.system()


Version 0.5.7 (2022-02-11)
--------------------------

Enhancements
    * Minor changes to support Landsat-9 from USGS

Version 0.5.6 (2021-10-19)
--------------------------

Enhancements
    * Cope with ESA's sudden inclusion of radiometric offsets in their
      Sentinel-2 reflectance imagery. Using earlier python-fmask versions 
      with the new ESA files will result in incorrect answers. 

Version 0.5.5 (2020-09-01)
--------------------------

Bug Fixes
    * Cope with the axis swapping behavour of GDAL 3.x. 
    * Display proper error messages if missing gdalwarp or gdal_merge.py.
    
Enhancements
    * Improve readability of code in masking step.
    * Allow sentinel2Stacked to be called directly from Python.


Version 0.5.4 (2019-08-06)
--------------------------

Bug Fixes
    * Removed the "darkness" test, suggested by Franz (2015), which had inadvertently been
      added from a test branch into version 0.5.0. This test aims to remove a certain 
      type of false positive cloud, but larger tests suggest that it removes more true 
      positives, and so should not have been included. 
    * Added new logic to cope with cases of nearly 100% cloud cover. In cases where only 
      very small amounts of land are visible, the dynamic threshold for the "clear land
      probability" is too contaminated to be usable. If less than 3% of the image is 
      clear land, then a fallback threshold is used instead. This avoids images which
      are entirely covered in cloud from being classed as almost entirely cloud-free. 
      This particularly affected Sentinel-2, where the thermal is not available to 
      catch these cases anyway. 

Documentation
    * Added a disclaimer to the front page, emphasizing that all errors are ours, and 
      the authors of the original papers bear no responsibility. 
    * Added a note to the front page that we have done some testing of the proposed
      Fmask4 changes from Qiu et al, 2019, and are unsure whether they help or not, 
      so have not implemented them. 

Version 0.5.3 (2019-01-15)
--------------------------

Bug Fixes
  * Fixed problem with new error checks which broke Landsat-7 case

Version 0.5.2 (2018-12-12)
--------------------------

Bug Fixes
  * Fix issue with entry points on Windows

Enhancements
  * Ensure temporary files are removed on Windows

Version 0.5.1 (2018-11-26)
--------------------------

Enhancements
  * Added better support for Conda packaging on Windows
  * Upgraded license to GPL v3

Version 0.5.0 (2018-10-18)
--------------------------

Enhancements
  * For Sentinel-2, added support for parallax test to remove false cloud in bright (typically
    urban) areas, as per Frantz (2018). Optional, defaults to off. Thanks to Vincent Schut
    for assistance in kicking that across the line. 
  * For Sentinel-2, added command line switch to directly use the SAFE directory as 
    provided by ESA, to avoid the need to manually stack up the individual bands externally. 
  * For Sentinel-2, turn off the (cloud prob > 99%) test, as suggested by Zhu (2015). It is 
    in effect for Landsat. 


Version 0.4.5 (2017-07-12)
--------------------------

Bug Fixes:
  * Handling old formats of USGS MTL files
  * Fixes for numpy 1.13


Version 0.4.4 (2017-04-04)
--------------------------

Enhancements
  * Added commandline options to set cloud probability threshold, and the two reflectance 
    thresholds used for the snow mask. 
  * Call gdal.UseExceptions(), so that when bad things happen, they are more likely
    to be reported accurately by the exception message, and reduce confusion for users. 


Version 0.4.3 (2016-11-26)
--------------------------

Bug Fixes:
  * Fix 32 bit builds
  * Fix help message for fmask_usgsLandsatStacked.py

Enhancements
  * Helper .bat file for Windows to expand wildcards
  * Changes to 'nodata' handling to make processing in parallel possible with RIOS


Version 0.4.2 (2016-09-01)
--------------------------

Bug Fixes
  * Fixed fall-back default values for Landsat brightness temperature equation constants, 
    as required when processing older USGS files which do not have these present in the MTL file. 
  * For Sentinel-2 only, added a work-around for the alarming random null pixels which
    ESA leave in the cirrus band. This avoids leaving corresponding null pixels in the 
    resulting output masks. 


Version 0.4.0 (2016-06-10)
--------------------------

Bug fixes
  * Proper null mask taken from all reflective bands combined, not just the blue band
  * Trap seg-faults in valueindexes C code
  * Use null value of 32767 for Landsat TOA image
  * Cope when Sentinel-2 metadata only has sensor angles for a subset of bands. 

Enhancements
  * Landsat angles code is now in a module, with a main program wrapper, consistent 
    with the rest of the package
  * Added :command:`--cloudbufferdistance`, :command:`--shadowbufferdistance` and 
    :command:`--mincloudsize` options to
    main program wrappers (both Landsat and Sentinel-2) to give user control over these
    parameters


Version 0.3.0 (2016-03-21)
--------------------------

Bug fixes
  * Added code for estimating per-pixel Landsat sun and sensor angles, to allow proper
    shadow tracking, as per original code
  * Full use of Sentinel-2 metadata XML, including per-pixel angles grid

