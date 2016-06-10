Release Notes
=============

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
  * Added --cloudbufferdistance, --shadowbufferdistance and --mincloudsize options to
    main program wrappers (both Landsat and Sentinel-2) to give user control over these
    parameters


Version 0.3.0 (2016-03-21)
--------------------------

Bug fixes
  * Added code for estimating per-pixel Landsat sun and sensor angles, to allow proper
    shadow tracking, as per original code
  * Full use of Sentinel-2 metadata XML, including per-pixel angles grid

