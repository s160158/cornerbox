# What does it do?

Extract a box from a raster file (.tif) and output boundary XYZ-file for mesh generation in MIKE. Compatible XYZ-file scatter data are also generated for mesh interpolation.

## MIKE boundary file

Extract a box by cornerpoint from a .tif raster image + output a .xyz-file of cornerpoints. The generated .xyz-file has the following form

722600 6184300 0.0 1 0 0

722600 6184550 0.0 0 0 0

722600 6184550 0.0 1 0 0

722850 6184550 0.0 0 0 0

722850 6184550 0.0 1 0 0

722850 6184300 0.0 0 0 0

722850 6184300 0.0 1 0 0

722600 6184300 0.0 0 0 0

This creates 4 arcs. One for each side of a rectangle.

## MIKE scatter data file


# Todo
restate purpose.. more singular

NO RESAMPLING.. its done in MIKE by interpolation


* Connectivity for nodes/vertices. What value
* Two ways: (1) separation into callable functions (2) main() function that generates all the needed files, i.e., .tif (full + reduced), .xyz (of reduced) and a .xyz of arcs connecting cornerpoints to be imported in MIKE
* Proper python notation
* Make a json(flat text, possibly without json) config file (tool paths, etc.) 
* GDAL instead of precompiled function