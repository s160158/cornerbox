# What does it do?

Generates two types of data used for testing meshes in MIKE

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

* Connectivity for nodes/vertices. What value
* Two ways: (1) separation into callable functions (2) main() function that generates all the needed files, i.e., .tif (full + reduced), .xyz (of reduced) and a .xyz of arcs connecting cornerpoints to be imported in MIKE
* Proper python notation