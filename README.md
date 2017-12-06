# cornerbox?

Names comes from a contraction of corners and box (Obvious, not?). A box is extracted from a smaller rectangular area on a source .tif DEM. A regular gridded rectangular interpolated mesh can then be created over the area.

## File formats
Listed are some important file formats this program works with.

### ESRI Shapefiles
Box boundary can be output in a shapefile to imported into a GIS software such as ArcGIS or QGIS.

### MIKE boundary file
Extract a box by cornerpoint from a .tif raster image + output a .xyz-file of cornerpoints. The generated .xyz-file has the following form

722600 6184300 0.0 1 0 0

722600 6184550 0.0 0 0 0

722600 6184550 0.0 1 0 0

722850 6184550 0.0 0 0 0

722850 6184550 0.0 1 0 0

722850 6184300 0.0 0 0 0

722850 6184300 0.0 1 0 0

722600 6184300 0.0 0 0 0

This creates 4 arcs. One for each side of a rectangle. Such a boundary file allows for creation of meshes within the MIKE Mesh Generator.

### MIKE scatter data file
In order to import a DEM into MIKE Mesh Generator it must have an ASCII type format. The Cornerbox class allows conversion from .tif to .xyz to be used in MIKE together with the .xyz boundary file.

### MIKE mesh file
MIKE 21 FM(Flexible Mesh) simulation software requires the DEM as a flexible mesh file .mesh particular to DHI products. It is usually created using the MIKE Mesh Generator. However, the MeshGenerator class in this project allows creation of regular gridded meshes in the .mesh fileformat. This makes it possible to omit mesh creations within Mesh Generator. A discription of the file format was obtained from <http://manuals.mikepoweredbydhi.help/2017/General/FM_FileSpecification.pdf>.



# TODO:
- give warning
* if resolution too crude for boundaries
* if resolution not divisible by 0.4
* if all crudes resolution not divisible by all other

- make parser

- better readme

- remove unused functionality

- make a test.py in ./tests that tests the combination of cornerbox and meshgenerator

- fix shp_boundary()
