#!/usr/bin/env python
from __future__ import division

import os  # filesystem
import subprocess
from osgeo import ogr, osr

global BIT  # Used to find OSGEO installation
BIT = 32

class Cornerbox:
    corners = []
    tif_s = ''; tif_e = ''; tif_r = ''
    borderx = 20.0; bordery = 20.0

    def __init__(self, tif_s = None):
        self.tif_s = tif_s

    @staticmethod
    def corners(ulx, uly, lx, ly):
        """
        Corner points from upper-left corner and side lengths of box
        :param ulx: Upper-left x-coordinate
        :param uly: Upper-left y-coordinate
        :param lx: length along x (positive: east)
        :param ly: length along y (positive: north)
        :return: List of corner coordinates
        """
        return [ulx, uly, ulx + lx, uly + ly]

    def set_corners(self, ulx, uly, lrx, lry):
        """
        Unnecessary function. Do cornerbox.corners = [...] instead
        :param ulx:
        :param uly:
        :param lrx:
        :param lry:
        :return:
        """
        self.corners = [ulx, uly, lrx, lry]

    @staticmethod
    def get_lengths(corners):
        """
        Side lengths of box from corners
        :param corners: List of corner coordinates
        :return: Pair of box side lengths
        """
        return abs(corners[0] - corners[2]), abs(corners[1] - corners[3])

    @staticmethod
    def num2str(num, type='resolution'):
        """
        Convert a number to a string with zero-padding. Used for informative naming of files
        :param num: Any number (float or integer)
        :param type: 'resolution'(default) or 'length'
        :return: Zero-padded number
        """
        if type == 'resolution':
            round(num, 2)
            num = num * 10
            num = int(num)
            num = '{:03d}'.format(num)
        elif type == 'length':
            num = str(num)

            num_nodot = ''
            for d in num:
                if d.isdigit():
                    num_nodot += d
            num = num_nodot

            num = int(num)
            num = '{:05d}'.format(num)
        else:
            return 1  # throw error!

        return num

    @staticmethod
    def folder_name(corners):
        """
        Folder name for box (saving the output files into)
        :param corners: List of corner coordinates
        :return: A string
        """
        # Folder name to store output file. An extract is referred to by the (ulx, uly) corner point.
        return 'box{}_{}'.format(str(corners[0]), str(corners[1]))  # Folder name "box_(ulx)_(uly)"

    def create_folder_structure(self):
        """
        Generates the used folder structure from a list of corners
        :param self:
        :return: None
        """
        foldername = './{}/{}'.format('extract', self.folder_name(self.corners))

        if not os.path.exists(foldername):
            os.mkdir(foldername)
            os.mkdir(foldername + '/e')  # Original (resolution equal to the original raster)
            os.mkdir(foldername + '/r')  # Resample (resolution different... )
            os.mkdir(foldername + '/s')  # Shape (containing shapefiles marking extent of extract)
            os.mkdir(foldername + '/m')  # Mesh (MIKE compatible mesh files)
            print 'created folders: {}(/e, /r, /s, /m)'.format(foldername)

    def box_xyz_boundary(self):
        """
        MIKE compatible .xyz boundary file marking extent of extract
        :param self:
        :return:
        """
        length = self.get_lengths(self.corners)
        fout = './extract/{}/s/box_{}x{}.xyz'.format(self.folder_name(self.corners),
                                                     self.num2str(length[0], 'length'),
                                                     self.num2str(length[1], 'length'))
        fopen = open(fout, 'w')

        z = '0.0'

        # Make MIKE arcs. 1 - nodes, 0 - vertices
        line = ''.join([str(self.corners[0]), ' ', str(self.corners[1]), ' ', z, ' 1 0 0\n'])  # (ulx, uly)
        fopen.write(line)
        line = ''.join([str(self.corners[0]), ' ', str(self.corners[3]), ' ', z, ' 0 0 0\n'])  # (ulx, lry)
        fopen.write(line)
        line = ''.join([str(self.corners[0]), ' ', str(self.corners[3]), ' ', z, ' 1 0 0\n'])  # (ulx, lry)
                                                                                               # End and beginning
        fopen.write(line)
        line = ''.join([str(self.corners[2]), ' ', str(self.corners[3]), ' ', z, ' 0 0 0\n'])  # (lrx, lry)
        fopen.write(line)
        line = ''.join([str(self.corners[2]), ' ', str(self.corners[3]), ' ', z, ' 1 0 0\n'])  # (lrx, lry)
        fopen.write(line)
        line = ''.join([str(self.corners[2]), ' ', str(self.corners[1]), ' ', z, ' 0 0 0\n'])  # (lrx, uly)
        fopen.write(line)
        line = ''.join([str(self.corners[2]), ' ', str(self.corners[1]), ' ', z, ' 1 0 0\n'])  # (lrx, uly)
        fopen.write(line)
        line = ''.join([str(self.corners[0]), ' ', str(self.corners[1]), ' ', z, ' 0 0 0\n'])  # (ulx, uly)
                                                                                               # Closes loop
        fopen.write(line)

        fopen.close()

    def box_shp_boundary(self):
        """
        Shapefile .shp marking extent of extract
        :param self:
        :return:
        """
        # set up the shapefile driver
        shp_driver = ogr.GetDriverByName("ESRI Shapefile")
        # create the data source
        length = self.get_lengths(self.corners)
        filename = './extract/{}/s/box_{}x{}.shp'.format(self.folder_name(self.corners),
                                                         self.num2str(length[0], 'length'),
                                                         self.num2str(length[1], 'length'))
        shp_ds = shp_driver.CreateDataSource(filename)

        # create the spatial reference, WGS84
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4647)

        # create and write the layer
        shp_lyr = shp_ds.CreateLayer('box', srs, ogr.wkbPolygon)
        print self.corners

        wkt = 'POLYGON (({} {}, {} {}, {} {}, {} {}, {} {}))'.format(
            self.corners[0], self.corners[1],
            self.corners[0], self.corners[3],
            self.corners[2], self.corners[3],
            self.corners[2], self.corners[1],
            self.corners[0], self.corners[1]  # Closes box (polygon)
        )

        geom = ogr.CreateGeometryFromWkt(wkt)

        feature = ogr.Feature(shp_lyr.GetLayerDefn())
        feature.SetGeometry(geom)
        shp_lyr.CreateFeature(feature)
        feature = None

        shp_ds = None  # Save and close the data source

    def extract_tif(self):
        """
        Extract a box from a larger raster file
        :param self:
        :return:
        """
        # define path to the used tool (installed through OSGEO4W)
        tool = 'gdal_translate.exe'
        if BIT == 64:
            path_tool = 'C:/OSGeo4W64/bin/%s' % tool
        else:
            path_tool = 'C:/OSGeo4W/bin/%s' % tool

        length = self.get_lengths(self.corners)

        print 'side lengths: {}'.format(length)
        print 'extracting square with corner points: (ulx, uly, lrx, lry) = ({},{},{},{}) \nand borders: ({},{}'.format(
                                                                                                 self.corners[0],
                                                                                                 self.corners[1],
                                                                                                 self.corners[2],
                                                                                                 self.corners[3],
                                                                                                 self.borderx,
                                                                                                 self.bordery)

        tif_e = './extract/{}/e/box_{}x{}_r004x004.tif'.format(self.folder_name(self.corners),
                                                               self.num2str(length[0], 'length'),
                                                               self.num2str(length[1], 'length'))  # hardcoded resolution!
                                                                                                   # new source (boundaries)
        try:
            subprocess.call('{} -projwin {} {} {} {} {} {}'.format(path_tool,
                                                                   self.corners[0] - self.borderx,
                                                                   self.corners[1] + self.bordery,
                                                                   self.corners[2] + self.borderx,
                                                                   self.corners[3] - self.bordery,
                                                                   self.tif_s,
                                                                   tif_e))
        except WindowsError:
            print 'Could not find tool! ({})'.format(path_tool)

        self.tif_e = tif_e

        tif_r = './extract/{}/r/box_{}x{}_r004x004.tif'.format(self.folder_name(self.corners),
                                                               self.num2str(length[0], 'length'),
                                                               self.num2str(length[1], 'length'))  # hardcoded resolution!

        try:
            subprocess.call('{} -projwin {} {} {} {} {} {}'.format(path_tool,
                                                                   self.corners[0],
                                                                   self.corners[1],
                                                                   self.corners[2],
                                                                   self.corners[3],
                                                                   tif_e,
                                                                   tif_r))
        except WindowsError:
            print 'Could not find tool! ({})'.format(path_tool)


    def resample(self, resolution_x, resolution_y):
        """
        Resample a tif at different resolution
        :param self:
        :param resolution_x:
        :param resolution_y:
        :return:
        """
        tool = 'gdalwarp.exe'
        if BIT == 64:
            path_tool = 'C:/OSGeo4W64/bin/%s' % tool
        else:
            path_tool = 'C:/OSGeo4W/bin/%s' % tool

        length = self.get_lengths(self.corners)

        tif_r = './extract/{}/r/box_{}x{}_r{}x{}.tif'.format(self.folder_name(self.corners),
                                                                   self.num2str(length[0], 'length'),
                                                                   self.num2str(length[1], 'length'),
                                                                   self.num2str(resolution_x, 'resolution'),
                                                                   self.num2str(resolution_y, 'resolution'))

        self.tif_r = tif_r

        tif = './extract/{}/r/box_{}x{}_r004x004.tif'.format(self.folder_name(self.corners),
                                                             self.num2str(length[0], 'length'),
                                                             self.num2str(length[1], 'length'))

        try:
            subprocess.call('%s -tr %s %s %s %s -overwrite' % (path_tool, resolution_x, resolution_y, tif, tif_r))
        except WindowsError:
            print 'Could not find tool! ({})'.format(path_tool)

    def box_tif2asc(self):
        """
        Convert the resampled .tif box to .asc compatible with MIKE grid conversion to .dfs2
        :param tif: Raster file (.tif) to convert to XYZ-file
        :return: None
        """
        tool = 'gdal_translate.exe'

        if BIT == 64:
            path_tool = 'C:/OSGeo4W64/bin/%s' % tool
        else:  # BIT == 32
            path_tool = 'C:/OSGeo4W/bin/%s' % tool

        basename = '.' + self.tif_r.split('.')[1]  # Remove extension

        try:
            subprocess.call('%s -a_srs EPSG:4647 -of AAIGrid %s.tif %s.asc' % (path_tool, basename, basename))
        except WindowsError:
            print 'Could not find tool! ({})'.format(path_tool)

    def box_tif2xyz(self):
        """
        Convert the resampled .tif box to .xyz
        :param tif: Raster file (.tif) to convert to XYZ-file
        :return: None
        """
        tool = 'gdal_translate.exe'

        if BIT == 64:
            path_tool = 'C:/OSGeo4W64/bin/%s' % tool
        else:  # BIT == 32
            path_tool = 'C:/OSGeo4W/bin/%s' % tool

        basename = '.' + self.tif_r.split('.')[1]  # Remove extension

        try:
            subprocess.call('%s -of XYZ %s.tif %s.xyz' % (path_tool, basename, basename))
        except WindowsError:
            print 'Could not find tool! ({})'.format(path_tool)

    def get_filename_mesh(self, resolution_x, resolution_y):
        """
        Not used!
        :param resolution_x:
        :param resolution_y:
        :return:
        """
        length = self.get_lengths(self.corners)

        return './extract/{}/m/box_{}x{}_r{}x{}.mesh'.format(self.folder_name(self.corners),
                                                                   self.num2str(length[0], 'length'),
                                                                   self.num2str(length[1], 'length'),
                                                                   self.num2str(resolution_x, 'resolution'),
                                                                   self.num2str(resolution_y, 'resolution'))

if __name__=='__main__':
    corners = Cornerbox.corners(722600, 6184300, 281.6, -281.6) # For small 0.4x0.4 m resolution area

    # .tif file to extract from:
    tif = "../01_DTM/DHYMRAIN.tif"

    box = Cornerbox(tif)
    box.corners = corners
    box.create_folder_structure()
    box.extract_tif()
    #box.resample(93.2, 93.2)
    box.resample(3.2, 3.2)
    #box.extract_tif()
    #for resolution in [0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2]:
    #    box.resample(resolution, resolution)
    box.box_shp_boundary()
    box.box_tif2xyz()
    box.box_tif2asc()
    box.box_xyz_boundary()



