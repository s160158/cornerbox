#!/usr/bin/env python


import os  # filesystem
import subprocess
import math  # sqrt,
from osgeo import ogr, osr


def create_folder_structure(cornerlist):
    """
    Generates the used folder structure
    :param cornerlist:
    :return: None
    """
    foldername = './{}/{}'.format('extract', folder_name(cornerlist))

    if not os.path.exists(foldername):
        print 'created folder: {}'.format(foldername)
        os.mkdir(foldername)
        os.mkdir(foldername + '/o')  # Original (resolution equal to the original raster)
        os.mkdir(foldername + '/r')  # Resample (resolution different... )
        os.mkdir(foldername + '/s')  # Shape (containing shapefiles marking extent of extract)
        os.mkdir(foldername + '/m')  # Mesh (MIKE compatible mesh files)


def corners(ulx, uly, lx, ly):
    """
    Corner points
    :param ulx: Upper-left x-coordinate
    :param uly: Upper-left y-coordinate
    :param lx: length along x (positive: east)
    :param ly: length along y (positive: north)
    :return: List of corner coordinates
    """
    return [ulx, uly, ulx + lx, uly + ly]


def box_xyz(cornerlist):
    """
    MIKE compatible .xyz shapefile marking extent of extract
    :param cornerlist: List of corner coordinates
    :return: None
    """
    length = get_lengths(cornerlist)
    fout = './extract/{}/s/box_{}x{}.xyz'.format(folder_name(cornerlist), num2str(length[0], 'length'), \
                                                 num2str(length[1], 'length'))
    fopen = open(fout, 'w')

    z = '0.0'

    # Make MIKE arcs. 1 - nodes, 0 - vertices
    line = ''.join([str(cornerlist[0]), ' ', str(cornerlist[1]), ' ', z, ' 1 0 0\n'])  # (ulx, uly)
    fopen.write(line)
    line = ''.join([str(cornerlist[0]), ' ', str(cornerlist[3]), ' ', z, ' 0 0 0\n'])  # (ulx, lry)
    fopen.write(line)
    line = ''.join([str(cornerlist[0]), ' ', str(cornerlist[3]), ' ', z, ' 1 0 0\n'])  # (ulx, lry) # End and beginning
    fopen.write(line)
    line = ''.join([str(cornerlist[2]), ' ', str(cornerlist[3]), ' ', z, ' 0 0 0\n'])  # (lrx, lry)
    fopen.write(line)
    line = ''.join([str(cornerlist[2]), ' ', str(cornerlist[3]), ' ', z, ' 1 0 0\n'])  # (lrx, lry)
    fopen.write(line)
    line = ''.join([str(cornerlist[2]), ' ', str(cornerlist[1]), ' ', z, ' 0 0 0\n'])  # (lrx, uly)
    fopen.write(line)
    line = ''.join([str(cornerlist[2]), ' ', str(cornerlist[1]), ' ', z, ' 1 0 0\n'])  # (lrx, uly)
    fopen.write(line)
    line = ''.join([str(cornerlist[0]), ' ', str(cornerlist[1]), ' ', z, ' 0 0 0\n'])  # (ulx, uly) # Closes loop
    fopen.write(line)

    fopen.close()


def box_shp(cornerlist):
    # set up the shapefile driver
    shp_driver = ogr.GetDriverByName("ESRI Shapefile")
    # create the data source
    length = get_lengths(cornerlist)
    filename = './extract/{}/s/box_{}x{}.shp'.format(folder_name(cornerlist),
                                                     num2str(length[0], 'length'),
                                                     num2str(length[1], 'length'))
    shp_ds = shp_driver.CreateDataSource(filename)

    # create the spatial reference, WGS84
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4647)

    # create and write the layer
    shp_lyr = shp_ds.CreateLayer('box', srs, ogr.wkbPolygon)
    print cornerlist

    feature = ogr.Feature(shp_lyr.GetLayerDefn())

    wkt = 'POLYGON (({} {}, {} {}, {} {}, {} {}, {} {}))'.format(
        cornerlist[0], cornerlist[1],
        cornerlist[0], cornerlist[3],
        cornerlist[2], cornerlist[3],
        cornerlist[2], cornerlist[1],
        cornerlist[0], cornerlist[1]  # Closes box (polygon)
    )

    box = ogr.CreateGeometryFromWkt(wkt)
    feature.SetGeometry(box)
    shp_lyr.CreateFeature(feature)

    shp_ds = None  # Save and close the data source



def get_lengths(cornerlist):
    """
    Side length of box
    :param cornerlist: List of corner coordinates
    :return: Pair of box side lengths
    """
    return abs(cornerlist[0] - cornerlist[2]), abs(cornerlist[1] - cornerlist[3])


def num2str(num, type='resolution'):
    """
    Convert a number to a string with zero-padding. Used for informative naming of files
    :param num: Any number (float or integer)
    :param type: 'resolution'(default) or 'length'
    :return: Zero-padded number
    """
    #
    if type == 'resolution':
        round(num, 2)
        num = num * 10
        num = int(num)
        num = '{:03d}'.format(num)
    elif type == 'length':
        num = int(num)
        num = '{:04d}'.format(num)
    else:
        return 1  # throw error!

    return num


def folder_name(cornerlist):
    """
    Folder name for box (saving the output files into)
    :param cornerlist: List of corner coordinates
    :return: A string
    """
    # Folder name to store output file. An extract is referred to by the (ulx, uly) corner point.
    return 'box{}_{}'.format(str(cornerlist[0]), str(cornerlist[1]))  # Folder name "box_(ulx)_(uly)"


def box_tif(tif, cornerlist):
    """
    Extract a box from a larger raster file
    :param tif: Large raster file to extract box from
    :param cornerlist: List of corner coordinates
    :return: None (Generates the files)
    """
    # define path to the used tool (installed through OSGEO4W)
    tool = 'gdal_translate.exe'
    path_tool = 'C:/OSGeo4W64/bin/%s' % tool

    length = get_lengths(cornerlist)

    print "side lengths: " + str(length)

    print 'extracting square with corner points: (ulx, uly, lrx, lry) = (%d, %d, %d, %d)' % (cornerlist[0], cornerlist[1], cornerlist[2], cornerlist[3])
    print tif


    tif_extract = './extract/{}/o/box_{}x{}_r04x04.tif'.format(folder_name(cornerlist), num2str(length[0], 'length'), num2str(length[1], 'length')) #find the resolution! not always 40cm
    print 'saving to :' + tif_extract

    try:
        subprocess.call('%s -projwin %s %s %s %s %s %s' % (path_tool, cornerlist[0], cornerlist[1], cornerlist[2], cornerlist[3], tif, tif_extract))
    except WindowsError:
        print 'Could not find tool! ({})'.format(path_tool)
        return 1


def resample(resolution_x, resolution_y, cornerlist):
    # Resample a tif at different resolution
    tool = 'gdalwarp.exe'
    path_tool = 'C:/OSGeo4W64/bin/%s' % tool

    length = get_lengths(cornerlist)

    tif = './extract/{}/o/box_{}x{}_r04x04.tif'.format(folder_name(cornerlist),
                                                       num2str(length[0], 'length'),
                                                       num2str(length[1], 'length'))


    tif_extract = './extract/{}/r/box_{}x{}_r{}x{}.tif'.format(folder_name(cornerlist),
                                                               num2str(length[0], 'length'),
                                                               num2str(length[1], 'length'),
                                                               num2str(resolution_x, 'resolution'),
                                                               num2str(resolution_y, 'resolution'))

    try:
        subprocess.call('%s -tr %s %s %s %s -overwrite' % (path_tool, resolution_x, resolution_y, tif, tif_extract))
    except WindowsError:
        print 'Could not find tool! ({})'.format(path_tool)
        return 1

def tif2xyz(tif):
    """
    Convert .tif to .xyz
    :param tif: Raster file (.tif) to convert to XYZ-file
    :return: None
    """
    tool = 'gdal_translate.exe'
    path_tool = 'C:/OSGeo4W64/bin/%s' % tool

    tif = '.' + tif.split('.')[1]  # Remove extension
    try:
        subprocess.call('%s -of XYZ %s.tif %s.xyz' % (path_tool, tif, tif))
    except WindowsError:
        print 'Could not find tool! ({})'.format(path_tool)
        return 1


if __name__=='__main__':
    cornerlist = corners(722600, 6184300, 280, -280) # For small 0.4x0.4 m resolution area
    #cornerlist = corners(722600, 6184300, 560, -560) # For larger resolution than 0.8x0.8 m
    create_folder_structure(cornerlist)

    box_xyz(cornerlist)

    # .tif file to extract from:
    tif = "../01_DTM/DHYMRAIN.tif"
    #box_tif(tif, cornerlist)
    resample(3.2, 3.2, cornerlist)
    #print rectangle_lengths(500000, 2, 0.8, 0.4)
    box_shp(cornerlist)
