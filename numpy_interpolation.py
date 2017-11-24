#!/usr/bin/env python


from osgeo import gdal
from cornerbox import *
import numpy as np

def read_raster(cornerlist):
    length = get_lengths(cornerlist)
    tif = './extract/{}/o/box_{}x{}_r04x04.tif'.format(folder_name(cornerlist),
                                                       num2str(length[0], 'length'),
                                                       num2str(length[1], 'length'))

    ds = gdal.Open(tif, gdal.GA_ReadOnly)
    dem = np.array(ds.GetRasterBand(1).ReadAsArray())  # only one band - DEM
    return dem

def create_grid(cornerlist, resolution_x = 0.4, resolution_y = 0.4):
    (length_x, length_y) = get_lengths(cornerlist)

    grid = np.zeros([length_x // resolution_x, length_y // resolution_y, 3])

    for i in xrange(int(length_x // resolution_x)):
        for j in xrange(int(length_y // resolution_y)):
            grid[i, j, 0] = cornerlist[0] + i * resolution_x
            grid[i, j, 1] = cornerlist[1] + j * resolution_y
            print grid[i, j]

    return grid

def interpolate_grid(cornerlist, tif, resolution_x = 0.4, resolution_y = 0.4):
    (length_x, length_y) = get_lengths(cornerlist)

    grid = np.zeros([length_x / resolution_x, length_y / resolution_y, 3])

    for i in xrange(int(length_x / resolution_x)):
        for j in xrange(int(length_y / resolution_y)):
            grid[i, j, 0] = cornerlist[0] + i * resolution_x
            grid[i, j, 1] = cornerlist[1] + j * resolution_y
            grid[i, j, 2] = tif[i, j]

    return grid

if __name__=='__main__':
    cornerlist = corners(722600, 6184300, 280, -280)  # For small 0.4x0.4 m resolution area
    dem = read_raster(cornerlist)
    interpolate_grid(cornerlist, dem, 3.2, 3.2)