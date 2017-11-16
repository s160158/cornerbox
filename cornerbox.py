#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os #filesystem
import re #regular expression
import gdal
import subprocess
import math # (sqrt, )


def corners(ulx, uly, lx, ly):
    return [ulx, uly, ulx + lx, uly + ly]

def box_xyz(cornerlist):
    fopen = open('./xyz/test.xyz', 'a')

    z = '0.0'

    line = ''.join([str(cornerlist[0]), ' ', str(cornerlist[1]), ' ', z, ' 1 0 0\n']) # first point 1 (node)
    fopen.write(line)
    line = ''.join([str(cornerlist[0]), ' ', str(cornerlist[3]), ' ', z, ' 0 0 0\n']) # following 0 points are vertices in the arc
    fopen.write(line)
    line = ''.join([str(cornerlist[0]), ' ', str(cornerlist[3]), ' ', z, ' 1 0 0\n']) # new arc (point 1)
    fopen.write(line)
    line = ''.join([str(cornerlist[2]), ' ', str(cornerlist[3]), ' ', z, ' 0 0 0\n'])
    fopen.write(line)
    line = ''.join([str(cornerlist[2]), ' ', str(cornerlist[3]), ' ', z, ' 1 0 0\n'])
    fopen.write(line)
    line = ''.join([str(cornerlist[2]), ' ', str(cornerlist[1]), ' ', z, ' 0 0 0\n'])
    fopen.write(line)
    line = ''.join([str(cornerlist[2]), ' ', str(cornerlist[1]), ' ', z, ' 1 0 0\n'])
    fopen.write(line)
    line = ''.join([str(cornerlist[0]), ' ', str(cornerlist[1]), ' ', z, ' 0 0 0\n'])
    fopen.write(line)

    fopen.close()


def square_length(npixels, res):
    return int(math.sqrt(npixels) * res)

def num2str(num, type='resolution'):
    '''
    convert a number to a string with zero-padding. Used for informative naming of files
    :param num: a float
    :return: zero-padded multiplied by 10
    '''
    if (type == 'resolution'):
        round(num, 2)
        num = num * 10
        num = int(num) # convert to integer
        num = '{:03d}'.format(num)
    elif (type == 'length'):
        num = int(num)
        num = '{:04d}'.format(num)

    return num


def box_tif(tif, tif_extract, cornerlist, npixels, res): #tif_extract is a basename
    # define path to the used tool (installed through OSGEO4W)
    tool = 'gdal_translate.exe'
    path_tool = 'C:/OSGeo4W64/bin/%s' % tool

    l = square_length(npixels, res) #square side length

    #tif_out = './t/resolution_full/test_%sx%s.tif' % (num2str(l, 'length'), num2str(l, 'length'))
    print 'extracting square with corner points: (ulx, uly, lrx, lry) = (%d, %d, %d, %d)' % (cornerlist[0], cornerlist[1], cornerlist[2], cornerlist[3])
    print tif
    tif_extract = tif_extract.join(['_%sx%s_r04x04.tif' % (num2str(l, 'length'), num2str(l, 'length'))])
    print tif_extract

    try:
        subprocess.call('%s -projwin %s %s %s %s %s %s' % (path_tool, cornerlist[0], cornerlist[1], cornerlist[2], cornerlist[3], tif, tif_extract))
    except WindowsError:
        print "Could not find tool!"
        return 1

    #resample file so that it has a maximumum number of pixels given by max_pixels
    resample(tif_extract, './resolution_reduced/test_%sx%s_r%sx%s.tif' % (num2str(l, 'length'), num2str(l, 'length'), num2str(res), num2str(res)), [res])
    return 0

def resample(tif, tif_reduced, cornerlist):
    tool = 'gdalwarp.exe'
    path_tool = 'C:/OSGeo4W64/bin/%s' % tool


    for res in cornerlist:
        try:
            subprocess.call('%s -tr %s %s %s %s -overwrite' % (path_tool, res, res, tif, tif_reduced))
        except WindowsError:
            print "Could not find tool!"
            return 1
        tif2xyz( '.' + tif_reduced.split('.')[1])

def tif2xyz(tif):
    tool = 'gdal_translate.exe'
    path_tool = 'C:/OSGeo4W64/bin/%s' % tool

    try:
        subprocess.call('%s -of XYZ %s.tif %s.xyz' % (path_tool, tif, tif))
    except WindowsError:
        print "Could not find tool!"
        return 1

if __name__=='__main__':
    cornerlist = corners(722600, 6184300, 250, -250)
    box_xyz(cornerlist)

    tif = "../01_DTM/DHYMRAIN.tif"
    tif_extract = "./resolution_full/test"
    box_tif(tif, tif_extract, cornerlist, 100000, 1.6)
