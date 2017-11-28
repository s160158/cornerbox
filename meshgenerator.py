#!/usr/bin/env python
from __future__ import division

"""
Generate regular rectangular interpolated grids from a .tif DEM file. Needs a two corners: upper-left and lower-right 
corners + the path to DEM .tif file extracted by that list using the Cornerbox class. 
"""


from osgeo import gdal
from cornerbox import *
from decimal import Decimal # Stupid floating point numbers and modulo
import numpy as np


class MeshGenerator():
    arr_e = None; arr_r = None  # To store NumPy arrays
    corners = []
    lorr = 0; uord = 0  # Direction of box edges from coordinate of first corner, i.e. (corners[0], corners[2])

    def __init__(self, tif):
        self.arr_e = self.tif2arr(tif)

    @staticmethod
    def tif2arr(tif):
        """
        Reads a TIF file into an NumPy array. (Remember there is no projection info in this array and thus alo no
        resolution info)
        :param tif: a TIF file
        :return: array holding TIF file values
        """
        ds = gdal.Open(tif, gdal.GA_ReadOnly)
        arr = np.array(ds.GetRasterBand(1).ReadAsArray())  # only one band - DEM
        return np.fliplr(np.transpose(arr))

    def get_resolution(self, arr):
        """
        Horizontal and vertical resolution of NumPy array
        :return:
        """
        # Horizontal resolution
        resolution_x = (self.corners[2] - self.corners[0]) / arr.shape[0]
        resolution_x = abs(resolution_x)  # Might not want to lose sign information
        self.lorr = np.sign(resolution_x)
        # Vertical resolution
        resolution_y = (self.corners[3] - self.corners[1]) / arr.shape[0]
        resolution_y = abs(resolution_y)
        self.uord = np.sign(resolution_y)

        return (resolution_x, resolution_y)

    def resample(self, resolution_x, resolution_y):
        """
        Resamples the input TIF (also ties the points geographically, i.e. creates mesh)
        :param resolution_x:
        :param resolution_y:
        :return:
        """
        # Resolution_x and resolution_y need both be divisible by the original tif resolutions in the respective
        # directions
        (resolution_x_e, resolution_y_e) = self.get_resolution(self.arr_e)

        if resolution_x != resolution_x_e or resolution_y != resolution_y_e: # Bad test!
            if Decimal(str(resolution_x)) % Decimal(str(resolution_x_e)) or \
                            Decimal(str(resolution_y)) % Decimal(str(resolution_y_e)):  # Really, this is a good solution?
                raise ValueError('Resolution not divisble!')                            # could do int(1000 * x) % int(...
                                                                                        # , or?

        jump_x = int(round(resolution_x / resolution_x_e))  # only integer jumps (again problems with float division)
        jump_y = int(round(resolution_y / resolution_y_e))

        # if the size of the image is not divisible by the new resolution there will be cutoffs (now warning displayed)
        self.arr_r = np.zeros([int(self.arr_e.shape[0] / jump_x), int(self.arr_e.shape[1] / jump_y), 3])

        for i in range(0, self.arr_r.shape[0]):
            for j in range(0, self.arr_r.shape[1]):
                self.arr_r[i, j, 0] = self.corners[0] + i * self.lorr * resolution_x
                self.arr_r[i, j, 1] = self.corners[1] + j * self.uord * resolution_y
                self.arr_r[i, j, 2] = self.arr_e[i * jump_x, j * jump_y]

    def write_header(self, filename, n_nodes):
        """
        Write header to a .mesh file (creates new .mesh file!)
        :param filename: path including filename of .mesh file to write
        :param n_nodes: number of nodes
        :return:
        """
        fh = open(filename, 'w')

        # eumIBathymetry = 100079, type meters = 1000, number of nodes, projections string
        fh.write('100079  1000  {}  PROJCS["WGS_1984_UTM_Zone_32N",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID'
                 '["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],'
                 'PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000],PARAMETER["False_Northing",0],'
                 'PARAMETER["Central_Meridian",9],PARAMETER["Scale_Factor",0.9996],'
                 'PARAMETER["Latitude_Of_Origin",0],UNIT["Meter",1]]'.format(int(n_nodes)))

        fh.close()

    def export_mesh(self):
        """
        Exports the resampled NumPy array to a .mesh file to be used with MIKE FM
        :return:
        """
        (resolution_x, resolution_y) = self.get_resolution(self.arr_r)
        (length_x, length_y) = Cornerbox.get_lengths(corners)
        filename = './extract/{}/m/box_{}x{}_r{}x{}.mesh'.format(Cornerbox.folder_name(self.corners), # take other (Cornerbox object)?
                                                                 Cornerbox.num2str(length_x, 'length'),
                                                                 Cornerbox.num2str(length_y, 'length'),
                                                                 Cornerbox.num2str(resolution_x, 'resolution'),
                                                                 Cornerbox.num2str(resolution_y, 'resolution'))
        print 'saving resample to: {}'.format(filename)
        self.write_header(filename, self.arr_r.size / 3)

        fh = open(filename, 'a')

        code = 0
        id = 1  # Node id
        for j in range(self.arr_r.shape[1]):
            for i in range(self.arr_r.shape[0]):

                # Node code: (land boundary is code 1)
                if i == 0 or j == 0 or i == self.arr_r.shape[0] - 1 or j == self.arr_r.shape[1] - 1:
                    code = 1

                fh.write('\n{} {} {} {} {}'.format(id,
                                               self.arr_r[i, j, 0],
                                               self.arr_r[i, j, 1],
                                               self.arr_r[i, j, 2],
                                               code))  # id X Y Z code
                code = 0
                id += 1

        fh.close()

        fh = open(filename, 'a')

        # number of elements, maximum number of vertices in an element, type/code
        element_header_line = '\n{} {} {}'.format((self.arr_r.shape[0] - 1) * (self.arr_r.shape[1] - 1), 4, 25)
        fh.write(element_header_line)
        print 'Wrote element header line: {}'.format(element_header_line)

        id_matrix = np.transpose(np.arange(1, self.arr_r.size / 3 + 1, dtype='uint16').reshape((self.arr_r.shape[0],
                                                                                                self.arr_r.shape[1])))
        id = 1  # Polygon id

        for j in range(id_matrix.shape[1] - 1):  # -1 as there are fewer elements than nodes/vertices
            for i in range(id_matrix.shape[0] - 1):
                element_line = '\n{} {} {} {} {}'.format(id,  # CCW connectivity direction
                                                         id_matrix[i,j],
                                                         id_matrix[i + 1,j],
                                                         id_matrix[i + 1,j + 1],
                                                         id_matrix[i, j + 1])  # id1, id2, id3, id4
                fh.write(element_line)
                id += 1

        print 'Wrote {} elements to file!'.format(id - 1)
        fh.close()


if __name__=='__main__':
    #corners = Cornerbox.corners(722600, 6184300, 280, -280) # For small 0.4x0.4 m resolution area
    corners = Cornerbox.corners(722600, 6184300, 560, -560)
    # .tif file to extract from:
    tif = "../01_DTM/DHYMRAIN.tif"

    box = Cornerbox(tif)
    box.corners = corners
    box.create_folder_structure()
    box.extract_tif()
    box.box_shp_boundary()

    mg = MeshGenerator(box.tif_e)
    mg.corners = corners
    #for resolution in [0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2]:
    #    mg.resample(resolution, resolution)
    #    mg.export_mesh()


    # 9x9 grid for debugging:
    mg.resample(400 * 0.4, 400 * 0.4)
    mg.export_mesh()
