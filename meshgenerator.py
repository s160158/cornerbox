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
    method = ''
    ds = None

    def __init__(self, tif):
        self.arr_e = self.tif2arr(tif)

    def tif2arr(self, tif):
        """
        Reads a TIF file into an NumPy array. (Remember there is no projection info in this array and thus alo no
        resolution info)
        :param tif: a TIF file
        :return: array holding TIF file values
        """
        ds = gdal.Open(tif, gdal.GA_ReadOnly)
        self.ds = ds
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

    def resample(self, resolution_x, resolution_y, method='near'):
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

        if method == 'near':
            for i in range(0, self.arr_r.shape[0]):
                for j in range(0, self.arr_r.shape[1]):
                    self.arr_r[i, j, 0] = self.corners[0] + (i + 0.5) * self.lorr * resolution_x
                    self.arr_r[i, j, 1] = self.corners[1] + (j + 0.5) * self.uord * resolution_y
                    self.arr_r[i, j, 2] = self.arr_e[int((i + 0.5) * jump_x), int((j + 0.5) * jump_y)]
        elif method == 'min':
            for i in range(0, self.arr_r.shape[0]):
                for j in range(0, self.arr_r.shape[1]):
                    self.arr_r[i, j, 0] = self.corners[0] + (i + 0.5) * self.lorr * resolution_x
                    self.arr_r[i, j, 1] = self.corners[1] + (j + 0.5) * self.uord * resolution_y
                    self.arr_r[i, j, 2] = np.min(self.arr_e[i * jump_x: (i + 1) * jump_x, j * jump_y: (j + 1) * jump_y])
        elif method == 'max':
            for i in range(0, self.arr_r.shape[0]):
                for j in range(0, self.arr_r.shape[1]):
                    self.arr_r[i, j, 0] = self.corners[0] + (i + 0.5) * self.lorr * resolution_x
                    self.arr_r[i, j, 1] = self.corners[1] + (j + 0.5) * self.uord * resolution_y
                    self.arr_r[i, j, 2] = np.max(self.arr_e[i * jump_x: (i + 1) * jump_x, j * jump_y: (j + 1) * jump_y])
        elif method == 'mean':
            for i in range(0, self.arr_r.shape[0]):
                for j in range(0, self.arr_r.shape[1]):
                    self.arr_r[i, j, 0] = self.corners[0] + (i + 0.5) * self.lorr * resolution_x
                    self.arr_r[i, j, 1] = self.corners[1] + (j + 0.5) * self.uord * resolution_y
                    self.arr_r[i, j, 2] = np.mean(self.arr_e[i * jump_x: (i + 1) * jump_x, j * jump_y: (j + 1) * jump_y])
        elif method == 'median':
            for i in range(0, self.arr_r.shape[0]):
                for j in range(0, self.arr_r.shape[1]):
                    self.arr_r[i, j, 0] = self.corners[0] + (i + 0.5) * self.lorr * resolution_x
                    self.arr_r[i, j, 1] = self.corners[1] + (j + 0.5) * self.uord * resolution_y
                    self.arr_r[i, j, 2] = np.median(self.arr_e[i * jump_x: (i + 1) * jump_x, j * jump_y: (j + 1) * jump_y])
        else:
            pass
        self.method = method

    def export_tif(self):
        # create the output image
        dvr = self.ds.GetDriver()
        # print driver
        (resolution_x, resolution_y) = self.get_resolution(self.arr_r)
        (length_x, length_y) = Cornerbox.get_lengths(corners)
        filename = './extract/{}/r/box_{}x{}_r{}x{}_{}.tif'.format(Cornerbox.folder_name(self.corners), # take other (Cornerbox object)?
                                                                 Cornerbox.num2str(length_x, 'length'),
                                                                 Cornerbox.num2str(length_y, 'length'),
                                                                 Cornerbox.num2str(resolution_x, 'resolution'),
                                                                 Cornerbox.num2str(resolution_y, 'resolution'),
                                                                 self.method)

        ds_out = dvr.Create(filename, self.arr_r.shape[0], self.arr_r.shape[1], 1, gdal.GDT_Float32)
        if ds_out is None:
            raise WindowsError('Could not create file')

        band = ds_out.GetRasterBand(1)

        # write the data
        band.WriteArray(np.transpose(np.fliplr(self.arr_r[:, :, 2])), 0, 0)

        # flush data to disk, set the NoData value and calculate stats
        band.FlushCache()
        band.SetNoDataValue(-99) # not a good value

        # georeference the image and set the projection
        transform = self.ds.GetGeoTransform()
        transform_out = (transform[0], resolution_x, transform[2], transform[3], transform[4], -resolution_y)  # set new pixel height/width
        ds_out.SetGeoTransform(transform_out)
        ds_out.SetProjection(self.ds.GetProjection())

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
        filename = './extract/{}/m/box_{}x{}_r{}x{}_{}.mesh'.format(Cornerbox.folder_name(self.corners), # take other (Cornerbox object)?
                                                                 Cornerbox.num2str(length_x, 'length'),
                                                                 Cornerbox.num2str(length_y, 'length'),
                                                                 Cornerbox.num2str(resolution_x, 'resolution'),
                                                                 Cornerbox.num2str(resolution_y, 'resolution'),
                                                                 self.method)
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

        id_matrix = np.transpose(np.arange(1, self.arr_r.size / 3 + 1, dtype='uint32').reshape((self.arr_r.shape[0],
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

    def export_mesh_new(self):
        """
        Exports the resampled NumPy array to a .mesh file to be used with MIKE FM. Node/vertex id ordering is done
        exactly as MIKE Mesh Generator does it (At least for regular rectangular meshes). Turns out not to be important.
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
        fh_nodes = open('./temp.mesh', 'w')

        id_matrix = np.zeros((self.arr_r.shape[0], self.arr_r.shape[1]), dtype='uint32')  # initiate matrix (shit need 32 bit values in many cases!)

        n_v = (self.arr_r.shape[0] - 2) * (self.arr_r.shape[1] - 2)  # number of vertices
        id_v = 1
        id_n = n_v + 1  # starting node id
        for j in range(self.arr_r.shape[1]):
            for i in range(self.arr_r.shape[0]):
                # check for special circumstances when we encounter a node:
                if i == 0 or j == 0 or i == self.arr_r.shape[0] - 1 or j == self.arr_r.shape[1] - 1:
                    # id to write to node list and increment id_n
                    if i == 0 and j == 1:
                        id_matrix[i, j] = (n_v + 1) + 2
                        id_matrix[i, j] = (n_v + 1) + 2
                    elif i == 0 and j == self.arr_r.shape[1] - 1:
                        id_matrix[i, j] = id_n + 1
                        id_matrix[i, j] = id_n + 1
                    elif i == 1 and j == self.arr_r.shape[1] - 1:
                        id_matrix[i, j] = id_n
                        id_n += 2
                    elif i == 1 and j == 0:
                        id_matrix[i, j] = id_n
                        id_n += 2
                    else:
                        id_matrix[i, j] = id_n
                        id_n += 1

                    # write to temporary node list
                    fh_nodes.write('\n{} {} {} {} {}'.format(id_matrix[i, j],
                                                       self.arr_r[i, j, 0],
                                                       self.arr_r[i, j, 1],
                                                       self.arr_r[i, j, 2],
                                                       1))  # id X Y Z code
                else:  # not a node; a vertex
                    id_matrix[i, j] = id_v
                    fh.write('\n{} {} {} {} {}'.format(id_matrix[i, j],
                                                       self.arr_r[i, j, 0],
                                                       self.arr_r[i, j, 1],
                                                       self.arr_r[i, j, 2],
                                                       0))  # id X Y Z code
                    id_v += 1
        fh_nodes.close()  # to get at beginning - better way?

        # append the nodes at the end of the vertices
        fh_nodes = open('./temp.mesh', 'r')
        for line in fh_nodes:
            fh.write(line)
        fh_nodes.close()

        fh.close() # close only to open with append immidieately after!

        fh = open(filename, 'a')

        # number of elements, maximum number of vertices in an element, type/code
        element_header_line = '\n{} {} {}'.format((self.arr_r.shape[0] - 1) * (self.arr_r.shape[1] - 1), 4, 25)
        fh.write(element_header_line)
        print 'Wrote element header line: {}'.format(element_header_line)

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
    #    mg.export_mesh_new()

    #
    #mg.resample(1.6, 1.6)
    #mg.export_mesh()

    # 9x9 grid for debugging:
    #mg.resample(400 * 0.4, 400 * 0.4, 'max')
    mg.resample(6.4, 6.4, 'mean')
    mg.export_mesh()
    mg.export_tif()
    #mg.resample(400 * 0.4, 400 * 0.4, 'min')
    #mg.export_mesh()
    #mg.resample(400 * 0.4, 400 * 0.4, 'near')
    #mg.export_mesh()
    #mg.resample(400 * 0.4, 400 * 0.4, 'mean')
    #mg.export_mesh()

    """
    for resolution in [0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2]:
        mg.resample(resolution, resolution, 'min')
        mg.export_tif()
        mg.resample(resolution, resolution, 'max')
        mg.export_tif()
        mg.resample(resolution, resolution, 'near')
        mg.export_tif()
        mg.resample(resolution, resolution, 'mean')
        mg.export_tif()
        mg.resample(resolution, resolution, 'median')
        mg.export_tif()
    """
