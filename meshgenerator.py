#!/usr/bin/env python
from __future__ import division

from osgeo import gdal
from cornerbox import *

import numpy as np


class MeshGenerator():
    arr_e = None; arr_r = None  # To store NumPy arrays (arr_e is the extracted tif with boundaries)
    corners = []
    method = ''
    ds = None
    resolution_x_e = 0.0; resolution_y_e = 0.0
    resolution_x = 0.0; resolution_y = 0.0

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
        # Get the resolution of the input raster
        transform = ds.GetGeoTransform()
        self.resolution_x_e = abs(transform[1])  # Resolution of original raster (x-direction)
        self.resolution_y_e = abs(transform[5])  # Resolution of original raster (y-direction)
        self.ds = ds
        arr = np.array(ds.GetRasterBand(1).ReadAsArray())  # only one band - DEM
        return np.fliplr(np.transpose(arr))

    def resample(self, resolution_x, resolution_y, method='near'):
        """
        Resamples the input TIF (also ties the points geographically, i.e. creates mesh)
        :param resolution_x:
        :param resolution_y:
        :return:
        """
        self.resolution_x = resolution_x
        self.resolution_y = resolution_y

        jump_x = int(round(resolution_x / self.resolution_x_e))  # only integer jumps (again problems with float division)
        jump_y = int(round(resolution_y / self.resolution_y_e))

        # if the size of the image is not divisible by the new resolution there will be cutoffs (now warning displayed)
        self.arr_r = np.zeros([int((self.arr_e.shape[0] - 2 * 20.0 / self.resolution_x_e) / jump_x) + 1,
                               int((self.arr_e.shape[1] - 2 * 20.0 / self.resolution_y_e) / jump_y) + 1,
                               3])

        corr_x = int(20.0 / self.resolution_x_e)
        corr_y = int(20.0 / self.resolution_y_e)

        if method == 'near':
            for i in range(0, self.arr_r.shape[0]):
                for j in range(0, self.arr_r.shape[1]):
                    self.arr_r[i, j, 0] = self.corners[0] + i * resolution_x
                    self.arr_r[i, j, 1] = self.corners[1] - j * resolution_y  # Go South
                    self.arr_r[i, j, 2] = self.arr_e[int(i * jump_x + corr_x), int(j * jump_y + corr_y)]
        elif method == 'min':
            for i in range(0, self.arr_r.shape[0]):
                for j in range(0, self.arr_r.shape[1]):
                    self.arr_r[i, j, 0] = self.corners[0] + i * resolution_x
                    self.arr_r[i, j, 1] = self.corners[1] - j * resolution_y
                    self.arr_r[i, j, 2] = np.min(self.arr_e[(i - 1) * jump_x + jump_x // 2 + corr_x: (i + 1) * jump_x + jump_x // 2 + corr_x,
                                                            (j - 1) * jump_y + jump_y // 2 + corr_y: (j + 1) * jump_y + jump_y // 2 + corr_y])
        elif method == 'max':
            for i in range(0, self.arr_r.shape[0]):
                for j in range(0, self.arr_r.shape[1]):
                    self.arr_r[i, j, 0] = self.corners[0] + i * resolution_x
                    self.arr_r[i, j, 1] = self.corners[1] - j * resolution_y
                    self.arr_r[i, j, 2] = np.max(self.arr_e[(i - 1) * jump_x + jump_x // 2 + corr_x: (i + 1) * jump_x + jump_x // 2 + corr_x,
                                                            (j - 1) * jump_y + jump_y // 2 + corr_y: (j + 1) * jump_y + jump_y // 2 + corr_y])
        elif method == 'mean':
            for i in range(0, self.arr_r.shape[0]):
                for j in range(0, self.arr_r.shape[1]):
                    self.arr_r[i, j, 0] = self.corners[0] + i * resolution_x
                    self.arr_r[i, j, 1] = self.corners[1] - j * resolution_y
                    self.arr_r[i, j, 2] = np.mean(self.arr_e[(i - 1) * jump_x + jump_x // 2 + corr_x: (i + 1) * jump_x + jump_x // 2 + corr_x,
                                                             (j - 1) * jump_y + jump_y // 2 + corr_y: (j + 1) * jump_y + jump_y // 2 + corr_y])
        elif method == 'median':
            for i in range(0, self.arr_r.shape[0]):
                for j in range(0, self.arr_r.shape[1]):
                    self.arr_r[i, j, 0] = self.corners[0] + i * resolution_x
                    self.arr_r[i, j, 1] = self.corners[1] - j * resolution_y
                    self.arr_r[i, j, 2] = np.median(self.arr_e[(i - 1) * jump_x + jump_x // 2 + corr_x: (i + 1) * jump_x + jump_x // 2 + corr_x,
                                                               (j - 1) * jump_y + jump_y // 2 + corr_y: (j + 1) * jump_y + jump_y // 2 + corr_y])
        else:
            raise ValueError('method: {} was not recognised.'.format(method))

        self.method = method

    def export_tif(self):
        # create the output image
        dvr = self.ds.GetDriver()
        # print driver
        (length_x, length_y) = Cornerbox.get_lengths(corners)
        filename = './extract/{}/r/box_{}x{}_r{}x{}_{}.tif'.format(Cornerbox.folder_name(self.corners), # take other (Cornerbox object)?
                                                                 Cornerbox.num2str(length_x, 'length'),
                                                                 Cornerbox.num2str(length_y, 'length'),
                                                                 Cornerbox.num2str(self.resolution_x, 'resolution'),
                                                                 Cornerbox.num2str(self.resolution_y, 'resolution'),
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
        transform_out = (self.corners[0] - 0.5 * self.resolution_x, self.resolution_x, transform[2],
                         self.corners[1] + 0.5 * self.resolution_y, transform[4], -self.resolution_y)  # set new pixel height/width

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

        (length_x, length_y) = Cornerbox.get_lengths(corners)
        filename = './extract/{}/m/box_{}x{}_r{}x{}_{}.mesh'.format(Cornerbox.folder_name(self.corners), # take other (Cornerbox object)?
                                                                 Cornerbox.num2str(length_x, 'length'),
                                                                 Cornerbox.num2str(length_y, 'length'),
                                                                 Cornerbox.num2str(self.resolution_x, 'resolution'),
                                                                 Cornerbox.num2str(self.resolution_y, 'resolution'),
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


if __name__=='__main__':
    corners = Cornerbox.corners(722600, 6184300, 281.6, -281.6) # For small 0.4x0.4 m resolution area

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
    #mg.resample(93.2, 93.2, 'near')


    mg.resample(0.8, 0.8, 'min')
    mg.export_mesh()
    mg.export_tif()

    mg.resample(0.8, 0.8, 'near')
    mg.export_mesh()
    mg.export_tif()


    #print mg.arr_e[0,0,1], mg.arr_e[0, 0, -1]
    #mg.resample(400 * 0.4, 400 * 0.4, 'min')
    #mg.export_mesh()
    #mg.resample(400 * 0.4, 400 * 0.4, 'near')
    #mg.export_mesh()
    #mg.resample(400 * 0.4, 400 * 0.4, 'mean')
    #mg.export_mesh()

"""
    mg.resample(140.8, 140.8, 'near')
    mg.export_mesh()
    mg.export_tif()

    mg.resample(140.8, 140.8, 'min')
    mg.export_mesh()
    mg.export_tif()

    mg.resample(70.4, 70.4, 'near')
    mg.export_mesh()
    mg.export_tif()


    for resolution in [0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2, 4.4, 6.4, 8.0]:
        mg.resample(resolution, resolution, 'min')
        mg.export_mesh()
        mg.export_tif()
        mg.resample(resolution, resolution, 'max')
        mg.export_mesh()
        mg.export_tif()
        mg.resample(resolution, resolution, 'near')
        mg.export_mesh()
        mg.export_tif()
        mg.resample(resolution, resolution, 'mean')
        mg.export_mesh()
        mg.export_tif()
        mg.resample(resolution, resolution, 'median')
        mg.export_mesh()
        mg.export_tif()

    mg.resample(0.4, 0.4, 'near')
    mg.export_mesh()
    mg.export_tif()
"""