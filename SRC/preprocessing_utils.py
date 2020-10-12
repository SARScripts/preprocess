#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 17:32:28 2020

@author: pacomoral
"""


import sys
import os
import subprocess
import pandas as pd
from datetime import datetime 
import glob
import xml.etree.ElementTree as ET
from zipfile import ZipFile
import shutil
import numpy as np
from scipy.interpolate import griddata
import gdal

def listar (a): #Genera una lista de un string separado por comas leído de un config file
    if pd.isnull(a):
        return ''
    else:
        return str(a).split(",") 

def unzipS1 (zipfile, dirouput):
    with ZipFile(zipfile, 'r') as zipObj:
        zipObj.extractall(dirouput)
        
def mountpairs (dflist, field): #Pairs combining the first element with the rest of the list elements
    pairs = []
    pairs = [((dflist[field][0]), (dflist[field][i+1])) for i in range(len(dflist)-1)]
    return pairs

class model():
    def __init__(self, process_type, imagelistfile, diroutorder, overwrite, datemaster, AOI, pathgpt, snappypath, DirProject):
        self.process_type = process_type
        self.imagelistfile = imagelistfile
        self.datemaster = datemaster
        if self.datemaster == '':
            self.datemaster = 'nodate'
        self.AOI = AOI
        if not self.AOI:
            self.AOI = ['', '', '', '']
        self.diroutorder = diroutorder
        self.overwrite = overwrite
        self.pathgpt = pathgpt
        self.snappypath = snappypath
        self.DirProject = DirProject
        self.processdf = pd.DataFrame()
        self.baselinelist = pd.DataFrame()
        self.maxbaseperp = ''
        self.maxbasetemp = ''
        self.baselinefiltered = pd.DataFrame()
        self.logfilename = 'process_log.csv'
        
        if not os.path.exists(self.diroutorder):
            os.makedirs(self.diroutorder)
            
        sys.path.append(self.snappypath)
        import snappy
        from snappy import ProductIO, GPF, HashMap

    def searchMetabytag(self, metafilepath, tag):
        meta = []
        tree = ET.parse(metafilepath)
        root = tree.getroot()
        for atr in root.iter('MDATTR'):
            if atr.get('name') == tag:
                meta.append(atr.text)
        return meta

    def unzipgetmeta (self, imagelist):
        #Unzip if zipfiles in imagelist and fill up ImageUnzip, ImageName, ImageDate columns
        ImageUnzip = []
        ImageName = []
        ImageDate = []
        diroutput = os.path.join(self.diroutorder, '00_datafiles')
        if not os.path.exists(diroutput):
            os.makedirs(diroutput)
        for i in range(len(imagelist)):
            if os.path.splitext(imagelist.Image[i])[1] == '.zip':
                #OVERWRITE CONDITION FOR UNZIP FILES
                if not os.path.isdir(os.path.join(diroutput, os.path.basename(os.path.dirname(imagelist.Image[i])))) or self.overwrite == '1':
                    unzipS1 (imagelist.Image[i], diroutput)
                ImageUnzip.append(os.path.join(diroutput, os.path.splitext(os.path.basename(imagelist.Image[i]))[0]+'.SAFE', 'manifest.safe'))
                ImageName.append(os.path.splitext(os.path.basename(imagelist.Image[i]))[0])
            elif os.path.splitext(imagelist.Image[i])[1] == '.safe':
                if not os.path.isdir(os.path.join(diroutput, os.path.basename(os.path.dirname(imagelist.Image[i])))) or self.overwrite == '1':
                    if os.path.isdir(os.path.join(diroutput, os.path.basename(os.path.dirname(imagelist.Image[i])))):
                        shutil.rmtree(os.path.join(diroutput, os.path.basename(os.path.dirname(imagelist.Image[i]))))
                    shutil.copytree(os.path.dirname(imagelist.Image[i]), os.path.join(diroutput, os.path.basename(os.path.dirname(imagelist.Image[i]))))
                ImageName.append(os.path.basename(os.path.dirname(imagelist.Image[i])).split('.')[0])
                ImageUnzip.append(os.path.join(diroutput, os.path.basename(os.path.dirname(imagelist.Image[i])), 'manifest.safe'))
            else:
                print('Image format not supported: ' + imagelist.Image[i])
                sys.exit()
            ImageDate.append(ImageName[i][17:25])
        imagelist['ImageUnzip'] = ImageUnzip
        imagelist['ImageName'] = ImageName
        imagelist['ImageDate'] = ImageDate
        return imagelist

    def run_alignmentlist(self, pairspath, pairsdate):
        gptxml_file = os.path.join(self.DirProject, 'RES', 'TOPSAR_Coreg_Interferogram_ESD_Polariz_args_wkt.xml')
        
        #CHANGE OUTPUT SUFIX
        output_sufix = '_Aligned'
        subdirout = '01_slc'
        extent = self.AOI
        # (lonmin, latmax, lonmax, latmin)
        lonmin = extent[0]
        latmax = extent[1]
        lonmax = extent[2]
        latmin = extent[3]
        outputfiles = []
        outputtimes = []
        datelist = []
        for i in range(len(pairspath)):
            time1 = datetime.now()
            input1 = pairspath[i][0]
            input2 = pairspath[i][1]
            if '_IW_SLC_' in input1:
                input3 = ['IW1', 'IW2', 'IW3']
            if '_1SDV_' in input1:
                input5 = ['VV', 'VH']
            if any (self.AOI):
                input4 = 'POLYGON (('+lonmin+' '+latmin+', '+lonmax+' '+latmin+', '+lonmax+' '+latmax+', '+lonmin+' '+latmax+', '+lonmin+' '+latmin+'))'
            else:
                input4 = ''
            #The process is executed for every subswath
            for polariz in input5:
                if os.path.splitext(input2)[1] == '.safe':
                    slavename = pairsdate[i][1] + '_SLC'
                else:
                    print('Image format not supported: ' + input2)
                    sys.exit()
                output1 = os.path.join(self.diroutorder, subdirout, self.datemaster + '_' + slavename + output_sufix + '_' + polariz)
                for k in range(len(input3)):
                    outputfile = os.path.join(output1 + '_' + input3[k] + '.dim')
                    #OVERWRITE CONDITION FOR SUBSWATH GENERATION
                    if not os.path.isfile(outputfile) or self.overwrite == '1':
                        try:
                            p = subprocess.check_output([self.pathgpt, gptxml_file, '-Pinput1='+input1, '-Pinput2='+input2, '-Pinput3='+input3[k], '-Pinput4='+input4, '-Pinput5='+polariz,'-Ptarget1='+outputfile])
                        except:
                            p = ''
                #After processing subswaths, they are put together with TOPS Merge Workflow
                #List of subswaths by checking the output files
                subswath_list = sorted(glob.glob(os.path.join(output1 + '*.dim')))
                #Function TOPS Merge is only called when more than 1 subswath is generated. Only IW products supported (max. 3 subswaths)
                if len(subswath_list)>1:
                    #OVERWRITE CONDITION FOR MERGE
                    if not os.path.isfile(output1 + '.dim') or self.overwrite == '1':
                        m = self.TOPS_Merge_subswaths(subswath_list, output1 + '.dim', self.pathgpt)
                    for filePath in subswath_list:
                        try:
                            os.remove(filePath)
                            shutil.rmtree(os.path.splitext(filePath)[0]+'.data')
                        except:
                            print("Error while deleting file : ", filePath)
                else:
                    #Change output name from IWx to 'merged' in case AOI covers only one subswath
                    shutil.move(subswath_list[0], output1+'.dim')
                    shutil.move(os.path.splitext(subswath_list[0])[0]+'.data', output1+'.data')
            if os.path.isfile(output1+'.dim'):
                outputfiles.append(output1+'.dim')
            else:
                outputfiles.append('Not generated')
            outputtimes.append((datetime.now()-time1).seconds/60)
            datelist.append(pairsdate[i][1])
        #Fill processing info in df
        self.processdf['Image_date'] = datelist
        self.processdf['Outputfiles_align'] = outputfiles
        self.processdf['Processing_minutes'] = outputtimes
        self.processdf.to_csv(os.path.join(self.diroutorder, self.logfilename), sep=';', index=False)

        
    def TOPS_Merge_subswaths (self, imagepathlist, outputfile, pathgpt):
        if(len(imagepathlist)==2):
            gptxml_file = os.path.join(self.DirProject, 'RES', 'merge_2subswaths.xml')
            input1 = imagepathlist[0]
            input2 = imagepathlist[1]
            try:
                m = subprocess.check_output([pathgpt, gptxml_file, '-Pinput1='+input1, '-Pinput2='+input2, '-Ptarget1='+outputfile])
            except:
                m = ''
        if(len(imagepathlist)==3):
            gptxml_file = os.path.join(self.DirProject, 'RES', 'merge_3subswaths.xml')
            input1 = imagepathlist[0]
            input2 = imagepathlist[1]
            input3 = imagepathlist[2]
            try:
                m = subprocess.check_output([pathgpt, gptxml_file, '-Pinput1='+input1, '-Pinput2='+input2, '-Pinput3='+input3, '-Ptarget1='+outputfile])
            except:
                m = ''
        return m
    
    def alignment_ifg (self):
        #Read the image file
        imagelist = pd.read_csv(self.imagelistfile, sep=';', decimal=',')
        #Unzip and get initial metadata
        imagelist = self.unzipgetmeta (imagelist)

        #IF MASTER IMAGE DATE MATCHES ANY DATE OF THE FILE LIST, RUN ALIGNMENT PROCESSING (OTHERWISE MASTER DATE IS COMPUTED AS THE 'MID POINT' DATE OF THE LIST)
        #Identify master image row and set up slave images dataframe
        if not any (imagelist['ImageDate'].str.contains(self.datemaster)):
            datelist = list(imagelist.ImageDate.sort_values().astype(int))
            self.datemaster = str(datelist[min(range(len(datelist)), key = lambda i: abs(datelist[i]-(datelist[-1]+datelist[0])/2))])
        #Master image is reindexed to the top of the df
        masterdf = imagelist[imagelist['ImageDate'].str.contains(self.datemaster)]
        slavedf = imagelist.copy()
        slavedf = slavedf.drop(slavedf[imagelist['ImageDate'].str.contains(self.datemaster)].index)
        pairspath = mountpairs(pd.concat([masterdf, slavedf], sort=False).reset_index(), 'ImageUnzip')
        pairsdate = mountpairs(pd.concat([masterdf, slavedf], sort=False).reset_index(), 'ImageDate')
        self.run_alignmentlist(pairspath, pairsdate)

    def generate_baselinelist (self):
        self.baselinelist = self.processdf.copy()
        self.baselinelist['Master'] = self.datemaster
        self.baselinelist['Slave'] = self.processdf['Image_date']
        self.baselinelist['Perp_Baseline'] = [self.searchMetabytag(self.processdf['Outputfiles'][i], 'Perp Baseline')[1] for i in range(len(self.processdf))]
        self.baselinelist['Temp_Baseline'] = [self.searchMetabytag(self.processdf['Outputfiles'][i], 'Temp Baseline')[1] for i in range(len(self.processdf))]
        self.baselinelist = self.baselinelist.drop(['Image_date', 'Outputfiles', 'Processing_minutes'], axis=1)
        waninglist = self.baselinelist.copy()
        while len(waninglist)>1:
            for i in range(len(waninglist)-1):
                perp = float(waninglist['Perp_Baseline'][0]) - float(waninglist['Perp_Baseline'][i+1])
                temp = float(waninglist['Temp_Baseline'][0]) - float(waninglist['Temp_Baseline'][i+1])
                self.baselinelist.loc[len(self.baselinelist)] = str(waninglist['Slave'][0]), str(waninglist['Slave'][i+1]), perp, temp
            waninglist = waninglist.drop(0)
            waninglist = waninglist.reset_index(drop=True)
        #Filter table with maxbasetemp and maxbaseperp
        self.baselinelist['Temp_Baseline'] = abs(pd.to_numeric(self.baselinelist['Temp_Baseline']))
        self.baselinelist['Perp_Baseline'] = abs(pd.to_numeric(self.baselinelist['Perp_Baseline']))
        if self.maxbasetemp != '':
            if self.maxbaseperp != '':
                self.baselinefiltered = self.baselinelist[(self.baselinelist.Temp_Baseline <= self.maxbasetemp) & (self.baselinelist.Perp_Baseline <= self.maxbaseperp)]
            else:
                self.baselinefiltered = self.baselinelist[(self.baselinelist.Temp_Baseline <= self.maxbasetemp)]
        else:
            if self.maxbaseperp != '':
                self.baselinefiltered = self.baselinelist[(self.baselinelist.Perp_Baseline <= self.maxbaseperp)]
            else:
                self.baselinefiltered = self.baselinelist.copy()
        self.baselinefiltered.to_csv(os.path.join(self.diroutorder, 'baselines.csv'), sep=';', index=False)

    def generate_interferogram(self):
        print('entra generate interf')
        gptxml_file = os.path.join(self.DirProject, 'RES', 'Interferogram_Multilook.xml')
        #for i in range(len(self.baselinefiltered)):
        #Fill output files if pair is already coregistered
        for i in range(len(self.processdf)):
            for j in range(len(self.baselinefiltered)):
                if self.processdf['Image_date'][i] == self.baselinefiltered['Slave'][j] and self.baselinefiltered['Master'][j] == self.datemaster:
                    self.baselinefiltered.loc[j, 'Outputfiles'] = self.processdf['Outputfiles'][i]

    def applyMultilook(self, inputfile, outputfile, azLooks, rgLooks):
        sys.path.append(self.snappypath)
        import snappy
        from snappy import ProductIO, GPF, HashMap
        HashMap = snappy.jpy.get_type('java.util.HashMap')
        #Get snappy Operators
        GPF.getDefaultInstance().getOperatorSpiRegistry().loadOperatorSpis() 
        if not os.path.isfile(inputfile):
            print('Input file ' + inputfile + ' does not exist')
            sys.exit()
        else:
            time1 = datetime.now()
            parameters = HashMap()
            parameters.put('grSquarePixel', True)
            parameters.put('nRgLooks', rgLooks)
            parameters.put('nAzLooks', azLooks)
            parameters.put('outputIntensity', False)
            multi_param = snappy.GPF.createProduct("Multilook", parameters, ProductIO.readProduct(inputfile))
            ProductIO.writeProduct(multi_param, outputfile, 'BEAM-DIMAP')
            return((datetime.now()-time1).seconds/60)
            
    def Multilook(self):
        if (self.processdf).empty:
            self.processdf = pd.read_csv(os.path.join(self.diroutorder, self.logfilename), sep=';', decimal=',')
        subdirout = '02_ml'
        output_sufix = '_ML_' + self.azLooks + '_' + self.rgLooks
        outputfiles = []
        outputtimes = []
        for i in range(len(self.processdf['Outputfiles'])):
            if os.path.isfile(self.processdf['Outputfiles'][i]):
                inputfile = self.processdf['Outputfiles'][i]
                outputfile = os.path.join(self.diroutorder, subdirout, os.path.splitext(os.path.basename(inputfile))[0] + output_sufix + '.dim')
            outputtimes.append(self.applyMultilook(inputfile, outputfile, self.azLooks, self.rgLooks))
            outputfiles.append(outputfile)
        self.processdf['Outputfiles_ML'] = outputfiles
        self.processdf['Processing_minutes_ML'] = outputtimes
        self.processdf.to_csv(os.path.join(self.diroutorder, 'process_log.csv'), sep=';', index=False)



#     def applyGeocoding(self, inputfile, outputfile):
#         inputfile = '/media/guadarrama/PROYECTOS/SAR2CUBE/OUTPUT_noESD_3secSRTM_clean/02_ml/20200327_20200303_SLC_Aligned_ML_4_19.dim'
#         prod = ProductIO.readProduct(inputfile)
#         lat = prod.getTiePointGrid('latitude')
#         lon = prod.getTiePointGrid('longitude')
#         w = prod.getSceneRasterWidth()
#         h = prod.getSceneRasterHeight()
#         array = np.zeros((w, h), dtype=np.float32)
#         latpixels = lat.readPixels(0, 0, w, h, array)
#         lat_arr = np.asarray(latpixels)
#         #lat_arr.shape = (h, w)
#         lonpixels = lon.readPixels(0, 0, w, h, array)
#         lon_arr = np.asarray(lonpixels)
#         #lon_arr.shape = (h, w)
#         band_names = prod.getBandNames()
#         list(band_names)
#         data_getband = prod.getBand('i_ifg_VV_27Mar2020_03Mar2020')
#         data_pixels = data_getband.readPixels(0, 0, w, h, array)
#         data_arr = np.asarray(data_pixels)
#         #Define reference grid
#         x = np.linspace(np.min(lon_arr), np.max(lon_arr), w)
#         y = np.linspace(np.min(lat_arr), np.max(lat_arr), h)
#         grid_x, grid_y = np.meshgrid(x,y)
        
#         grid_data = griddata((lon_arr, lat_arr), data_arr.flatten(), (grid_x, grid_y), method='nearest')
        
#         import matplotlib.pyplot as plt
#         plt.imshow(grid_data)
        
        
#         (x,y) = data.shape
# 	format = "GTiff"
# 	noDataValue = -9999
# 	driver = gdal.GetDriverByName(format)
# 	# you can change the dataformat but be sure to be able to store negative values including -9999
# 	dst_datatype = gdal.GDT_Float32

# 	#print(data)

# 	dst_ds = driver.Create(filename,y,x,1,dst_datatype)
# 	dst_ds.GetRasterBand(1).WriteArray(data)
# 	dst_ds.GetRasterBand(1).SetNoDataValue( noDataValue )
# 	dst_ds.SetGeoTransform(geotransform)
# 	dst_ds.SetProjection(geoprojection)
        
        
        
        
        
        
        
        
#     def Geocoding(self):
#         subdirout = '03_gc'
#         output_sufix = '_GC'
#         for i in range(len(self.processdf['Multilook file'])):
#             if os.path.isfile(self.processdf['Multilook file'][i]):
#                 inputfile = self.processdf['Outputfiles'][i]
#                 outputfile = os.path.join(self.diroutorder, subdirout, os.path.splitext(os.path.basename(inputfile))[0] + output_sufix + '.dim')
#                 self.applyGeocoding(inputfile, outputfile, self.azLooks, self.rgLooks)
            
#     def array2raster(newRasterfn, rasterOrigin, pixelWidth, pixelHeight, array):
#         cols = array.shape[1]
#         rows = array.shape[0]
#         originX = rasterOrigin[0]
#         originY = rasterOrigin[1]
#         driver = gdal.GetDriverByName('ENVI')
#         outRaster = driver.Create(newRasterfn, cols, rows, 1, gdal.GDT_Byte)
#         outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
#         outband = outRaster.GetRasterBand(1)
#         outband.WriteArray(array)
#         outRasterSRS = osr.SpatialReference()
#         outRasterSRS.ImportFromEPSG(4326)
#         outRaster.SetProjection(outRasterSRS.ExportToWkt())
#         outband.FlushCache()            
            
#             indir='/media/guadarrama/PROYECTOS/SAR2CUBE/OUTPUT_SNAP/S1A_IW_SLC__1SDV_20200502T060119_20200502T060146_032382_03BFB9_57CF_Orb_Stack_Ifg_Deb.dim'
# prod = ProductIO.readProduct(indir)
# list(prod.getTiePointGridNames())
# lat = prod.getTiePointGrid('latitude')
# lon = prod.getTiePointGrid('longitude')
 
# lat_data = lat.getGridData()
# #latnode=prod.getRasterDataNode('latitude')

# w = prod.getSceneRasterWidth()
# h = prod.getSceneRasterHeight()
# array = np.zeros((w, h), dtype=np.float32)

# array.shape
# latpixels = lat.readPixels(0, 0, w, h, array)
# lat_arr = np.asarray(latpixels)
# lat_arr.shape = (h, w)

# lonpixels = lon.readPixels(0, 0, w, h, array)
# lon_arr = np.asarray(lonpixels)
# lon_arr.shape = (h, w)

# #Read data from product
# band_names = prod.getBandNames()
# list(band_names)
# data_getband = prod.getBand('i_ifg_IW1_VV_02May2020_14May2020')
# data_pixels = data_getband.readPixels(0, 0, w, h, array)
# data_pixels.shape = h, w

# #Subset for testing
# lat_sub = lat_arr[4900:5000, 4900:5000]
# lon_sub = lon_arr[4900:5000, 4900:5000]
# data_sub = data_pixels[4900:5000, 4900:5000]


# lat_sub_flatten = lat_sub.flatten()
# lon_sub_flatten = lon_sub.flatten()
# data_sub_flatten = data_sub.flatten()




# #Define reference grid
# x = np.linspace(np.min(lon_sub_flatten), np.max(lon_sub_flatten),100)
# y = np.linspace(np.min(lat_sub_flatten), np.max(lat_sub_flatten),100)
# grid_x, grid_y = np.meshgrid(x,y)

# from scipy.interpolate import griddata
# #grid_z0 = griddata(points, values, (grid_x, grid_y), method='nearest')


# grid_z0 = griddata((lon_sub_flatten, lat_sub_flatten), data_sub_flatten, (grid_x, grid_y), method='nearest')
# grid_z1 = griddata((lon_sub_flatten, lat_sub_flatten), data_sub_flatten, (grid_x, grid_y), method='linear')
# grid_z2 = griddata((lon_sub_flatten, lat_sub_flatten), data_sub_flatten, (grid_x, grid_y), method='cubic')
