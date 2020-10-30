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
from pyproj import Proj, transform
#import snappy
from snappy import ProductIO, GPF, HashMap
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

    def run_alignmentlist(self, pairspath, pairsdate, sufix):
        gptxml_file = os.path.join(self.DirProject, 'RES', 'TOPSAR_Coreg_Interferogram_ESD_Polariz_args_wkt.xml')
        #gptconvert_file = os.path.join(self.DirProject, 'RES', 'convert_DIMAP.xml')
        gptsubset_master = os.path.join(self.DirProject, 'RES', 'applyAOI_master_convertDIMAP.xml')

        #CHANGE OUTPUT SUFIX
        output_sufix = '_Aligned'
        subdirout = '01_slc' + sufix
        # (lonmin, latmax, lonmax, latmin)
        lonmin = self.AOI[0]
        latmax = self.AOI[1]
        lonmax = self.AOI[2]
        latmin = self.AOI[3]
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
            
            if os.path.splitext(input2)[1] == '.safe':
                slavename = pairsdate[i][1] + '_SLC'
            else:
                print('Image format not supported: ' + input2)
                sys.exit()
            imagebasename = os.path.join(self.diroutorder, subdirout, self.datemaster + '_' + slavename + output_sufix + sufix)
            if not os.path.isfile(imagebasename + '.dim') or self.overwrite == '1':
                for polariz in input5:
                    output1 = imagebasename + '_' + polariz
                    for k in range(len(input3)):
                        outputfile = os.path.join(output1 + '_' + input3[k] + '.dim')
                        #OVERWRITE CONDITION FOR SUBSWATH GENERATION
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
                polariz_list = sorted(glob.glob(imagebasename + '*V*.dim'))
                if len(polariz_list)>1:
                    #Stack polarizations (dual)
                    m = self.stack_polariz(polariz_list, imagebasename + '.dim', self.pathgpt)
                    for filePath in polariz_list:
                        try:
                            os.remove(filePath)
                            shutil.rmtree(os.path.splitext(filePath)[0]+'.data')
                        except:
                            print("Error while deleting file : ", filePath)
                else:
                    #Change output name from IWx to 'merged' in case AOI covers only one subswath
                    shutil.move(polariz_list[0], imagebasename+'.dim')
                    shutil.move(os.path.splitext(polariz_list[0])[0]+'.data', imagebasename+'.data')
            if os.path.isfile(imagebasename+'.dim'):
                outputfiles.append(imagebasename+'.dim')
            else:
                outputfiles.append('Not generated')
            outputtimes.append((datetime.now()-time1).seconds/60)
            datelist.append(pairsdate[i][1])
        #Subset master image by AOI and convert to DIMAP format for next processing, fill processing table.
        imagebasename = os.path.join(self.diroutorder, subdirout, self.datemaster + '_' + self.datemaster + '_SLC' + output_sufix)
        if not os.path.isfile(imagebasename+'.dim'):
            time1 = datetime.now()
            for k in range(len(input3)):
                outputtemp = os.path.join(imagebasename + '_' + input3[k] + '.dim')
                #OVERWRITE CONDITION FOR SUBSWATH GENERATION
                try:
                    p = subprocess.check_output([self.pathgpt, gptsubset_master, '-Pinput1='+input1, '-Pinput2='+input3[k], '-Pinput3='+input4, '-Ptarget1='+outputtemp])
                except:
                    p = ''
            #After processing subswaths, they are put together with TOPS Merge Workflow
            subswath_list = sorted(glob.glob(os.path.join(imagebasename + '*.dim')))
            if len(subswath_list)>1:
                #OVERWRITE CONDITION FOR MERGE
                m = self.TOPS_Merge_subswaths(subswath_list, imagebasename + '.dim', self.pathgpt)
                for filePath in subswath_list:
                    try:
                        os.remove(filePath)
                        shutil.rmtree(os.path.splitext(filePath)[0]+'.data')
                    except:
                        print("Error while deleting file : ", filePath)
            else:
                #Change output name from IWx to 'merged' in case AOI covers only one subswath
                shutil.move(subswath_list[0], imagebasename+'.dim')
                shutil.move(os.path.splitext(subswath_list[0])[0]+'.data', imagebasename+'.data')
            outputtimes.append((datetime.now()-time1).seconds/60)
        else:
            outputtimes.append('0')
        datelist.append(self.datemaster)
        outputfiles.append(imagebasename + '.dim')
        #Fill processing info in df
        self.processdf['Image_date'] = datelist
        self.processdf['Outputfiles_align'+sufix] = outputfiles
        self.processdf['Processing_minutes'+sufix] = outputtimes
        self.processdf.to_csv(os.path.join(self.diroutorder, self.logfilename), sep=';', index=False)
        
    def stack_polariz (self, imagepathlist, outputfile, pathgpt):
        if len(imagepathlist)==2:
            gptxml_file = os.path.join(self.DirProject, 'RES', 'stack_dualpolariz.xml')
            input1 = imagepathlist[0]
            input2 = imagepathlist[1]
            try:
                m = subprocess.check_output([pathgpt, gptxml_file, '-Pinput1='+input1, '-Pinput2='+input2, '-Ptarget1='+outputfile])
            except:
                m = ''
        return m
    
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
        gpt_calibration = os.path.join(self.DirProject, 'RES', 'calibration_outputcomplex.xml')
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
        sufix=''
        self.run_alignmentlist(pairspath, pairsdate, sufix)
        #Calibration of original products
        outputfiles = []
        outputtimes = []
        for i in range(len(imagelist)):
            time1 = datetime.now()
            outfile = os.path.join(self.diroutorder, '00_calib', imagelist['ImageDate'][i]+'_calib.dim')
            try:
                p = subprocess.check_output([self.pathgpt, gpt_calibration, '-Pinput1='+imagelist['Image'][i], '-Ptarget1='+outfile])
            except:
                p = ''
            outputfiles.append(outfile)
            outputtimes.append((datetime.now()-time1).seconds/60)
        self.processdf['Outputfiles_calib'] = outputfiles
        self.processdf['Processing_minutes_calib'] = outputtimes
        #TODO LANZAR PROCESO COREGISTR SOBRE CALIB, MONTAR PARES, ALMACENAR RUTAS Y TIEMPOS
        pairspath_calib = mountpairs(pd.concat([self.processdf[imagelist['ImageDate'].str.contains(self.datemaster)], self.processdf.drop(self.processdf[imagelist['ImageDate'].str.contains(self.datemaster)].index)], sort=False).reset_index(), 'Outputfiles_calib')
        sufix = '_calib'
        self.run_alignmentlist(pairspath_calib, pairsdate, sufix)

    def generate_baselinelist (self):
        self.baselinelist = self.processdf.copy()
        #Master record removed (no baseline between the same image as master and slave)
        self.baselinelist = self.baselinelist[:-1]
        self.baselinelist['Master'] = self.datemaster
        self.baselinelist['Slave'] = self.processdf['Image_date']
        self.baselinelist['Perp_Baseline'] = [self.searchMetabytag(self.processdf['Outputfiles_align'][i], 'Perp Baseline')[1] for i in range(len(self.baselinelist))]
        self.baselinelist['Temp_Baseline'] = [self.searchMetabytag(self.processdf['Outputfiles_align'][i], 'Temp Baseline')[1] for i in range(len(self.baselinelist))]
        self.baselinelist = self.baselinelist.drop(['Image_date', 'Outputfiles_align', 'Processing_minutes'], axis=1)
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

    # def generate_interferogram(self):
    #     gptxml_file = os.path.join(self.DirProject, 'RES', 'Interferogram_Multilook.xml')
    #     #for i in range(len(self.baselinefiltered)):
    #     #Fill output files if pair is already coregistered
    #     for i in range(len(self.processdf)):
    #         for j in range(len(self.baselinefiltered)):
    #             if self.processdf['Image_date'][i] == self.baselinefiltered['Slave'][j] and self.baselinefiltered['Master'][j] == self.datemaster:
    #                 self.baselinefiltered.loc[j, 'Outputfiles_align'] = self.processdf['Outputfiles_align'][i]

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
        subdirout = '02_ml'
        output_sufix = '_ML_' + self.azLooks + '_' + self.rgLooks
        outputfiles = []
        outputtimes = []
        for i in range(len(self.processdf['Outputfiles_align'])):
            if os.path.isfile(self.processdf['Outputfiles_align'][i]):
                inputfile = self.processdf['Outputfiles_align'][i]
                outputfile = os.path.join(self.diroutorder, subdirout, os.path.splitext(os.path.basename(inputfile))[0] + output_sufix + '.dim')
                if not os.path.isfile(outputfile) or self.overwrite == '1':
                    outputtimes.append(self.applyMultilook(inputfile, outputfile, self.azLooks, self.rgLooks))
                else:
                    outputtimes.append('0')
                outputfiles.append(outputfile)
            else:
                continue
        self.processdf['Outputfiles_ML'] = outputfiles
        self.processdf['Processing_minutes_ML'] = outputtimes
        self.processdf.to_csv(os.path.join(self.diroutorder, 'process_log.csv'), sep=';', index=False)

    def applyGeocoding(self, inputfile, outputdir, epsg, taglist):
        time1 = datetime.now()
        prod = ProductIO.readProduct(inputfile)
        lat = prod.getTiePointGrid('latitude')
        lon = prod.getTiePointGrid('longitude')
        w = prod.getSceneRasterWidth()
        h = prod.getSceneRasterHeight()
        
        array = np.zeros((w, h), dtype=np.float32)
        latpixels = lat.readPixels(0, 0, w, h, array)
        lat_arr = np.asarray(latpixels)
        lat_arr.shape = (h, w)
        lonpixels = lon.readPixels(0, 0, w, h, array)
        lon_arr = np.asarray(lonpixels)
        lon_arr.shape = (h, w)
        
        if not os.path.exists(outputdir):
            os.makedirs(outputdir)
        np.save(os.path.join(outputdir, 'lat'), lat_arr)
        np.save(os.path.join(outputdir, 'lon'), lon_arr)
        #np.savetxt(os.path.join(outputdir, 'lat'), lat_arr, delimiter=';')
        #np.savetxt(os.path.join(outputdir, 'lon'), lon_arr, delimiter=';')
        
        inProj = Proj('epsg:4326')
        outProj = Proj('epsg:'+epsg)
        x, y = transform(inProj,outProj,lat_arr,lon_arr)
        np.save(os.path.join(outputdir, 'x'), x)
        np.save(os.path.join(outputdir, 'y'), y)
        #Check for bands for geocoding and resampling
        bandlist = []
        try:
            for tag in taglist:
                for band in list(prod.getBandNames()):
                    if tag in band:
                        bandlist.append(band)
        except:
            print('No tags')
        if bandlist is not None:
            x1 = np.linspace(np.min(x), np.max(x), w)
            y1 = np.linspace(np.min(y), np.max(y), h)
            grid_x, grid_y = np.meshgrid(x1,y1)
            for band in bandlist:
                data_getband = prod.getBand(band)
                data_pixels = data_getband.readPixels(0, 0, w, h, array)
                data_arr = np.asarray(data_pixels)
                grid_data = griddata((x.flatten(), y.flatten()), data_arr.flatten(), (grid_x, grid_y), method='nearest')
                np.save(os.path.join(outputdir, band), grid_data)
        return((datetime.now()-time1).seconds/60)
        
    def Geocoding(self):
        subdirout = '03_gc'
        output_sufix = '_GC'
        outputfiles = []
        outputtimes = []
        for i in range(len(self.processdf['Outputfiles_ML'])):
            if os.path.isfile(self.processdf['Outputfiles_ML'][i]):
                inputfile = self.processdf['Outputfiles_ML'][i]
                outputdir = os.path.join(self.diroutorder, subdirout, os.path.splitext(os.path.basename(inputfile))[0])
                if not os.path.isdir(outputdir) or self.overwrite == '1':
                    outputtimes.append(self.applyGeocoding(inputfile, outputdir, self.epsg, self.taglist))
                else:
                    outputtimes.append('0')
                outputfiles.append(outputdir)
            else:
                continue
        self.processdf['Outputfiles_GC'] = outputfiles
        self.processdf['Processing_minutes_GC'] = outputtimes
        self.processdf.to_csv(os.path.join(self.diroutorder, 'process_log.csv'), sep=';', index=False)



    # def do_calibration(source, polarization, pols):
    #     print('\tCalibration...')
    #     parameters = HashMap()
    #     parameters.put('outputSigmaBand', True)
    #     if polarization == 'DH':
    #         parameters.put('sourceBands', 'Intensity_HH,Intensity_HV')
    #     elif polarization == 'DV':
    #         parameters.put('sourceBands', 'Intensity_VH,Intensity_VV')
    #     elif polarization == 'SH' or polarization == 'HH':
    #         parameters.put('sourceBands', 'Intensity_HH')
    #     elif polarization == 'SV':
    #         parameters.put('sourceBands', 'Intensity_VV')
    #     else:
    #         print("different polarization!")
    #     parameters.put('selectedPolarisations', pols)
    #     parameters.put('outputImageScaleInDb', False)
    #     output = GPF.createProduct("Calibration", parameters, source)
    #     return output