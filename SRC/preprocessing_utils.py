#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 17:32:28 2020

"""


import sys
import os
import subprocess
import pandas as pd
from datetime import datetime 
from glob import glob
import xml.etree.ElementTree as ET
from zipfile import ZipFile
import shutil
import numpy as np
from scipy.interpolate import griddata
import gdal, osr
from pyproj import Proj, transform

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
    def __init__(self, process_type, imagelistfile, diroutorder, overwrite, datemaster, AOI, calib, pathgpt, snappypath, DirProject):
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
        self.calib = calib
        self.pathgpt = pathgpt
        self.snappypath = snappypath
        self.DirProject = DirProject
        self.processdf = pd.DataFrame()
        self.baselinelist = pd.DataFrame()
        self.maxbaseperp = ''
        self.maxbasetemp = ''
        self.baselinefiltered = pd.DataFrame()
        self.logfilename = 'process_log.csv'
        self.epsg = '4326'
        
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
        diroutput = os.path.join(self.diroutorder, '00_data')
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
        gptxml_file = os.path.join(self.DirProject, 'RES', 'TOPSAR_Coreg_Interferogram_ESD_Polariz_wkt_outimage_flatearth_topophase.xml')
        gptsubset_master = os.path.join(self.DirProject, 'RES', 'applyAOI_master_convertDIMAP.xml')

        #OUTPUT SUFIX
        output_sufix = '_Coregistered'
        subdiroutimag = '01_slc' + sufix
        subdiroutifg = '01_ifg' + sufix
        nocorrection_tag = '_nocorrect'
        tempdir = os.path.join(self.diroutorder, subdiroutimag, 'temp')
        if not os.path.exists(tempdir):
            os.makedirs(tempdir)
        lonmin = self.AOI[0]
        latmax = self.AOI[1]
        lonmax = self.AOI[2]
        latmin = self.AOI[3]
        outputfiles = []
        outputtimes = []
        outputifg1 = []
        outputifg2 = []
        datelist = []
        for i in range(len(pairspath)):
            time1 = datetime.now()
            input1 = pairspath[i][0]
            input2 = pairspath[i][1]
            input3 = ['IW1', 'IW2', 'IW3']
            if any (self.AOI):
                input4 = 'POLYGON (('+lonmin+' '+latmin+', '+lonmax+' '+latmin+', '+lonmax+' '+latmax+', '+lonmin+' '+latmax+', '+lonmin+' '+latmin+'))'
            else:
                input4 = ''
            input5 = ['VV', 'VH']
            #The process is executed for every subswath
            if os.path.splitext(input2)[1] == '.safe' or os.path.splitext(input2)[1] == '.dim':
                slavename = pairsdate[i][1] + '_SLC' + sufix
            else:
                print('Image format not supported: ' + input2)
                sys.exit()
            imagebasename = self.datemaster + '_' + slavename + output_sufix
            if not os.path.isfile(os.path.join(self.diroutorder, subdiroutimag, imagebasename + '.dim')) or self.overwrite == '1':
                for polariz in input5:
                    output1 = imagebasename + '_' + polariz
                    for k in range(len(input3)):
                        outputfile1 = os.path.join(self.diroutorder, subdiroutifg, output1 + '_' + input3[k] + '.dim')
                        outputfile2 = os.path.join(self.diroutorder, subdiroutifg, output1 + '_' + input3[k] + nocorrection_tag + '.dim')
                        outputfile3 = os.path.join(self.diroutorder, subdiroutimag, output1 + '_' + input3[k] + 'im.dim')
                        #OVERWRITE CONDITION FOR SUBSWATH GENERATION
                        try:
                            p = subprocess.check_output([self.pathgpt, gptxml_file, '-Pinput1='+input1, '-Pinput2='+input2, '-Pinput3='+input3[k], '-Pinput4='+input4, '-Pinput5='+polariz,'-Ptarget1='+outputfile1,'-Ptarget2='+outputfile2 ,'-Ptarget3='+outputfile3])
                        except:
                            p = ''
                    #After processing subswaths, they are put together with TOPS Merge Workflow
                    #List of subswaths by checking the output files. The output is added at the end of the list
                    subswath_list1 = sorted(glob(os.path.join(self.diroutorder, subdiroutifg, output1 + '*' + nocorrection_tag + '*.dim')))
                    subswath_list2 = sorted(glob(os.path.join(self.diroutorder, subdiroutifg, output1 + '*.dim')))
                    subswath_list3 = sorted(glob(os.path.join(self.diroutorder, subdiroutimag, output1 + '*.dim')))

                    for f in subswath_list1:
                        subswath_list2.remove(f)
                    #Append output filename to the list for merging
                    subswath_list1.append(os.path.join(self.diroutorder, subdiroutifg, output1 + nocorrection_tag))
                    subswath_list2.append(os.path.join(self.diroutorder, subdiroutifg, output1))
                    subswath_list3.append(os.path.join(self.diroutorder, subdiroutimag, output1+'im'))

                    #Function TOPS Merge is only called when more than 1 subswath is generated. Only IW products supported (max. 3 subswaths)
                    for listmerge in [subswath_list1, subswath_list2, subswath_list3]:
                        if len(listmerge)==1:
                            #Change output name from IWx to 'merged' in case AOI covers only one subswath
                            shutil.move(listmerge[0], listmerge[-1])
                            shutil.move(os.path.splitext(listmerge[0])[0]+'.data', os.path.join(tempdir, os.path.splitext(os.path.basename(filePath))[0]+'.data'))
                        else:
                            #OVERWRITE CONDITION FOR MERGE
                            m = self.TOPS_Merge_subswaths(listmerge[:-1], listmerge[-1] + '.dim', self.pathgpt)
                            #Delete/Move temp files
                            del listmerge[-1]    #Remove last element of list (output filename)
                            for filePath in listmerge:
                                shutil.move(filePath, os.path.join(tempdir, os.path.basename(filePath)))
                                shutil.move(os.path.splitext(filePath)[0]+'.data', os.path.join(tempdir, os.path.splitext(os.path.basename(filePath))[0]+'.data'))
                polariz_list1 = sorted(glob(os.path.join(self.diroutorder, subdiroutifg, imagebasename) + '*V*nocorrect*.dim'))
                polariz_list2 = sorted(glob(os.path.join(self.diroutorder, subdiroutifg, imagebasename) + '*V*.dim'))
                polariz_list3 = sorted(glob(os.path.join(self.diroutorder, subdiroutimag, imagebasename) + '*V*.dim'))
                for f in polariz_list1:
                    polariz_list2.remove(f) 
                for listpol in [polariz_list1, polariz_list2, polariz_list3]:
                    if listpol==polariz_list1:
                        outputfile = os.path.join(self.diroutorder, subdiroutifg, imagebasename) + nocorrection_tag + '.dim'
                    if listpol == polariz_list2:
                        outputfile = os.path.join(self.diroutorder, subdiroutifg, imagebasename) + '.dim'
                    if listpol == polariz_list3:
                        outputfile = os.path.join(self.diroutorder, subdiroutimag, imagebasename) + '.dim'
                    if len(listpol)>1:
                        #Stack polarizations (dual)
                        m = self.stack_polariz(listpol, outputfile, self.pathgpt)
                        #Delete/Move temp files
                        for filePath in listpol:
                            try:
                                shutil.move(filePath, os.path.join(tempdir, os.path.basename(filePath)))
                                shutil.move(os.path.splitext(filePath)[0]+'.data', os.path.join(tempdir, os.path.splitext(os.path.basename(filePath))[0]+'.data'))
                            except:
                                print("Error while deleting file : ", filePath)
                    else:
                        #Change output name from IWx to 'merged' in case AOI covers only one subswath
                        shutil.move(listpol[0], outputfile)
                        shutil.move(os.path.splitext(listpol[0])[0]+'.data', os.path.splitext(outputfile)+'.data')
            if os.path.isfile(os.path.join(self.diroutorder, subdiroutimag, imagebasename + '.dim')):
                outputfiles.append(os.path.join(self.diroutorder, subdiroutimag, imagebasename + '.dim'))
                outputifg1.append(os.path.join(self.diroutorder, subdiroutifg, imagebasename + '.dim'))
                outputifg2.append(os.path.join(self.diroutorder, subdiroutifg, imagebasename + '_nocorrect.dim'))

            else:
                outputfiles.append('Not generated')
            outputtimes.append((datetime.now()-time1).seconds/60)
            datelist.append(pairsdate[i][1])
        #Subset master image by AOI and convert to DIMAP format for next processing, fill processing table.
        imagebasename = os.path.join(self.diroutorder, subdiroutimag, self.datemaster + '_' + self.datemaster + '_SLC' + sufix + output_sufix)
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
            subswath_list = sorted(glob(os.path.join(imagebasename + '*.dim')))
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
            timemaster = (datetime.now()-time1).seconds/60
        else:
            timemaster = '0'
        
        #Process differential phase (interferograms with and without flat Earth and topographic corrections)
        phaselist = [os.path.splitext(f)[0] + '.data' for f in outputifg1]
        phasedifflist = [os.path.splitext(f)[0] + '.data' for f in outputifg2]
        ifgpairs = [(phaselist[j], phasedifflist[j]) for j in range(len(phaselist))]
        self.differential_phase(ifgpairs)
        
        
        #Get position of master
        masterid = self.processdf[self.processdf['ImageDate_str']==self.datemaster].index[0]
        datelist.insert(masterid, self.datemaster)
        outputfiles.insert(masterid, imagebasename + '.dim')
        outputtimes.insert(masterid, timemaster)

        #Fill processing info into df
        self.processdf['Outputfiles_align'+sufix] = outputfiles
        self.processdf['Processing_minutes_align'] = outputtimes
        self.processdf.to_csv(os.path.join(self.diroutorder, self.logfilename), sep=';', index=False)
        #shutil.rmtree(tempdir)
        
        #Generate reference matrices (lat/lon) for further geocoding
        self.Geocoding(self.processdf['Outputfiles_align'+sufix], '')
        
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
        
        #Calibration of original products
        outputfiles = []
        outputtimes = []
        if self.calib == '1':
            sufix = '_calib'
            for i in range(len(imagelist)):
                outfile = os.path.join(self.diroutorder, '00'+sufix, imagelist['ImageDate'][i]+sufix+'.dim')
                time1 = datetime.now()
                if not os.path.isfile(outfile) or self.overwrite == '1':
                    try:
                        p = subprocess.check_output([self.pathgpt, gpt_calibration, '-Pinput1='+imagelist['Image'][i], '-Ptarget1='+outfile])
                    except:
                        p = ''
                outputfiles.append(outfile)
                outputtimes.append(str((datetime.now()-time1).seconds/60))
            self.processdf['Outputfiles_calib'] = outputfiles
            self.processdf['Processing_minutes_calib'] = outputtimes
        else:
            sufix =''
            outputfiles = imagelist['ImageUnzip']
            outputtimes = ['0'] * len(imagelist)
        self.processdf['Outputfiles'] = outputfiles
        self.processdf['Processing_minutes'] = outputtimes
        self.processdf['ImageDate'] = pd.to_datetime(imagelist['ImageDate'])
        self.processdf = self.processdf.sort_values(by=['ImageDate'])
        #IF MASTER IMAGE DATE MATCHES ANY DATE OF THE FILE LIST, RUN ALIGNMENT PROCESSING (OTHERWISE MASTER DATE IS COMPUTED AS THE 'MID POINT' DATE OF THE LIST)
        #Identify master image row and set up slave images dataframe
        if not any (imagelist['ImageDate'].str.contains(self.datemaster)):
            datelist = list(imagelist.ImageDate.sort_values().astype(int))
            self.datemaster = str(datelist[min(range(len(datelist)), key = lambda i: abs(datelist[i]-(datelist[-1]+datelist[0])/2))])
        self.processdf['ImageDate_str'] = [date_obj.strftime('%Y%m%d') for date_obj in self.processdf['ImageDate']]
        #Master image is reindexed to the top of the df
        masterdf = self.processdf[self.processdf['ImageDate'] == datetime.strptime(self.datemaster, '%Y%m%d')]
        slavedf = self.processdf.copy()
        slavedf = slavedf.drop(slavedf[self.processdf['ImageDate'] == datetime.strptime(self.datemaster, '%Y%m%d')].index)
        
        pairspath = mountpairs(pd.concat([masterdf['Outputfiles'], slavedf['Outputfiles']], sort=False).reset_index(), 'Outputfiles')
        pairsdate = mountpairs(pd.concat([masterdf, slavedf], sort=False).reset_index(), 'ImageDate_str')
        self.run_alignmentlist(pairspath, pairsdate, sufix)
        msg = 'Coregistering generated'
        return msg

    def generate_baselinelist (self):
        self.baselinelist = self.processdf.copy()
        #Master record removed (no baseline between the same image as master and slave)
        self.baselinelist = self.baselinelist.drop(self.baselinelist[self.processdf['ImageDate'] == datetime.strptime(self.datemaster, '%Y%m%d')].index)
        self.baselinelist.reset_index(drop=True, inplace=True) 
        self.baselinelist['Master'] = self.datemaster
        self.baselinelist['Slave'] = self.baselinelist['ImageDate_str']
        self.baselinelist['Perp_Baseline'] = [self.searchMetabytag(self.baselinelist['Outputfiles_align'][i], 'Perp Baseline')[1] for i in range(len(self.baselinelist))]
        self.baselinelist['Temp_Baseline'] = [self.searchMetabytag(self.baselinelist['Outputfiles_align'][i], 'Temp Baseline')[1] for i in range(len(self.baselinelist))]
        fieldlist = []
        #Delete fields not needed
        for tag in ['ImageDate', 'Outputfiles', 'Processing_minutes']:
            for band in list(self.baselinelist):
                if tag in band:
                    fieldlist.append(band)
        self.baselinelist = self.baselinelist.drop(fieldlist, axis=1)
        waninglist = self.baselinelist.copy()
        while len(waninglist)>1:
            for i in range(len(waninglist)-1):
                perp = float(waninglist['Perp_Baseline'][0]) - float(waninglist['Perp_Baseline'][i+1])
                temp = float(waninglist['Temp_Baseline'][0]) - float(waninglist['Temp_Baseline'][i+1])
                self.baselinelist.loc[len(self.baselinelist)] = [str(waninglist['Slave'][0]), str(waninglist['Slave'][i+1]), perp, temp]
            waninglist = waninglist.drop(0)
            waninglist.reset_index(drop=True, inplace=True)
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
        msg = []
        #If calibration has been applied, ML over calibrated products, if not, ML over non-calibrated products
        if 'Outputfiles_align_calib' in list(self.processdf):
            filesML = 'Outputfiles_align_calib'
        else:
            filesML = 'Outputfiles_align'
        for i in range(len(self.processdf[filesML])):
            if os.path.isfile(self.processdf[filesML][i]):
                inputfile = self.processdf[filesML][i]
                outputfile = os.path.join(self.diroutorder, subdirout, os.path.splitext(os.path.basename(inputfile))[0] + output_sufix + '.dim')
                if not os.path.isfile(outputfile) or self.overwrite == '1':
                    outputtimes.append(self.applyMultilook(inputfile, outputfile, self.azLooks, self.rgLooks))
                else:
                    outputtimes.append('0')
                outputfiles.append(outputfile)
            else:
                msg.append = 'Multilook not generated file ' + self.processdf[filesML][i] + ' not found'
                continue
        self.processdf['Outputfiles_ML'] = outputfiles
        self.processdf['Processing_minutes_ML'] = outputtimes
        self.processdf.to_csv(os.path.join(self.diroutorder, 'process_log.csv'), sep=';', index=False)
        msg = 'Multilook generated'
        return msg
    
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
        
    def applyTerrainCorrection(self, inputfile, outputfile, pixelspacing, epsg):
        from snappy import ProductIO, GPF, HashMap
        time1 = datetime.now()
        parameters = HashMap()
        prod = ProductIO.readProduct(inputfile)
        parameters.put('demResamplingMethod', 'NEAREST_NEIGHBOUR') 
        parameters.put('imgResamplingMethod', 'NEAREST_NEIGHBOUR') 
        parameters.put('demName', 'SRTM 3Sec') 
        parameters.put('mapProjection', epsg)
        parameters.put('saveDEM', True)
        parameters.put('nodataValueAtSea', False)
        parameters.put('saveLatLon', True)
        parameters.put('saveIncidenceAngleFromEllipsoid', True)
        parameters.put('saveLocalIncidenceAngle', True)
        parameters.put('saveProjectedLocalIncidenceAngle', True)
        if pixelspacing == '':
            prodTC = GPF.createProduct("Terrain-Correction", parameters, prod) 
        else:
            parameters.put('pixelSpacingInMeter', spatialres) 
            prodTC = GPF.createProduct("Terrain-Correction", parameters, prod) 
        ProductIO.writeProduct(prodTC, outputfile, 'BEAM-DIMAP')
        return datetime.now() - time1

    def Geocoding(self, listinputfiles, pixelspacing):
        #Generate the baseline list file if it is not generated yet.
        #if not os.path.isfile(os.path.join(self.diroutorder, 'baselines.csv')):
        self.generate_baselinelist()
        subdirout = '03_gc'
        outputfiles = [''] * len(self.processdf)
        outputtimes = ['0'] * len(self.processdf)
        indexGC = int(self.baselinefiltered[['Perp_Baseline']].idxmax())
        inputfile = listinputfiles[indexGC]
        if os.path.isfile(inputfile):
            outputfile = os.path.join(self.diroutorder, subdirout, os.path.splitext(os.path.basename(inputfile))[0])
            if not os.path.isdir(outputfile) or self.overwrite == '1':
                outputtimes[indexGC] = self.applyTerrainCorrection(inputfile, outputfile, pixelspacing, self.epsg)
            outputfiles[indexGC] = outputfile
        else:
            msg = 'Geocoding not generated in file ' + listinputfiles[indexGC] + ' not found'
            return msg
        self.processdf['Outputfiles_GC'] = outputfiles
        self.processdf['Processing_minutes_GC'] = outputtimes
        self.processdf.to_csv(os.path.join(self.diroutorder, 'process_log.csv'), sep=';', index=False)
        msg = 'Geocoding generated'
        return msg

    def compute_phase(self, inputdir, polariz, outputfile):
            realband = gdal.Open(glob(os.path.join(inputdir, '*i_ifg' + '*' + polariz + '*.img'))[0])
            imagband = gdal.Open(glob(os.path.join(inputdir, 'q_ifg' + '*' + polariz + '*.img'))[0])
            realarr = realband.ReadAsArray()
            imagarr = imagband.ReadAsArray()
            complex_mat = realarr + imagarr * 1j
            phase = np.angle(complex_mat)
            self.array2raster(outputfile + '_'+ polariz + '.img', (0, 0), phase.shape[1], phase.shape[0], phase)
            return phase
    
    def differential_phase(self, ifgpairs):
        for k in range(len(ifgpairs)):
            polarizations = ['VV', 'VH']
            for polariz in polarizations:
                phasecomplete = self.compute_phase(ifgpairs[k][0], polariz, os.path.join(ifgpairs[k][0], 'ifg_compl'))
                phasediff = self.compute_phase(ifgpairs[k][1], polariz, os.path.join(ifgpairs[k][0], 'ifg_diff'))
                phasesubtr = phasediff-phasecomplete
                #eit= cost+isint     https://www.math.wisc.edu/~angenent/Free-Lecture-Notes/freecomplexnumbers.pdf
                #cmath.exp(a))    https://www.askpython.com/python/python-complex-numbers
                re_result = np.cos(phasesubtr)
                im_result = np.sin(phasesubtr)
                # phase_result = 
                self.array2raster(os.path.join(ifgpairs[k][0], 'ifg_result' + polariz) + '.img', (0, 0), phasesubtr.shape[1], phasesubtr.shape[0], phasesubtr)
        
    def array2raster(self, newRasterfn, rasterOrigin, pixelWidth, pixelHeight, array):
        cols = array.shape[1]
        rows = array.shape[0]
        originX = rasterOrigin[0]
        originY = rasterOrigin[1]
        driver = gdal.GetDriverByName('ENVI')
        outRaster = driver.Create(newRasterfn, cols, rows, 1, gdal.GDT_Float32)
        outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
        outband = outRaster.GetRasterBand(1)
        outband.WriteArray(array)
        outRasterSRS = osr.SpatialReference()
        outRasterSRS.ImportFromEPSG(4326)
        outRaster.SetProjection(outRasterSRS.ExportToWkt())
        outband.FlushCache()

