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
    def __init__(self, process_type, imagelistfile, diroutorder, overwrite, datemaster, AOI, pathgpt, DirProject):
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
        self.DirProject = DirProject
        self.processdf = pd.DataFrame()
        self.baselinelist = pd.DataFrame()
        self.maxbaseperp = ''
        self.maxbasetemp = ''
        self.baselinefiltered = pd.DataFrame()
        
        if not os.path.exists(self.diroutorder):
            os.makedirs(self.diroutorder)

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
                if not os.path.isdir(os.path.join(diroutput, os.path.basename(os.path.dirname(imagelist.Image[i])))):
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
        gptxml_file = os.path.join(self.DirProject, 'RES', 'TOPSAR_Coreg_Interferogram_ESD_args_wkt.xml')
        #gptxml_file = os.path.join('/media/guadarrama/PROYECTOS/SAR2CUBE/Git/RES', 'TOPSAR_Coreg_Interferogram_ESD_args_wkt.xml')
        
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
            if any (self.AOI):
                input4 = 'POLYGON (('+lonmin+' '+latmin+', '+lonmax+' '+latmin+', '+lonmax+' '+latmax+', '+lonmin+' '+latmax+', '+lonmin+' '+latmin+'))'
            else:
                input4 = ''
            #The process is executed for every subswath
            for k in range(len(input3)):
                if os.path.splitext(input2)[1] == '.safe':
                    slavename = pairsdate[i][1] + '_SLC'
                else:
                    print('Image format not supported: ' + input2)
                    sys.exit()
                output1 = os.path.join(self.diroutorder, subdirout, self.datemaster + '_' + slavename + output_sufix)
                outputfile = os.path.join(self.diroutorder, subdirout, self.datemaster + '_' + slavename + output_sufix + '_' + input3[k] + '.dim')
                #OVERWRITE CONDITION FOR SUBSWATH GENERATION
                if not os.path.isfile(outputfile) or self.overwrite == '1':
                    try:
                        p = subprocess.check_output([self.pathgpt, gptxml_file, '-Pinput1='+input1, '-Pinput2='+input2, '-Pinput3='+input3[k], '-Pinput4='+input4,'-Ptarget1='+outputfile])
                    except:
                        p = ''
            #After processing subswaths, they are put together with TOPS Merge Workflow
            #List of subswaths by checking the output files
            subswath_list = sorted(glob.glob(os.path.join(output1 + '*.dim')))
            #Function TOPS Merge is only called when more than 1 subswath is generated. Only IW products supported (max. 3 subswaths)
            if len(subswath_list)>1:
                #OVERWRITE CONDITION FOR MERGE
                if not os.path.isfile(output1 + '.dim') or self.overwrite == '1':
                    m = self.TOPS_Merge_subswaths(subswath_list, output1, self.pathgpt)
                    for filePath in subswath_list:
                        try:
                            os.remove(filePath)
                            shutil.rmtree(os.path.splitext(filePath)[0]+'.data')
                        except:
                            print("Error while deleting file : ", filePath)
            else:
                #Change output name from IWx to 'merged' in case AOI touches only one subswath
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
        self.processdf['Outputfiles'] = outputfiles
        self.processdf['Processing_minutes'] = outputtimes
        
    def TOPS_Merge_subswaths (self, imagepathlist, outputfile, pathgpt):
        #gptxml_file = '/media/guadarrama/PROYECTOS/SAR2CUBE/WP2000/GPT_command_line/merge_3subswaths.xml'
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

        #TO-DO CHECK IF IMAGES ARE TYPE SLC

        #IF MASTER IMAGE DATE MATCHES ANY DATE OF THE FILE LIST, ONE ROUND PROCESSING (OTHERWISE MASTER DATE IS COMPUTED AS THE 'MEDIAN' DATE OF THE LIST)
        #Identify master image row and set up slave images dataframe
        if any (imagelist['ImageDate'].str.contains(self.datemaster)):
            #Master image is reindexed to the top of the df
            masterdf = imagelist[imagelist['ImageDate'].str.contains(self.datemaster)]
            slavedf = imagelist.copy()
            slavedf = slavedf.drop(slavedf[imagelist['ImageDate'].str.contains(self.datemaster)].index)
            pairspath = mountpairs(pd.concat([masterdf, slavedf], sort=False).reset_index(), 'ImageUnzip')
            pairsdate = mountpairs(pd.concat([masterdf, slavedf], sort=False).reset_index(), 'ImageDate')
            self.run_alignmentlist(pairspath, pairsdate)
        else:
            datelist = list(imagelist.ImageDate.sort_values().astype(int))
            self.datemaster = str(datelist[min(range(len(datelist)), key = lambda i: abs(datelist[i]-(datelist[-1]+datelist[0])/2))])
            masterdf = imagelist[imagelist['ImageDate'].str.contains(self.datemaster)]
            slavedf = imagelist.copy()
            slavedf = slavedf.drop(slavedf[imagelist['ImageDate'].str.contains(self.datemaster)].index)
            pairspath = mountpairs(pd.concat([masterdf, slavedf], sort=False).reset_index(), 'ImageUnzip')
            pairsdate = mountpairs(pd.concat([masterdf, slavedf], sort=False).reset_index(), 'ImageDate')
            self.run_alignmentlist(pairspath, pairsdate)
        self.processdf.to_csv(os.path.join(self.diroutorder, 'process.csv'), sep=';', index=False)
        
    def generate_baselinelist (self):
        #TO-DO CHECK IF PROCESSDF AND COREGISTERED DATA EXIST
        if self.processdf is None:
            print('Process database not generated, run alignment over a stack of images first')
            sys.exit()
            #TODO Check si existen las imágenes alineadas
        else:
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
                    #self.baselinelist.append(pd.Series([str(self.baselinelist['Slave'][0]), str(self.baselinelist['Slave'][i+1]), perp, temp]), ignore_index=True)
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
            #Fill output files if pair is already coregistered
            for i in range(len(self.processdf)):
                for j in range(len(self.baselinefiltered)):
                    if self.processdf['Image_date'][i] == self.baselinefiltered['Slave'][j] and self.baselinefiltered['Master'][j] == self.datemaster:
                        self.baselinefiltered.loc[j, 'Outputfiles'] = self.processdf['Outputfiles'][i]
            self.baselinefiltered.to_csv(os.path.join(self.diroutorder, 'baselines.csv'), sep=';', index=False)
    
    def generate_interferogram(self):
        print('entra generate interf')
        gptxml_file = os.path.join(self.DirProject, 'RES', 'Interferogram_Multilook.xml')
        for i in range(len(self.baselinefiltered)):
            pass
            

        # output_sufix = '_Aligned'
        # subdirout = '01_slc'
        # extent = self.AOI
        # # (lonmin, latmax, lonmax, latmin)
        # lonmin = extent[0]
        # latmax = extent[1]
        # lonmax = extent[2]
        # latmin = extent[3]
        # outputfiles = []
        # outputtimes = []
        # datelist = []
        # for i in range(len(pairspath)):
        #     time1 = datetime.now()
        #     input1 = pairspath[i][0]
        #     input2 = pairspath[i][1]
        #     if '_IW_SLC_' in input1:
        #         input3 = ['IW1', 'IW2', 'IW3']
        #     if any (self.AOI):
        #         input4 = 'POLYGON (('+lonmin+' '+latmin+', '+lonmax+' '+latmin+', '+lonmax+' '+latmax+', '+lonmin+' '+latmax+', '+lonmin+' '+latmin+'))'
        #     else:
        #         input4 = ''
        #     #The process is executed for every subswath
        #     for k in range(len(input3)):
        #         #Check file type input for assign filename output (valid formats: .safe, .zip)
        #         if os.path.splitext(input2)[1] == '.safe':
        #             slavename = pairsdate[i][1] + '_SLC'
        #         else:
        #             print('Image format not supported: ' + input2)
        #             sys.exit()
        #         output1 = os.path.join(self.diroutorder, subdirout, slavename + output_sufix + '_' + input3[k])
        #         try:
        #             p = subprocess.check_output([self.pathgpt, gptxml_file, '-Pinput1='+input1, '-Pinput2='+input2, '-Pinput3='+input3[k], '-Pinput4='+input4,'-Ptarget1='+output1])
        #         except:
        #             p = ''
        #     #After processing subswaths, they are put together with TOPS Merge Workflow
        #     #List of subswaths by checking the output files
        #     subswath_list = sorted(glob.glob(os.path.join(self.diroutorder, subdirout, slavename + output_sufix + '*.dim')))
        #     #Function TOPS Merge is only called when more than 1 subswath is generated. Only IW products supported (max. 3 subswaths)
        #     if (len(subswath_list)>1):
        #         output_Merge = os.path.join(self.diroutorder, subdirout, slavename + output_sufix + '_MrgIW')
        #         m = self.TOPS_Merge_subswaths(subswath_list, output_Merge, self.pathgpt)
        #     if os.path.isfile(output_Merge+'.dim'):
        #         outputfiles.append(output_Merge+'.dim')
        #     else:
        #         outputfiles.append('Not generated')
        #     outputtimes.append((datetime.now()-time1).seconds/60)
        #     datelist.append(pairsdate[i][1])
        #     #TO-DO REMOVE TEMP FILES? (IW1, IW2...)
        # #Fill processing info in df
        # self.processdf['Image_date'] = datelist
        # self.processdf['Outputfiles'] = outputfiles
        # self.processdf['Processing_minutes'] = outputtimes
