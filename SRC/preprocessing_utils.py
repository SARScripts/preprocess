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
from osgeo import gdal, osr
from pyproj import Proj, transform, Transformer, CRS
from joblib import Parallel, delayed
import multiprocessing

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
        self.calib = str(calib)
        if self.calib == '1':
            self.sufix = '_calib'
        else:
            self.sufix = ''
        self.subdiroutslc = '01_slc' + self.sufix
        self.subdiroutifg = '01_ifg' + self.sufix
        self.subdiroutml = '02_ml'
        self.subdiroutgc = '03_gc'
        self.pathgpt = pathgpt
        self.snappypath = snappypath
        self.DirProject = DirProject
        self.processdf = pd.DataFrame()
        self.baselinelist = pd.DataFrame()
        self.maxbaseperp = ''
        self.maxbasetemp = ''
        self.baselinefiltered = pd.DataFrame()
        self.logfilename = 'process_log.csv'
        self.epsgin = '4326'
        #self.epsgout = 
        #self.spatialres = 10
        #Number of cores fits the ratio between RAM memory and 23 Gb/core (need for coregistration)
        self.num_cores = 7

        if not os.path.exists(self.diroutorder):
            os.makedirs(self.diroutorder)            
        sys.path.append(self.snappypath)
        import snappy
        from snappy import ProductIO, GPF, HashMap

    def searchMetabytag(self, metafilepath, tag):
        if pd.isnull(metafilepath):
            meta = '0'
        else:
            meta = []
            tree = ET.parse(metafilepath)
            root = tree.getroot()
            for atr in root.iter('MDATTR'):
                if atr.get('name') == tag:
                    meta.append(atr.text)
        return meta


    def mpUnzipAssembly(self, imagelist, ImageDate, date):
        time1 = datetime.now()
        if ImageDate.count(date) == 1:
            imagefile = imagelist.iloc[imagelist[imagelist['Date']==date].index[0], imagelist.columns.get_loc('Image')]
            outputpath = os.path.join(self.diroutputzip, os.path.splitext(os.path.basename(imagefile))[0]+'.SAFE', 'manifest.safe')
            if not os.path.isfile(outputpath) or self.overwrite == '1':
                if os.path.splitext(imagefile)[1] == '.zip':
                    #OVERWRITE CONDITION FOR UNZIP FILES
                    unzipS1 (imagefile, self.diroutputzip)
                    inputimage = imagefile
                    timeproc = (datetime.now()-time1).seconds/60
                elif os.path.splitext(os.path.basename(imagefile))[1] == '.safe':
                    if os.path.isdir(os.path.join(self.diroutputzip, os.path.basename(os.path.dirname(imagefile)))):
                        shutil.copytree(os.path.dirname(imagefile), os.path.join(self.diroutputzip, os.path.basename(os.path.dirname(imagefile))))
                    outputpath = os.path.join(self.diroutputzip, os.path.basename(os.path.dirname(imagefile)), 'manifest.safe')
                    inputimage = imagefile
                    timeproc = (datetime.now()-time1).seconds/60
                else:
                    print('Image format not supported: ' + imagelist.Image[i])
                    sys.exit()
            else:
                inputimage = imagefile
                timeproc = '0'
                
        elif ImageDate.count(date) == 2:
            #Slice-assembly for joining products
            assemblist = imagelist[imagelist.Image.str.contains(date)]['Image'].tolist()
            outputpath = os.path.join(self.diroutputzip, date, date+'_assembl.dim')
            if not os.path.isfile(outputpath) or self.overwrite == '1':
                subprocess.check_output([self.pathgpt, self.gptxml_unzip, '-Pinput1='+assemblist[0], '-Pinput2='+assemblist[1], '-Ptarget1='+outputpath])
            timeproc = (datetime.now()-time1).seconds/60
            inputimage = assemblist
        else:
            print('More than 2 images (not allowed) with date:', date)
        return (inputimage, outputpath, timeproc)

    # def unzipassemblygetmeta (self, imagelist):
    #     #Unzip if zipfiles, slice assembly if two products for the same date in imagelist and fill up ImagePreprocess, ImageName, ImageDate columns
    #     self.gptxml_unzip = os.path.join(self.DirProject, 'RES', 'sliceassembly.xml')
    #     ImagePreprocess = []
    #     ImageName = []
    #     ImageDate = []
    #     ImageOrigin = []
    #     TimeZipAssem = []
    #     self.diroutputzip = os.path.join(self.diroutorder, '00_data')
    #     if not os.path.exists(self.diroutputzip):
    #         os.makedirs(self.diroutputzip)
    #     #Extract dates from filenames
    #     for i in range(len(imagelist)):
    #         if os.path.splitext(imagelist.Image[i])[1] == '.zip':
    #             name = os.path.splitext(os.path.basename(imagelist.Image[i]))[0]
    #         elif os.path.splitext(imagelist.Image[i])[1] == '.safe':
    #             name = os.path.basename(os.path.dirname(imagelist.Image[i])).split('.')[0]
    #         else:
    #             print('Image format not supported: ' + imagelist.Image[i])
    #             sys.exit()
    #         ImageDate.append(name[name.find('T')-8:name.find('T')])
    #     imagelist['Date'] = ImageDate
    #     uniquedates = list(set(ImageDate))
    #     num_cores = multiprocessing.cpu_count()-1
    #     df_unzipassemb = pd.DataFrame(Parallel(n_jobs=num_cores)(delayed(self.mp_unzip_assembly)(imagelist, ImageDate, uniquedates[i]) for i in range(len(uniquedates))))
    #     imagelist['ImageOrigin'] = df_unzipassemb[0]
    #     imagelist['ImagePreprocess'] = df_unzipassemb[1]
    #     imagelist['TimeunzipAsseml'] = df_unzipassemb[2]
    #     imagelist['ImageDate'] = uniquedates
    #     return imagelist

    def pair_coregistration(self, paths, dates, ortholatlon=0):
        self.gptxml_coreg = os.path.join(self.DirProject, 'RES', 'TOPSAR_Coreg_Interferogram_ESD_Polariz_wkt_outimage_flatearth_topophase.xml')
        gptsubset_master = os.path.join(self.DirProject, 'RES', 'applyAOI_master_convertDIMAP.xml')
        output_sufix = '_Coregistered'
        nocorrection_tag = '_nocorrect'
        if ortholatlon == 0:
            self.tempdir = os.path.join(self.diroutorder, self.subdiroutslc, 'temp')
            self.subdiroutslc = '01_slc' + self.sufix
            self.subdiroutifg = '01_ifg' + self.sufix
        else:
            self.tempdir = os.path.join(self.diroutorder, self.subdiroutgc, 'temp')
            self.subdiroutslc, self.subdiroutifg = self.subdiroutgc, self.subdiroutgc              
        if not os.path.isdir(self.tempdir):
            os.makedirs(self.tempdir)
        lonmin = self.AOI[0]
        latmax = self.AOI[1]
        lonmax = self.AOI[2]
        latmin = self.AOI[3]
        outputfiles = []
        outputtimes = []
        outputifg1 = []
        outputifg2 = []
        datelist = []
        time1 = datetime.now()
        input1 = paths[0]
        input2 = paths[1]
        input3 = ['IW1', 'IW2', 'IW3']
        if any (self.AOI):
            input4 = 'POLYGON (('+lonmin+' '+latmin+', '+lonmax+' '+latmin+', '+lonmax+' '+latmax+', '+lonmin+' '+latmax+', '+lonmin+' '+latmin+'))'
        else:
            input4 = ''
        input5 = ['VV', 'VH']
        if ortholatlon == 0:
            input6 = 'false'
        else:
            input6 = 'true'
        #The process is executed for every subswath
        if os.path.splitext(input2)[1] == '.safe' or os.path.splitext(input2)[1] == '.dim':
            slavename = dates[1] + '_SLC' + self.sufix
        else:
            print('Image format not supported: ' + input2)
            sys.exit()
        imagebasename = dates[0] + '_' + slavename + output_sufix
        if not os.path.isfile(os.path.join(self.diroutorder, self.subdiroutslc, imagebasename + '.dim')) or self.overwrite == '1':
            for polariz in input5:
                output1 = imagebasename + '_' + polariz
                for k in range(len(input3)):
                    outputfile1 = os.path.join(self.diroutorder, self.subdiroutifg, output1 + '_' + input3[k] + '.dim')
                    outputfile2 = os.path.join(self.diroutorder, self.subdiroutifg, output1 + '_' + input3[k] + nocorrection_tag + '.dim')
                    outputfile3 = os.path.join(self.diroutorder, self.subdiroutslc, output1 + '_' + input3[k] + 'im.dim')
                    #OVERWRITE CONDITION FOR SUBSWATH GENERATION
                    try:
                        p = subprocess.check_output([self.pathgpt, self.gptxml_coreg, '-Pinput1='+input1, '-Pinput2='+input2, '-Pinput3='+input3[k], '-Pinput4='+input4, '-Pinput5='+polariz,'-Pinput6='+input6,'-Ptarget1='+outputfile1,'-Ptarget2='+outputfile2 ,'-Ptarget3='+outputfile3])
                    except:
                        p = ''
                #After processing subswaths, they are put together with TOPS Merge Workflow
                #List of subswaths by checking the output files. The output is added at the end of the list
                subswath_list1 = sorted(glob(os.path.join(self.diroutorder, self.subdiroutifg, output1 + '*' + nocorrection_tag + '*.dim')))
                subswath_list2 = sorted(glob(os.path.join(self.diroutorder, self.subdiroutifg, output1 + '*.dim')))
                subswath_list3 = sorted(glob(os.path.join(self.diroutorder, self.subdiroutslc, output1 + '*.dim')))
                for f in subswath_list1:
                    subswath_list2.remove(f)
                #Append output filename to the list for merging
                subswath_list1.append(os.path.join(self.diroutorder, self.subdiroutifg, output1 + nocorrection_tag))
                subswath_list2.append(os.path.join(self.diroutorder, self.subdiroutifg, output1))
                subswath_list3.append(os.path.join(self.diroutorder, self.subdiroutslc, output1+'im'))
                #Function TOPS Merge is only called when more than 1 subswath is generated. Only IW products supported (max. 3 subswaths)
                for listmerge in [subswath_list1, subswath_list2, subswath_list3]:
                    if len(listmerge)<2:
                        print('AOI does not match scene extents')
                        sys.exit()
                    elif len(listmerge)==2:
                        #Change output name from IWx to 'merged' in case AOI covers only one subswath
                        shutil.move(listmerge[0], listmerge[-1]+'.dim')
                        shutil.move(os.path.splitext(listmerge[0])[0]+'.data', os.path.splitext(listmerge[-1])[0]+'.data')
                    else:
                        m = self.TOPS_Merge_subswaths(listmerge[:-1], listmerge[-1] + '.dim', self.pathgpt)
                        #Delete/Move temp files
                        del listmerge[-1]    #Remove last element of list (output filename)
                        for filePath in listmerge:
                            if os.path.isfile(os.path.join(self.tempdir, os.path.basename(filePath))):
                                os.remove(os.path.join(self.tempdir, os.path.basename(filePath)))    
                            shutil.move(filePath, os.path.join(self.tempdir, os.path.basename(filePath)))
                            if os.path.isdir(os.path.join(self.tempdir, os.path.splitext(os.path.basename(filePath))[0]+'.data')):
                                shutil.rmtree(os.path.join(self.tempdir, os.path.splitext(os.path.basename(filePath))[0]+'.data'))    
                            shutil.move(os.path.splitext(filePath)[0]+'.data', os.path.join(self.tempdir, os.path.splitext(os.path.basename(filePath))[0]+'.data'))
            polariz_list1 = sorted(glob(os.path.join(self.diroutorder, self.subdiroutifg, imagebasename) + '*V*nocorrect*.dim'))
            polariz_list2 = sorted(glob(os.path.join(self.diroutorder, self.subdiroutifg, imagebasename) + '*V*.dim'))
            polariz_list3 = sorted(glob(os.path.join(self.diroutorder, self.subdiroutslc, imagebasename) + '*V*.dim'))
            for f in polariz_list1:
                polariz_list2.remove(f) 
            for listpol in [polariz_list1, polariz_list2, polariz_list3]:
                if listpol==polariz_list1:
                    outputfile = os.path.join(self.diroutorder, self.subdiroutifg, imagebasename) + nocorrection_tag + '.dim'
                if listpol == polariz_list2:
                    outputfile = os.path.join(self.diroutorder, self.subdiroutifg, imagebasename) + '.dim'
                if listpol == polariz_list3:
                    outputfile = os.path.join(self.diroutorder, self.subdiroutslc, imagebasename) + '.dim'
                if len(listpol)>1:
                    #Stack polarizations (dual)
                    m = self.stack_polariz(listpol, outputfile, self.pathgpt)
                    #Delete/Move temp files
                    for filePath in listpol:
                        try:
                            shutil.move(filePath, os.path.join(self.tempdir, os.path.basename(filePath)))
                            shutil.move(os.path.splitext(filePath)[0]+'.data', os.path.join(self.tempdir, os.path.splitext(os.path.basename(filePath))[0]+'.data'))
                        except:
                            print("Error while deleting file : ", filePath)
                else:
                    #Change output name in case the product has only one polarization
                    shutil.move(listpol[0], outputfile)
                    shutil.move(os.path.splitext(listpol[0])[0]+'.data', os.path.splitext(outputfile)+'.data')
        if os.path.isfile(os.path.join(self.diroutorder, self.subdiroutslc, imagebasename + '.dim')):
            outputfiles.append(os.path.join(self.diroutorder, self.subdiroutslc, imagebasename + '.dim'))
            outputifg1.append(os.path.join(self.diroutorder, self.subdiroutifg, imagebasename + '.dim'))
            outputifg2.append(os.path.join(self.diroutorder, self.subdiroutifg, imagebasename + '_nocorrect.dim'))
        else:
            outputfiles.append('Not generated')
        outputtimes.append((datetime.now()-time1).seconds/60)
        datelist.append(dates[1])
        return [outputfiles, outputifg1, outputifg2, datelist, outputtimes]
        
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
    
    def coregistration_ifg (self):
        #Read the image file
        imagelist = pd.read_csv(self.imagelistfile, sep=';', decimal=',')
        #Unzip if zipfiles, slice assembly if two products for the same date in imagelist and fill up ImagePreprocess, ImageName, ImageDate columns
        self.gptxml_unzip = os.path.join(self.DirProject, 'RES', 'sliceassembly.xml')
        ImagePreprocess = []
        ImageName = []
        ImageDate = []
        ImageOrigin = []
        TimeZipAssem = []
        self.diroutputzip = os.path.join(self.diroutorder, '00_data')
        if not os.path.exists(self.diroutputzip):
            os.makedirs(self.diroutputzip)
        #Extract dates from filenames
        for i in range(len(imagelist)):
            if os.path.splitext(imagelist.Image[i])[1] == '.zip':
                name = os.path.splitext(os.path.basename(imagelist.Image[i]))[0]
            elif os.path.splitext(imagelist.Image[i])[1] == '.safe':
                name = os.path.basename(os.path.dirname(imagelist.Image[i])).split('.')[0]
            else:
                print('Image format not supported: ' + imagelist.Image[i])
                sys.exit()
            ImageDate.append(name[name.find('T')-8:name.find('T')])
        imagelist['Date'] = ImageDate
        uniquedates = list(set(ImageDate))
        df_unzipassemb = pd.DataFrame(Parallel(n_jobs=self.num_cores)(delayed(self.mpUnzipAssembly)(imagelist, ImageDate, uniquedates[i]) for i in range(len(uniquedates))))
        self.processdf['ImageOrigin'] = df_unzipassemb[0]
        self.processdf['ImagePreprocess'] = df_unzipassemb[1]
        self.processdf['TimeunzipAsseml'] = df_unzipassemb[2]
        self.processdf['ImageDate'] = uniquedates
        #Calibration of original products
        gpt_calibration = os.path.join(self.DirProject, 'RES', 'calibration_outputcomplex.xml')
        outputfiles = []
        outputtimes = []
        if self.calib == '1':
            for i in range(len(self.processdf)):
                outfile = os.path.join(self.diroutorder, '00'+sufix, self.processdf['ImageDate'][i]+sufix+'.dim')
                time1 = datetime.now()
                if not os.path.isfile(outfile) or self.overwrite == '1':
                    try:
                        p = subself.check_output([self.pathgpt, gpt_calibration, '-Pinput1='+self.processdf['ImagePreprocess'][i], '-Ptarget1='+outfile])
                    except:
                        p = ''
                outputfiles.append(outfile)
                outputtimes.append(str((datetime.now()-time1).seconds/60))
            self.processdf['Outputfiles_calib'] = outputfiles
            self.processdf['Processing_minutes_calib'] = outputtimes
        else:
            self.sufix =''
            outputfiles = self.processdf['ImagePreprocess']
            outputtimes = ['0'] * len(self.processdf)
        self.processdf['Outputfiles'] = outputfiles
        self.processdf['Processing_minutes'] = outputtimes
        self.processdf['ImageDate'] = pd.to_datetime(self.processdf['ImageDate'])
        self.processdf = self.processdf.sort_values(by=['ImageDate'])
        #IF MASTER IMAGE DATE MATCHES ANY DATE OF THE FILE LIST, RUN ALIGNMENT PROCESSING (OTHERWISE MASTER DATE IS COMPUTED AS THE 'MID POINT' DATE OF THE LIST)
        #Identify master image row and set up slave images dataframe
        #if not any (self.processdf['ImageDate'].str.contains(self.datemaster)):
        if not any (self.processdf['ImageDate'] == self.datemaster):
            datelist = list(self.processdf.ImageDate.sort_values().astype(int)/10**9)
            self.datemaster = pd.to_datetime(datelist[min(range(len(datelist)), key = lambda i: abs(datelist[i]-(datelist[-1]+datelist[0])/2))], unit='s', origin='unix').strftime('%Y%m%d')
        self.processdf['ImageDate_str'] = [date_obj.strftime('%Y%m%d') for date_obj in self.processdf['ImageDate']]
        #Master image is reindexed to the top of the df
        masterdf = self.processdf[self.processdf['ImageDate_str'] == self.datemaster]
        slavedf = self.processdf.copy()
        #slavedf = slavedf.drop(slavedf[self.processdf['ImageDate'] == datetime.strptime(self.datemaster, '%Y%m%d')].index)
        pairspath = mountpairs(pd.concat([masterdf['Outputfiles'], slavedf['Outputfiles']], sort=False).reset_index(), 'Outputfiles')
        pairsdate = mountpairs(pd.concat([masterdf, slavedf], sort=False).reset_index(), 'ImageDate_str')
#        [outputfiles, outputifg1, outputifg2, datelist, outputtimes] = Parallel(n_jobs=num_cores)(delayed(self.pair_coregistration)(pairspath[i], pairsdate[i]) for i in range(len(pairspath)))
        df_coreg = pd.DataFrame(Parallel(n_jobs=self.num_cores)(delayed(self.pair_coregistration)(pairspath[i], pairsdate[i]) for i in range(len(pairspath))))             
        outputfiles = [item[0] for item in df_coreg[0]]
        datelist = [item[0] for item in df_coreg[3]]
        outputtimes = [item[0] for item in df_coreg[4]]
        #Get position of master and add master record
        self.processdf.reset_index(drop=True, inplace=True) 
        masterid = self.processdf[self.processdf['ImageDate_str']==self.datemaster].index[0]
        #datelist.insert(masterid, self.datemaster)
        #outputfiles.insert(masterid, os.path.join(self.diroutorder, '01_slc/'+self.datemaster+'_'+self.datemaster+'_SLC_Coregistered.dim'))
        #outputtimes.insert(masterid, '0')
        self.processdf['Outputfiles_align'+self.sufix] = outputfiles
        self.cleaning_coreg()
        #Generate reference matrices (lat/lon) for further geocoding
        self.generateGeocoding()
        #Process differential phase (interferograms with and without flat Earth and topographic corrections)
        #phaselist = [os.path.splitext(f)[0] + '.data' for f in df_coreg[1]]
        #phasedifflist = [os.path.splitext(f)[0] + '.data' for f in df_coreg[2]]
        #ifgpairs = [(phaselist[j], phasedifflist[j]) for j in range(len(phaselist))]
        #self.differential_phase(ifgpairs)
        #Fill processing info into df
        self.processdf['Outputfiles_align'+self.sufix] = outputfiles
        self.processdf['Processing_minutes_align'] = outputtimes
        self.processdf.to_csv(os.path.join(self.diroutorder, self.logfilename), sep=';', index=False)
        msg = 'Coregistration generated'
        return msg

    def cleaning_coreg(self):
        #Go through the coregistered images and change file names, remove master image...
        for imagecoreg in self.processdf['Outputfiles_align'+self.sufix]:
            imageref = os.path.splitext(os.path.basename(imagecoreg))[0]
            imagedir = os.path.splitext(imagecoreg)[0] + '.data/'
            ifgdir = os.path.join(self.diroutorder, self.subdiroutifg, imageref+'.data')
            ifgnocorrectdir = os.path.join(self.diroutorder, self.subdiroutifg, imageref+'_nocorrect.data')
            deletelist = glob(os.path.join(imagedir, '*mst*'))
            for file in deletelist:
                os.remove(file)
            changenamelist = glob(os.path.join(imagedir, '*slv*'))
            for file in changenamelist:
                os.rename(file, file.replace('_slv1', ''))
                filetemp = file.replace('_slv1', '')
                os.rename(filetemp, filetemp.replace(filetemp[filetemp.find('_slv3'):filetemp.find('_slv3')+15], ''))
            for fold in [ifgdir, ifgnocorrectdir]:
                for file in glob(os.path.join(fold, '*.*')):
                    filetemp = file
                    if 'slv' in file:
                        filetemp = file.replace(file[file.find('_slv'):file.find('_slv')+15], '')
                    filedest = filetemp.replace(filetemp[filetemp.find('_ifg_')+8:filetemp.find('_ifg_')+18], '')
                    os.rename(file, filedest)
                    
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
        output_sufix = '_ML_' + self.azLooks + '_' + self.rgLooks
        outputfiles = []
        outputtimes = []
        msg = []
        #If calibration has been applied, ML over calibrated products, if not, ML over non-calibrated products
        if 'Outputfiles_align_calib' in list(self.processdf):
            filesML = 'Outputfiles_align' + self.sufix
        else:
            filesML = 'Outputfiles_align'
        for i in range(len(self.processdf[filesML])):
            if os.path.isfile(self.processdf[filesML][i]):
                inputfile = self.processdf[filesML][i]
                outputfile = os.path.join(self.diroutorder, self.subdiroutml, os.path.splitext(os.path.basename(inputfile))[0] + output_sufix + '.dim')
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
    
    def applydatelistGeocoding(self, gcdatelist):
        #Check if coregistration geocoding has been generated
        if not (os.path.isfile(os.path.join(self.diroutorder, self.subdiroutgc, 'orthorectifiedLat_ML.img')) and os.path.isfile(os.path.join(self.diroutorder, self.subdiroutgc, 'orthorectifiedLon_ML.img'))):
            self.generateGeocoding()
        outputdir = os.path.join(self.diroutorder, self.subdiroutgc)
        for date in gcdatelist:
            inputfilegc = self.processdf[self.processdf['ImageDate_str'] == date]['Outputfiles_ML'][0]
            #outputdirgc = os.path.join(outputdir, gcdatelist[i]+'_ifg_gc.dim')
            timegc = self.applyGeocoding(inputfilegc, outputdir, date)
        
    def applyGeocoding(self, inputfile, outputdir, date):
        time1 = datetime.now()
        #Load data
        #prod = ProductIO.readProduct(inputfile)
        inputfolder = os.path.splitext(inputfile)[0] + '.data'
        lat = gdal.Open(os.path.join(self.diroutorder, self.subdiroutgc, 'orthorectifiedLat_ML.img'))
        lon = gdal.Open(os.path.join(self.diroutorder, self.subdiroutgc, 'orthorectifiedLon_ML.img'))
        #lon = prod.getBand('orthorectifiedLon')
        #w = prod.getSceneRasterWidth()
        #h = prod.getSceneRasterHeight()
        #array = np.zeros((w, h), dtype=np.float32)
        lat_arr = lat.ReadAsArray()
        #lat_arr = np.asarray(latpixels)
        #lat_arr.shape = (h, w)
        #lonpixels = lon.readPixels(0, 0, w, h, array)
        lon_arr = lon.ReadAsArray()
        #lon_arr.shape = (h, w)


        #Check for bands for geocoding and resampling
        bandlist = glob(inputfolder+'/*.img')
        # try:
        #     for tag in taglist:
        #         for band in list(prod.getBandNames()):
        #             if tag in band:
        #                 bandlist.append(band)
        # except:
        #     print('No tags')
        if bandlist is not None:
            #Adjust to the Sentinel 2 grid
            ulxgrid, ulygrid, lrxgrid, lrygrid = self.check_S2grid(np.min(self.x), np.max(self.y), np.max(self.x), np.min(self.y))
            wout = int((lrxgrid - ulxgrid) / self.spatialres)
            hout = int((ulygrid - lrygrid) / self.spatialres)
            x1 = np.linspace(lrxgrid, ulxgrid, wout)
            y1 = np.linspace(lrygrid, ulygrid, hout)
            grid_x, grid_y = np.meshgrid(x1,y1)
            for band in bandlist:
                data_getband = gdal.Open(band)
                #data_pixels = data_getband.readPixels(0, 0, w, h, array)
                data_arr = data_getband.ReadAsArray()
                grid_data = griddata((self.x.flatten(), self.y.flatten()), data_arr.flatten(), (grid_x, grid_y), method='nearest')
                np.save(os.path.join(outputdir, date, os.path.basename(band)), grid_data)
        return((datetime.now()-time1).seconds/60)
        
    def check_S2grid(self, ulx, uly, lrx, lry):
        #Extract x,y reference coordinates for each zone
        if not (hasattr(self, 's2gridrefx') and hasattr(self, 's2gridrefy')):
            s2gridpath = os.path.join(self.DirProject, 'RES', 'tabularize_s2_footprint.csv')
            if not os.path.isfile(s2gridpath):
                print('Sentinel grid file not found in /RES')
                sys.exit()
            #load s2gridpath
            s2grid = pd.read_csv(s2gridpath, sep=',')
            epsglist = set(s2grid['epsg'])
            #Extract unique zones
            s2grid['zone_ref'] = [Proj(s2grid['utm_zone'][i]).crs.utm_zone for i in range(len(s2grid))]
            if all(s2grid.groupby('zone_ref')['ul_x'].apply(lambda x: len(set(x%float(self.spatialres)))==1)):
                self.s2gridrefx = {ep:s2grid[s2grid['epsg']==ep].iloc[0]['ul_x'] for ep in epsglist}
                self.s2gridrefy = {ep:s2grid[s2grid['epsg']==ep].iloc[0]['ul_y'] for ep in epsglist}
            else:
                print('Sentinel grid file provided not coherent with coordinates ')
                sys.exit()
        if self.s2gridrefx[int(self.epsgout)] > ulx:
            ulxgrid = self.s2gridrefx[int(self.epsgout)] - (int((self.s2gridrefx[float(self.epsgout)] - ulx) / self.spatialres)*self.spatialres) + self.spatialres
        else:
            ulxgrid = self.s2gridrefx[int(self.epsgout)] - int((self.s2gridrefx[float(self.epsgout)] - ulx) / self.spatialres)*self.spatialres
        if self.s2gridrefy[int(self.epsgout)] > uly:
            ulygrid = self.s2gridrefy[int(self.epsgout)] - int((self.s2gridrefy[float(self.epsgout)] - uly) / self.spatialres)*self.spatialres
        else:
            ulygrid = self.s2gridrefy[int(self.epsgout)] - (int((self.s2gridrefy[float(self.epsgout)] - uly) / self.spatialres)*self.spatialres - self.spatialres)
    
        if self.s2gridrefx[int(self.epsgout)] >= lrx:
            lrxgrid = self.s2gridrefx[int(self.epsgout)] - int((self.s2gridrefx[float(self.epsgout)] - lrx) / self.spatialres)*self.spatialres
        else:
            lrxgrid = self.s2gridrefx[int(self.epsgout)] - (int((self.s2gridrefx[float(self.epsgout)] - lrx) / self.spatialres)*self.spatialres - self.spatialres)
        if self.s2gridrefy[int(self.epsgout)] > lry:
            lrygrid = self.s2gridrefy[int(self.epsgout)] - (int((self.s2gridrefy[float(self.epsgout)] - lry) / self.spatialres)*self.spatialres) + self.spatialres
        else:
            lrygrid = self.s2gridrefy[int(self.epsgout)] - int((self.s2gridrefy[float(self.epsgout)] - lry) / self.spatialres)*self.spatialres
        return (ulxgrid, ulygrid, lrxgrid, lrygrid)

    def backup_filenames(file_list):
        for file in file_list:
            shutil.move(file, os.path.join(os.path.splitext(file)[0]+'_bckp'+os.path.splitext(file)[1]))
    
    def mergeortholatlon(self, inputfile, outputdir):
        subswath_list = sorted(glob(os.path.join(self.tempdir, os.path.splitext(inputfile)[0] + '*VH_IW[1-3].dim')))
        #shutil.copytree(subswath_list, outputdir)
        refdir = []
        for subswath in subswath_list:
            refdir.append(os.path.join(outputdir, os.path.splitext(os.path.basename(subswath))[0]+'.data'))
            if os.path.exists(refdir[-1]):
                shutil.rmtree(refdir[-1])
            shutil.copytree(os.path.splitext(subswath)[0]+'.data', refdir[-1])
            shutil.copy(subswath, os.path.join(outputdir, os.path.basename(subswath)))
            i_list = glob(os.path.join(refdir[-1], '*i_ifg*'))
            q_list = glob(os.path.join(refdir[-1], '*q_ifg*'))
            backup_filenames(i_list)
            backup_filenames(q_list)
            for el in i_list:
                shutil.copy(os.path.join(refdir[-1], 'orthorectifiedLat'+el[-4:]), el)
            for el in q_list:
                shutil.copy(os.path.join(refdir[-1], 'orthorectifiedLon'+el[-4:]), el)
        self.TOPS_Merge_subswaths (refdir+'.dim', os.path.join(outputdir, 'merge'), self.pathgpt)
            
    def generateGeocoding(self):
        from snappy import ProductIO, GPF, HashMap
        #Generate the baseline list file if it is not generated yet.
        self.generate_baselinelist()
        gptxml_file = os.path.join(self.DirProject, 'RES', 'ifg_ortholatlon.xml')
        indexGC = int(self.baselinefiltered[['Perp_Baseline']].idxmax())
        inputfile1 = self.processdf[self.processdf['ImageDate_str']==self.baselinefiltered.loc[indexGC]['Master']]['Outputfiles_ML'].tolist()[0]
        inputfile2 = self.processdf[self.processdf['ImageDate_str']==self.baselinefiltered.loc[indexGC]['Slave']]['Outputfiles_ML'].tolist()[0]
        date1 = self.baselinefiltered.loc[indexGC]['Master']
        date2 = self.baselinefiltered.loc[indexGC]['Slave']
        outputdir = os.path.join(self.diroutorder, self.subdiroutgc)
        outputfile = os.path.join(outputdir, 'ifg_gc_ML.dim')
        subprocess.check_output([self.pathgpt, gptxml_file, '-Pinput1='+inputfile1, '-Pinput2='+inputfile2, '-Ptarget1='+outputfile])
        shutil.copy(os.path.join(os.path.splitext(outputfile)[0]+'.data', 'orthorectifiedLon.img'), os.path.join(outputdir, 'orthorectifiedLon_ML.img'))
        shutil.copy(os.path.join(os.path.splitext(outputfile)[0]+'.data', 'orthorectifiedLon.hdr'), os.path.join(outputdir, 'orthorectifiedLon_ML.hdr'))
        shutil.copy(os.path.join(os.path.splitext(outputfile)[0]+'.data', 'orthorectifiedLat.img'), os.path.join(outputdir, 'orthorectifiedLat_ML.img'))
        shutil.copy(os.path.join(os.path.splitext(outputfile)[0]+'.data', 'orthorectifiedLat.hdr'), os.path.join(outputdir, 'orthorectifiedLat_ML.hdr'))
        #If ortho latlon numpy matrices needed
        prod = ProductIO.readProduct(outputfile)
        ortho_lat = prod.getBand('orthorectifiedLat')
        ortho_lon = prod.getBand('orthorectifiedLon')
        w = prod.getSceneRasterWidth()
        h = prod.getSceneRasterHeight()
        array = np.zeros((w, h), dtype=np.float32)
        latpix = ortho_lat.readPixels(0, 0, w, h, array)
        lonpix = ortho_lon.readPixels(0, 0, w, h, array)
        np.save(os.path.join(outputdir, 'lat_ML'), latpix)
        np.save(os.path.join(outputdir, 'lon_ML'), lonpix)
        #Management of CRS, by default UTM, zone according to longitude
        try:
            print(self.epsgout)
        except AttributeError:
            zone = int((np.mean(lonpix) + 180) / 6) + 1
            if np.mean(latpix)>0:
                south = False
            else:
                south = True
            crs = CRS.from_dict({'proj': 'utm', 'zone': zone, 'south': south})
            self.epsgout = crs.to_authority()[1]
        #Resampling
        inProj = 'epsg:'+self.epsgin
        outProj = 'epsg:'+self.epsgout
        transformer = Transformer.from_crs(inProj, outProj)
        self.x, self.y = transformer.transform(lonpix,latpix)
        np.save(os.path.join(outputdir, 'x_ML'), self.x)
        np.save(os.path.join(outputdir, 'y_ML'), self.y)
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
                complex_mat = re_result + im_result * 1j
                phase_result = np.angle(complex_mat)
                self.array2raster(os.path.join(ifgpairs[k][0], 'ifg_result' + polariz) + '.img', (0, 0), phase_result.shape[1], phase_result.shape[0], phase_result)
        
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

