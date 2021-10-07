#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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
import math


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
    def __init__(self, process_type, imagelistfile, diroutorder, overwrite, datemaster, AOI, calib, pathgpt, snappypath, snaphupath, snaphuconf, DirProject, num_cores):
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
        if self.calib != '0':
            self.sufix = '_calib'
        else:
            self.sufix = ''
        self.subdiroutslc = '01_slc' + self.sufix
        self.subdiroutifg = '01_ifg' + self.sufix
        self.subdiroutml = '02_ml'
        self.subdiroutgc = '03_gc'
        self.pathgpt = pathgpt
        self.snappypath = snappypath
        self.snaphupath = snaphupath
        self.snaphuconf = snaphuconf
        self.DirProject = DirProject
        self.processdf = pd.DataFrame()
        self.baselinelist = pd.DataFrame()
        self.maxbaseperp = ''
        self.maxbasetemp = ''
        self.baselinefiltered = pd.DataFrame()
        self.logfilename = 'process_log.csv'
        self.epsgin = '4326'
        self.polarizphaseunwrap = 'VV'
        #Number of cores should fit the ratio between RAM memory and 23 Gb/core (need for coregistration)
        self.num_cores = num_cores
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
        print("[*] Running mpUnzipAssembly...")
        time1 = datetime.now()
        if ImageDate.count(date) == 1:
            imagefile = imagelist.iloc[imagelist[imagelist['Date']==date].index[0], imagelist.columns.get_loc('Image')]
            if os.path.splitext(imagefile)[1] == '.zip':
                print("[*] Unzipping Sentinel data from file " + imagefile)
                outputpath = os.path.join(self.diroutputzip, os.path.splitext(os.path.basename(imagefile))[0]+'.SAFE', 'manifest.safe')
                if not os.path.isfile(outputpath) or self.overwrite == '1':
                    unzipS1 (imagefile, self.diroutputzip)
                    inputimage = imagefile
                    timeproc = (datetime.now()-time1).seconds/60
                else:
                    inputimage = imagefile
                    timeproc = '0'
            elif os.path.splitext(imagefile)[1] == '.safe':
                outputpath = os.path.join(self.diroutputzip, os.path.basename(os.path.normpath(os.path.dirname(imagefile))), 'manifest.safe')
                if not os.path.isfile(outputpath) or self.overwrite == '1':
                    shutil.copytree(os.path.dirname(imagefile), os.path.join(self.diroutputzip, os.path.basename(os.path.dirname(imagefile))))
                    timeproc = (datetime.now()-time1).seconds/60
                else:
                    timeproc = '0'
                #outputpath = os.path.join(self.diroutputzip, os.path.basename(os.path.dirname(imagefile)), 'manifest.safe')
                inputimage = imagefile
            else:
                print('Image format not supported: ' + imagelist.Image[i])
                sys.exit()
        elif ImageDate.count(date) == 2:
            #Slice-assembly for joining products
            assemblist = imagelist[imagelist.Image.str.contains(date)]['Image'].tolist()
            outputpath = os.path.join(self.diroutputzip, date, date+'_assembl.dim')
            if not os.path.isfile(outputpath) or self.overwrite == '1':
                try:
                    p = subprocess.check_output([self.pathgpt, self.gptxml_unzip, '-Pinput1='+assemblist[0], '-Pinput2='+assemblist[1], '-Ptarget1='+outputpath])
                except Exception as e:
                    raise e
#                 finally:
#                     p.terminate() # send sigterm, or ...
#                     p.kill()      # send sigkill    
            timeproc = (datetime.now()-time1).seconds/60
            inputimage = assemblist
        else:
            print('More than 2 images (not allowed) with date:', date)
        return (inputimage, outputpath, timeproc)

    def pair_coregistration(self, paths, dates, ortholatlon=0):
        self.gptxml_coreg = os.path.join(self.DirProject, 'RES', 'TOPSAR_Coreg_ifg_ESD_Polariz_wkt.xml')
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
            #Check and delete temporary files 
            searchfolders = [os.path.join(self.diroutorder, self.subdiroutslc), os.path.join(self.diroutorder, self.subdiroutifg), self.tempdir] 
            for sfold in searchfolders:
                trashfiles = sorted(glob(sfold + '/*' + imagebasename + '*'))
                for f in trashfiles:
                    if os.path.isfile(f):
                        os.remove(f)
                    else:
                        shutil.rmtree(f)
            #Go through polarizations
            for polariz in input5:
                output1 = imagebasename + '_' + polariz
                #Generate coregistration for every subswath
                for k in range(len(input3)):
                    outputfile1 = os.path.join(self.diroutorder, self.subdiroutifg, output1 + '_' + input3[k] + '.dim')
                    outputfile2 = os.path.join(self.diroutorder, self.subdiroutifg, output1 + '_' + input3[k] + nocorrection_tag + '.dim')
                    outputfile3 = os.path.join(self.diroutorder, self.subdiroutslc, output1 + '_' + input3[k] + 'im.dim')
                    try:
                        p = subprocess.check_output([self.pathgpt, self.gptxml_coreg, '-Pinput1='+input1, '-Pinput2='+input2, '-Pinput3='+input3[k], '-Pinput4='+input4, '-Pinput5='+polariz,'-Pinput6='+input6,'-Ptarget1='+outputfile1,'-Ptarget2='+outputfile2 ,'-Ptarget3='+outputfile3])
                    except Exception as e:
                        print(e)
#                     finally:
#                         p.terminate() # send sigterm, or ...
#                         p.kill()      # send sigkill    
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
            #Compute differential phase
            self.differential_phase([os.path.join(self.diroutorder, self.subdiroutifg, imagebasename + '.data'), os.path.join(self.diroutorder, self.subdiroutifg, imagebasename + '_nocorrect.data')]) 
        else:
            #Compute differential phase if not exists
            if not os.path.isfile(os.path.join(self.diroutorder, self.subdiroutifg, imagebasename + '.data', 'ifg_result_unwrap').replace('01_ifg', '01_slc') + '.img') or self.overwrite == '1':
                 self.differential_phase([os.path.join(self.diroutorder, self.subdiroutifg, imagebasename + '.data'), os.path.join(self.diroutorder, self.subdiroutifg, imagebasename + '_nocorrect.data')]) 
        #Clean/rename temp files
        #Store production data in table
        if os.path.isfile(os.path.join(self.diroutorder, self.subdiroutslc, imagebasename + '.dim')):
            outputfiles.append(os.path.join(self.diroutorder, self.subdiroutslc, imagebasename + '.dim'))
            outputifg1.append(os.path.join(self.diroutorder, self.subdiroutifg, imagebasename + '.dim'))
            outputifg2.append(os.path.join(self.diroutorder, self.subdiroutifg, imagebasename + '_nocorrect.dim'))
        else:
            outputfiles.append('Not generated')
        outputtimes.append((datetime.now()-time1).seconds/60)
        datelist.append(dates[1])        
        return [outputfiles, outputifg1, outputifg2, datelist, outputtimes,imagebasename]
        
    def stack_polariz (self, imagepathlist, outputfile, pathgpt):
        if len(imagepathlist)==2:
            gptxml_file = os.path.join(self.DirProject, 'RES', 'stack_dualpolariz.xml')
            input1 = imagepathlist[0]
            input2 = imagepathlist[1]
            try:
                p = subprocess.check_output([pathgpt, gptxml_file, '-Pinput1='+input1, '-Pinput2='+input2, '-Ptarget1='+outputfile])
            except Exception as e:
                raise e
#             finally:
#                 p.terminate() # send sigterm, or ...
#                 p.kill()      # send sigkill    
        return p
    
    def TOPS_Merge_subswaths (self, imagepathlist, outputfile, pathgpt):
        if(len(imagepathlist)==2):
            gptxml_file = os.path.join(self.DirProject, 'RES', 'merge_2subswaths.xml')
            input1 = imagepathlist[0]
            input2 = imagepathlist[1]
            try:
                p = subprocess.check_output([pathgpt, gptxml_file, '-Pinput1='+input1, '-Pinput2='+input2, '-Ptarget1='+outputfile])
            except Exception as e:
                raise e
#             finally:
#                 p.terminate() # send sigterm, or ...
#                 p.kill()      # send sigkill    
        if(len(imagepathlist)==3):
            gptxml_file = os.path.join(self.DirProject, 'RES', 'merge_3subswaths.xml')
            input1 = imagepathlist[0]
            input2 = imagepathlist[1]
            input3 = imagepathlist[2]
            try:
                p = subprocess.check_output([pathgpt, gptxml_file, '-Pinput1='+input1, '-Pinput2='+input2, '-Pinput3='+input3, '-Ptarget1='+outputfile])
            except Exception as e:
                raise e
#             finally:
#                 p.terminate() # send sigterm, or ...
#                 p.kill()      # send sigkill    
        return p
    
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
        print("[*] Unzipping the data...")
        
        df_unzipassemb = pd.DataFrame(Parallel(n_jobs=self.num_cores)(delayed(self.mpUnzipAssembly)(imagelist, ImageDate, uniquedates[i]) for i in range(len(uniquedates))))
#         out_0_list = []
#         out_1_list = []
#         out_2_list = []
#         for i in range(len(uniquedates)):
#             print("Inside the for loop")
#             print(uniquedates[i])
#             out = self.mpUnzipAssembly(imagelist, ImageDate, uniquedates[i])
#             out_0_list.append(out[0])
#             out_1_list.append(out[1])
#             out_2_list.append(out[2])
#         print(out_0_list)
#         print(out_1_list)
#         print(out_2_list)
#         df_unzipassemb = [out_0_list,out_1_list,out_2_list]
        self.processdf['ImageOrigin'] = df_unzipassemb[0]
        self.processdf['ImagePreprocess'] = df_unzipassemb[1]
        self.processdf['TimeunzipAsseml'] = df_unzipassemb[2]
        self.processdf['ImageDate'] = uniquedates
        #Calibration of original products
        outputfiles = []
        outputtimes = []
        if self.calib != '0':
            if self.calib == 'sigma' or self.calib == 'Sigma':
                sigmaparam, gammaparam, betaparam = True, False, False
            elif self.calib == 'gamma' or self.calib == 'Gamma':
                sigmaparam, gammaparam, betaparam = False, True, False
            elif self.calib == 'beta' or self.calib == 'Beta':
                sigmaparam, gammaparam, betaparam = False, False, True  
            calib_bool = [sigmaparam, gammaparam, betaparam]
            outputfiles = [os.path.join(self.diroutorder, '00'+self.sufix, self.processdf['ImageDate'][i]+self.sufix+'.dim') for i in range(len(self.processdf))]
            outputtimes = list(Parallel(n_jobs=self.num_cores)(delayed(self.mpCalibration)(self.processdf['ImagePreprocess'][i], outputfiles[i], calib_bool) for i in range(len(self.processdf))))
            self.processdf['Outputfiles_calib'] = outputfiles
            self.processdf['Processing_minutes_calib'] = outputtimes
        else:
            self.sufix =''
            outputfiles = self.processdf['ImagePreprocess']
            outputtimes = ['0'] * len(self.processdf)
        self.processdf['Outputfiles'] = outputfiles
        self.processdf['Processing_minutes_calib'] = outputtimes
        self.processdf['ImageDate'] = pd.to_datetime(self.processdf['ImageDate'])
        self.processdf = self.processdf.sort_values(by=['ImageDate'])
        #IF MASTER IMAGE DATE MATCHES ANY DATE OF THE FILE LIST, RUN ALIGNMENT PROCESSING (OTHERWISE MASTER DATE IS COMPUTED AS THE 'MID POINT' DATE OF THE LIST)
        #Identify master image row and set up slave images dataframe
        if not any (self.processdf['ImageDate'] == self.datemaster):
            datelist = list(self.processdf.ImageDate.sort_values().astype(int)/10**9)
            self.datemaster = pd.to_datetime(datelist[min(range(len(datelist)), key = lambda i: abs(datelist[i]-(datelist[-1]+datelist[0])/2))], unit='s', origin='unix').strftime('%Y%m%d')
        print("[*] Master date: ",self.datemaster)
        self.processdf['ImageDate_str'] = [date_obj.strftime('%Y%m%d') for date_obj in self.processdf['ImageDate']]
        #Master image is reindexed to the top of the df
        masterdf = self.processdf[self.processdf['ImageDate_str'] == self.datemaster]
        slavedf = self.processdf.copy()
        pairspath = mountpairs(pd.concat([masterdf['Outputfiles'], slavedf['Outputfiles']], sort=False).reset_index(), 'Outputfiles')
        pairsdate = mountpairs(pd.concat([masterdf, slavedf], sort=False).reset_index(), 'ImageDate_str')
        #Coregistration process        
        df_coreg = pd.DataFrame(Parallel(n_jobs=self.num_cores)(delayed(self.pair_coregistration)(pairspath[i], pairsdate[i]) for i in range(len(pairspath))))             
        outputfiles = [item[0] for item in df_coreg[0]]
        datelist = [item[0] for item in df_coreg[3]]
        outputtimes = [item[0] for item in df_coreg[4]]
        imagebasenames = [item for item in df_coreg[5]]
        #Get position of master and add master record
        self.processdf.reset_index(drop=True, inplace=True) 
        masterid = self.processdf[self.processdf['ImageDate_str']==self.datemaster].index[0]
        #Fill processing info into df
        self.processdf['Outputfiles_align'+self.sufix] = outputfiles
        self.processdf['Processing_minutes_align'] = outputtimes
        self.processdf.to_csv(os.path.join(self.diroutorder, self.logfilename), sep=';', index=False)
        #Generate reference matrices (lat/lon) for further geocoding. Baselines file is also generated
        self.generateGeocoding()
        #Remove temp dir
        for im in imagebasenames:
            self.cleaning_coreg(os.path.join(self.diroutorder, self.subdiroutslc, im + '.dim'))        
        if os.path.isdir(self.tempdir):
            shutil.rmtree(self.tempdir)
        msg = 'Coregistration generated'
        return msg
    
    def mpCalibration(self, inputfile, outputfile, calib_bool):
        gpt_calibration = os.path.join(self.DirProject, 'RES', 'calibration_outputcomplex.xml')
        time1 = datetime.now()
        if not os.path.isfile(outputfile) or self.overwrite == '1':
            try:
                p = subprocess.check_output([self.pathgpt, gpt_calibration, '-Pinput1='+inputfile, '-Psigmaparam='+str(calib_bool[0]), '-Pgammaparam='+str(calib_bool[1]), '-Pbetaparam='+str(calib_bool[2]), '-Ptarget1='+outputfile])
            except Exception as e:
                raise e
#             finally:
#                 p.terminate() # send sigterm, or ...
#                 p.kill()      # send sigkill
        return((datetime.now()-time1).seconds/60)

    def cleaning_coreg(self, imagecoreg):
        print("[*] Cleaning unnecessary files...")
        #Go through the coregistered images and change file names, remove master image...
        imageref = os.path.splitext(os.path.basename(imagecoreg))[0]
        imagedir = os.path.splitext(imagecoreg)[0] + '.data/'
        ifgdir = os.path.join(self.diroutorder, self.subdiroutifg, imageref+'.data')
        ifgnocorrectdir = os.path.join(self.diroutorder, self.subdiroutifg, imageref+'_nocorrect.data')
        deletelist = glob(os.path.join(imagedir, '*mst*'))
        if os.path.isfile(os.path.join(imagedir, 'ifg_result.img')):
            deletelist.append(os.path.join(imagedir, 'ifg_result.img'))
            deletelist.append(os.path.join(imagedir, 'ifg_result.hdr'))
        for file in deletelist:
            os.remove(file)
        changenamelist = glob(os.path.join(imagedir, '*slv*'))
        for file in changenamelist:
            os.rename(file, file.replace('_slv1', ''))
            filetemp = file.replace('_slv1', '')
            os.rename(filetemp, filetemp.replace(filetemp[filetemp.find('_slv3'):filetemp.find('_slv3')+15], ''))
        for fold in [ifgdir, ifgnocorrectdir]:
            for file in glob(os.path.join(fold, '*.*')):
                if 'coh' not in file:
                    filetemp = file
                    if 'slv' in file:
                        filetemp = file.replace(file[file.find('_slv'):file.find('_slv')+15], '')
                    filedest = filetemp.replace(filetemp[filetemp.find('_ifg_V')+8:filetemp.find('_ifg_')+18], '')
                    os.rename(file, filedest)
                
    def generate_baselinelist (self):
        self.baselinelist = self.processdf.copy()
        #Master record removed (no baseline between the same image as master and slave)
        self.baselinelist = self.baselinelist.drop(self.baselinelist[self.processdf['ImageDate'] == datetime.strptime(self.datemaster, '%Y%m%d')].index)
        self.baselinelist.reset_index(drop=True, inplace=True) 
        self.baselinelist['Master_date'] = self.datemaster
        self.baselinelist['Slave_date'] = self.baselinelist['ImageDate_str']
        self.baselinelist['Perp_Baseline'] = [round(-1*float(self.searchMetabytag(self.baselinelist['Outputfiles_align'+self.sufix][i], 'Perp Baseline')[1]), 2) for i in range(len(self.baselinelist))]
        self.baselinelist['Temp_Baseline'] = [int(-1*float(self.searchMetabytag(self.baselinelist['Outputfiles_align'+self.sufix][i], 'Temp Baseline')[1])) for i in range(len(self.baselinelist))]
        self.baselinelist['Coregistered_SLCs'] = self.baselinelist['Outputfiles_align'+self.sufix]
        fieldlist = []
        #Delete fields not needed
        for tag in ['Timeunzip', 'ImageDate', 'ImageOrigin', 'Outputfiles', 'ImagePreprocess', 'Processing_minutes']:
            for band in list(self.baselinelist):
                if tag in band:
                    fieldlist.append(band)
        self.baselinelist = self.baselinelist.drop(fieldlist, axis=1)
        waninglist = self.baselinelist.copy()
        while len(waninglist)>1:
            for i in range(len(waninglist)-1):
                perp = round(float(waninglist['Perp_Baseline'][i+1]) - float(waninglist['Perp_Baseline'][0]), 2)
                temp = int(round(float(waninglist['Temp_Baseline'][i+1]) - float(waninglist['Temp_Baseline'][0]), 0))
                self.baselinelist = self.baselinelist.append(pd.Series([str(waninglist['Slave_date'][0]), str(waninglist['Slave_date'][i+1]), perp, temp, None], index=self.baselinelist.columns), ignore_index=True)
            waninglist = waninglist.drop(0)
            waninglist.reset_index(drop=True, inplace=True)
        #Filter table with maxbasetemp and maxbaseperp
        self.baselinelist['Temp_Baseline'] = pd.to_numeric(self.baselinelist['Temp_Baseline'])
        self.baselinelist['Perp_Baseline'] = pd.to_numeric(self.baselinelist['Perp_Baseline'])
        if self.maxbasetemp != '':
            if self.maxbaseperp != '':
                self.baselinefiltered = self.baselinelist[(abs(self.baselinelist.Temp_Baseline) <= self.maxbasetemp) & (abs(self.baselinelist.Perp_Baseline) <= self.maxbaseperp)]
            else:
                self.baselinefiltered = self.baselinelist[(abs(self.baselinelist.Temp_Baseline) <= self.maxbasetemp)]
        else:
            if self.maxbaseperp != '':
                self.baselinefiltered = self.baselinelist[(abs(self.baselinelist.Perp_Baseline) <= self.maxbaseperp)]
            else:
                self.baselinefiltered = self.baselinelist.copy()
        self.baselinefiltered.to_csv(os.path.join(self.diroutorder, 'baselines_filtered.csv'), sep=';', index=False)
        self.baselinelist.to_csv(os.path.join(self.diroutorder, 'baselines.csv'), sep=';', index=False)
            
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
                print('Sentinel 2 grid file provided not coherent with coordinates ')
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
     
    def generateGeocoding(self):
        if self.overwrite == '1':
            from snappy import ProductIO, GPF, HashMap
            self.generate_baselinelist()
            gptxml_file = os.path.join(self.DirProject, 'RES', 'ifg_ortholatlon.xml')
            indexGC = int(self.baselinefiltered[['Perp_Baseline']].idxmax())
            inputfile1 = self.processdf[self.processdf['ImageDate_str']==self.baselinefiltered.loc[indexGC]['Master_date']]['Outputfiles_align'+self.sufix].tolist()[0]
            inputfile2 = self.processdf[self.processdf['ImageDate_str']==self.baselinefiltered.loc[indexGC]['Slave_date']]['Outputfiles_align'+self.sufix].tolist()[0]
            date1 = self.baselinefiltered.loc[indexGC]['Master_date']
            date2 = self.baselinefiltered.loc[indexGC]['Slave_date']
            outputdir = os.path.join(self.diroutorder, self.subdiroutgc)
            outputfile = os.path.join(outputdir, 'ifg_gc_'+date1+'_'+date2+'.dim')
            print("[*] Geocoding step")
            print([self.pathgpt, gptxml_file, '-Pinput1='+inputfile1, '-Pinput2='+inputfile2, '-Ptarget1='+outputfile])
            try:
                p = subprocess.check_output([self.pathgpt, gptxml_file, '-Pinput1='+inputfile1, '-Pinput2='+inputfile2, '-Ptarget1='+outputfile])
            except Exception as e:
                raise e
    #         finally:
    #             p.terminate() # send sigterm, or ...
    #             p.kill()      # send sigkill        
            shutil.copy(os.path.join(os.path.splitext(outputfile)[0]+'.data', 'orthorectifiedLon.img'), os.path.join(outputdir, 'orthorectifiedLon.img'))
            shutil.copy(os.path.join(os.path.splitext(outputfile)[0]+'.data', 'orthorectifiedLon.hdr'), os.path.join(outputdir, 'orthorectifiedLon.hdr'))
            shutil.copy(os.path.join(os.path.splitext(outputfile)[0]+'.data', 'orthorectifiedLat.img'), os.path.join(outputdir, 'orthorectifiedLat.img'))
            shutil.copy(os.path.join(os.path.splitext(outputfile)[0]+'.data', 'orthorectifiedLat.hdr'), os.path.join(outputdir, 'orthorectifiedLat.hdr'))
            shutil.copy(os.path.join(os.path.splitext(outputfile)[0]+'.data', 'tie_point_grids', 'incident_angle.img'), os.path.join(outputdir, 'incident_angle.img'))
            shutil.copy(os.path.join(os.path.splitext(outputfile)[0]+'.data', 'tie_point_grids', 'incident_angle.hdr'), os.path.join(outputdir, 'incident_angle.hdr'))
            shutil.copy(os.path.join(os.path.splitext(outputfile)[0]+'.data', 'elevation.img'), os.path.join(outputdir, 'elevation.img'))
            shutil.copy(os.path.join(os.path.splitext(outputfile)[0]+'.data', 'elevation.hdr'), os.path.join(outputdir, 'elevation.hdr')) 
            #If ortho latlon numpy matrices needed
            prod = ProductIO.readProduct(outputfile)
    #         ortho_lat = prod.getBand('orthorectifiedLat')
    #         ortho_lon = prod.getBand('orthorectifiedLon')
    #         w = prod.getSceneRasterWidth()
    #         h = prod.getSceneRasterHeight()
    #         array = np.zeros((w, h), dtype=np.float32)
    #         latpix = ortho_lat.readPixels(0, 0, w, h, array)
    #         lonpix = ortho_lon.readPixels(0, 0, w, h, array)
    #         np.save(os.path.join(outputdir, 'lat'), latpix)
    #         np.save(os.path.join(outputdir, 'lon'), lonpix)
            #Management of CRS, by default UTM, zone according to longitude
    #         try:
    #             print(self.epsgout)
    #         except AttributeError:
    #             zone = int((np.mean(lonpix) + 180) / 6) + 1
    #             if np.mean(latpix)>0:
    #                 south = False
    #             else:
    #                 south = True
    #             crs = CRS.from_dict({'proj': 'utm', 'zone': zone, 'south': south})
    #             self.epsgout = crs.to_authority()[1]
    #         if not hasattr(self, 'spatialres'):
    #             self.spatialres = 10
            #Resampling
    #         inProj = 'epsg:'+self.epsgin
    #         outProj = 'epsg:'+self.epsgout
    #         transformer = Transformer.from_crs(inProj, outProj)
    #         self.x, self.y = transformer.transform(lonpix,latpix)
    #         np.save(os.path.join(outputdir, 'x'), self.x)
    #         np.save(os.path.join(outputdir, 'y'), self.y)
            #Incident angle
            inc = prod.getRasterDataNode('incident_angle')
            theta = np.zeros((h, w), dtype=np.float32)
            theta = inc.readPixels(0, 0, w, h, theta)
    #         np.save(os.path.join(outputdir, 'incid_angle'), theta)
            self.array2raster(os.path.join(outputdir, 'incid_angle.img'), (0, 0), 1, -1, theta)
            #ProductIO.writeProduct(theta, os.path.join(outputdir, 'incident_angle.img'), None)
            #Adjust to the Sentinel 2 grid
    #         ulxgrid, ulygrid, lrxgrid, lrygrid = self.check_S2grid(np.min(self.x), np.max(self.y), np.max(self.x), np.min(self.y))
    #         wout = int((lrxgrid - ulxgrid) / self.spatialres)
    #         hout = int((ulygrid - lrygrid) / self.spatialres)
    #         x1 = np.linspace(lrxgrid, ulxgrid, wout+1)
    #         y1 = np.linspace(lrygrid, ulygrid, hout+1)
    #         grid_x, grid_y = np.meshgrid(x1,y1)
            # inci_grid = griddata((self.x.flatten(), self.y.flatten()), theta.flatten(), (grid_x, grid_y), method='nearest')
            # np.save(os.path.join(outputdir, 'incid_angle_EPSG'+self.epsgout, inci_grid))   
        msg = 'Geocoding generated'
        return msg

    def compute_phase(self, inputdir, polariz, outputfile):
            realband = gdal.Open(glob(os.path.join(inputdir, 'i_ifg_' + polariz + '*.img'))[0])
            imagband = gdal.Open(glob(os.path.join(inputdir, 'q_ifg_' + polariz + '*.img'))[0])
            realarr = realband.ReadAsArray()
            imagarr = imagband.ReadAsArray()
            complex_mat = realarr + imagarr * 1j
            phase = np.angle(complex_mat)
            self.array2raster(outputfile + '_'+ polariz + '.img', (0, 0), 1, -1, phase)
            return phase
    
    def differential_phase(self, ifgpair):
        polariz = self.polarizphaseunwrap
        output = os.path.join(ifgpair[0], 'ifg_result').replace('01_ifg', '01_slc') + '.img'
        time1 = datetime.now()
        if not os.path.isfile(output) or self.overwrite == '1':
            phasecomplete = self.compute_phase(ifgpair[0], polariz, os.path.join(os.path.splitext(ifgpair[0])[0] + '_ifg_compl.data'))
            phasediff = self.compute_phase(ifgpair[1], polariz, os.path.join(os.path.splitext(ifgpair[0])[0] + '_ifg_diff.data'))
            phasesubtr = phasediff - phasecomplete
            re_result = np.cos(phasesubtr)
            im_result = np.sin(phasesubtr)
            complex_mat = re_result + im_result * 1j
            phase_result = np.angle(complex_mat)
            self.array2raster(output, (0, 0), 1, -1, phase_result)
        else:
            phase_result =  gdal.Open(output).ReadAsArray()
        outputphasefile = os.path.join(ifgpair[0], 'ifg_result_unwrap').replace('01_ifg', '01_slc') + '.img'
        if not os.path.isfile(outputphasefile) or self.overwrite == '1':
            self.unwrapphase(phase_result, outputphasefile)
        return(outputphasefile ,(datetime.now()-time1).seconds/60)
    
    def blocksplitter(self, width, height, block_w, block_h, overlapw, overlaph):
        nblockx = int(width/block_w) + (width%block_w>0)
        nblocky = int(height/block_h) + (height%block_h>0)
        splitblock = []
        for bly in range(nblocky):
            for blx in range(nblockx):
                ulx = (block_w - overlapw) * blx
                uly = (block_h - overlaph) * bly
                if bly == 0:
                    uly = 0
                if blx == 0:
                    ulx = 0
                lrx = ulx + block_w
                lry = uly + block_h
                if blx == nblockx - 1:
                    lrx = width
                if bly == nblocky - 1:
                    lry = height   
                splitblock.append((bly, blx, ulx, uly, lrx, lry))
        return(np.array(splitblock))
    
    def unwrapphase(self, phase_result, outputphasefile):
        #Call to SNAPHU, which has a limited size for phase unwrapping
        #Slicing parameters block size, overlap between blocks
        block_w = 3000
        block_h = 3000
        overlapw = 10
        overlaph = 10
        phase_w = phase_result.shape[1]
        phase_h = phase_result.shape[0]
        splitblock = self.blocksplitter(phase_w, phase_h, block_w, block_h, overlapw, overlaph)
        nblockx = max(splitblock[:,1])+1
        nblocky = max(splitblock[:,2])+1
        #Dummy coherence (one's). Filtered (where phase ==0)=0
        coh = np.where(phase_result==0, 0, 1)
        outcoh = os.path.join(os.path.dirname(outputphasefile), 'coh_dummy.img')
        #self.array2raster(outcoh, (0, 0), 1, -1, coh)
        unwrapfilelist = []
        if not os.path.exists(os.path.join(os.path.dirname(outputphasefile), 'unwrap')):
            os.makedirs(os.path.join(os.path.dirname(outputphasefile), 'unwrap'))
        for row in splitblock:
            outblockfile = os.path.join(os.path.dirname(outputphasefile), 'unwrap', 'block_'+str(row[0])+'_'+str(row[1])+'.img')
            self.array2raster(outblockfile, (row[2], -row[3]), 1, -1, phase_result[row[3]:row[5], row[2]:row[4]])
            cohblockfile = os.path.join(os.path.dirname(outputphasefile), 'unwrap', 'cohblock_'+str(row[0])+'_'+str(row[1])+'.img')
            self.array2raster(cohblockfile, (row[2], -row[3]), 1, -1, coh[row[3]:row[5], row[2]:row[4]])
            outblockunwrap = os.path.join(os.path.dirname(outputphasefile), 'unwrap', 'unwrapblock_'+str(row[0])+'_'+str(row[1])+'.img')
            try:
                p = subprocess.check_output([self.snaphupath, '-s', outblockfile, str(row[4]-row[2]), '-o', outblockunwrap, '-c', cohblockfile, '-f', self.snaphuconf])
            except Exception as e:
                raise e
#             finally:
#                 p.terminate() # send sigterm, or ...
#                 p.kill()      # send sigkill               
            shutil.copyfile(outblockfile.replace('img', 'hdr'), outblockunwrap.replace('img', 'hdr'))
            unwrapfilelist.append(outblockunwrap)
        #Setup of complete matrix
        unwrap_complete = np.zeros((phase_h, phase_w))
        #Merge blocks
        for bl in range(len(splitblock)-1):   
            rowref = splitblock[bl, 0]
            colref = splitblock[bl, 1]
            rowtar = splitblock[bl+1, 0]
            coltar = splitblock[bl+1, 1]
            #Relative positions of overlap to every image
            ref_ov_coords = [block_w-overlapw, 0, block_w, (splitblock[(splitblock[:, 0]==rowref) & (splitblock[:, 1]==colref), 5].item() - splitblock[(splitblock[:, 0]==rowref) & (splitblock[:, 1]==colref), 3].item())]
            tar_ov_coords = [0, 0, overlapw, (splitblock[(splitblock[:, 0]==rowref) & (splitblock[:, 1]==colref), 5].item() - splitblock[(splitblock[:, 0]==rowref) & (splitblock[:, 1]==colref), 3].item())]
            im_ref_index = bl
            im_tar_index = bl+1
            if rowref != rowtar:
                rowref = rowtar-1
                colref = coltar
                ref_ov_coords = [0, block_h-overlaph, block_w, block_h]
                tar_ov_coords = [0, 0, block_w, overlaph]
                im_ref_index = bl-(nblockx-1)
                im_tar_index = bl+1                 
            if bl==0:
                refblock = gdal.Open(unwrapfilelist[im_ref_index]).ReadAsArray()
            else:
                if rowref != rowtar:
                    refblock = unwrap_complete[splitblock[im_ref_index, 3]:splitblock[im_ref_index, 5], splitblock[im_ref_index, 2]:splitblock[im_ref_index, 4]]
                else:
                    refblock = tarblock
            tarblock = gdal.Open(unwrapfilelist[im_tar_index]).ReadAsArray()
            ref_overlap = refblock[ref_ov_coords[1]:ref_ov_coords[3], ref_ov_coords[0]:ref_ov_coords[2]]
            tar_overlap = tarblock[tar_ov_coords[1]:tar_ov_coords[3], tar_ov_coords[0]:tar_ov_coords[2]]
            #Overlap values may differ, most frequent value used for offsetting the unwrapped blocks
            diff_overlap = ref_overlap - tar_overlap
            offvalue = np.percentile(sorted(diff_overlap.flatten()), [10, 90])[-1]
            print('row:'+str(rowref), 'col:'+str(colref), 'offs:'+str(offvalue))
            tarblock = tarblock + offvalue
            #Fill unwrap matrix
            if bl == 0:
                unwrap_complete[splitblock[bl, 3]:splitblock[bl, 5], splitblock[bl, 2]:splitblock[bl, 4]] = refblock
            unwrap_complete[splitblock[bl+1, 3]:splitblock[bl+1, 5], splitblock[bl+1, 2]:splitblock[bl+1, 4]] = tarblock
        self.array2raster(outputphasefile, (0, 0), 1, -1, unwrap_complete)
        #Delete temp unwrap dir
        shutil.rmtree(os.path.join(os.path.dirname(outputphasefile), 'unwrap'))        
        
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

