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
import logging
import stsa
import shapely.wkt
from shapely.geometry import Polygon
import geopandas as gpd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("sar2cube_debug.log"),
    ]
)

def unzipS1(zipfile, dirouput):
    with ZipFile(zipfile, 'r') as zipObj:
        zipObj.extractall(dirouput)

def mountpairs(dflist, field): #Pairs combining the first element with the rest of the list elements
    pairs = []
    pairs = [((dflist[field][0]), (dflist[field][i+1])) for i in range(len(dflist)-1)]
    return pairs

class model():
    def __init__(self,
                 process_type,
                 imagelistfile,
                 diroutorder,
                 overwrite,
                 datemaster,
                 AOI,
                 calib,
                 pathgpt,
                 snappypath,
                 DirProject,
                 num_cores,
                 do_interferogram,
                 do_phase):

        self.process_type  = process_type
        self.imagelistfile = imagelistfile
        self.datemaster    = datemaster
        if self.datemaster == '':
            self.datemaster = 'nodate'
        self.AOI = AOI
        if not self.AOI:
            self.AOI = ['', '', '', '']
        self.diroutorder = diroutorder
        self.overwrite   = overwrite
        self.calib       = str(calib)
        if self.calib != '0':
            self.sufix = '_calib'
        else:
            self.sufix = ''
        self.subdiroutslc = '01_slc' + self.sufix
        self.subdiroutifg = '01_ifg' + self.sufix
        self.subdiroutml  = '02_ml'
        self.subdiroutgc  = '03_gc'
        self.pathgpt      = pathgpt
        self.snappypath   = snappypath
        self.DirProject   = DirProject
        self.processdf    = pd.DataFrame()
        self.baselinelist = pd.DataFrame()
        self.maxbaseperp  = ''
        self.maxbasetemp  = ''
        self.baselinefiltered = pd.DataFrame()
        self.logfilename = 'process_log.csv'
        self.epsgin = '4326'
        self.polarizphaseunwrap = 'VH'
        # Number of cores should fit the ratio between RAM memory and 23 Gb/core (need for coregistration).
        # We can try to use more for the other steps
        self.num_cores = num_cores
        self.tempdir = None
        self.do_interferogram = int(do_interferogram)
        self.do_phase = int(do_phase)
        self.iw = None
        
        self.gptxml_unzip         = os.path.join(self.DirProject, 'RES', 'sliceassembly.xml')
        self.gpt_calibration      = os.path.join(self.DirProject, 'RES', 'calibration_outputcomplex.xml')
        self.gptxml_ESD           = os.path.join(self.DirProject, 'RES', 'TOPSAR_Coreg_ESD.xml')
        self.gptxml_Deburst       = os.path.join(self.DirProject, 'RES', 'Deburst_Write.xml')
        self.gptxml_IFG           = os.path.join(self.DirProject, 'RES', 'interferogram_deburst.xml')
        self.gptxml_stackdual     = os.path.join(self.DirProject, 'RES', 'stack_dualpolariz.xml')
        self.gptxml_merge2IW      = os.path.join(self.DirProject, 'RES', 'merge_2subswaths.xml')
        self.gptxml_merge3IW      = os.path.join(self.DirProject, 'RES', 'merge_3subswaths.xml')
        self.gptxml_geocod        = os.path.join(self.DirProject, 'RES', 'ifg_ortholatlon.xml')

        if not os.path.exists(self.diroutorder):
            os.makedirs(self.diroutorder)            
        sys.path.append(self.snappypath)
        # import esa_snappy
        from esa_snappy import ProductIO


    def coregistration_ifg (self):
        #Read the image file
        imagelist = pd.read_csv(self.imagelistfile, sep=';', decimal=',')
        #Unzip if zipfiles, slice assembly if two products for the same date in imagelist and fill up ImagePreprocess, ImageName, ImageDate columns
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
            elif os.path.splitext(imagelist.Image[i])[1].lower() == '.safe':
                name = os.path.basename(imagelist.Image[i]).split('.')[0]
            else:
                logging.error("[!] Image format not supported: " + imagelist.Image[i])
                sys.exit()
            ImageDate.append(name[name.find('T')-8:name.find('T')])
        imagelist['Date'] = ImageDate
        uniquedates = list(set(ImageDate))
        
        #######################################
        # Unzipping SLC data & Slice Assembly # 
        #######################################
        
        logging.info("Unzipping the data...")
        df_unzipassemb = pd.DataFrame(Parallel(n_jobs=self.num_cores)(delayed(self.mpUnzipAssembly)(imagelist, ImageDate, uniquedates[i]) for i in range(len(uniquedates))))
        self.processdf['ImageOrigin'] = df_unzipassemb[0]
        self.processdf['ImagePreprocess'] = df_unzipassemb[1]
        self.processdf['TimeunzipAsseml'] = df_unzipassemb[2]
        self.processdf['ImageDate'] = uniquedates
        
        #######################
        # Calibration process # 
        #######################
        
        logging.info("Computing calibration...")
        #Calibration of original products
        outputfiles = []
        outputtimes = []
        if self.calib != '0':
            sigmaparam, gammaparam, betaparam = False, False, False
            if self.calib.lower() == 'sigma':
                sigmaparam = True
            elif self.calib.lower() == 'gamma':
                gammaparam = True
            elif self.calib.lower() == 'beta':
                betaparam  = True  
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
            
        logging.info("Master date: " + str(self.datemaster))
        
        self.processdf['ImageDate_str'] = [date_obj.strftime('%Y%m%d') for date_obj in self.processdf['ImageDate']]
        #Master image is reindexed to the top of the df
        masterdf = self.processdf[self.processdf['ImageDate_str'] == self.datemaster]
        slavedf = self.processdf.copy()
        pairspath = mountpairs(pd.concat([masterdf['Outputfiles'], slavedf['Outputfiles']], sort=False).reset_index(), 'Outputfiles')
        pairsdate = mountpairs(pd.concat([masterdf, slavedf], sort=False).reset_index(), 'ImageDate_str')
        
        ################################
        # AOI - Subswaths intersection #
        ################################
        
        logging.info("Computing AOI - Subswaths intersection...")

        intersected_subswaths = ['IW1', 'IW2', 'IW3']
        if any (self.AOI):
            lonmin = self.AOI[0]
            latmax = self.AOI[1]
            lonmax = self.AOI[2]
            latmin = self.AOI[3]
            input4 = 'POLYGON (('+lonmin+' '+latmin+', '+lonmax+' '+latmin+', '+lonmax+' '+latmax+', '+lonmin+' '+latmax+', '+lonmin+' '+latmin+'))'
            crs = 'epsg:4326'
            poly = shapely.wkt.loads(input4)
            polygon = gpd.GeoDataFrame(index=[0], crs=crs, geometry=[poly])   
            # polygon.to_file('wktAoi.shp')
            s1 = stsa.TopsSplitAnalyzer(target_subswaths=['iw1', 'iw2', 'iw3'], polarization='vh')
            s1.load_data(zip_path=imagelist.Image[0])
            s1._create_subswath_geometry()
            poly1 = s1.df.reset_index(drop=True)
            poly1_intersection = poly1.copy()
            for index, orig in poly1.iterrows():
                for index2, ref in polygon.iterrows():
                    if ref['geometry'].intersects(orig['geometry']):
                        continue
                    else:
                        poly1_intersection = poly1_intersection.drop(index=index)
            subswaths_list = list(poly1_intersection['subswath'])
            intersected_subswaths = list(set(subswaths_list))
        
        self.iw = intersected_subswaths
        logging.info(f"The subswaths intersecting the provided AOI are: {self.iw}")
        
        ##########################
        # Coregistration process # 
        ##########################
        
        logging.info("Starting the coregistration process...")
        
        self.tempdir = os.path.join(self.diroutorder, self.subdiroutslc, 'temp')
        if not os.path.isdir(self.tempdir):
            os.makedirs(self.tempdir)
        
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
        
        #####################
        # Geocoding process #
        #####################
        
        #Generate reference matrices (lat/lon) for further geocoding. Baselines file is also generated
        self.generateGeocoding()
        #Remove temp dir
        for im in imagebasenames:
            self.cleaning_coreg(os.path.join(self.diroutorder, self.subdiroutslc, im + '.dim'))
        if self.tempdir is not None:
            if os.path.isdir(self.tempdir):
                shutil.rmtree(self.tempdir)
        
        return
        
    def pair_coregistration(self, paths, dates):
    
        output_sufix = '_Coregistered'
        
        self.subdiroutslc = '01_slc' + self.sufix
        self.subdiroutifg = '01_ifg' + self.sufix  
        
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
        input4 = ''
        if any (self.AOI):
            input4 = 'POLYGON (('+lonmin+' '+latmin+', '+lonmax+' '+latmin+', '+lonmax+' '+latmax+', '+lonmin+' '+latmax+', '+lonmin+' '+latmin+'))'

        input3 = self.iw            
        input5 = ['VV', 'VH']
        input6 = 'false' #ortholatlon output flag
        
        #The process is executed for every subswath
        if os.path.splitext(input2)[1] == '.safe' or os.path.splitext(input2)[1] == '.dim':
            slavename = dates[1] + '_SLC' + self.sufix
        else:
            logging.error("[!] Image format not supported: " + input2)
            sys.exit()
        imagebasename = dates[0] + '_' + slavename + output_sufix
        if not os.path.isfile(os.path.join(self.diroutorder, self.subdiroutslc, imagebasename + '.dim')) \
        or not os.path.isfile(os.path.join(self.diroutorder, self.subdiroutifg, imagebasename + '.dim')) \
        or self.overwrite == '1':
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
                logging.info("Processing " + output1)
                #Generate coregistration for every subswath
                for k in range(len(input3)):
                    ################################################################################################
                    # First step: Read, TOPSAR-Split, Apply orbit file, back-geocoding-enhanced spectral diversity #
                    ################################################################################################
                    
                    time_start = datetime.now()

                    outputfile_esd = os.path.join(self.tempdir, output1 + '_' + input3[k] + '_ESD.dim')
                    logging.info("Read, TOPSAR-Split, Apply orbit file, back-geocoding-enhanced spectral diversity.")
                    logging.info(f"Output file: {outputfile_esd}")
                    if not os.path.isfile(outputfile_esd) or self.overwrite == '1':
                        try:
                            p = subprocess.check_output([self.pathgpt, self.gptxml_ESD, '-Pinput1='+input1, '-Pinput2='+input2, '-Pinput3='+input3[k], '-Pinput4='+input4, '-Pinput5='+polariz,'-Pinput6='+input6,'-Ptarget1='+outputfile_esd])
                        except Exception as e:
                            print(e)    
                    
                    elapsed_time = (datetime.now()-time_start).seconds/60
                    logging.info("Elapsed time: " + str(elapsed_time))

                    # Define output file paths for the three outputs:
                    outputfile1 = os.path.join(self.diroutorder, self.subdiroutifg, output1 + '_' + input3[k] + '.dim')
                    outputfile3 = os.path.join(self.diroutorder, self.subdiroutslc, output1 + '_' + input3[k] + 'im.dim')
                    
                    #####################################
                    # Second step: Read, Deburst, Write #
                    #####################################
                    time_start = datetime.now()
                    logging.info("Read, Deburst, Write.")
                    if not os.path.isfile(outputfile3) or self.overwrite == '1':
                        try:
                            p = subprocess.check_output([self.pathgpt, self.gptxml_Deburst, '-Pinput1='+outputfile_esd,'-Ptarget1='+outputfile3])
                        except Exception as e:
                            print(e)

                    elapsed_time = (datetime.now()-time_start).seconds/60
                    logging.info("Elapsed time: " + str(elapsed_time))
                    
                    ##############################################
                    # Second step: Interferogram, Deburst, Write #
                    ##############################################
                    if self.do_interferogram:
                        time_start = datetime.now()
                        logging.info("Interferogram, Deburst, Write.")

                        if not os.path.isfile(outputfile1) or self.overwrite == '1':
                            try:
                                p = subprocess.check_output([self.pathgpt, self.gptxml_IFG, '-Pinput1='+outputfile_esd,'-Ptarget1='+outputfile1])
                            except Exception as e:
                                print(e)     

                        elapsed_time = (datetime.now()-time_start).seconds/60
                        logging.info("Elapsed time: " + str(elapsed_time))
                    
                #After processing subswaths, they are put together with TOPS Merge Workflow
                #List of subswaths by checking the output files. The output is added at the end of the list
                subswath_list2 = sorted(glob(os.path.join(self.diroutorder, self.subdiroutifg, output1 + '*IW*.dim')))
                subswath_list3 = sorted(glob(os.path.join(self.diroutorder, self.subdiroutslc, output1 + '*.dim')))

                subswath_list2.append(os.path.join(self.diroutorder, self.subdiroutifg, output1))
                subswath_list3.append(os.path.join(self.diroutorder, self.subdiroutslc, output1+'im'))
                #Function TOPS Merge is only called when more than 1 subswath is generated. Only IW products supported (max. 3 subswaths)
                
                ##############
                # TOPS Merge #
                ############## 
                
                time_start = datetime.now()
                logging.info("TOPS Merge.")
                
                ## subswath_list2 will not exist if we do not generate the interferogram
                merged_list = [subswath_list2, subswath_list3]
                if not self.do_interferogram:
                    merged_list = [subswath_list3]
                for listmerge in merged_list:
                    if len(listmerge)<2:
                        logging.error("AOI does not match scene extents")
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
            
                elapsed_time = (datetime.now()-time_start).seconds/60
                logging.info("Elapsed time: " + str(elapsed_time))
                        
            #######################
            # Polarizations Stack #
            ####################### 
            
            time_start = datetime.now()
            logging.info("Polarizations Stack.")
            polariz_list2 = sorted(glob(os.path.join(self.diroutorder, self.subdiroutifg, imagebasename) + '*V*.dim'))
            polariz_list3 = sorted(glob(os.path.join(self.diroutorder, self.subdiroutslc, imagebasename) + '*V*.dim'))

            merged_pol_list = [polariz_list2, polariz_list3]

            if not self.do_interferogram:
                merged_pol_list = [polariz_list3]
            for listpol in merged_pol_list:
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
                            logging.error(f"Error while deleting file: {str(filePath)}")
                else:
                    #Change output name in case the product has only one polarization
                    shutil.move(listpol[0], outputfile)
                    shutil.move(os.path.splitext(listpol[0])[0]+'.data', os.path.splitext(outputfile)+'.data')
                    
            elapsed_time = (datetime.now()-time_start).seconds/60
            logging.debug("Elapsed time: " + str(elapsed_time))

        ###################
        # Synthetic Phase #
        ################### 
              
        if self.do_phase:
            time_start = datetime.now()
            logging.info("Computing synthetic phase.")

            self.synt_phase(os.path.join(self.diroutorder, self.subdiroutifg, imagebasename + '.data'))

            elapsed_time = (datetime.now()-time_start).seconds/60
            logging.debug("Elapsed time: " + str(elapsed_time))

        #Store production data in table
        if os.path.isfile(os.path.join(self.diroutorder, self.subdiroutslc, imagebasename + '.dim')):
            outputfiles.append(os.path.join(self.diroutorder, self.subdiroutslc, imagebasename + '.dim'))
            outputifg1.append(os.path.join(self.diroutorder, self.subdiroutifg, imagebasename + '.dim'))
        else:
            outputfiles.append('Not generated')
        outputtimes.append((datetime.now()-time1).seconds/60)
        datelist.append(dates[1])        
        return [outputfiles, outputifg1, outputifg2, datelist, outputtimes,imagebasename]
     
    def generateGeocoding(self):
        #############
        # Geocoding #
        ############# 
        if not os.path.exists(os.path.join(self.diroutorder, self.subdiroutgc)) or self.overwrite == '1':
            time_start = datetime.now()
            logging.info("Geocoding step.")
            try:
                os.makedirs(os.path.join(self.diroutorder, self.subdiroutgc))
            except:
                pass
            self.generate_baselinelist()
            indexGC = int(self.baselinefiltered[['Perp_Baseline']].idxmax())
            inputfile1 = self.processdf[self.processdf['ImageDate_str']==self.baselinefiltered.loc[indexGC]['Master_date']]['Outputfiles_align'+self.sufix].tolist()[0]
            inputfile2 = self.processdf[self.processdf['ImageDate_str']==self.baselinefiltered.loc[indexGC]['Slave_date']]['Outputfiles_align'+self.sufix].tolist()[0]
            date1 = self.baselinefiltered.loc[indexGC]['Master_date']
            date2 = self.baselinefiltered.loc[indexGC]['Slave_date']
            outputdir = os.path.join(self.diroutorder, self.subdiroutgc)
            outputfile = os.path.join(outputdir, 'ifg_gc_'+date1+'_'+date2+'.dim')
            try:
                p = subprocess.check_output([self.pathgpt, self.gptxml_geocod, '-Pinput1='+inputfile1, '-Pinput2='+inputfile2, '-Ptarget1='+outputfile])
            except Exception as e:
                raise e        
            shutil.copy(os.path.join(os.path.splitext(outputfile)[0]+'.data', 'orthorectifiedLon.img'), os.path.join(outputdir, 'orthorectifiedLon.img'))
            shutil.copy(os.path.join(os.path.splitext(outputfile)[0]+'.data', 'orthorectifiedLon.hdr'), os.path.join(outputdir, 'orthorectifiedLon.hdr'))
            shutil.copy(os.path.join(os.path.splitext(outputfile)[0]+'.data', 'orthorectifiedLat.img'), os.path.join(outputdir, 'orthorectifiedLat.img'))
            shutil.copy(os.path.join(os.path.splitext(outputfile)[0]+'.data', 'orthorectifiedLat.hdr'), os.path.join(outputdir, 'orthorectifiedLat.hdr'))
            shutil.copy(os.path.join(os.path.splitext(outputfile)[0]+'.data', 'tie_point_grids', 'incident_angle.img'), os.path.join(outputdir, 'incident_angle.img'))
            shutil.copy(os.path.join(os.path.splitext(outputfile)[0]+'.data', 'tie_point_grids', 'incident_angle.hdr'), os.path.join(outputdir, 'incident_angle.hdr'))
            shutil.copy(os.path.join(os.path.splitext(outputfile)[0]+'.data', 'elevation.img'), os.path.join(outputdir, 'elevation.img'))
            shutil.copy(os.path.join(os.path.splitext(outputfile)[0]+'.data', 'elevation.hdr'), os.path.join(outputdir, 'elevation.hdr')) 
            #If ortho latlon numpy matrices needed
            from esa_snappy import ProductIO
            prod = ProductIO.readProduct(outputfile)
            w = prod.getSceneRasterWidth()
            h = prod.getSceneRasterHeight()
            #Incident angle
            inc = prod.getRasterDataNode('incident_angle')
            theta = np.zeros((h, w), dtype=np.float32)
            theta = inc.readPixels(0, 0, w, h, theta)
            self.array2raster(os.path.join(outputdir, 'incid_angle.img'), (0, 0), 1, -1, theta)

            elapsed_time = (datetime.now()-time_start).seconds/60
            logging.debug("Elapsed time: " + str(elapsed_time))
        
        return  
    
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
            if os.path.splitext(imagefile)[1] == '.zip':
                outputpath = os.path.join(self.diroutputzip, os.path.splitext(os.path.basename(imagefile))[0]+'.SAFE', 'manifest.safe')
                if not os.path.isfile(outputpath) or self.overwrite == '1':
                    logging.info(f"Unzipping Sentinel data from file {imagefile}")
                    unzipS1(imagefile, self.diroutputzip)
                    inputimage = imagefile
                    timeproc = (datetime.now()-time1).seconds/60
                else:
                    inputimage = imagefile
                    timeproc = '0'
            elif os.path.splitext(imagefile)[1].lower() == '.safe':
                logging.error("Input are SAFE folders, currently not supported by the AOI Subswaths check")
                sys.exit()
                outputpath = imagefile + '/manifest.safe'
                if not os.path.isfile(outputpath) or self.overwrite == '1':
                    shutil.copytree(os.path.dirname(imagefile), os.path.join(self.diroutputzip, os.path.basename(os.path.dirname(imagefile))))
                    timeproc = (datetime.now()-time1).seconds/60
                else:
                    timeproc = '0'
                #outputpath = os.path.join(self.diroutputzip, os.path.basename(os.path.dirname(imagefile)), 'manifest.safe')
                inputimage = imagefile
            else:
                logging.error(f"Image format not supported: {imagelist.Image[i]}")
                sys.exit()
        elif ImageDate.count(date) == 2:
            logging.info("Running slice assembly...")
            #Slice-assembly for joining products
            assemblist = imagelist[imagelist.Image.str.contains(date)]['Image'].tolist()
            outputpath = os.path.join(self.diroutputzip, date, date+'_assembl.dim')
            if not os.path.isfile(outputpath) or self.overwrite == '1':
                try:
                    p = subprocess.check_output([self.pathgpt, self.gptxml_unzip, '-Pinput1='+assemblist[0], '-Pinput2='+assemblist[1], '-Ptarget1='+outputpath])
                except Exception as e:
                    raise e

            timeproc = (datetime.now()-time1).seconds/60
            inputimage = assemblist
        else:
            logging.error(f"More than 2 images (not allowed) with date: {str(date)}")
        return (inputimage, outputpath, timeproc)
    
    def stack_polariz(self, imagepathlist, outputfile, pathgpt):
        if len(imagepathlist)==2:
            input1 = imagepathlist[0]
            input2 = imagepathlist[1]
            try:
                p = subprocess.check_output([pathgpt, self.gptxml_stackdual, '-Pinput1='+input1, '-Pinput2='+input2, '-Ptarget1='+outputfile])
            except Exception as e:
                raise e  
        return p
    
    def TOPS_Merge_subswaths (self, imagepathlist, outputfile, pathgpt):
        if(len(imagepathlist)==2):
            input1 = imagepathlist[0]
            input2 = imagepathlist[1]
            try:
                p = subprocess.check_output([pathgpt, self.gptxml_merge2IW, '-Pinput1='+input1, '-Pinput2='+input2, '-Ptarget1='+outputfile])
            except Exception as e:
                raise e
        if(len(imagepathlist)==3):
            input1 = imagepathlist[0]
            input2 = imagepathlist[1]
            input3 = imagepathlist[2]
            try:
                p = subprocess.check_output([pathgpt, self.gptxml_merge3IW, '-Pinput1='+input1, '-Pinput2='+input2, '-Pinput3='+input3, '-Ptarget1='+outputfile])
            except Exception as e:
                raise e  
        return p
    
    def mpCalibration(self, inputfile, outputfile, calib_bool):
        logging.info("Running Calibration...")
        time1 = datetime.now()
        if not os.path.isfile(outputfile) or self.overwrite == '1':
            try:
                p = subprocess.check_output([self.pathgpt, self.gpt_calibration, '-Pinput1='+inputfile, '-Psigmaparam='+str(calib_bool[0]), '-Pgammaparam='+str(calib_bool[1]), '-Pbetaparam='+str(calib_bool[2]), '-Ptarget1='+outputfile])
            except Exception as e:
                raise e
                
        return((datetime.now()-time1).seconds/60) 

    def cleaning_coreg(self, imagecoreg):
        logging.info("Cleaning unnecessary files...")
        #Go through the coregistered images and change file names, remove master image...
        imageref = os.path.splitext(os.path.basename(imagecoreg))[0]
        imagedir = os.path.splitext(imagecoreg)[0] + '.data/'
        ifgdir = os.path.join(self.diroutorder, self.subdiroutifg, imageref+'.data')
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
        for file in glob(os.path.join(ifgdir, '*.*')):
            if 'coh' not in file and 'fep' not in file and 'tgp' not in file:
                os.remove(file)

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
                data_dict = dict(zip(list(self.baselinelist.columns),[str(waninglist['Slave_date'][0]), str(waninglist['Slave_date'][i+1]), perp, temp, None]))
                dataframe_to_append = pd.DataFrame(data=data_dict,index=range(len(data_dict)-1))
                self.baselinelist = pd.concat([self.baselinelist, dataframe_to_append], ignore_index=True)
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
            s2gridpath = os.path.join(self.DirProject, 'AUXDATA', 'tabularize_s2_footprint.csv')
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
    
    def synt_phase(self, ifg_folder):
        polariz = self.polarizphaseunwrap
                
        output = os.path.join(ifg_folder, 'phase').replace('01_ifg', '01_slc') + '.img'
        time1 = datetime.now()
        
        if not os.path.isfile(output) or self.overwrite == '1':
            if len(glob(os.path.join(ifg_folder, 'fep_' + polariz + '*.img'))) == 0:
                logging.error(f"Flat Earth Phase Component image (fep_*.img) not found in: {ifg_folder}. Did you recompile the Sentinel-1 toolbox before running the preprocessing? Please check this dicsussion:\nhttps://forum.step.esa.int/t/snap-compiled-from-source-to-set-output-phase-true/40096")
                return

            if len(glob(os.path.join(ifg_folder, 'tgp_' + polariz + '*.img'))) == 0:
                logging.error(f"Topographic Phase Component image (tgp_*.img) not found in: {ifg_folder}. Did you recompile the Sentinel-1 toolbox before running the preprocessing? Please check this dicsussion:\nhttps://forum.step.esa.int/t/snap-compiled-from-source-to-set-output-phase-true/40096")
                return
        
            flat_earth_phase = gdal.Open(glob(os.path.join(ifg_folder, 'fep_' + polariz + '*.img'))[0]).ReadAsArray()
            topo_phase = gdal.Open(glob(os.path.join(ifg_folder, 'tgp_' + polariz + '*.img'))[0]).ReadAsArray()
            synt_phase_result = flat_earth_phase + topo_phase
            self.array2raster(output, (0, 0), 1, -1, synt_phase_result)

        return

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

