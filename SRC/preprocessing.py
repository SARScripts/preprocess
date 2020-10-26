#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 23:21:09 2020

@author: paco.moral
"""

import sys
import os
import subprocess
import pandas as pd
from datetime import datetime 
import glob

def controlslash(path): #Control OS path generation
    return path.replace('\\', '/')

if __name__ == "__main__":
    try:
        arg1 = sys.argv[1]
    except IndexError:
        print ("Usage: python " + os.path.basename(__file__) + " <Path/to/project/config/file>")
        sys.exit(1)

    # Configuration variables from config file
    try:
        in_file = open(arg1, 'r')
        for line in in_file.readlines():
            if "PROJECT_DIRECTORY" in line:
                DirProj = controlslash(line.split('=')[1].strip())
                print ("PROJECT_DIRECTORY:" + DirProj)
            if "ORDER_FILE" in line:
                OrderFile = controlslash(line.split('=')[1].strip())
                print ("ORDER_FILE:" + OrderFile)
            if "GPTBIN_PATH" in line:
                pathgpt = controlslash(line.split('=')[1].strip())
                print ("GPTBIN_PATH:" + pathgpt)
            if "SNAPPY_FOLDER" in line:
                snappypath = controlslash(line.split('=')[1].strip())
                print ("SNAPPY_FOLDER:" + snappypath)

                #ASSESS WHETHER CACHE AND CPU ARE VARIABLES TO TAKE INTO ACCOUNT
            # if "CACHE" in line:
            #     CACHE = line.split('=')[1].strip()
            # if "CPU" in line:
            #     CPU = line.split('=')[1].strip()
    finally:
            in_file.close()
            
    sys.path.append(os.path.join(DirProj, 'SRC'))
    import preprocessing_utils as preproc
    
    orders = pd.read_csv(OrderFile, sep=';', decimal=',')
    
    #PARAMETERS NEEDED FOR EVERY PROCESSING
    #'alignment_interferogram': list_of_images (essential), output directory (essential), overwrite=[0, 1], masterimage_date (if not specified, master date taken as the closest date to mid-point between first and last date), AOI (if not specified whole scene considered)
    #'generate_pairs_list': list_of_images (essential)
    
    #A PROCESSING OBJECT IS CREATED FOR THE WHOLE PROCESSING CHAIN
    for i in range(len(orders)):
        #First process will always be alignment of images (coregister)
        if orders.Process_type[0] != 'alignment_interferogram':
            sys.exit('Preprocessing must start with alignment of a list of images')
            
        if orders.Process_type[i] == 'alignment_interferogram':
            Process = preproc.model(orders.Process_type[i], 
                            orders.Parameter_list.values[i].split(',')[0], 
                            orders.Parameter_list.values[i].split(',')[1].replace(' ', ''), 
                            orders.Parameter_list.values[i].split(',')[2].replace(' ', '').split('overwrite=')[1], 
                            orders.Parameter_list.values[i].split(',')[3].replace(' ', ''), 
                            orders.Parameter_list.values[i].split(',')[4].split(), 
                            pathgpt,
                            snappypath,
                            DirProj)
            #Parameters check
            if (not os.path.isfile(orders.Parameter_list.values[0].split(',')[0])):
                sys.exit('Introduce a file with a list of image paths')
            Process.alignment_ifg ()
        if orders.Process_type[i] == 'generate_list':
            if Process.processdf is not None:
                for im in Process.processdf['Outputfiles_align']:
                    if not os.path.isfile(im):
                        sys.exit('Outputs from coregistration processing not found (' + im + ')')
                Process.maxbasetemp = float(orders.Parameter_list.values[i].split(',')[0].replace(' ', ''))
                Process.maxbaseperp = float(orders.Parameter_list.values[i].split(',')[1].replace(' ', ''))
                Process.generate_baselinelist()
            else:
                sys.exit('Image list does not match with alignment preprocess list')
                
        if orders.Process_type[i] == 'generate_interferogram':
            if Process.baselinefiltered is not None:
                Process.generate_interferogram()
            else:
                sys.exit('Generate list not found')
                
        if orders.Process_type[i] == 'Multilook':
            if Process.processdf is not None:
                Process.azLooks = str(orders.Parameter_list.values[i].split(',')[0].replace(' ', ''))
                Process.rgLooks = str(orders.Parameter_list.values[i].split(',')[1].replace(' ', ''))
                Process.Multilook()
        
        if orders.Process_type[i] == 'Geocoding':
            if Process.processdf is not None:
                Process.epsg = str(orders.Parameter_list.values[i].split(',')[0].replace(' ', ''))
                Process.taglist = list(orders.Parameter_list.values[i].split(',')[1].replace("'", "").split())
                Process.Geocoding()
    
