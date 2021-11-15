#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import subprocess
import pandas as pd
from datetime import datetime 
import glob

def controlslash(path): # Control OS path generation
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
            if "SNAPHU_FOLDER" in line:
                snaphupath = controlslash(line.split('=')[1].strip())
                print ("SNAPHU_FOLDER:" + snaphupath)
            if "SNAPHU_CONFFILE" in line:
                snaphuconf = controlslash(line.split('=')[1].strip())
                print ("SNAPHU_CONFFILE:" + snaphuconf)    
            if "NUMBER_CORES" in line:
                num_cores = controlslash(line.split('=')[1].strip())
                print ("NUMBER_CORES:" + num_cores)    
            if "COMPUTE_INTERFEROGRAM" in line:
                compute_interferogram = controlslash(line.split('=')[1].strip())
                print ("COMPUTE_INTERFEROGRAM:" + compute_interferogram)    
            if "COMPUTE_PHASE" in line:
                compute_phase = controlslash(line.split('=')[1].strip())
                print ("COMPUTE_PHASE:" + compute_phase)
    finally:
            in_file.close()
            
    sys.path.append(os.path.join(DirProj, 'SRC'))
    import preprocessing_utils as preproc
    orders = pd.read_csv(OrderFile, sep=';', decimal=',')
    
    # PARAMETERS NEEDED FOR EACH TYPE OF PROCESSING
    # * 'alignment_interferogram': list_of_images (essential), 
    # output directory (essential), overwrite=[0, 1], 
    # masterimage_date (if not specified, master date taken as 
    # the closest date to mid-point between first and last date), 
    # AOI ULx, ULy, LRx, LRy (if not specified whole scene considered),
    # calib=[0, ]
    # * 'generate_pairs_list': list_of_images (essential)
    
    time1 = datetime.now()
    # A PROCESSING OBJECT IS CREATED FOR THE WHOLE PROCESSING CHAIN
    for i in range(len(orders)):
        # First process will always be coregistration of images
        if orders.Process_type[0] != 'alignment_interferogram':
            sys.exit('Preprocessing must start with coregistration of a list of images')
        if orders.Process_type[i] == 'alignment_interferogram':
            Process = preproc.model(orders.Process_type[i], 
                            orders.Parameter_list.values[i].split(',')[0], 
                            orders.Parameter_list.values[i].split(',')[1].replace(' ', ''), 
                            orders.Parameter_list.values[i].split(',')[2].replace(' ', '').split('overwrite=')[1], 
                            orders.Parameter_list.values[i].split(',')[3].replace(' ', ''), 
                            orders.Parameter_list.values[i].split(',')[4].split(), 
                            orders.Parameter_list.values[i].split(',')[5].replace(' ', '').split('calib=')[1], 
                            pathgpt,
                            snappypath,
                            snaphupath,
                            snaphuconf,
                            DirProj,
                            int(num_cores),
                            int(compute_interferogram),
                            int(compute_phase))
            # Parameters check
            if (not os.path.isfile(orders.Parameter_list.values[0].split(',')[0])):
                sys.exit('Introduce a file with a list of image paths')
            Process.coregistration_ifg()

    print('Processing time: ', (datetime.now()-time1).seconds/60)
