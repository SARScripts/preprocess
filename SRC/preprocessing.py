#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 23:21:09 2020

@author: paco
"""

import sys
import os
import subprocess
import pandas as pd
from datetime import datetime 
import glob

#configfile = sys.argv[1]
configfile = '/media/guadarrama/PROYECTOS/SAR2CUBE/Git/CONF/project.conf'
#FUNCIÓN NO UTILIZADA
def getpath (path): #Controls separator for paths in different OS
    splitpath = path.split('/')
    splitlist = [i.split('\\')[0] for i in splitpath]
    return os.path.join(*[i for i in splitlist])

def controlslash(path): #Control OS path generation
    return path.replace('\\', '/')

# Configuration variables from config file
try:
    in_file = open(configfile, 'r')
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
        # if "OUTPUT_DIRECTORY" in line:
        #     OutDir = controlslash(line.split('=')[1].strip())
        #     print ("OUTPUT_DIRECTORY:" + OutDir)
            #VALORAR SI CACHE Y CPU SON VARIABLES A TENER EN CUENTA
        if "CACHE" in line:
            CACHE = line.split('=')[1].strip()
        if "CPU" in line:
            CPU = line.split('=')[1].strip()
finally:
        in_file.close()
        
sys.path.append(os.path.join(DirProj, 'SRC'))
import preprocessing_utils as preproc

orders = pd.read_csv(OrderFile, sep=';', decimal=',')

#TO-DO CHECK THE PARAMETERS INTRODUCED CORRESPOND TO EACH PROCESSING NEEDS
#'alignment_interferogram' list_of_images (essential), masterimage_date (if not specified, all combinations processed), AOI (if not specified)
#'generate_pairs_list' list_of_images (essential)


#AHORA SE CREA UN OBJETO POR CADA ORDEN DE PROCESAMIENTO, PERO ESOS OBJETOS NO PUEDEN SER IGUALES PQ LOS PROCESAMIENTOS SON DIFERENTES CON ARGUMENTOS DISTINTOS
for i in range(len(orders)):
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
                        DirProj)
        #Parameters check
        if (not os.path.isfile(orders.Parameter_list.values[0].split(',')[0])):
            sys.exit('Introduce a file with a list of image paths')
        Process.alignment_ifg ()

    if orders.Process_type[i] == 'generate_list':
        if Process.processdf is not None:
        #if Process.imagelistfile == orders.Parameter_list.values[0].split(',')[0]:
            #TO-DO Parameters check Process.processdf must exist and aligned images from list
            for im in Process.processdf['Outputfiles']:
                if not os.path.isfile(im):
                    sys.exit('Outputs from coregistration processing not found (' + im + ')')
            Process.maxbasetemp = float(orders.Parameter_list.values[i].split(',')[1].replace(' ', ''))
            Process.maxbaseperp = float(orders.Parameter_list.values[i].split(',')[2].replace(' ', ''))
            Process.generate_baselinelist()
        else:
            sys.exit('Image list does not match with alignment preprocess list')
            
    if orders.Process_type[i] == 'generate_interferogram':
        if Process.baselinefiltered is not None:
        #if Process.imagelistfile == orders.Parameter_list.values[0].split(',')[0]:
            Process.generate_interferogram()
        else:
            sys.exit('Generate list not found')

