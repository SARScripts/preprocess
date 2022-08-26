"""Collection of utilities to work with Sentinel-1 SLC data."""

# !/usr/bin/env python
# -*- coding: utf-8 -*-

from glob import glob
import numpy as np
import zipfile
import os
from shapely.geometry import Polygon

def get_orbit_from_zip(s1_zip_path):
    """
    Returns an integer value indicating the orbit number of the
    provided Sentinel-1 SLC zip file.

    Parameters
    ----------
    s1_zip_path : str
        A string containing the path to the Sentinel-1 SLC zip file.
    Returns
    -------
    int :
        The orbit number extracted from the manifest.safe file in the zip file.
    """
    if not os.path.exists(s1_zip_path):
        raise('The specified zip file does not exist: {}'.format(s1_zip_path))
    try:
        archive = zipfile.ZipFile(s1_zip_path, 'r')
    except:
        print('Error reading the archive: {}'.format(s1_zip_path))
        return None
    for name in archive.namelist():
        if 'manifest.safe' in name:
            manifest = archive.read(name)
            index = str(manifest).find('relativeOrbitNumber')
            orbit_str = str(manifest)[index:index+70]
            orbit_number = orbit_str.split('>')[1].split('<')[0]
            return int(orbit_number)
        
    print('Orbit not found! Please check the zip file and if it contains a valid manifest.safe file.')
    return None

def check_zip_orbit(s1_zip_path,orbit):
    """
    Checks if the provided zip file matches with the orbit number provided.
    Returns again the filepath if the orbit provided and the orbit of the zip file match.
    Returns None otherwise.

    Parameters
    ----------
    s1_zip_path : str
        A string containing the path to the Sentinel-1 SLC zip file.
    orbit : int
        An integer value expressing the Sentinel-1 orbit
        
    Returns
    -------
    str or None :
        Returns the filepath provided if the orbits match, None otherwise.
    """
    try:
        orbit_number = get_orbit_from_zip(s1_zip_path)
    except Exception as e:
        raise e
    if orbit_number == orbit:
        return s1_zip_path
    else:
        return None


def get_s1_footprint(s1_zip_path):
    """
    Returns the Sentinel-1 footprint as a shapely.geometry.Polygon object.
    The code reads the manifest.safe file and looks for the <gml:coordinates> field.

    Parameters
    ----------
    s1_zip_path : str
        A string containing the path to the Sentinel-1 SLC zip file.
        
    Returns
    -------
    shapely.geometry.Polygon :
        A shapely polygon with the Sentinel-1 data footprint.
    """
    if not os.path.exists(s1_zip_path):
        raise('The specified zip file does not exist: {}'.format(s1_zip_path))
    try:
        archive = zipfile.ZipFile(s1_zip_path, 'r')
    except:
        print('Error reading the archive: {}'.format(s1_zip_path))
        return None

    for name in archive.namelist():
        if 'manifest.safe' in name:
            manifest = archive.read(name)
            index = str(manifest).find('<gml:coordinates>')
            coords = (str(manifest)[index+17:index+94])
            splitted_coords = coords.split(' ')
            S1_coverage = Polygon([(float(splitted_coords[0].split(',')[1]),float(splitted_coords[0].split(',')[0])),
                                   (float(splitted_coords[1].split(',')[1]),float(splitted_coords[1].split(',')[0])),
                                   (float(splitted_coords[2].split(',')[1]),float(splitted_coords[2].split(',')[0])),
                                   (float(splitted_coords[3].split(',')[1]),float(splitted_coords[3].split(',')[0])),
                                   (float(splitted_coords[0].split(',')[1]),float(splitted_coords[0].split(',')[0]))])
            return S1_coverage
        
    print('Coordinates not found! Please check the zip file and if it contains a valid manifest.safe file.')
    return None
            
def get_s1_processing_date(s1_zip_path):
    """
    Sentinel-1 data could be reprocessed by ESA multiple times. This could result in
    having to choose between zip files which apparently are the same, since they have
    exactly the same filename except the last 4 digits.
    See https://sentinel.esa.int/web/sentinel/technical-guides/sentinel-1-sar/products-algorithms/level-1-product-formatting
    for more information about the Sentinel-1 SLC filename formatting.
    
    Returns the Sentinel-1 processing date as a numpy.datetime64 object.
    The code reads the manifest.safe file and looks for the
    <safe:processing name="SLC Post Processing" field.

    Parameters
    ----------
    s1_zip_path : str
        A string containing the path to the Sentinel-1 SLC zip file.

    Returns
    -------
    numpy.datetime64 :
        A numpy.datetime64 object expressing the Sentinel-1 data processing time.
    """
    if not os.path.exists(s1_zip_path):
        raise('The specified zip file does not exist: {}'.format(s1_zip_path))
    try:
        archive = zipfile.ZipFile(s1_zip_path, 'r')
    except:
        print('Error reading the archive: {}'.format(s1_zip_path))
        return None
        
    for name in archive.namelist():
        if 'manifest.safe' in name:
            manifest = archive.read(name)
            index = str(manifest).find('<safe:processing name="SLC Post Processing"')
            substring = (str(manifest)[index+50:index+120])
            index = substring.find('stop=\"')
            substring = (substring)[index+6:index+32]
            date = np.datetime64(substring)
            return date
        
    print('Processing date not found! Please check the zip file and if it contains a valid manifest.safe file.')
    return None
            
def get_latest_processed_zip(s1_zips_path):
    """
    Sentinel-1 data could be reprocessed by ESA multiple times. This could result in
    having to choose between zip files which apparently are the same, since they have
    exactly the same filename except the last 4 digits.
    See https://sentinel.esa.int/web/sentinel/technical-guides/sentinel-1-sar/products-algorithms/level-1-product-formatting
    for more information about the Sentinel-1 SLC filename formatting.
    
    This function returns the most up-to-date Sentinel-1 filepath (the one with the latest processing date)
    comparing the processing dates of the zips provided.

    Parameters
    ----------
    s1_zips_path : list
        List of strings containing the paths to the Sentinel-1 SLC zip file.

    Returns
    -------
    str :
        The filepath of the latest processed zip file of the provided list.
    """
    if not isinstance(s1_zips_path,list):
        print('Error with the input data, please provide a list.')
        return None
        
    first_date = get_s1_processing_date(s1_zips_path[0])
    first_zip = s1_zips_path[0]
    if len(s1_zips_path)==1:
        return first_zip
    
    for z in s1_zips_path[1:]:
        if first_date > get_s1_processing_date(z):
            continue
        else:
            first_date = get_s1_processing_date(z)
            first_zip = z
            
    return first_zip

