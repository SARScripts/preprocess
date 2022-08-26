"""Sample script for generating the Sentinel-1 zip list that will be processed in SAR2Cube."""

# !/usr/bin/env python
# -*- coding: utf-8 -*-

from s1_slc_utils import *
from joblib import Parallel
from joblib import delayed as joblibDelayed
import geopandas as gpd
import logging
import sys

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("sar2cube_imagelist.log"),
        logging.StreamHandler(sys.stdout)
    ]
)

def main():
    s1_zips_path = '/mnt/CEPH_BASEDATA/SATELLITE/SENTINEL/SENTINEL-1/IW/SLC/'
    
    logging.info('[*] Looking for S1 zip files in {}'.format(s1_zips_path))
    
    s1_files_list = glob(s1_zips_path + '*IW_SLC__1SDV*.zip')
    
    logging.info('[*] Found {} S1 SLC zip files.'.format(len(s1_files_list)))
    
    orbit = 117
    output_list = []

    # Get the list of zips for the given orbit.
    output_list = Parallel(n_jobs=64)(joblibDelayed(check_zip_orbit)(zip_path,orbit) for zip_path in s1_files_list)
    
    filtered_output_list = [i for i in output_list if i is not None]
    # If you want to get only data for specific years, you can use and modify the following line:
    # filtered_output_list = [i for i in filtered_output_list if ('1SDV_2018' in i) | ('1SDV_2019' in i) | ('1SDV_2020' in i) | ('1SDV_2021' in i) | ('1SDV_2022' in i)]
    
    logging.info('[*] Found {} S1 SLC zip files for the orbit {}'.format(len(filtered_output_list),orbit))

    # Order the list based on the sensing date

    sensing_dates = [i.split('SDV_')[1].split('_')[0] for i in filtered_output_list]
    sorted_zipped_lists = sorted(zip(filtered_output_list,sensing_dates), key = lambda x: x[1])
    sorted_output_list = list(list(zip(*sorted_zipped_lists))[0])

    # Filter the zips and keep only the ones covering South Tyrol

    st_gdf = gpd.read_file('./AUXDATA/ST.shp')
    st_gdf_4326 = st_gdf.to_crs('EPSG:4326')
    st_bounds = st_gdf_4326.bounds
    st = st_gdf_4326.geometry.values[0]

    st_list = []

    for s1_zip_path in sorted_output_list:
        s1_coverage = get_s1_footprint(s1_zip_path)
        if st.intersects(s1_coverage):
            st_list.append(s1_zip_path)

    # Check if there are duplicate dates and keep only the most up to date

    for i,f in enumerate(st_list):
        current_date = f.split('IW_SLC__1SDV_')[1][:46]
        # same_date_zips = glob(S1_zips_path + '*{}*.zip'.format(current_date))
        same_date_zips = [x for x in st_list if current_date in x]
        if len(same_date_zips)>1:
            zip_to_keep = get_latest_processed_zip(same_date_zips)
            zips_to_remove = [x for x in same_date_zips if zip_to_keep not in x]
            for z in zips_to_remove:
                st_list.remove(z)
        else:
            continue

    # The DSC orbit over South Tyrol is covered by two zip files. If only one is present, skip this date.

    st_list_filtered = st_list.copy()

    if orbit==168:
        for i,f in enumerate(st_list):
            current_date = f.split('IW_SLC__1SDV_')[1].split('_')[0][:8]
            count = 0
            for f2 in st_list_filtered:
                other_date = f2.split('IW_SLC__1SDV_')[1].split('_')[0][:8]
                if current_date == other_date:
                    count += 1
                if count == 2:
                    continue
                if count > 2:
                    print('This date is problematic: ', current_date)
                    continue
            if count < 2: # We need to remove this zip
                st_list_filtered.remove(f)
                
    logging.info('[*] Found {} S1 SLC zip files for the orbit {} over South Tyrol'.format(len(st_list_filtered),orbit))

    output_lines = [str(i+1) + ';' + f + '\n' for i,f in enumerate(st_list_filtered)]
    output_lines = ['Index;Image\n'] + output_lines
    
    filename = 'images_ASC_2016_2022_SouthTyrol.csv'
    
    with open(filename,'w') as file:
        file.writelines(output_lines)

    logging.info('[*] Zip list successfully written to {}.\nThe file can now be used for the SAR2Cube preprocessing pipeline.'.format(filename))

if __name__=='__main__':
    main()