from osgeo import gdal, osr
import numpy as np
import xarray as xr
from glob import glob
import os
import shutil
from tqdm import tqdm
import fileinput
from pathlib import Path
import logging
import sys
import rasterio
import rioxarray

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout)
    ]
)

# INPUTS:
data_folder = "/mnt/large_drive/work_spaces/mclaus/OUTPUT_ASC_STyrol_4images_PHASE_2/01_slc_calib/"
geocoding_folder = "/mnt/large_drive/work_spaces/mclaus/OUTPUT_ASC_STyrol_4images_PHASE_2/03_gc/"
baselines_folder = "/mnt/large_drive/work_spaces/mclaus/OUTPUT_ASC_STyrol_4images_PHASE_2/"
vrts_folder = "/mnt/large_drive/work_spaces/mclaus/OUTPUT_ASC_STyrol_4images_PHASE_2/04_vrts"
# vrt template with added sample projection
vrt_template = "vrt_template.vrt"

# Check that folders and files exist:
def check_exists(list_to_check):
    for f in list_to_check:
        if not Path(f).exists():
            raise Exception(f"{f} does not exist!")

check_exists([data_folder,
              geocoding_folder,
              baselines_folder,
              vrt_template
             ])

data_folder_pattern = "*.data"
phase_filename = "phase.img"
grid_lat_filename = "orthorectifiedLat.img"
grid_lon_filename = "orthorectifiedLon.img"
dem_filename = "elevation.img"
lia_filename = "incident_angle.img"
baselines_filename = "baselines.csv"

def get_date_string(path):
    date = path.split('/')[-1].split('.')[0]
    return date

def get_date_string_snap(path):
    date = path.split('/')[-1].split('_')[1]
    return date

# Create VRTs
vrt_folders_path = Path(vrts_folder) 
if not vrt_folders_path.exists():
    try:
        vrt_folders_path.mkdir()
    except Exception as e:
        logging.error(e)

for f in Path(data_folder).glob(data_folder_pattern):
    date = get_date_string_snap(f.as_posix())

    i_VV_file = list(f.glob("i_VV*.img"))[0]
    q_VV_file = list(f.glob("q_VV*.img"))[0]
    i_VH_file = list(f.glob("i_VH*.img"))[0]
    q_VH_file = list(f.glob("q_VH*.img"))[0]

    if (f / phase_filename).exists() and i_VV_file.exists() and q_VV_file.exists() and i_VH_file.exists() and q_VH_file.exists():
        vrt_out_folder = vrt_folders_path / date
        print(vrt_out_folder)
        if not vrt_out_folder.exists():
            try:
                vrt_out_folder.mkdir()
            except Exception as e:
                print(e)
                logging.error(e)
        # Copy the template vrts to the destination folder
        shutil.copyfile(vrt_template,vrt_out_folder / "i_VV.vrt")
        shutil.copyfile(vrt_template,vrt_out_folder / "q_VV.vrt")
        shutil.copyfile(vrt_template,vrt_out_folder / "i_VH.vrt")
        shutil.copyfile(vrt_template,vrt_out_folder / "q_VH.vrt")
        shutil.copyfile(vrt_template,vrt_out_folder / "grid_lat.vrt")
        shutil.copyfile(vrt_template,vrt_out_folder / "grid_lon.vrt")
        shutil.copyfile(vrt_template,vrt_out_folder / "phase.vrt")
        shutil.copyfile(vrt_template,vrt_out_folder / "DEM.vrt")
        shutil.copyfile(vrt_template,vrt_out_folder / "LIA.vrt")

        ds = rioxarray.open_rasterio(f / i_VV_file, chunks={})
        x_size = str(len(ds.x))
        y_size = str(len(ds.y))

        vrt_general_files = ["i_VV.vrt",
                             "q_VV.vrt",
                             "i_VH.vrt",
                             "q_VH.vrt",
                             "grid_lat.vrt",
                             "grid_lon.vrt",
                             "phase.vrt",
                             "DEM.vrt",
                             "LIA.vrt"]

        for gen_file in vrt_general_files:          
            with fileinput.FileInput(vrt_out_folder / gen_file, inplace=True) as file:
                for line in file:
                    print(line.replace("X_SIZE", x_size), end='')

            with fileinput.FileInput(vrt_out_folder / gen_file, inplace=True) as file:
                for line in file:
                    print(line.replace("Y_SIZE", y_size), end='')

            with fileinput.FileInput(vrt_out_folder / gen_file, inplace=True) as file:
                for line in file:
                    if "grid_lat" in gen_file:
                        filepath = Path(geocoding_folder) / grid_lat_filename
                    elif "grid_lon" in gen_file:
                        filepath = Path(geocoding_folder) / grid_lon_filename
                    elif "phase" in gen_file:
                        filepath = f / phase_filename
                    elif "DEM" in gen_file:
                        filepath = Path(geocoding_folder) / dem_filename
                    elif "LIA" in gen_file:
                        filepath = Path(geocoding_folder) / lia_filename
                    elif "i_VV" in gen_file:
                        filepath = i_VV_file
                    elif "q_VV" in gen_file:
                        filepath = q_VV_file
                    elif "i_VH" in gen_file:
                        filepath = i_VH_file
                    elif "q_VH" in gen_file:
                        filepath = q_VH_file
                    print(line.replace("FILEPATH", filepath.as_posix()), end='')