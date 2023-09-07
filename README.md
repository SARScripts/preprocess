# Introduction

This repository contains the necessary code to replicate the SAR2CUBE infrastructure.

A paper was presented at IGARSS 2023, you can find it here: https://hdl.handle.net/10863/36007

The main idea consist in coregistering a stack of Sentinel-1 SLC data, which will be used to generate different products on-the-fly like:
- Temporal subset
- Spatial subset
- Intensity/Amplitude
- Multilook
- Box-car filter
- Interferometry
- Pixel Selection for PSI
- Geocoding

# Installation

## ESA SNAP

The preprocessing was tested with ESA SNAP version 9 and the Sentinel toolbox (not tested with the latest version 10).
 - You can install SNAP following the official instructions here: https://step.esa.int/main/download/snap-download/ 
 - To enable the storage of the phase components, you need to recompile the Sentinel-1 toolbox, please read the discussion here:
https://forum.step.esa.int/t/snap-compiled-from-source-to-set-output-phase-true/40096/8
 - You have to configure python snappy after creating the provided Python environment:
 https://senbox.atlassian.net/wiki/spaces/SNAP/pages/50855941/Configure+Python+to+use+the+SNAP-Python+snappy+interface

## Python Environment

Since ESA Snappy, the Python bindings for SNAP are compatible only with Python 3.6, we are limited to use this old Python version.
Once SNAP 10 will be officially released, newer Python versions will be supported and these instructions will be updated.

1. Install Anaconda to manage virtual environments. You can follow the instructions [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
2. Clone the repository and get into the repo folder:
 ```
git clone https://github.com/SARScripts/preprocess
cd preprocess
```
3. Create a new conda environment with the following command:
```
conda env create -f sar2cube.yml
```
4. Once the process is complete, you can activate the environment:
```
conda activate sar2cube
```
5. If necessary, install the required libgfortran library by typing:
```
sudo apt-get install libgfortran5
```

# SAR2Cube Preprocessing Workflow
Preprocessing operations over SLC images for SAR2CUBE:

## Generate the coregistered data stack

- Update paths in CONF/config.project

- Update paths in CONF/order_file_test.csv

- Update CONF/list_images.csv


Navigate to the path of preprocessing.py and run it with the config file path:
```
python preprocessing.py path/to/project.config
```

After a successfull preprocessing, this is how the output directory should look like this (with more date pairs):

<details>
<summary>Output directory content:</summary>
    
```
├── 00_calib
│   ├── 20181010_calib.data
│   │   ├── i_IW1_VH.hdr
│   │   ├── i_IW1_VH.img
│   │   ├── i_IW1_VV.hdr
│   │   ├── i_IW1_VV.img
│   │   ├── i_IW2_VH.hdr
│   │   ├── i_IW2_VH.img
│   │   ├── i_IW2_VV.hdr
│   │   ├── i_IW2_VV.img
│   │   ├── i_IW3_VH.hdr
│   │   ├── i_IW3_VH.img
│   │   ├── i_IW3_VV.hdr
│   │   ├── i_IW3_VV.img
│   │   ├── q_IW1_VH.hdr
│   │   ├── q_IW1_VH.img
│   │   ├── q_IW1_VV.hdr
│   │   ├── q_IW1_VV.img
│   │   ├── q_IW2_VH.hdr
│   │   ├── q_IW2_VH.img
│   │   ├── q_IW2_VV.hdr
│   │   ├── q_IW2_VV.img
│   │   ├── q_IW3_VH.hdr
│   │   ├── q_IW3_VH.img
│   │   ├── q_IW3_VV.hdr
│   │   ├── q_IW3_VV.img
│   │   ├── tie_point_grids
│   │   └── vector_data
│   ├── 20181010_calib.dim
├── 00_data
│   ├── S1A_IW_SLC__1SDV_20230710T170721_20230710T170748_049364_05EFA1_25AD.SAFE
│   │   ├── S1A_IW_SLC__1SDV_20230710T170721_20230710T170748_049364_05EFA1_25AD.SAFE-report-20230710T175811.pdf
│   │   ├── annotation
│   │   ├── manifest.safe
│   │   ├── measurement
│   │   │   ├── s1a-iw1-slc-vh-20230710t170722-20230710t170747-049364-05efa1-001.tiff
│   │   │   ├── s1a-iw1-slc-vv-20230710t170722-20230710t170747-049364-05efa1-004.tiff
│   │   │   ├── s1a-iw2-slc-vh-20230710t170723-20230710t170748-049364-05efa1-002.tiff
│   │   │   ├── s1a-iw2-slc-vv-20230710t170723-20230710t170748-049364-05efa1-005.tiff
│   │   │   ├── s1a-iw3-slc-vh-20230710t170721-20230710t170746-049364-05efa1-003.tiff
│   │   │   └── s1a-iw3-slc-vv-20230710t170721-20230710t170746-049364-05efa1-006.tiff
│   │   ├── preview
│   │   └── support
├── 01_ifg_calib
│   ├── 20181010_20181010_SLC_calib_Coregistered.data
│   │   ├── fep_VH_10Oct2018_10Oct2018.hdr
│   │   ├── fep_VH_10Oct2018_10Oct2018.img
│   │   ├── tgp_VH_10Oct2018_10Oct2018.hdr
│   │   ├── tgp_VH_10Oct2018_10Oct2018.img
│   │   ├── tie_point_grids
│   │   └── vector_data
│   ├── 20181010_20181010_SLC_calib_Coregistered.dim
├── 01_slc_calib
│   ├── 01_slc_calib/20181010_20181010_SLC_calib_Coregistered.data/
│   │   ├── i_VH_10Oct2018.hdr
│   │   ├── i_VH_10Oct2018.img
│   │   ├── i_VV_10Oct2018.hdr
│   │   ├── i_VV_10Oct2018.img
│   │   ├── phase.hdr
│   │   ├── phase.img
│   │   ├── q_VH_10Oct2018.hdr
│   │   ├── q_VH_10Oct2018.img
│   │   ├── q_VV_10Oct2018.hdr
│   │   ├── q_VV_10Oct2018.img
│   │   ├── tie_point_grids
│   │   ├── vector_data
│   ├── 20181010_20181010_SLC_calib_Coregistered.dim
├── 03_gc
│   ├── elevation.hdr
│   ├── elevation.img
│   ├── ifg_gc_20181010_20230710.data
│   │   ├── coh_VH_10Oct2018_10Oct2018.hdr
│   │   ├── coh_VH_10Oct2018_10Oct2018.img
│   │   ├── coh_VV_10Oct2018_10Oct2018.hdr
│   │   ├── coh_VV_10Oct2018_10Oct2018.img
│   │   ├── elevation.hdr
│   │   ├── elevation.img
│   │   ├── fep_VH_10Oct2018_10Oct2018.hdr
│   │   ├── fep_VH_10Oct2018_10Oct2018.img
│   │   ├── fep_VV_10Oct2018_10Oct2018.hdr
│   │   ├── fep_VV_10Oct2018_10Oct2018.img
│   │   ├── i_ifg_VH_10Oct2018_10Oct2018.hdr
│   │   ├── i_ifg_VH_10Oct2018_10Oct2018.img
│   │   ├── i_ifg_VV_10Oct2018_10Oct2018.hdr
│   │   ├── i_ifg_VV_10Oct2018_10Oct2018.img
│   │   ├── orthorectifiedLat.hdr
│   │   ├── orthorectifiedLat.img
│   │   ├── orthorectifiedLon.hdr
│   │   ├── orthorectifiedLon.img
│   │   ├── q_ifg_VH_10Oct2018_10Oct2018.hdr
│   │   ├── q_ifg_VH_10Oct2018_10Oct2018.img
│   │   ├── q_ifg_VV_10Oct2018_10Oct2018.hdr
│   │   ├── q_ifg_VV_10Oct2018_10Oct2018.img
│   │   ├── tgp_VH_10Oct2018_10Oct2018.hdr
│   │   ├── tgp_VH_10Oct2018_10Oct2018.img
│   │   ├── tgp_VV_10Oct2018_10Oct2018.hdr
│   │   ├── tgp_VV_10Oct2018_10Oct2018.img
│   │   ├── tie_point_grids
│   │   └── vector_data
│   ├── ifg_gc_20181010_20230710.dim
│   ├── incid_angle.hdr
│   ├── incid_angle.img
│   ├── incident_angle.hdr
│   ├── incident_angle.img
│   ├── orthorectifiedLat.hdr
│   ├── orthorectifiedLat.img
│   ├── orthorectifiedLon.hdr
│   ├── orthorectifiedLon.img
├── baselines.csv
├── baselines_filtered.csv
├── process_log.csv
```

</details>

# OpenDataCube data indexing

The generated data stack can be indexed in OpenDataCube with the following instructions. Run the following commands in the ODC folder.

1. Since ODC doesn't allow to index not-georeferenced data, or data on an irregular grid, we generate VRTs https://gdal.org/drivers/raster/vrt.html , assigning to each file a fake projection. Please modify the script with your paths and then run it:
```
python prepare_sar2cube_vrts.py
```

2. With the generated VRTs, we can create the corresponding ODC dataset files. You can do so using the script: `prepare_sar2cube_datasets.py`:
```
python prepare_sar2cube_datasets.py --help
usage: prepare_sar2cube_datasets.py [-h] [-t TEMPLATE] [-p PRODUCT] [-b BASELINES] [-o OVERWRITE] source target

Create the eo3 ODC datasets files given a folder with the VRTs. Each VRTs must have a single band.

positional arguments:
  source                               Path to the input folder with VRTs files.
  target                               Path to write the output datasets files.

optional arguments:
  -h, --help                           show this help message and exit
  -t TEMPLATE, --template TEMPLATE     Path to the template yaml dataset file.
  -p PRODUCT, --product PRODUCT        Path to the yaml product file.
  -b BASELINES, --baselines BASELINES  Path to the baselines file.
  -o OVERWRITE, --overwrite OVERWRITE  Overwrite the destination file if already exists. (default: 0).

Author: Michele Claus, michele.claus@eurac.edu
```

3. You can finally add the product to ODC:
```
datacube product add SAR2Cube_ASC_117_datacube.yaml
```

4. Add the datasets to the product:
```
find ./datasets/ -type f -name "*.yaml"| xargs --max-procs=8 -n1 datacube dataset add -p SAR2Cube_ASC_117_datacube
```