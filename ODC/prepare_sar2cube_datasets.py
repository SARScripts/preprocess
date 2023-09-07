import argparse
import logging
import sys
import pathlib
import datetime
from glob import glob
import uuid
import os
import rasterio
from ruamel import yaml
import numpy as np
from joblib import Parallel
from joblib import delayed as joblibDelayed
import warnings

warnings.filterwarnings(action='ignore') #Ignore warnings from ODC when a dataset is already in the datacube
yaml.default_flow_style = None

# epilogue to display at the end of each parser
EPILOGUE = 'Author: Michele Claus, michele.claus@eurac.edu'

# default values
DEFAULT = '(default: %(default)s)'


def args_parser():

    # define command line argument parser
    parser = argparse.ArgumentParser(
        description='Create the eo3 ODC datasets files given a folder with the VRTs. Each VRTs must have a single band.',
        epilog=EPILOGUE,
        formatter_class=lambda prog: argparse.RawDescriptionHelpFormatter(
            prog, max_help_position=50, indent_increment=2))
    
    # positional arguments

    parser.add_argument('source', type=pathlib.Path,
                        help='Path to the input folder with VRTs files.')

    parser.add_argument('target', type=pathlib.Path,
                        help='Path to write the output datasets files.')
    
    parser.add_argument('-t','--template', type=pathlib.Path,
                        help='Path to the template yaml dataset file.')
    
    parser.add_argument('-p','--product', type=pathlib.Path,
                        help='Path to the yaml product file.')
        
    parser.add_argument('-b','--baselines', type=pathlib.Path,
                        help='Path to the baselines file.')
    
    parser.add_argument('-o', '--overwrite', type=int,
                        help='Overwrite the destination file if already exists. {}.'.format(DEFAULT),
                        default=0)
    return parser



def get_date(path):
    date = path.split('/')[-1]
    date = date[:4] + '-' + date[4:6] + '-' + date[6:]
    return date

def get_title(yaml_product_path):
    with open(yaml_product_path) as file:
        try:
            product = yaml.safe_load(file)   
        except Exception as e:
            raise e 
    return product['metadata']['title']

def get_name(yaml_product_path):
    with open(yaml_product_path) as file:
        try:
            product = yaml.safe_load(file)   
        except Exception as e:
            raise e
    return product['name']

def prepare_SAR2CUBE_dataset(folder,name,title,template,bands,output_path,baselines,overwrite):

    images = list(folder.glob('*.vrt'))
    if images == []:
        print(f'No .vrt images found in {f}!')
        return

    with open(template) as file:
        try:
            dataset = yaml.safe_load(file)   
        except Exception as e:
            raise e

    ds = rasterio.open(images[0])
    transform = list(ds.transform)
    shape = list(ds.shape)
    crs = str(ds.crs)

    id = str(uuid.uuid4())
    
    try:
        date = get_date(folder.as_posix())
    except Exception as e:
        raise e

    output_file = output_path + '/' + date + '.yaml'
    if os.path.exists(output_file):
        if not overwrite:
            return

    dataset['id'] = str(uuid.uuid4())
    dataset['product'] = {'name':name}
    dataset['title'] = title
    dataset['crs'] = crs
    dataset['grids']['default']['shape'] = shape
    dataset['grids']['default']['transform'] = transform

    bands_dict = {}
    for b in bands:
        bands_dict[b] = {'path':[x.as_posix() for x in images if b in x.as_posix()][0]}
    dataset['measurements'] = bands_dict
    dataset['properties']['datetime'] = datetime.datetime(int(date[:4]),int(date[5:7]),int(date[8:10]),tzinfo=datetime.timezone.utc)
    dataset['accessories']['metadata:baselines']['path'] = baselines
    with open(output_file, 'w') as file:
        outputs = yaml.dump(dataset, file)
    return

def main():
    
    # -------------------------------------------------------------------------
    # INITIALIZE LOGGING ------------------------------------------------------
    # -------------------------------------------------------------------------
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )

    # define command line argument parser
    parser = args_parser()

    # parse command line arguments
    args = sys.argv[1:]
    if not args:
        parser.print_help()
        sys.exit()
    else:
        args = parser.parse_args(args)

    # check if source file exists
    if not args.source.exists():
        logging.error(f'Source folder {args.source} does not exist!')
        sys.exit()
        
    # check if target folder exists
    if not args.target.exists():
        logging.info(f'Target folder {args.target} does not exist, creating it...')
        args.target.mkdir()
        
    # check if template file exists
    if not args.template.exists():
        logging.error(f'Template file {args.template} does not exist! You must provide one.')
        sys.exit(0)
        
    # check if template file exists
    if not args.baselines.exists():
        logging.error(f'Baselines file {args.baselines} does not exist! You must provide one.')
        sys.exit(0)
        
    INPUT_PATH = args.source.as_posix()
    OUTPUT_PATH = args.target.as_posix()
    PRODUCT = args.product.as_posix()
    BASELINES = args.baselines.as_posix()
    
    title = get_title(PRODUCT)
    name = get_name(PRODUCT)

    
    BANDS = ['i_VV', 'q_VV', 'i_VH', 'q_VH', 'DEM', 'LIA', 'phase' ,'grid_lon', 'grid_lat']

    folders = list(args.source.glob('2*'))
    for f in folders:
        prepare_SAR2CUBE_dataset(f,name,title,args.template,BANDS,OUTPUT_PATH,BASELINES,args.overwrite)


            
    # Parallel(n_jobs=-1)(joblibDelayed(prepare_SAR2CUBE_dataset)(i,im,transform,shape,crs,title) for i,im in enumerate(images_list))
    # for i,im in enumerate(images_list):
        # prepare_ADO_dataset(i,im,transform,shape,crs)

                    
if __name__ == '__main__':
    main()