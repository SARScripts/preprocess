# Instructions

1. Install Anaconda to manage virtual environments. You can follow the instructions [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
2. Clone the repository and get into the repo folder:
 ```
        git clone https://github.com/SARScripts/preprocess
        cd preprocess
```
3. Create a new conda environment with the following command:
```
        conda env create -f environment.yml
```
4. Once the process is complete, you can activate the environment:
```
        conda activate s2cprod
```
5. Install the required libgfortran library by typing:
```
sudo apt-get install libgfortran5
```

# preprocess
Preprocessing operations over SLC images for SAR2CUBE

# Usage
Update CONF/config.project with local paths

Update CONF/order_file_test.csv

Update CONF/list_images.csv


Navigate to the path of preprocessing.py and run it with the config file path:

python preprocessing.py path/to/project.config
