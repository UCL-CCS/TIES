The easiest way to start is to use the latest 
published conda-forge package. Create a new 
environment and use: 

```shell
conda install ties
```

## Development version

Whereas most of the dependancies can be installed with pip, 
ambertools has to be either compiled. 
The easiest way is to use the conda-forge environment:   

```shell
mamba env create -f environment.yml
conda activate ties 
pip install --no-deps . 
```