# Pangeo-Enabled ESM Pattern Scaling (PEEPS) for CMIP6 Earth System Models 

This repository houses the Jupyter notebook for Pangeo-Enabled ESM Pattern Scaling (PEEPS): a flexible tool for accessing CMIP6 Earth System Model (ESM) data via Pangeo without housing the entire archive in-house, and to process and perform linear pattern scaling on the data. Due to Pangeo's capabilities, this Jupyter notebook is effectively the dataset of patterns described in Kravitz and Snyder. The data set itself is housed on zenodo (DOI:nnn) due to size.

## Relevant citation
Ben Kravitz and Abigail Snyder, Pangeo-Enabled ESM Pattern Scaling (PEEPS): A customizable dataset of emulated Earth System Model output. Submitted. 

Patterns described in this data paper are located on zenodo: DOI nnn

## setup



### packages necessary to access Pangeo data and pattern scale

`pip install` the following:
```

matplotlib

xarray

numpy

pandas

sklearn

intake

intake-esm

fsspec

seaborn

nc_time_axis

gcsfs
```


## One way to launch and work with a Jupyter notebook
```
unset XDG_RUNTIME_DIR
grep `hostname` /etc/hosts
export MYPORT=`expr $CUDA_VISIBLE_DEVICES + 8888`
export LAST_DEV=`echo $CUDA_VISIBLE_DEVICES | awk -F, '{print $NF}'`
export MYPORT=`expr $LAST_DEV + 8888`
echo $MYPORT
CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES jupyter-notebook --ip=0.0.0.0 --no-browser --port=$MYPORT
```

From the above set of commands, you will get some output that looks like a URL. Copy/paste this URL into your browser.

## Pattern scaling
The `pattern_main` Jupyter notebook contains the calls to create the set of CMIP6 linear patterns for specified variables,
the dataset described in Kravitz and Snyder. This notebook pulls data directly from Pangeo and sources the `helpers.py`
script in this repository that contains the functions for reshaping CMIP6-style data and performing linear pattern scaling.


