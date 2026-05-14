# Agent-based model of Hepatitis B in Thai-Myanmar border region

This repository contains the code used for our study of Hepatitis B in Thai-Myanmar border region. A preprint version of the report can be found here: https://arxiv.org/abs/2511.16096

The model is built upon an existing framework of dynamic population & disease transmission described in https://www.sciencedirect.com/science/article/pii/S175543651500081X . 

## Repository structure

- `hepb/`: the python package of the model 

- `main/`: 
  - `data/`: input data used in the model, mainly growth rates and life tables
  - `output/`: model output in hd5 format

- `rscripts/`: R scripts used to analyse and plot figures from model output

## Generating outputs

### Installing required packages

Python version 3.10 or later is needed. Additionally, please install required packages using

```sh
python3 -m pip install -r requirement.txt
```

Moreover, `population` package from https://bitbucket.org/ngeard/simodd-pop/src/master/ and `disease` package from https://bitbucket.org/ngeard/simodd-dis/src/master/ must be downloaded and put in the top level directory. Additionally, the local path of the repository in the first line of `main/main-combo.py` needs to be updated.

### Setting the parameters

Parameters are set in `main/params.py` and `main-combo.py`. For single value parameters, simply edit them in `params.py` file.

Parameters with multiple values can be specified in `main-combo.py`.

### Running simulations

Simulations are run by simply running the `main-combo.py`

```sh
python main-combo.py
```

## Analysing output

R scripts used to read the hd5 files, analyse the output and plot the results are provided. The provided scripts were specifically written for the scenarios investigated in our study and thus might need modifying if the parameters in `main-combo.py` are significantly changed.