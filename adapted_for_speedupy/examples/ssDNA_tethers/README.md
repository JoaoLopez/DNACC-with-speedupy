# ssDNA Tethers Experiment

## Experiment dependencies
This experiment has the following dependencies: **numpy, scipy, cython**

To install them, execute: **pip install numpy, scipy, cython**

## Experiment Setup
Before executing the experiment, it is necessary to set it up, executing:
```bash
mkdir dnacc/
cp -r ../../dnacc_speedupy/* dnacc/* 
cd dnacc
cython generic.pyx
python setup.py build_ext --inplace
cd ..
```

## Trials Used in the Article
The five trials used on the memoization techniques article are:

- ssDNA_tethers_1.py
- ssDNA_tethers_2.py
- ssDNA_tethers_3.py
- ssDNA_tethers_4.py
- ssDNA_tethers.py

To execute a trial, type:

```bash
mkdir dnacc/
cp -r ../../dnacc_speedupy/* dnacc/* 
python speedupy/setup_exp/setup.py SCRIPT_NAME.py
python SCRIPT_NAME.py --exec-mode manual
```
