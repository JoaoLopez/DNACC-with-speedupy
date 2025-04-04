# Sphere Potentials Experiment

## Experiment dependencies
This experiment has the following dependencies: **numpy, scipy, cython**

To install them, execute: **pip install numpy, scipy, cython**

To run the experiment, make sure you have the appropriate release of SpeeduPy installed. Instructions on which release to use and how to install it can be found [here](https://github.com/dew-uff/memoization/blob/main/README.md#reproducing-the-article-analyses). Ensure that the **speedupy/** folder is located in the same directory as the main experiment script.

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

- sphere_potentials_5.py
- sphere_potentials_10.py
- sphere_potentials_15.py
- sphere_potentials_25.py
- sphere_potentials.py

To execute a trial, type:

```bash
python speedupy/setup_exp/setup.py SCRIPT_NAME.py
python SCRIPT_NAME.py --exec-mode manual
```
