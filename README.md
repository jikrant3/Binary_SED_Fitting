# Binary_SED_Fitting v3.4.0
Fitting stellar SEDs for single and multi-component systems.

## Documentation at [jikrant3.github.io/Binary_SED_Fitting](https://jikrant3.github.io/Binary_SED_Fitting)

## Requirements
- `python` >= 3.12.2

- `astropy` >= 6.0.1

- `matplotlib` >= 3.8.4

- `numpy` >= 1.26.4

- `pandas` >= 2.2.1

- `scipy` >= 1.13.0

- `xarray` >= 2024.3.0
   - Installation: `python3 -m pip install xarray`
   - Install optional dependencies for handling I/O: `python3 -m pip install "xarray[io]"`

## Installation
- Add the `binary_sed_fitting.py` to your working directory.
   - ⚠️ Edit the `DIR_MODELS` in the `binary_sed_fitting.py` file according to where your model files are located

- Download the required models from [Google Drive](https://drive.google.com/drive/folders/1UdpMiPVj-q91IpcmcLBmSXmYh_iDgwTq?usp=sharing). 
  More information about how to make these model files is given in [models_and_tools](https://github.com/jikrant3/models_and_tools/).
   
   - Currently supported spectral models:

      - `kurucz_synthetic_photometry.nc` (Kurucz 2003IAUS..210P.A20C)

      - `koester_synthetic_photometry.nc` (Koester 2010MmSAI..81..921K) 

      - `kurucz_uvblue_synthetic_photometry.nc` (Hybrid Kurucz_UVBLUE 2005ApJ...626..411R)

   - Isochrones: `master_isochrone.csv` ([Parsec isochrones](http://stev.oapd.inaf.it/cmd))

   - WD cooling curves: `master_Bergeron_WD.csv` ([Bergeron models](https://www.astro.umontreal.ca/~bergeron/CoolingModels/))

## Caution
⚠️ The software is under active development. 

The code has been designed for 3 component SEDs. 
However, only double component fitting is thoroughly tested.

## Reference:
This code is developed by Vikrant V. Jadhav (Universität Bonn). 

If you use the code, please refer to Jadhav et al. (2021), UOCS. IV. Discovery of diverse hot companions to blue stragglers in the old open cluster King 2, Journal of Astrophysics and Astronomy, Volume 42, Issue 2, article id.89 (https://doi.org/10.1007/s12036-021-09746-y)