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
   - Install the dependencies for handling I/O: `python3 -m pip install "xarray[io]"`

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

## Reference [![DOI](https://zenodo.org/badge/304547682.svg)](https://doi.org/10.5281/zenodo.13928317):
This code is developed by Vikrant V. Jadhav (Universität Bonn). 


If you use the code, please refer to the software repository and Jadhav et al. JApA, 42, 89, (2021) 
```text
@software{jadhav_2024_13928317,
  author       = {Jadhav, Vikrant V.},
  title        = {Binary\_SED\_Fitting: v3.4.2},
  month        = oct,
  year         = 2024,
  publisher    = {Zenodo},
  version      = {v3.4.2},
  doi          = {10.5281/zenodo.13928317},
  url          = {https://doi.org/10.5281/zenodo.13928317}
}
@ARTICLE{2021JApA...42...89J,
       author = {{Jadhav}, Vikrant V. and {Pandey}, Sindhu and {Subramaniam}, Annapurni and {Sagar}, Ram},
        title = "{UOCS. IV. Discovery of diverse hot companions to blue stragglers in the old open cluster King 2}",
      journal = {Journal of Astrophysics and Astronomy},
     keywords = {Open star clusters (1160), blue straggler stars (168), extreme horizontal branch stars (513), B subdwarf stars (129), ultraviolet astronomy (1736), spectral energy distribution (2129), binary stars (154), Astrophysics - Solar and Stellar Astrophysics, Astrophysics - Astrophysics of Galaxies},
         year = 2021,
        month = oct,
       volume = {42},
       number = {2},
          eid = {89},
        pages = {89},
          doi = {10.1007/s12036-021-09746-y},
archivePrefix = {arXiv},
       eprint = {2102.13375},
 primaryClass = {astro-ph.SR},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2021JApA...42...89J},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

### Used in the following publications
1. Saketh et al. AJ, 168, 97, (2024) 

1. Pal et al. ApJL, 970, L39, (2024) 

1. Jadhav et al. A&A, 688, A152, (2024) 

1. Panthi and Vaidya et al. MNRAS, 527, 10335, (2024) 

1. Panthi et al. MNRAS, 527, 8325, (2024) 

1. Panthi et al. MNRAS, 525, 1311, (2023) 

1. Jadhav et al. A&A, 676, A47, (2023) 

1. Panthi et al. MNRAS, 516, 5318, (2022) 

1. Rao et al. MNRAS, 516, 2444, (2022) 

1. Vaidya et al. MNRAS, 511, 2274, (2022) 

1. Pandey et al. MNRAS, 507, 2373, (2021) 

1. Jadhav et al. JApA, 42, 89, (2021) 

1. Subramaniam et al. JApA, 41, 45, (2020) 

1. Jadhav et al. ApJ, 886, 13, (2019) 

1. Sindhu et al. ApJ, 882, 43, (2019)
