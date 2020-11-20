# Binary_SED_Fitting
Fitting SEDs of binary stars. The code is optimised for chi2 SED fitting of a double component star. It is assumed that one component has already been fit (using VOSA for our case).


File structure:
``` bash
|   Binary_SED_fitting.ipynb                # Fitting binary SEDs and plotting
│   Binary_SED_fitting.py                   
│   Creating_VOSA_and_model_files.ipynb     # Creating files compatible with VOSA 
│   Creating_VOSA_and_model_files.py
│   LICENSE
│   log_file.csv                            # log file
│   README.md
│
├───data
│   │   example_isochrone.txt               
│   │   example_photomety_file.csv          
│   │   example_VOSA_input_file.txt
│   │
│   └───vosa_results_38604
│       ├───objects
│       │   └───WOCS2002
│       │       ├───bestfitp
│       │       │       WOCS2002.bfit.phot.dat
│
├───models
│       Koe_logg7.0.csv
│       Koe_logg8.0.csv
│       Koe_logg9.0.csv
│       Kr_logg3.0_Zm05.csv
│       Kr_logg3.0_Zp00.csv
│       ...
└───outputs
    ├───chi_files
    │       WOCS2002_ChiSqur_logg_B7.0_00_finer_Kr_Koe_Koe_1.csv
    │       WOCS2002_ChiSqur_logg_B7.0_00_rough_Kr_Koe_Koe_1.csv
    │
    ├───finer_Kr_Koe
    │       WOCS2002_14500_logg7.0_Z00_Koe_1.png
    │
    └───rough_Kr_Koe
            WOCS2002_35000_logg7.0_Z00_Koe_1.png
```
