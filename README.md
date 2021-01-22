# Binary_SED_Fitting
Fitting SEDs of binary stars. The code is optimised for chi2 SED fitting of a double component star. It is assumed that one component has already been fit (using VOSA for our case).

## Disclaimer
Please make sure to familiarise with the inputs and working of the algotithm. The accuracy of the results depends solely on the user. <br />
The contents of this repository are part of an ongoing research project. <br />
The form and function are generally limited to what is required for that project and may not be appropriate for your project. <br />
Feedback is welcomed, but feature requests will probably not be honored. <br />
The contents of this repository might change dramatically without notice. <br />

## Acknowledgement:
This code is developed by Vikrant V. Jadhav and Sindhu Pandey. 
If you use the code, please refer to Jadhav et al. (2021, under review), UOCS. IV. Discovery of diverse hot companions to blue stragglers in the old open cluster King 2.

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
