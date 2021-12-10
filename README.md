# Binary_SED_Fitting
Fitting SEDs of binary stars. The code is optimised for chi2 SED fitting of a double component star. It is assumed that one component has already been fit (using VOSA for our case).

## Installation
The code is entirely made of Jupyter notebooks.
The code is tested on `python 3.9.5`, `numpy 1.20.3`, `pandas 1.3.4`, `scipy 1.7.1` and `matplotlib 3.4.3`
#### Synthetic photometry model files
Currently, the code has Koester models and Kurucz models (in compressed form). 
If user requires more model files, they can create them using `Creating_VOSA_and_model_files_v1.1.ipynb`

## Reference:
This code is developed by Vikrant V. Jadhav. 
If you use the code, please refer to Jadhav et al. (2021), UOCS. IV. Discovery of diverse hot companions to blue stragglers in the old open cluster King 2, Journal of Astrophysics and Astronomy, Volume 42, Issue 2, article id.89 (https://doi.org/10.1007/s12036-021-09746-y)

## Disclaimer
Please make sure to familiarise with the inputs and working of the algotithm. The accuracy of the results depends solely on the user. <br />
The contents of this repository are part of an ongoing research project. <br />
The form and function are generally limited to what is required for that project and may not be appropriate for your project. <br />
Feedback is welcomed, but feature requests will probably not be honored. <br />
The contents of this repository might change dramatically without notice. <br />

File structure:
``` bash
|   Binary_SED_fitting_v2.5.ipynb                # Fitting binary SEDs and plotting
│   Creating_VOSA_and_model_files_v1.1.ipynb     # Creating files compatible with VOSA 
│   LICENSE
│   log_file_binary.csv                          # log file
│   README.md
│
├───data
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
│       readme_models.txt                        
│       isochrone.txt                            # Parsec isochrones of logAge = 8,9 and 10
│       Table_Mass_0.2                           # WD cooling curve of mass 0.2
│       Table_Mass_1.3                           # WD cooling curve of mass 1.3
│       Koe_logg7.0.csv                          # Koester models of various logg
│       Koe_logg8.0.csv
│       Koe_logg9.0.csv
│       Kr_logg3.0_Zm05.csv                      # Kurucz models of various logg and Z
│       Kr_logg3.0_Zp00.csv
│       Kurucz_models_logg0_to_3.zip             # Compressed file including Kurucz models
│       Kurucz_models_logg3.5_to_5.zip           # Compressed file including Kurucz models
│       ...
└───outputs
    ├───chi_files
    │       WOCS2002_chi2_binary_Koe_logg7.0_1.csv                 # Initial chi2 file for WOCS2002 star
    │       WOCS2002_chi2_iterations100_binary_Koe_logg7.0_1.csv   # chi2 file after iterations
    │
    ├───binary_SEDs
        └───WOCS2002
                A_SED_WOCS2002_5250.0_logg3.5.png                 # SED fitted with only A component
                WOCS2002_Koe_logg7.0_14750_1.png                  # Initial binary SED
                WOCS2002_Koe_logg7.0_14750_1_iterations.png       # binary SED after iterations        
```
