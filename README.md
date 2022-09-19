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
│   1_Creating_VOSA_files_v2.1.0.ipynb                       # Creating files compatible with VOSA 
│   2_example_single_star.ipynb                              # Example of single SED fitting and plotting
│   3_example_binary_star.ipynb                              # Example of double SED fitting and plotting
│   4_example_star_cluster.ipynb                             # Example of double SED fitting of multiple stars at once
│   log_binary_fitting.csv                                   # log file for binary fits
│   log_single_fitting.csv                                   # log file for single fits
│   README.md
│   util_functions.py
│
├───data
│   │   example_photomety_file.csv
│   │   example_VOSA_input_file.txt
│   │
│   ├───extinction_corrected_flux_files                      # extinction corrected flux for stars
│   │       WOCS1007.csv
│   │       WOCS2002.csv
│   │       Y1168.csv
│   │
│   └───vosa_results_53985                                   # Results of single fitting by VOSA
│       ├───objects
│       │   ├───WOCS1007
│       │   │   ├───bestfitp
│       │   │   │       WOCS1007.bfit.phot.dat
│       │   │
│       │   ├───WOCS2002
│       │   │   ├───bestfitp
│       │   │   │       WOCS2002.bfit.phot.dat
│       │   │
│       │   └───Y1168
│       │       ├───bestfitp
│       │       │       Y1168.bfit.phot.dat
│
├───outputs
│   ├───chi_files                                            # chi2 files of SEDs  
│   │       binary_WOCS1007_Koester_noisy_13.csv
│   │       ...
│   │
│   └───pickels                                              # pickels of SEDs fitted class objects   
│           BinaryStar_WOCS1007_13.pkl
│           ...
│
└───plots
    │   single_SEDsWOCS2002_Kurucz_1_noisy.jpg
    │
    ├───binary_SEDs                                          # binary SED fits
    │       WOCS1007_Koester_11750_13.jpg
    │       ...
    │
    ├───single_SEDs                                          # single SED fits
    │       blackbody_SED_WOCS2002_6181.457794845878_1.jpg
    │       ...
    │
    └───test                                                 # test SED fits
```
