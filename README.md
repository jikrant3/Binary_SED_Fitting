# Binary_SED_Fitting
Fitting SEDs of binary stars. The code is optimised for chi2 SED fitting of a double component star. It is assumed that one component has already been fit (using VOSA for our case).


File structure:

Working_dir
    |- Binary_SED_fitting.ipynb                     # Fitting binary SEDs and plotting  
    |- Creating_VOSA_and_model_files.ipynb          # Creating files compatible with VOSA  
    |                                                 (http://svo2.cab.inta-csic.es/theory/vosa/index.php)   
    |                                                 and creating model files using VOSA synthetic photometry  
    |                                                 (http://svo2.cab.inta-csic.es/theory/newov2/syph.php)  
    |
    |- data  
    |    |-photometry_file.csv                      # photometric information in .csv format  
    |    |                                                Must include "name, ra, dec, magnitude, magnitude_errors"  
    |    |-vosa_input.txt                           # VOSA upload file  
    |    |-vosa_results_38604                       # Keeping single fits from VOSA  
    |         |-objects  
    |               |-WOCS2002                      # Star name  
    |                    |-bestfit  
    |                         |-WOCS2002.bfit.phot.dat  # File with single fit parameters and observed flux  
    |  
    |- outputs  
    |    |- chi_files                               # Saving chi2 files  
    |    |- <mode1>                                 # Individual folders for SED plots of each mode  
    |    |    |- SED plots  
    |    |  
    |    |- <mode2>                                 # for example, finer_Kr_Kr, rough_Kr_Koe ...  
    |         |- SED plots  
    |  
    |- models                                       # Keeping synthetic photometry model files  
        |- Koe_logg7.0.csv                          # Koester model with logg = 7  
        |- Kr_logg5.0_Zm05.csv                      # Kurucz model with logg=5 and metallicity = -0.5  
        |  
        |- Raw_synthetic_files                      # If you want to create other models using Creating_VOSA_and_model_files.ipynb  
             |-koester2_da05000_700.dk.phot.dat     # example for Koester model  
             .  
             .  
