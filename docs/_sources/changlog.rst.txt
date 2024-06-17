.. currentmodule:: src.binary_sed_fitting

Changlog
========

3.4.1
-----

Fixed calculations of fractional errors in multi-component Systems.

3.4.0
-----

Added support for single SED fitting.

Added support for a hybrid Kurucz_UVBLUE model. 
UVBLUE (2005ApJ...626..411R) spectra for 850-4700 Angstroms and Kurucz spectra outside the wavelength 

Refactored the code to create seperate classes to deal with single stars and systems. 
:class:`Star`, :class:`System` :class:`Model` and :class:`Fitter`.

All code moved into ``src.binary_sed_fitting.py`` which works as a module.

Added logging. The raw log is saved in ``data\log_<data>.txt``. 

3.3.0
-----

The model files were updated and reformatted into ``xarray`` DataArrays.

The binary fitting code modified to work with xarrays.

All code moved into ``util_functions.py`` which works as a module.

Supported models: Kurucz (2003IAUS..210P.A20C), Koester (2010MmSAI..81..921K) 
