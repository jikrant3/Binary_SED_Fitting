import copy
import linecache
import os
import pickle
from itertools import product
from time import gmtime, process_time, strftime

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import xarray as xr
from scipy.optimize import curve_fit
from scipy.stats import zscore

DIR_MODELS  = 'C:/Users/user/Documents/GitHub/models_and_tools/models/'
DIR_OUTPUTS = 'outputs/'


def get_blackbody_spectrum_loglog(log_wavelength, Te, log_scaling_factor):
    '''
    Adapted from the IUE RDAF 1989
    
    Parameters
    ----------
    log_wavelength : float or list
        log(Wavelength [A]) to calculater blackbody spectrum.
    Te : float
        Temperature [K] of the source.
    scaling_factor : float
        scaling factor to namalize flux (for a star: `(radius/distance^2)`).
    Returns
    -------
    bbflux : float or list
        Flux [ergs/cm2/s/a] for given temperature and wavelength range.
    '''
    c1 = 3.7417749e-5                    # =2*!dpi*h*c*c
    c2 = 1.4387687                       # =h*c/k

    wave = (10**log_wavelength)/1.e8     # angstroms to cm

    bbflux = c1 / (wave**5 * (np.exp(c2/wave/Te)-1.))
    bbflux = bbflux*1.e-8                # ergs/cm2/s/a

    log_bbflux = np.log10(bbflux)+log_scaling_factor
    
    return log_bbflux


def calc_radius(sf, distance):
    '''
    Radius for given scaling factor and distance
    
    Parameters
    ----------
    sf : float or list
        Scaling factor
    distance : float
        distance [pc]
    Returns
    -------
    radius : float or list
        radius [Rsun]
    '''
    return sf**0.5 * (distance*44353566.0)


def calc_sf(radius, distance):
    '''
    Calculate scaling factor for given radius and distance

    Parameters
    ----------
    radius : float or list
        Radius [Rsun]
    distance : float
        distance [pc]
    Returns
    -------
    sf : float or list
        Scaling factor
    '''
    return (radius/(distance*44353566.0))**2   


def calc_luminosity(radius,Te):
    '''
    Calculate luminosity from radius and temperature

    Parameters
    ----------
    radius : float or list
        radius [Rsun]
    Te : float
        Temperature [K]
    Returns
    -------
    luminosity : float or list
        luminosity [Lsun]
    '''
    sigma = 5.67e-8  #W m−2 K−4
    return (sigma * 4 * 3.141592 * (radius*6.957e+8)**2 * Te**4)/3.828e+26


def addition_in_quadrature(a,b):
    '''
    Square root of sum of squares
    '''
    return (a**2+b**2)**0.5


def find_nearest_with_index(array, value):
    '''
    Finds the nearest value available in an array
    '''
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]


def print_progress(current_iteration, total_iterations=100, step=10, message='Progress'): 
    '''
    Print progress
    '''
    if (current_iteration%step==0):
        print ('\r%s:  %d/%d  (%d%%)' %(message, current_iteration+1, total_iterations, 100*(current_iteration+1)/total_iterations), end='')


def get_number_from_string(line): 
    '''
    Get a float number from string line.
    '''
    for t in line.split():
        try: l = float(t)
        except ValueError: pass
    return l


def get_parameter_from_VOSA_file(param, file_name, verbose=False):
    '''
    Get value of parameter from VOSA file
    
    Parameters
    ----------
    param :str
        Name of the parameter
    file_name : str
        Name of the file to seach the ``param``
    Returns
    -------
    radius : float or list
        radius [Rsun]
    '''
    with open(file_name,'r') as f:
        for (i, line) in enumerate(f):
            if param in line:
                value = get_number_from_string(line)
                if verbose: print(line, i)
                return value
    raise Exception('    ERROR: %s not found in %s'%(param, file_name))

    
def get_model_file_name(model_to_fit):
    '''
    Get the model file name and path.
    
    Returns
    -------
    model_name : str
        Path of the file
    '''
    if model_to_fit == 'Kurucz':
        return DIR_MODELS + 'kurucz_synthetic_photometry.nc'
    if model_to_fit == 'Koester':
        return DIR_MODELS + 'koester_synthetic_photometry.nc'

    
def get_realistic_errors_from_iterations(parameter_values,grid_steps,para_type='parameter'):
    '''
    Estimates errors using spread in the noisy fits and edge cases.
    
    - If the best fit parameter is near boundry         --> Exaggerate errors and issue warning
    - If spread in iterations is less than step size    --> Keep errors similar to the step size
    - Otherwise:
    
        - Calculates the 32, 50 and 68th percentiles of a given array. 
        - In case the percentile lies between two grid values, a lower, nearest and higher value is chosen as the 23th, 50th and 68th percentile, respectively. 

    Parameters
    ----------
    parameter_values : list
        List of parameter values
    grid_steps : list
        List of grid steps available
    para_type : str
        Name of the parameter
    Returns
    -------
    para_50 : float
        Median of the ``parameter_values``
    error_lower : float
        Lower bound error in the parameter
    error_upper : float
        Upper bound error in the parameter
        
    '''
    if len(grid_steps)==1:
        para_50, error_lower, error_upper = np.percentile(parameter_values,50,interpolation='nearest'), np.nan, np.nan
        return para_50, error_lower, error_upper
    
    para_32, para_50, para_68 = np.percentile(parameter_values,31.7,interpolation='lower'), np.percentile(parameter_values,50,interpolation='nearest'), np.percentile(parameter_values,68.3,interpolation='higher')
    error_lower               = para_50 - para_32
    error_upper               = para_68 - para_50

    median_index, _           = find_nearest_with_index(grid_steps, para_50)
        
    if median_index==0:
        print('    WARNING: Best fit %s is at lower limit of the model.' %para_type)
        if error_upper  > 0:
        # Lower error is kept as 3x of upper limit errors
            error_lower = 3 * error_upper               
        else:
            error_upper = grid_steps[median_index+1]-grid_steps[median_index]
            error_lower = 3 * error_upper               

    if median_index==len(grid_steps)-1:
        print('    WARNING: Best fit %s is at upper limit of the model.' %para_type)
        if error_lower  > 0:
        # upper erro is kept as 3x of lower limit errors
            error_upper = 3 * error_lower              
        else:
            error_lower = grid_steps[median_index]-grid_steps[median_index-1]
            error_upper = 3 * error_lower              

    if (median_index>0) & (median_index<len(grid_steps)-1):
        if error_lower == 0:
            error_lower = grid_steps[median_index]-grid_steps[median_index-1]
        if error_upper == 0:
            error_upper = grid_steps[median_index+1]-grid_steps[median_index] 
            
    return para_50, error_lower, error_upper


def load_SingleStar(name, cycle, verbose=False):
    '''
    Loads pickeled SingleStar object
    
    Returns
    -------
    self: SingleStar
        SingleStar object initialised using previous fitting 

    '''
    file_name = DIR_OUTPUTS + 'pickels/SingleStar_%s_%d.pkl'%(name,cycle)
    with open(file_name, 'rb') as f:
        self = pickle.load(f)
        
    self.read_model_file(verbose=verbose)
    self.create_star_dataarrays(verbose=verbose)
    self.calculate_chi2(cycle, refit=False,verbose=verbose)
    self.calculate_chi2_noisy(cycle, refit=False,verbose=verbose)

    return self


def load_BinaryStar(name, cycle,verbose=False):
    '''
    Loads pickeled BinaryStar object
    
    Returns
    -------
    self: BinaryStar
        BinaryStar object initialised using previous fitting 
    '''
    file_name = DIR_OUTPUTS+'pickels/BinaryStar_%s_%d.pkl'%(name,cycle)
    with open(file_name, 'rb') as f:
        self = pickle.load(f)

    self.read_model_file(verbose=verbose)
    self.create_star_dataarrays(verbose=verbose)
    self.calculate_chi2(cycle, refit=False,verbose=verbose)
    self.calculate_chi2_noisy(cycle, refit=False,verbose=verbose)

    return self


class SingleStar:
    '''
    Class for SED fitting of a single star
    '''
    def __init__(self, name, model_to_fit, dir_obs, distance, distance_err, free_para, filters_to_drop=[], wavelength_range= [0,100_000_000],verbose=False):
        '''
        Initialises the star using name, distance, model_to_fit etc.
        
        Parameters
        ----------
        name : str
            Name of the source
        model_to_fit : str
            Model name of the component
        dir_obs : str
            Path to the observed SED '.csv' file  
        distance : float
            Distance [pc]
        distance_err : float
            Error in distance [pc]    
        free_para : int
            Number of free parameters
        filters_to_drop : list
            List of filters to be dropped before fitting SED. Default value is ``[]``.
        wavelength_range : list
            minimum and maximum value of filters to consider fitting
        '''
        self.name               = name
        self.model_to_fit       = model_to_fit
        self.dir_obs            = dir_obs
        self.distance           = distance
        self.distance_err       = distance_err
        self.free_para          = free_para
        self.filters_to_drop    = filters_to_drop
        self.wavelength_range   = wavelength_range
        self.type               = 'single'
        
        if verbose: print('================\n%s\n----------------' %self.name)

    def read_observed_SED(self, verbose=False):
        '''
        Reading single component observed SED.
        
        - Input file should be extinction corrected .csv with columns = [filter, wavelength, flux, error]
        
        Returns
        -------
        flux : DataFrame
            DataFrame with filter names, wavelengths, observed flux and errors (including log format)
        flux_all : DataFrame
            Copy of ``flux`` for backup. Used later for plotting.
        N_points : int
            Number of data points
        N_Np : int
            Degrees of freedom (``N_points - free_para``)
        '''
        file = self.dir_obs + self.name + '.csv'
        flux = pd.read_csv(file, engine='python', header=0)
        flux = flux.set_index('FilterID')

        if len(flux.index) > len(flux.index.unique()):
            raise Exception('    ERROR! Some filters are repeated. Remove repeated filters.')
            
        # Replacing zeros in errors with 110% of max error
        flux['error_fraction']             = flux['error']/flux['flux']
        flux.loc[flux.error == 0, 'error'] = flux['flux']*(flux['error_fraction'].max()*1.1)
        # Recalculating errors_fraction 
        flux['error_fraction']             = flux['error']/flux['flux']

        # error modification for calculating vgf (minimum error = 2%) and vgfb (minimum error = 10%)
        flux['error_2percent']  = np.where(flux['error_fraction']<0.02, 0.02, flux['error_fraction'])*flux['flux']
        flux['error_10percent'] = np.where(flux['error_fraction']<0.10, 0.10, flux['error_fraction'])*flux['flux']
        
        flux['log_wavelength']  = np.log10(flux['wavelength'])
        flux['log_flux']        = np.log10(flux['flux'])
        flux['log_error']       = 0.434 * flux['error']/flux['flux']
 
        self.flux     = flux
        self.flux_all = flux
        self.N_points = len(self.flux)
        self.N_Np     = self.N_points-self.free_para
        
        # Remove some filters from the above file if necessary
        self.drop_filters(verbose=verbose)

        if verbose:
            print('\n    RUNNING: read_observed_SED')
            print('    Total filters: %d' %len(flux))

    def drop_filters(self, verbose=False):
        '''
        Dropping filters before further SED fitting.
        
        - Some filters will have to be removed while fitting.
        - You can remove individual filters based on science case. For example:
        
            - In case of IR excess, you can remove Wise filters
            - For bad chi2, you can remove specific filters
            
        Returns
        -------
        flux : DataFrame
            DataFrame containing only required filters
        not_fitted : DataFrame
            DataFrame containing only removed filters           
        N_points : int
            Updated number of data points
        N_Np : int
            Updated degrees of freedom (``N_points - free_para``)
        '''    
        self.not_fitted     = self.flux[(self.flux['wavelength']<0)]

        for filter_name in self.filters_to_drop:
            self.not_fitted = pd.concat([self.not_fitted, self.flux[(self.flux.index==filter_name)]])
            self.flux       = self.flux.drop(index=filter_name)
            
        _excess_filters = self.flux[(self.flux.wavelength<self.wavelength_range[0]) | (self.flux.wavelength>self.wavelength_range[1])]        
        for filter_name in _excess_filters.index:
            self.not_fitted = pd.concat([self.not_fitted, self.flux[(self.flux.index==filter_name)]])
            self.flux       = self.flux.drop(index=filter_name)
         
        self.N_points       = len(self.flux)
        self.N_Np           = self.N_points-self.free_para
        
        
        # printing filters to be fitted and not_fitted
        if verbose:
            _t1, _t2            = pd.DataFrame(),pd.DataFrame()    
            _t1['to_be_fitted'] = self.flux.index.values
            _t1['wavelength']   = self.flux.wavelength.values
            _t1                 = _t1.set_index('wavelength')
            _t2['not_fitted']   = self.not_fitted.index.values
            _t2['wavelength']   = self.not_fitted.wavelength.values
            _t2                 = _t2.set_index('wavelength')
            _filter_table       = pd.concat([_t1,_t2],sort=True)    
            print('\n    RUNNING: drop_filters')
            print(_filter_table.sort_index().fillna(''))
            print('    Filters to fit: %d' %len(self.flux))

    def fit_blackbody(self, cycle, p0=[5000.,-20], plot=False, save_plot=False, verbose=False,show_plot=True, folder_path='plots/single_SEDs/'):
        '''
        Fitting and plotting a blackbody SED to the observed SED.
        
        Parameters
        ----------
        cycle : int
            Fitting cycle
        p0 : list
            List of initial guesses for Temperature [K] and log(scaling factor)
        plot : boolean
            If True, will plot the blackbody fitted SED
        save_plot : boolean
            if True, will save the blackbody fitted SED
            
        Return 
        ------
        Te_blackbody : float
            Temperature from blackbody fit
        log_sf_blackbody : float
            log(scaling factor) from the blackbody fit
        '''
        popt, pcov            = curve_fit(get_blackbody_spectrum_loglog, self.flux['log_wavelength'], self.flux['log_flux'], p0=p0)
        
        self.Te_blackbody     = popt[0]
        self.log_sf_blackbody = popt[1]
        
        if verbose:
            print('\n    RUNNING: fit_blackbody ')
            print('    Fit parameters: T=%d K, log_sf=%.2f' % tuple(popt))
        
        if plot:
            f, axes = plt.subplots(figsize=(8,6),nrows = 3, ncols = 1)
            [axi.set_axis_off() for axi in axes.ravel()]
            axes[0] = f.add_axes([0.1, 0.40, 0.85, 0.55])
            axes[1] = f.add_axes([0.1, 0.25, 0.85, 0.15])
            axes[2] = f.add_axes([0.1, 0.10, 0.85, 0.15])

            ####################### SED
            # observed data
            axes[0].scatter(self.not_fitted['log_wavelength'], self.not_fitted['log_flux'], color='k', marker='o',label ='No Fit', s=30, facecolors='none',zorder=2)
            axes[0].errorbar(self.flux['log_wavelength'], self.flux['log_flux'], yerr=self.flux['log_error'],color='k', label='Obs',fmt='none',lw=2, zorder=3, capsize=4)
            # BB spectrum
            xvalues = np.linspace(self.flux_all['log_wavelength'].min(),self.flux_all['log_wavelength'].max(), 100)
            axes[0].plot(xvalues, get_blackbody_spectrum_loglog(xvalues, *popt), color='pink',ls='--', label='Blackbody spectrum', zorder=0)
            # BB SED
            axes[0].plot(self.flux_all['log_wavelength'], get_blackbody_spectrum_loglog(self.flux_all['log_wavelength'], *popt), 'r--', label='Blackbody SED (T=%d, log_sf=%.2f)' % tuple(popt), zorder=1)
            ####################### Residual
            fractional_residual = (self.flux_all['flux'] - 10**get_blackbody_spectrum_loglog(self.flux_all['log_wavelength'], *popt))/self.flux_all['flux']
            axes[1].plot(self.flux_all['log_wavelength'], fractional_residual, 'r--', label='')
            ####################### chi2
            chi2_i = (self.flux['flux']-10**get_blackbody_spectrum_loglog(self.flux['log_wavelength'], *popt))**2 / self.flux['error']**2
            axes[2].plot(self.flux['log_wavelength'], chi2_i, 'r--', label='$\chi^2$ = %.2f\n$\chi^2_r$ = %.2f'%(chi2_i.sum(),chi2_i.sum()/self.N_Np))
            ####################### Titles and labels
            axes[0].set_title(self.name + ' (blackbody fit)', x=0, y=1, ha='left')
            axes[0].set_ylabel('log(Flux) (erg s$^{-1}$ cm$^{-2}$ $\AA$$^{-1}$)')
            axes[1].set_ylabel('Fractional\nResidual')
            axes[1].set_xlabel('log(Wavelength) ($\AA$)')
            ####################### decoration    
            plt.setp(axes[0].get_xticklabels(),visible=False)
            # taking xlim and ylims from one axis to another
            axes[2].set_xlim(axes[0].get_xlim())
            axes[0].grid()
            axes[1].grid()
            axes[2].grid()
            axes[0].tick_params(which='both', direction='out', length=4)
            axes[1].tick_params(which='both', direction='out', length=4)
            axes[0].legend()
            axes[2].legend()
            ####################### Saving file    
            if save_plot:
                if not os.path.exists(folder_path): 
                    os.makedirs(folder_path)
                plt.savefig (folder_path + 'blackbody_SED_' +self.name+'_'+str(self.Te_blackbody)+'_'+str(cycle)+'.jpg', format='jpg', dpi=300)
            if not show_plot:
                plt.close()
                
    def read_model_file(self,verbose=False):
        '''
        Reads model file from the 'DIR_MODELS' folder as a DataArray.
        
        Returns
        -------
        da_model : DataArray
            nD array of all available SED models. The dimentions are FilterID, Te, logg (and more if applicable).
        '''
        self.model_file_name = get_model_file_name(self.model_to_fit)

        with xr.open_dataarray(self.model_file_name) as da_model:
            # editing the "Wavelengths" attribute to present filters
            df1                           = pd.DataFrame()
            df1['FilterID']               = da_model['FilterID']
            df1                           = df1.set_index('FilterID')
            df1['Wavelength']             = da_model.attrs['Wavelengths']
            df2                           = pd.DataFrame()
            df2['FilterID']               = self.flux.index
            df2                           = df2.set_index('FilterID')
            df2['Wavelength']             = df1['Wavelength']
            da_model.attrs['Wavelengths'] = df2['Wavelength'].values

            da_model = da_model.sel(FilterID=self.flux.index)

        if verbose:
            print('\n    RUNNING: read_model_file ')
            print(da_model.coords)
            print('\n    Provide contrains on ', da_model.dims[1:])

        self.da_model = da_model

    def constrain_fitting_parameters(self, limits, verbose=False):
        '''
        Cropping the model DataArray according to givel limits
        
        Returns
        -------
        da_model : DataArray
            Cropped DataArray
        '''
        if self.model_to_fit[:6] == 'Kurucz':
            self.da_model = self.da_model.sel(Te   =slice(limits['Te'][0],   limits['Te'][1]))
            self.da_model = self.da_model.sel(logg =slice(limits['logg'][0], limits['logg'][1]))
            self.da_model = self.da_model.sel(MH   =slice(limits['MH'][0],   limits['MH'][1]))
            self.da_model = self.da_model.sel(alpha=slice(limits['alpha'][0],limits['alpha'][1]))
        if self.model_to_fit == 'Koester':
            self.da_model = self.da_model.sel(Te   =slice(limits['Te'][0],  limits['Te'][1]))
            self.da_model = self.da_model.sel(logg =slice(limits['logg'][0],limits['logg'][1]))
        if verbose: print('    da model:\n',self.da_model.coords)
        
        self.create_star_dataarrays(verbose=verbose)
    
    def create_star_dataarrays(self, verbose=False):
        '''
        Create observed flux and flux error DataArrays with dimentions of ``da_model`` for vectorised arithmatic later.
        
        Returns
        -------
        da_obs : DataArray
            nD array in the shape of ``da_model``. All filters have same value: observed SED
        da_obs_error : DataArray
            nD array in the shape of ``da_model``. All filters have same value: observed error
        da_obs_error_2percent : DataArray
            nD array in the shape of ``da_model``. All filters have same value: observed error (min 2%)
        da_obs_error_10percent : DataArray
            nD array in the shape of ``da_model``. All filters have same value: observed error (min 10%)
        '''    

        da_obs                 = self.da_model.copy()
        da_obs_error           = self.da_model.copy()
        da_obs_error_2percent  = self.da_model.copy()
        da_obs_error_10percent = self.da_model.copy()

        for filter_name in self.flux.index:
            da_obs                  = da_obs.where(da_obs.FilterID!=filter_name,self.flux['flux'][filter_name])
            da_obs_error            = da_obs_error.where(da_obs_error.FilterID!=filter_name,self.flux['error'][filter_name])
            da_obs_error_2percent   = da_obs_error_2percent.where(da_obs_error_2percent.FilterID!=filter_name,self.flux['error_2percent'][filter_name])
            da_obs_error_10percent  = da_obs_error_10percent.where(da_obs_error_10percent.FilterID!=filter_name,self.flux['error_10percent'][filter_name])

        self.da_obs                 = da_obs
        self.da_obs_error           = da_obs_error
        self.da_obs_error_2percent  = da_obs_error_2percent
        self.da_obs_error_10percent = da_obs_error_10percent

        self.da_obs_error.attrs['long_name']           = 'Error'
        self.da_obs_error_2percent.attrs['long_name']  = 'Error (2\% min)'
        self.da_obs_error_10percent.attrs['long_name'] = 'Error (10\% min)'
        if verbose: 
            print('    da observed flux:\n',da_obs.head(2))
            print('    da observed error:\n',da_obs_error.head(2))
            print('    da observed error (2\% min):\n',da_obs_error_2percent.head(2))
            print('    da observed error (10\% min):\n',da_obs_error_10percent.head(2))

    def create_sf_list(self, log_sf_flexibility=3., log_sf_stepsize=0.01):
        '''
        Create list of scaling factors
        
        Parameters
        ----------
        log_sf_flexibility : float
            Allowed flexibility in the ``log(sf)`` around mean of ``log(sf_blackbody)``
        log_sf_stepsize : float
            Stepsize for ``sf`` in log space
        
        Returns
        -------
        sf_list : list
            List of scaling factor (NOT in log space) within ``log_sf_blackbody +- log_sf_flexibility`` with stepsize of ``log_sf_stepsize``.
        '''
        self.sf_list = 10**np.arange(self.log_sf_blackbody-log_sf_flexibility,self.log_sf_blackbody+log_sf_flexibility, log_sf_stepsize) 

    def calculate_chi2(self, cycle, refit=True, verbose=False, trim=True, minimising_param='chi2'):
        '''
        Calculate chi2
        
        - Takes ``da_model`` and scales it for each sf in ``sf_list``.
        - Subtracts the ``da_model`` from ``da_obs`` and calculates corresponding chi2_i for each datapoint.
        - Sums the chi2_i for each model and saves the chi2 values into DataArray (``da_chi2_stacked``). 
        - Simultaneously calculates vgfb2
        - Converts ``da_chi2_stacked`` into a dataframe (``df_chi2``)
        - Calculates R [Rsun] and L [Lsun] for each fit
        - Saves ``df_chi2`` into a .csv file
        
        minimising_param : str
            The parameter over which to determine the best fit. Select within 'chi2', 'vgf2', 'vgfb2'
        Returns
        -------
        da_chi2 : DataFrame
            chi2 sorted DataFrame with Te [K], logg, (MH), sf, R [Rsun], L [Lsun], vgfb2.
        '''
        chi2_file_name = DIR_OUTPUTS + 'chi_files/'+self.type +'_'+self.name+'_'+self.model_to_fit+'_'+str(cycle)+'.csv'

        # If the chi2 file exists, it will be read (depends on refit = True or False)
        if os.path.isfile(chi2_file_name): 
            if not refit:
                df_chi2      = pd.read_csv(chi2_file_name)
                if verbose:     
                    print('    Reading %s'%chi2_file_name)
                    print('    WARNING! give "refit=True" if you want to rerun the fitting process.')
                    print(df_chi2.head())
                self.df_chi2 = df_chi2
                return 
            if refit: print('    WARNING: '+chi2_file_name+' file will be overwritten.')

        da_model = self.da_model
        if self.model_to_fit[:6]=='Kurucz':
            nd_arr           = np.empty((len(da_model.Te), len(da_model.logg), len(da_model.MH), len(da_model.alpha), len(self.sf_list)))
            nd_arr.fill(np.nan)
            da_chi2_stacked  = xr.DataArray(nd_arr, coords={'Te'   :da_model.Te,
                                                            'logg' :da_model.logg,
                                                            'MH'   :da_model.MH,
                                                            'alpha':da_model.alpha,
                                                            'sf'   :self.sf_list})
            da_vgf2_stacked  = da_chi2_stacked.copy()
            da_vgfb2_stacked = da_chi2_stacked.copy()
        if self.model_to_fit=='Koester':
            nd_arr           = np.empty((len(da_model.Te), len(da_model.logg), len(self.sf_list)))
            nd_arr.fill(np.nan)
            da_chi2_stacked  = xr.DataArray(nd_arr, coords={'Te'  :da_model.Te,
                                                            'logg':da_model.logg,
                                                            'sf'  :self.sf_list})
            da_vgf2_stacked  = da_chi2_stacked.copy()
            da_vgfb2_stacked = da_chi2_stacked.copy()

        for idx in range (len(self.sf_list)):
            print_progress(idx, len(self.sf_list), 10,  message='   Calculating chi2')

            sf          = self.sf_list[idx]
            da_residual = self.da_obs - (da_model * sf)

            da_chi2_i   = (da_residual/self.da_obs_error)**2
            da_chi2     = da_chi2_i.sum('FilterID',min_count=1)
            da_chi2_stacked.loc[dict(sf=sf)] = da_chi2

            da_vgf2_i  = (da_residual/self.da_obs_error_2percent)**2
            da_vgf2    = da_vgf2_i.sum('FilterID',min_count=1)
            da_vgf2_stacked.loc[dict(sf=sf)] = da_vgf2

            da_vgfb2_i  = (da_residual/self.da_obs_error_10percent)**2
            da_vgfb2    = da_vgfb2_i.sum('FilterID',min_count=1)
            da_vgfb2_stacked.loc[dict(sf=sf)] = da_vgfb2

        df_chi2         = da_chi2_stacked.to_dataframe(name='chi2').reset_index()
        df_vgf2         = da_vgf2_stacked.to_dataframe(name='vgf2').reset_index()
        df_vgfb2        = da_vgfb2_stacked.to_dataframe(name='vgfb2').reset_index()

        df_chi2['R']    = calc_radius(df_chi2['sf'], self.distance)
        df_chi2['L']    = calc_luminosity(df_chi2['R'], df_chi2['Te'])
        df_chi2['vgf2'] = df_vgf2['vgf2'].values
        df_chi2['vgfb2']= df_vgfb2['vgfb2'].values

        if minimising_param not in ['chi2', 'vgf2', 'vgfb2']:
            raise Exception('    Error: The minimising_param should be one of the following: chi2, vgf2, vgfb2')
        df_chi2 = df_chi2.sort_values(minimising_param)
        
        df_chi2.reset_index(drop=True, inplace=True)

        if verbose: 
            print('\n    Calculated chi2 for %d models. \n    Best 5 models:'%len(df_chi2))
            print(df_chi2.head())
            
        if trim: df_chi2 = df_chi2.head(5000)
        self.df_chi2 = df_chi2

        if not os.path.exists(DIR_OUTPUTS + 'chi_files/'): os.makedirs(DIR_OUTPUTS + 'chi_files/')
        if verbose: print('    Saving %s'%chi2_file_name)
        df_chi2.to_csv(chi2_file_name, index=False, header=True, sep=',')
        
    def get_parameters_from_chi2_minimization(self, verbose=False):
        '''
        Estimates best fit parameters from least chi2 fit
        
        - model_flux, Te, sf, R, L, chi2, vgf, vgfb
        - Updates these parameter to parent class object
        - Raises warning if any datapoints have significantly (3-sigma) higher chi2_i than average
        '''
        self.Te                   = self.df_chi2.Te[0]
        self.logg                 = self.df_chi2.logg[0]

        if self.model_to_fit[:6] =='Kurucz':
            self.MH               = self.df_chi2.MH[0]
            self.alpha            = self.df_chi2.alpha[0]
            best_fit_flux         = self.da_model.sel(Te=self.Te).sel(logg=self.logg).sel(MH=self.MH).sel(alpha=self.alpha)

        if self.model_to_fit     =='Koester':
            self.MH               = np.nan
            self.alpha            = np.nan
            best_fit_flux         = self.da_model.sel(Te=self.Te).sel(logg=self.logg)

        self.sf                   = self.df_chi2.sf[0]
        self.R                    = self.df_chi2.R[0]
        self.L                    = self.df_chi2.L[0]

        self.flux['model_flux'] = self.sf * best_fit_flux
        _residual_flux          = self.flux['flux'] - self.flux['model_flux']
        self.flux['chi2_i']     = _residual_flux**2 / self.flux['error']**2
        self.flux['vgf2_i']     = _residual_flux**2 / self.flux['error_2percent']**2
        self.flux['vgfb2_i']    = _residual_flux**2 / self.flux['error_10percent']**2
        self.chi2               = self.flux['chi2_i'].sum()
        self.chi2_r             = self.chi2/self.N_Np
        self.vgf2               = self.flux['vgf2_i'].sum()/self.N_Np
        self.vgfb2              = self.flux['vgfb2_i'].sum()/self.N_Np

        with xr.open_dataarray(self.model_file_name) as da_model_all:
            da_model_all = da_model_all.sel(FilterID=self.flux_all.index)
            if self.model_to_fit[:6]=='Kurucz':
                da_model_all = da_model_all.sel(Te=self.Te).sel(logg=self.logg).sel(MH=self.MH).sel(alpha=self.alpha)
            if self.model_to_fit=='Koester':
                da_model_all = da_model_all.sel(Te=self.Te).sel(logg=self.logg)
            df_model_all = da_model_all.to_dataframe(name='model_flux')
        self.flux_all['model_flux'] = df_model_all['model_flux']*self.sf    

        if self.model_to_fit[:6]=='Kurucz':
            if verbose: print('\n    Fitting parameters: T=%d,logg=%.2f,MH=%.2f,alpha=%.1f,sf=%.2e,R=%f,L=%f with chi2=%.2f'%(self.Te,self.logg, self.MH, self.alpha, self.sf, self.R, self.L,self.chi2))
        if self.model_to_fit       =='Koester':
            if verbose: print('\n    Fitting parameters: T=%d,logg=%.2f,sf=%.2e,R=%f,L=%f with chi2=%.2f'%(self.Te,self.logg, self.sf, self.R, self.L,self.chi2))

        # Printing the filtes with too large chi2 values i.e. 3sigma away from other chi2 values
        abs_zscore = np.abs(zscore(self.flux['chi2_i']))
        outliers = self.flux[abs_zscore>=3].index.values
        if len(outliers)>0: print('\n    Based on chi2, I recommend removal of following filters: ', outliers)
            
    def plot_skeletal_SED(self):
        '''
        Builds the skeletal plots (along with labels, limits etc) for binary SED.
        
        Returns
        -------
        f : Figure
        axes : axes.Axes or array of axes
        '''
        ###################### initialising
        f, axes = plt.subplots(figsize=(12,6),nrows = 3, ncols = 3)
        [axi.set_axis_off() for axi in axes.ravel()]

        axes[0,0] = f.add_axes([0.06, 0.44, 0.49, 0.50])
        axes[1,0] = f.add_axes([0.06, 0.27, 0.49, 0.17])
        axes[2,0] = f.add_axes([0.06, 0.10, 0.49, 0.17])

        axes[0,1] = f.add_axes([0.63, 0.66, 0.30, 0.28])
        axes[1,1] = f.add_axes([0.63, 0.38, 0.30, 0.28])
        axes[2,1] = f.add_axes([0.63, 0.10, 0.30, 0.28])

        axes[0,2] = f.add_axes([0.91, 0.10, 0.02, 0.28])
        ####################### SED
        axes[0,0].scatter(self.not_fitted['wavelength'], self.not_fitted['flux'], color='k', marker='o',label ='No Fit', s=30, facecolors='none',zorder=2)
        axes[0,0].errorbar(self.flux['wavelength'], self.flux['flux'], yerr=self.flux['error'],color='k', label='Obs',fmt='none',lw=2, capsize=4,zorder=3)
        ####################### residual
        axes[1,0].errorbar(self.flux['wavelength'], self.flux['flux']-self.flux['flux'], yerr=self.flux['error_fraction'],color='k', label='Obs',fmt='none',lw=2, capsize=4) 
        axes[1,0].errorbar(self.not_fitted['wavelength'], self.not_fitted['flux']-self.not_fitted['flux'], yerr=self.not_fitted['error_fraction'],color='0.5', label='',fmt='none',lw=2, capsize=4) 
        ####################### isochrones and WD cooling curves
        iso        = pd.read_csv(DIR_MODELS + 'master_isochrone.csv')
        iso_8      = iso[iso.logAge==8]
        iso_9      = iso[iso.logAge==9]
        iso_10     = iso[iso.logAge==10]
        Bergeron_WD= pd.read_csv(DIR_MODELS + 'master_Bergeron_WD.csv')
        WD_02      = Bergeron_WD[(Bergeron_WD.mass==0.2) & (Bergeron_WD.spectral_type=='DA')]
        WD_13      = Bergeron_WD[(Bergeron_WD.mass==1.3) & (Bergeron_WD.spectral_type=='DA')]

        axes[0,1].plot(10**(iso_8.logTe),10**(iso_8.logL), label='',c='0',lw=0.5, rasterized = True, zorder=1)
        axes[0,1].plot(10**(iso_9.logTe),10**(iso_9.logL), label='',c='0',lw=0.5, rasterized = True, zorder=1)
        axes[0,1].plot(10**(iso_10.logTe),10**(iso_10.logL), label='',c='0',lw=0.5, rasterized = True, zorder=1)
        axes[0,1].plot(WD_02.Teff,10**WD_02.logL, label='',c='0', ls=(0,(5,10)),lw=0.5, rasterized = True, zorder=1)
        axes[0,1].plot(WD_13.Teff,10**WD_13.logL, label='',c='0', ls=(0,(5,10)),lw=0.5, rasterized = True, zorder=1)
        ####################### Labels
        axes[0,0].set_ylabel('Flux (erg s$^{-1}$ cm$^{-2}$ $\AA$$^{-1}$)')
        axes[1,0].set_ylabel('Fractional\nResidual')
        axes[2,0].set_ylabel('$\chi^2_i$')
        axes[2,0].set_xlabel('Wavelength ($\AA$)')
        axes[0,1].set_ylabel('L ($L_{\odot}$)')
        axes[1,1].set_ylabel('R ($R_{\odot}$)')
        axes[2,1].set_ylabel('$\chi^2_r$')
        axes[2,1].set_xlabel('Te (K)')
        ####################### axes range and scales
        axes[0,0].set_xscale('log')
        axes[0,0].set_yscale('log')
        axes[1,0].set_xscale('log')
        axes[2,0].set_xscale('log')
        axes[0,1].set_yscale('log')
        axes[1,1].set_yscale('log')
        axes[0,1].set_xscale('log')
        axes[1,1].set_xscale('log')
        axes[2,1].set_xscale('log')
        axes[2,1].set_yscale('log')

        wave_min = self.flux_all['wavelength'].min()
        wave_max = self.flux_all['wavelength'].max()

        axes[0,0].set_xlim([wave_min/1.2,wave_max*1.2])
        axes[1,0].set_xlim([wave_min/1.2,wave_max*1.2])
        axes[2,0].set_xlim([wave_min/1.2,wave_max*1.2])

        axes[1,0].set_ylim([-0.2, 1.1])

        T_min = min(self.df_chi2.Te.min(), 4500)
        T_max = max(self.df_chi2.Te.max(), 40000)

        for ax in [axes[0,1],axes[1,1],axes[2,1]]:
            ax.set_xlim(T_max*1.1,T_min/1.4)

        flux_min = self.flux_all['flux'].min()
        flux_max = self.flux_all['flux'].max()
        axes[0,0].set_ylim([flux_min/5,flux_max*5])
        axes[0,1].set_ylim(self.L/1e3,self.L*1e3)

        plt.setp(axes[0,0].get_xticklabels(),visible=False)
        plt.setp(axes[1,0].get_xticklabels(),visible=False)
        plt.setp(axes[0,1].get_xticklabels(),visible=False)
        plt.setp(axes[1,1].get_xticklabels(),visible=False)

        ####################### decoration    
        for i,j in product(range(3),range(2)):
            axes[i,j].tick_params(which='both', direction='out', length=4)
            axes[i,j].grid(which='both', axis='x', zorder=0)
            axes[i,j].grid(which='major', axis='y', zorder=0)
        
        axes[0,2].tick_params(which='both', direction='out', length=4)

        return f, axes

    def plot_fitted_SED(self,cycle,save_plot=True, show_plot=True,excess_cutoff=None, folder_path='plots/single_SEDs/'):
        '''
        Plotting single SED
        
        - Creates SED plot after observed flux is fitted.
        - Includes HR diagram with isochrone of logAge 8/9/10, WD models for mass 0.2/1.3.
        - Includes radius vs Te plot
        - Includes chi2 vs Te plot. This plot should have explicit global minima in a good fit.
        '''
        f, axes = self.plot_skeletal_SED()
 
        axes[0,0].plot(self.flux_all['wavelength'], self.flux_all['model_flux'], color='green', linestyle='-',label ='Model', lw=1)
            ########## Fractional residual
        _fractional_residual = (self.flux_all['flux']-self.flux_all['model_flux'])/self.flux_all['flux']
        axes[1,0].plot(self.flux_all['wavelength'], _fractional_residual,label='',marker='',color='green',lw=1, linestyle='-')
        if excess_cutoff!=None:
            _excess = self.flux_all[_fractional_residual>excess_cutoff]
            _fractional_residual = (_excess['flux']-_excess['model_flux'])/_excess['flux']
            axes[1,0].scatter(_excess['wavelength'], _fractional_residual,label='Excess',color='r')
            axes[1,0].axhline(excess_cutoff, c='pink')
        ########## chi2_i
        axes[2,0].plot(self.flux['wavelength'], self.flux['chi2_i'],label='',marker='o',color='green', linestyle='-',lw=1)
        ########## L vs T (HR diagram)
        if len(self.df_chi2)>10_000: 
            top_chi2        = np.nanpercentile(self.df_chi2.chi2, 1)
        else: 
            top_chi2        = self.df_chi2.chi2.max()
        top_df_chi2     = self.df_chi2[self.df_chi2.chi2<=top_chi2]
        Te_data_top     = top_df_chi2.Te
        L_data_top      = top_df_chi2.L
        vgfb2_data_top  = top_df_chi2.vgfb2
        R_data_top      = top_df_chi2.R
        chi2_r_data_top = top_df_chi2.chi2/self.N_Np

        axes[0,1].scatter(Te_data_top,L_data_top, marker='.',c=vgfb2_data_top ,label='',s=5,rasterized = True, zorder=2, cmap='summer')
        axes[0,1].scatter(self.Te,self.L, marker='s', label='Best fit',c='r', s=40,rasterized = True, zorder=3)
        ########## R vs T    
        # axes[1,1].scatter(x_data, y_data, marker='.', c='0.8', s=0.2,rasterized = True, zorder=0)
        cs1 = axes[1,1].scatter(Te_data_top,R_data_top,c=vgfb2_data_top, cmap='summer',s=10,rasterized = True, zorder=2, norm=matplotlib.colors.LogNorm())
        axes[1,1].scatter(self.Te,self.R, marker='s', label='',c='r', s=40,rasterized = True, zorder=3)
        ########## chi2 vs T
        cs2 = axes[2,1].scatter(Te_data_top,chi2_r_data_top,c=vgfb2_data_top, cmap='summer',s=10,rasterized = True,zorder=2, norm=matplotlib.colors.LogNorm())  
        axes[2,1].scatter(self.Te,self.chi2_r, marker='s', label='',c='r', s=40,rasterized = True, zorder=3)
        ########## colorbar
        f.colorbar(cs2, cax=axes[0,2])
        cs1.set_clim(vgfb2_data_top.min(),vgfb2_data_top.max())
        cs2.set_clim(vgfb2_data_top.min(),vgfb2_data_top.max())
        axes[0,2].set_ylabel('$vgf_b^2$')
        axes[0,2].yaxis.set_label_position("right")
        ####################### Titles and labels
        if self.model_to_fit[:6]=='Kurucz':
            label = '%d $K$, %.3f $R_{\odot}$, %.4f $L_{\odot}$, logg=%.2f, MH=%.2f, alpha=%.1f'%(self.Te, self.R, self.L, self.logg, self.MH, self.alpha)
        if self.model_to_fit=='Koester':
            label = '%d $K$, %.3f $R_{\odot}$, %.4f $L_{\odot}$, logg=%.2f'%(self.Te, self.R, self.L, self.logg)
        axes[0,0].set_title(self.name+'       ' + self.model_to_fit + '       ' + label, x=0, y=1, ha='left')
        axes[2,0].set_title('$\chi^2$ = %.2f\n$\chi_r^2$ = %.2f\n$vgf_b^2$ = %.2f'%(self.chi2, self.chi2_r, self.vgfb2),
                             x=0.98,y=0.9, ha='right', va='top')
        ####################### decoration and saving
        axes[1,1].set_ylim(R_data_top.min()/1.1,R_data_top.max()*1.1)

        axes[0,0].legend(scatterpoints=1, loc='upper center', ncol=5,frameon=False,handletextpad=0.3, borderpad=0.1)
        axes[0,1].legend(scatterpoints=1)
        if save_plot:   
            if not os.path.exists(folder_path):
                os.makedirs(folder_path)
            plt.savefig (folder_path +self.name+'_' + self.model_to_fit +'_'+str(self.Te)+'_'+str(cycle)+'.jpg', format='jpg', dpi=300)
        if not show_plot:
                plt.close()        
                
    def calculate_chi2_noisy(self, cycle, total_iterations=100, verbose=False, refit=True, percentile_threshold=5, plot=False, minimising_param='chi2'):
        '''
        Calculate chi2 after adding noise to the observed data
        
        - Takes ``da_model`` and scales it for each sf in ``sf_list_noisy``.
        - Creates ``total_iterations`` versions of ``da_obs`` as ``da_obs_noisy``.
        - Subtracts the ``da_model`` from ``da_obs_noisy`` and calculates corresponding chi2_i for each datapoint.
        - Sums the chi2_i for each model and saves the chi2 values into DataArray (``da_chi2_stacked``). 
        - Simultaneously calculates vgfb2
        - Converts ``da_chi2_stacked`` into a dataframe (``df_chi2_noisy``)
        - Saves ``df_chi2_noisy`` into a .csv file
        
        Returns
        -------
        df_chi2_noisy : DataFrame
            chi2 sorted DataFrame with Te [K], logg, (MH), sf.            
        '''
        top_chi2           = np.nanpercentile(self.df_chi2.chi2, percentile_threshold)
        top_df_chi2        = self.df_chi2[self.df_chi2.chi2<top_chi2]
        sf_min             = top_df_chi2.sf.min()
        sf_max             = top_df_chi2.sf.max()
        self.sf_list_noisy = 10**np.arange(np.log10(sf_min),np.log10(sf_max), 0.01) 

        chi2_file_name = DIR_OUTPUTS + 'chi_files/'+self.type +'_'+self.name+'_'+self.model_to_fit+'_noisy_'+str(cycle)+'.csv'

        # If the chi2 file exists, it will be read (depends on refit = True or False)
        if os.path.isfile(chi2_file_name): 
            if not refit:
                df_chi2_noisy      = pd.read_csv(chi2_file_name)
                if verbose:     
                    print('    Reading %s'%chi2_file_name)
                    print('    WARNING! give "refit=True" if you want to rerun the fitting process.')
                    print(df_chi2_noisy.head())
                if plot: 
                    size_arr    = 20 * df_chi2_noisy['logg']**0
                    size_arr[0] = 100
                    sns.pairplot(df_chi2_noisy, vars=['Te','logg','sf','chi2'],corner=True, diag_kind="hist", plot_kws={"s": size_arr})
                self.df_chi2_noisy = df_chi2_noisy
                return 
            if refit: print('    WARNING: '+chi2_file_name+' file will be overwritten.')

        da_model = self.da_model
        if self.model_to_fit[:6]=='Kurucz':
            nd_arr = np.empty((len(da_model.Te), len(da_model.logg), len(da_model.MH), len(da_model.alpha), len(self.sf_list_noisy)))
            nd_arr.fill(np.nan)
            da_chi2_stacked = xr.DataArray(nd_arr, coords={'Te'   :da_model.Te,
                                                           'logg' :da_model.logg,
                                                           'MH'   :da_model.MH,
                                                           'alpha':da_model.alpha,
                                                           'sf'   :self.sf_list_noisy})
        if self.model_to_fit=='Koester':
            nd_arr = np.empty((len(da_model.Te), len(da_model.logg), len(self.sf_list_noisy)))
            nd_arr.fill(np.nan)
            da_chi2_stacked = xr.DataArray(nd_arr, coords={'Te'  :da_model.Te,
                                                           'logg':da_model.logg,
                                                           'sf'  :self.sf_list_noisy})

        da_obs_noisy                    = self.da_model.copy()
        da_obs_noisy.attrs['long_name'] = 'Noisy flux'

        df_chi2_noisy = pd.DataFrame()

        for seed in range (total_iterations):
            print_progress(seed, total_iterations=total_iterations, step=2, message='    Calculating noisy chi2')

            np.random.seed(seed)
            noise = np.random.normal(0, 1, [self.N_points]) 
            if seed == 0:    # No noise added for 0th iteration
                noisy_flux = self.flux['flux']
            else:            # Adding Gaussian noise to the flux 
                noisy_flux = self.flux['flux'] + self.flux['error']*noise

            for filter_name in self.flux.index:
                da_obs_noisy = da_obs_noisy.where(da_obs_noisy.FilterID!=filter_name,noisy_flux[filter_name])

            for idx in range (len(self.sf_list_noisy)):
                sf = self.sf_list_noisy[idx]
                da_residual = da_obs_noisy - (da_model * sf)

                da_chi2_i   = (da_residual/self.da_obs_error)**2
                da_chi2     = da_chi2_i.sum('FilterID',min_count=1)
                da_chi2_stacked.loc[dict(sf=sf)] = da_chi2

                if minimising_param != 'chi2':
                    da_vgf2_i   = (da_residual/self.da_obs_error_2percent)**2
                    da_vgf2     = da_vgf2_i.sum('FilterID',min_count=1)
                    da_vgf2_stacked.loc[dict(sf=sf)] = da_vgf2

                    da_vgfb2_i  = (da_residual/self.da_obs_error_10percent)**2
                    da_vgfb2    = da_vgfb2_i.sum('FilterID',min_count=1)
                    da_vgfb2_stacked.loc[dict(sf=sf)] = da_vgfb2
                    
            df_chi2         = da_chi2_stacked.to_dataframe(name='chi2').reset_index()

            if minimising_param not in ['chi2', 'vgf2', 'vgfb2']:
                raise Exception('    Error: The minimising_param should be one of the following: chi2, vgf2, vgfb2')
            df_chi2 = df_chi2.sort_values(minimising_param)
            df_chi2.reset_index(drop=True, inplace=True)

            df_chi2_noisy = pd.concat([df_chi2_noisy, df_chi2.head(1)], ignore_index = True)

        if verbose: print('\n    ',df_chi2_noisy.head())
        if plot: 
            size_arr    = 20 * df_chi2_noisy['logg']**0
            size_arr[0] = 100
            sns.pairplot(df_chi2_noisy, vars=['Te','logg','sf','chi2'],corner=True, diag_kind="hist",plot_kws={"s": size_arr})

        self.df_chi2_noisy = df_chi2_noisy

        if not os.path.exists(DIR_OUTPUTS + 'chi_files/'): os.makedirs(DIR_OUTPUTS + 'chi_files/')
        if verbose: print('    Saving %s'%chi2_file_name)
        df_chi2_noisy.to_csv(chi2_file_name, index=False, header=True, sep=',')

    def get_parameters_from_noisy_chi2(self, verbose=False):
        '''
        Estimates best fit parameters from noisy least chi2 fit
        
        - model_flux, Te, sf, R, L, chi2, vgf, vgfb
        - Errors in Te, sf, logg, MH, R, L based on their noisy distribution
        - Updates these parameter to parent class object
        '''
        self.Te_median, self.Te_error_lower, self.Te_error_upper = get_realistic_errors_from_iterations(self.df_chi2_noisy['Te'],self.da_model['Te'].values,para_type='Te')
        self.logg_median, self.logg_error_lower, self.logg_error_upper = get_realistic_errors_from_iterations(self.df_chi2_noisy['logg'],self.da_model['logg'].values,para_type='logg')
        self.sf_median, self.sf_error_lower, self.sf_error_upper = get_realistic_errors_from_iterations(self.df_chi2_noisy['sf'],self.sf_list_noisy,para_type='sf')

        self.R_median      = calc_radius(self.sf_median, self.distance)
        # error contribution from noisy fitting
        self.R_error_lower = self.R_median - calc_radius(self.sf_median-self.sf_error_lower, self.distance)
        self.R_error_upper = calc_radius(self.sf_median+self.sf_error_upper, self.distance) - self.R_median
        # error after including distance error,     deltaR = R * deltaD/D
        self.R_error_lower = addition_in_quadrature(self.R_error_lower, self.R_median * self.distance_err/self.distance)
        self.R_error_upper = addition_in_quadrature(self.R_error_upper, self.R_median * self.distance_err/self.distance)

        self.L_median      = calc_luminosity(self.R_median, self.Te_median)
        # error contribution from noisy fitting
        self.L_error_lower = self.L_median - calc_luminosity(self.R_median-self.R_error_lower, self.Te_median-self.Te_error_lower)
        self.L_error_upper = calc_luminosity(self.R_median+self.R_error_upper, self.Te_median+self.Te_error_upper) - self.L_median
        # error after including distance error,     deltaL = L * (2 * deltaD/D)
        self.L_error_lower = addition_in_quadrature(self.L_error_lower, self.L_median * 2. * self.distance_err/self.distance)
        self.L_error_upper = addition_in_quadrature(self.L_error_upper, self.L_median * 2. * self.distance_err/self.distance)

        if self.Te!=self.Te_median:
            print ('    WARNING! Median of noisy "Te" is not same as "Te".\n        THE FIT MAY NOT BE RELIABLE.')

        if self.model_to_fit[:6]  =='Kurucz':
            self.MH_median, self.MH_error_lower, self.MH_error_upper = get_realistic_errors_from_iterations(self.df_chi2_noisy['MH'],self.da_model['MH'].values,para_type='MH')
            self.alpha_median, self.alpha_error_lower, self.alpha_error_upper = get_realistic_errors_from_iterations(self.df_chi2_noisy['alpha'],self.da_model['alpha'].values,para_type='alpha')
            best_fit_flux_median         = self.da_model.sel(Te=self.Te_median).sel(logg=self.logg_median).sel(MH=self.MH_median).sel(alpha=self.alpha_median)

        if self.model_to_fit       =='Koester':
            self.MH_median, self.MH_error_lower, self.MH_error_upper  = np.nan,np.nan,np.nan
            self.alpha_median, self.alpha_error_lower, self.alpha_error_upper  = np.nan,np.nan,np.nan
            best_fit_flux_median           = self.da_model.sel(Te=self.Te_median).sel(logg=self.logg_median)

        self.flux['model_flux_median'] = self.sf_median * best_fit_flux_median
        _residual_flux                 = self.flux['flux'] - self.flux['model_flux_median']
        self.flux['chi2_i_median']     = _residual_flux**2 / self.flux['error']**2
        self.flux['vgf2_i_median']     = _residual_flux**2 / self.flux['error_2percent']**2
        self.flux['vgfb2_i_median']    = _residual_flux**2 / self.flux['error_10percent']**2
        self.chi2_median               = self.flux['chi2_i_median'].sum()
        self.chi2_r_median             = self.chi2_median/self.N_Np
        self.vgf2_median               = self.flux['vgf2_i_median'].sum()/self.N_Np
        self.vgfb2_median              = self.flux['vgfb2_i_median'].sum()/self.N_Np

        with xr.open_dataarray(self.model_file_name) as da_model_all:
            da_model_all = da_model_all.sel(FilterID=self.flux_all.index)
            if self.model_to_fit[:6]=='Kurucz':
                da_model_all = da_model_all.sel(Te=self.Te_median).sel(logg=self.logg_median).sel(MH=self.MH_median).sel(alpha=self.alpha_median)
            if self.model_to_fit=='Koester':
                da_model_all = da_model_all.sel(Te=self.Te_median).sel(logg=self.logg_median)
            df_model_all = da_model_all.to_dataframe(name='model_flux_median')
        self.flux_all['model_flux_median'] = df_model_all['model_flux_median']*self.sf_median

        if verbose: 
            print('\n    Fitting parameters from noisy iterations:')
            print('        Te   = %d (+%d-%d)'%(self.Te_median, self.Te_error_lower, self.Te_error_upper))
            print('        logg = %.2f (+%.2f-%.2f)'%(self.logg_median, self.logg_error_lower, self.logg_error_upper))
            print('        MH   = %.2f (+%.2f-%.2f)'%(self.MH_median, self.MH_error_lower, self.MH_error_upper))
            print('        alpha= %.2f (+%.2f-%.2f)'%(self.alpha_median, self.alpha_error_lower, self.alpha_error_upper))
            print('        sf   = %.2e (+%.2e-%.2e)'%(self.sf_median, self.sf_error_lower, self.sf_error_upper))
            print('        R    = %.4f (+%.4f-%.4f)'%(self.R_median, self.R_error_lower, self.R_error_upper))
            print('        L    = %.4f (+%.4f-%.4f)'%(self.L_median, self.L_error_lower, self.L_error_upper))
            print('        chi2 = %.1f, vgf2 = %.1f, vgfb2 = %.1f'%(self.chi2_median, self.vgf2_median, self.vgfb2_median))

    def plot_fitted_noisy_SED(self,cycle,plot_noisy_SEDs=True, save_plot=True,show_plot=True,folder_path='plots/single_SEDs'):
        '''
        Creates noisy SED
        
        - Includes HR diagram with isochrone of logAge 8/9/10, WD models for mass 0.2/1.3.
        - Includes radius vs Te plot
        - Includes chi2 vs Te plot. This plot should have explicit global minima in a good fit.
        '''
        f, axes = self.plot_skeletal_SED()

        axes[0,0].plot(self.flux_all['wavelength'], self.flux_all['model_flux'], color='green', linestyle='-',label ='Model', lw=1)

        ########## plotting noisy fits
        if plot_noisy_SEDs:
            for idx, sf_noisy in enumerate(self.df_chi2_noisy['sf']):
                print_progress(idx, total_iterations=100, step=5, message='    Plotting noisy SEDs')
                Te_noisy   = self.df_chi2_noisy['Te'][idx]
                logg_noisy   = self.df_chi2_noisy['logg'][idx]
                if self.model_to_fit[:6]  =='Kurucz':
                    MH_noisy     = self.df_chi2_noisy['MH'][idx]
                    alpha_noisy  = self.df_chi2_noisy['alpha'][idx]
                    flux_B_noisy = self.da_model.sel(Te=Te_noisy).sel(logg=logg_noisy).sel(MH=MH_noisy).sel(alpha=alpha_noisy).values * sf_noisy 
                if self.model_to_fit  =='Koester':
                    flux_B_noisy = self.da_model.sel(Te=Te_noisy).sel(logg=logg_noisy).values * sf_noisy
                axes[0,0].plot(self.flux['wavelength'], flux_B_noisy, color='cyan', linestyle='-',label='', lw=0.2, alpha=0.2,zorder=0)
                residual_noisy = (self.flux['flux'] - flux_B_noisy)/self.flux['flux']
                axes[1,0].plot(self.flux['wavelength'], residual_noisy, color='cyan', linestyle='-',label='', lw=0.2, alpha=0.2,zorder=0)   
        ########## Fractional residual
        _fractional_residual = (self.flux_all['flux']-self.flux_all['model_flux'])/self.flux_all['flux']
        axes[1,0].plot(self.flux_all['wavelength'], _fractional_residual,label='',marker='',color='green',lw=1, linestyle='-')
        ########## chi2_i
        axes[2,0].plot(self.flux['wavelength'], self.flux['chi2_i'],label='',marker='o',color='green', linestyle='-',lw=1)
        ########## L vs T (HR diagram)
        self.df_chi2_noisy['R']    = calc_radius(self.df_chi2_noisy['sf'], self.distance)
        self.df_chi2_noisy['L']    = calc_luminosity(self.df_chi2_noisy['R'], self.df_chi2_noisy['Te'])
        Te_data     = self.df_chi2_noisy.Te
        L_data      = self.df_chi2_noisy.L
        R_data      = self.df_chi2_noisy.R
        chi2_r_data = self.df_chi2_noisy.chi2/self.N_Np
        axes[0,1].scatter(Te_data,L_data, marker='.',label='',s=10,rasterized = True, zorder=2)
        axes[0,1].scatter(self.Te_median,self.L_median, marker='s', label='Best fit',c='r', s=40,rasterized = True, zorder=3)
        ########## R vs T    
        axes[1,1].scatter(Te_data,R_data,s=10,rasterized = True, zorder=2)
        axes[1,1].scatter(self.Te_median,self.R_median, marker='s', label='',c='r', s=40,rasterized = True, zorder=3)
        ########## chi2 vs T
        cs2 = axes[2,1].scatter(Te_data,chi2_r_data,c=R_data, cmap='summer',s=10,rasterized = True,zorder=2)  
        axes[2,1].scatter(self.Te_median,self.chi2_r_median, marker='s', label='',c='r', s=40,rasterized = True, zorder=3)
        ########## colorbar
        f.colorbar(cs2, cax=axes[0,2])
        axes[0,2].set_ylabel('R/$R_{\odot}$')
        axes[0,2].yaxis.set_label_position("right")
        ####################### errorbars
        axes[0,1].axvspan(self.Te-self.Te_error_lower, self.Te+self.Te_error_upper, alpha=0.2, color='dodgerblue', zorder=1)
        axes[0,1].axhspan(self.L-self.L_error_lower, self.L+self.L_error_upper, alpha=0.2, color='dodgerblue', zorder=1)        
        axes[1,1].axvspan(self.Te-self.Te_error_lower, self.Te+self.Te_error_upper, alpha=0.2, color='dodgerblue', zorder=1)
        axes[1,1].axhspan(self.R-self.R_error_lower, self.R+self.R_error_upper, alpha=0.2, color='dodgerblue', zorder=1)        
        axes[2,1].axvspan(self.Te-self.Te_error_lower, self.Te+self.Te_error_upper, alpha=0.2, color='dodgerblue', zorder=1)
        ####################### Titles and labels
        if self.model_to_fit[:6]=='Kurucz':
            label = '%d$^{+%d}_{-%d}$ $K$, '%(self.Te, self.Te_error_lower, self.Te_error_upper) +\
                    '%.4f$^{+%.4f}_{-%.4f}$ $R_{\odot}$, '%(self.R, self.R_error_lower, self.R_error_upper) +\
                    '%.4f$^{+%.4f}_{-%.4f}$ $L_{\odot}$, '%(self.L, self.L_error_lower, self.L_error_upper) +\
                    'logg=%.2f, MH=%.2f, alpha=%.1f'%(self.logg, self.MH, self.alpha)
        if self.model_to_fit=='Koester':
            label = '%d$^{+%d}_{-%d}$ $K$, '%(self.Te, self.Te_error_lower, self.Te_error_upper) +\
                    '%.4f$^{+%.4f}_{-%.4f}$ $R_{\odot}$, '%(self.R, self.R_error_lower, self.R_error_upper) +\
                    '%.4f$^{+%.4f}_{-%.4f}$ $L_{\odot}$, '%(self.L, self.L_error_lower, self.L_error_upper) +\
                    'logg=%.2f'%(self.logg)
        axes[0,0].set_title(self.name+'       ' + self.model_to_fit + '       ' + label, x=0, y=1, ha='left')
        axes[2,0].set_title('$\chi^2$ = %.2f\n$\chi_r^2$ = %.2f\n$vgf_b^2$ = %.2f'%(self.chi2, self.chi2_r,self.vgfb2),
                            x=0.98,y=0.9, ha='right', va='top')
        ####################### decoration and saving
        axes[1,1].set_ylim(self.R-3*self.R_error_lower, self.R+3*self.R_error_upper)

        axes[1,1].set_yscale('linear')
        axes[2,1].set_yscale('linear')

        axes[0,0].legend(scatterpoints=1, loc='upper center', ncol=5,frameon=False,handletextpad=0.3, borderpad=0.1)
        axes[0,1].legend(scatterpoints=1)
        if save_plot:   
            if not os.path.exists(folder_path):
                os.makedirs(folder_path)
            plt.savefig (folder_path + '%s_%s_%d_noisy.jpg'%(self.name, self.model_to_fit, cycle), format='jpg', dpi=300)#,bbox_inches='tight')
        if not show_plot:
                plt.close()   
                
    def save_object(self, cycle, save_log=True):
        '''
        Saving class object and log

        - Saves the SingleStar object as a pickel.
        - Saves class attributes in a .csv log file 
        '''
        _self_copy = copy.deepcopy(self)
        _self_copy.time  = strftime("%Y-%m-%d %H:%M:%S", gmtime())
        _self_copy.cycle = cycle
        try:
            for attrName in ['da_obs', 'da_obs_error', 'da_obs_error_2percent', 'da_obs_error_10percent','df_chi2', 'df_chi2_noisy']:
                delattr(_self_copy, attrName)
        except: pass

        file_name = DIR_OUTPUTS + 'pickels/SingleStar_%s_%d.pkl'%(self.name,cycle)
        if not os.path.exists(DIR_OUTPUTS + 'pickels/'): os.makedirs(DIR_OUTPUTS + 'pickels/')
        with open(file_name, 'wb') as outp:
            pickle.dump(_self_copy, outp, pickle.HIGHEST_PROTOCOL)

        if save_log: 
            log_file_name              = 'log_single_fitting.csv'

            # remove arrays and dataframes from the log file
            try:
                for attrName in ['flux', 'flux_all','not_fitted', 'da_model', 'sf_list', 'sf_list_noisy']:
                    delattr(_self_copy, attrName)
            except: pass
        
            _self_copy.filters_to_drop = ' '.join(_self_copy.__dict__['filters_to_drop'])
            _self_copy.wavelength_range = ' '.join(map(str, _self_copy.__dict__['wavelength_range']))

            for idx, key in enumerate(_self_copy.__dict__.keys()):
                if type(_self_copy.__dict__[key])==list:
                    if len(_self_copy.__dict__[key])>1:
                        raise Exception('    ERROR: Some parameter (e.g. %s) in the object is list-like. Convert the list into string for saving (similar to "filter_to_drop" in above lines).'%key)
                                      
            df_log = pd.DataFrame(_self_copy.__dict__, index=[0])
                
            if not os.path.isfile(log_file_name):
                print('    Creating %s and saving log'%log_file_name)
            else:
                print('    Saving log in %s'%log_file_name)
                df_past_log = pd.read_csv(log_file_name)
                df_log = pd.concat([df_past_log,df_log])

            df_log.to_csv(log_file_name, index=False)

    def to_BinaryStar(self, model_to_fit_B, free_para, filters_to_drop=[], verbose=False):
        '''
        Initialises BinaryStar object from a SingleStar object
        '''
        bi_star = BinaryStar(name            = self.name,
                             model_to_fit_A  = self.model_to_fit,
                             model_to_fit_B  = model_to_fit_B,
                             dir_obs         = self.dir_obs,
                             distance        = self.distance,
                             distance_err    = self.distance_err,
                             free_para       = free_para,
                             filters_to_drop = filters_to_drop)

        bi_star.flux                = self.flux_all
        bi_star.flux                = bi_star.flux.drop(columns=['model_flux'])
        bi_star.flux.rename({'model_flux_median': 'model_flux_A'}, axis=1, inplace=True)
        bi_star.flux_all            = bi_star.flux.copy()
        bi_star.drop_filters(verbose=verbose)
        bi_star.model_file_name_A   =  self.model_file_name
        bi_star.Te_A                =  self.Te_median  
        bi_star.Te_error_lower_A    = self.Te_error_lower  
        bi_star.Te_error_upper_A    = self.Te_error_upper  
        bi_star.logg_A              = self.logg_median  
        bi_star.logg_error_lower_A  = self.logg_error_lower  
        bi_star.logg_error_upper_A  = self.logg_error_upper  
        bi_star.sf_A                = self.sf_median  
        bi_star.sf_error_lower_A    = self.sf_error_lower  
        bi_star.sf_error_upper_A    = self.sf_error_upper  
        bi_star.MH_A                = self.MH_median  
        bi_star.MH_error_lower_A    = self.MH_error_lower  
        bi_star.MH_error_upper_A    = self.MH_error_upper
        bi_star.R_A                 = self.R_median  
        bi_star.R_error_lower_A     = self.R_error_lower  
        bi_star.R_error_upper_A     = self.R_error_upper
        bi_star.L_A                 = self.L_median  
        bi_star.L_error_lower_A     = self.L_error_lower  
        bi_star.L_error_upper_A     = self.L_error_upper

        if verbose:
            print('    Fitting parameters of A component:')
            print('        Te   = %d (+%d-%d)'%(self.Te, self.Te_error_lower, self.Te_error_upper))
            print('        logg = %.2f (+%.2f-%.2f)'%(self.logg, self.logg_error_lower, self.logg_error_upper))
            print('        MH   = %.2f (+%.2f-%.2f)'%(self.MH, self.MH_error_lower, self.MH_error_upper))
            print('        sf   = %.2e (+%.2e-%.2e)'%(self.sf, self.sf_error_lower, self.sf_error_upper))
        return bi_star

    
def initialise_BinaryStar_from_VOSA(name, model_to_fit_A, model_to_fit_B, dir_obs, distance, distance_err, free_para, filters_to_drop=[], verbose=False):
    '''
    Initialises BinaryStar object using a file downloded from VOSA.
    
    Parameters
    ----------
    name : str
        Name of the source
    model_to_fit_A : str
        Model name of the A component
    model_to_fit_B : str
        Model name of the B component
    dir_obs : str
        Path to the VOSA folder (till ``.../object/``)
    distance : float
        Distance [pc]
    distance_err : float
        Error in distance [pc]    
    free_para : int
        Number of free parameters
    filters_to_drop : list
        List of filters to be dropped before fitting SED. Default value is ``[]``.
    
    Returns
    -------
    bi_star : BinaryStar
        BinaryStar object
    '''
    bi_star = BinaryStar(name            = name,
                         model_to_fit_A  = model_to_fit_A,
                         model_to_fit_B  = model_to_fit_B,
                         dir_obs         = dir_obs,
                         distance        = distance,
                         distance_err    = distance_err,
                         free_para       = free_para,
                         filters_to_drop = filters_to_drop)

    # Reading parameters derived by VOSA
    file_name    = dir_obs + name +'/bestfitp/'+ name +'.bfit.phot.dat'
    flux         = pd.read_csv(file_name, engine='python', comment='#', delim_whitespace= True, skipinitialspace=True, header=None)
    flux.columns = ['FilterID','wavelength','Obs.Flux','Obs.Error','flux','error','model_flux_A','Fitted','Excess','FitExc','UpLim']
    
    # Removing filters noted as "upper limit"
    if verbose: print('    WARNING: %s removed due to upper limit'%(flux[flux['UpLim']=='1']['FilterID'].values))
    flux         = flux[flux['UpLim'] == '---']
    flux         = flux.drop(columns=['Obs.Flux','Obs.Error','Fitted','Excess','FitExc','UpLim'])

    # Replacing zeros in errors with 110% of max error
    flux['error_fraction']             = flux['error']/flux['flux']
    flux.loc[flux.error == 0, 'error'] = flux['flux']*(flux['error_fraction'].max()*1.1)
    # Recalculating errors_fraction 
    flux['error_fraction']             = flux['error']/flux['flux']

    # error modification for calculating vgf (minimum error = 2%) and vgfb (minimum error = 10%)
    flux['error_2percent']  = np.where(flux['error_fraction']<0.02, 0.02, flux['error_fraction'])*flux['flux']
    flux['error_10percent'] = np.where(flux['error_fraction']<0.10, 0.10, flux['error_fraction'])*flux['flux']

    flux['log_wavelength']  = np.log10(flux['wavelength'])
    flux['log_flux']        = np.log10(flux['flux'])
    flux['log_error']       = 0.434 * flux['error']/flux['flux']
    
    bi_star.flux            = flux.set_index('FilterID')
    # Saving the original flux file for plotting purpose
    bi_star.flux_all        = bi_star.flux.copy()
    
    bi_star.drop_filters(verbose=verbose)
    
    bi_star.Te_A                   = get_parameter_from_VOSA_file(param='Teff', file_name=file_name)
    bi_star.Te_error_lower_A       = bi_star.Te_A - get_parameter_from_VOSA_file(param='Teff_min', file_name=file_name)
    bi_star.Te_error_upper_A       = get_parameter_from_VOSA_file(param='Teff_max', file_name=file_name) - bi_star.Te_A

    bi_star.logg_A                 = get_parameter_from_VOSA_file(param='logg', file_name=file_name)
    bi_star.logg_error_lower_A     = bi_star.logg_A - get_parameter_from_VOSA_file(param='logg_min', file_name=file_name)
    bi_star.logg_error_upper_A     = get_parameter_from_VOSA_file(param='logg_max', file_name=file_name) - bi_star.logg_A

    if bi_star.model_to_fit_A[:6] == 'Kurucz':
        bi_star.MH_A               = get_parameter_from_VOSA_file(param='Meta', file_name=file_name)
        bi_star.MH_error_lower_A   = bi_star.MH_A - get_parameter_from_VOSA_file(param='Meta_min', file_name=file_name)
        bi_star.MH_error_upper_A   = get_parameter_from_VOSA_file(param='Meta_max', file_name=file_name) - bi_star.MH_A
    else:
        bi_star.MH_A               = np.nan
        bi_star.MH_error_lower_A   = np.nan
        bi_star.MH_error_upper_A   = np.nan

    bi_star.sf_A                   = get_parameter_from_VOSA_file(param='Md', file_name=file_name)

    bi_star.R_A                    = get_parameter_from_VOSA_file(param='Rad1', file_name=file_name)
    bi_star.R_error_lower_A        = get_parameter_from_VOSA_file(param='eRad1', file_name=file_name)
    bi_star.R_error_upper_A        = bi_star.R_error_lower_A
    
    bi_star.L_A                    = get_parameter_from_VOSA_file(param='Lbol', file_name=file_name)
    bi_star.L_error_lower_A        = get_parameter_from_VOSA_file(param='Lberr', file_name=file_name)
    bi_star.L_error_upper_A        = bi_star.L_error_lower_A

    if verbose:
        print('    Fitting parameters of A component:')
        print('        Te   = %d (+%d-%d)'%(bi_star.Te_A, bi_star.Te_error_lower_A, bi_star.Te_error_upper_A))
        print('        logg = %.2f (+%.2f-%.2f)'%(bi_star.logg_A, bi_star.logg_error_lower_A, bi_star.logg_error_upper_A))
        print('        MH   = %.2f (+%.2f-%.2f)'%(bi_star.MH_A, bi_star.MH_error_lower_A, bi_star.MH_error_upper_A))
        print('        sf   = %.2e'%(bi_star.sf_A))
        print('        Total filters: %d' %len(bi_star.flux))
    return bi_star


class BinaryStar:
    '''
    Class for SED fitting of a binary star
    '''
    def __init__(self, name, model_to_fit_A,model_to_fit_B, dir_obs, distance, distance_err, free_para, filters_to_drop=[], verbose=False):
        '''
        Initialises the binary star using name, distance, model_to_fit etc.
        
        Parameters
        ----------
        name : str
            Name of the source
        model_to_fit_A : str
            Model name of the A component
        model_to_fit_B : str
            Model name of the B component
        dir_obs : str
            Path to the observed SED '.csv' file  
        distance : float
            Distance [pc]
        distance_err : float
            Error in distance [pc]    
        free_para : int
            Number of free parameters
        filters_to_drop : list
            List of filters to be dropped before fitting SED. Default value is ``[]``.
        '''
        self.name               = name
        self.model_to_fit_A     = model_to_fit_A
        self.model_to_fit_B     = model_to_fit_B
        self.dir_obs            = dir_obs
        self.distance           = distance
        self.distance_err       = distance_err
        self.free_para          = free_para
        self.filters_to_drop    = filters_to_drop
        self.type               = 'binary'
    
        if verbose: print('================\n%s\n----------------' %self.name)
        
    def drop_filters(self, verbose=False):
        '''
        Dropping filters before further SED fitting.
        
        - Some filters will have to be removed while fitting.
        - You can remove individual filters based on science case. For example:
        
            - In case of IR excess, you can remove Wise filters
            - For bad chi2, you can remove specific filters
            
        Returns
        -------
        flux : DataFrame
            DataFrame containing only required filters
        not_fitted : DataFrame
            DataFrame containing only removed filters           
        N_points : int
            Updated number of data points
        N_Np : int
            Updated degrees of freedom (``N_points - free_para``)
        '''    
        self.not_fitted     = self.flux[(self.flux['wavelength']<0)]

        for filter_name in self.filters_to_drop:
            self.not_fitted = pd.concat([self.not_fitted, self.flux[(self.flux.index==filter_name)]])
            self.flux       = self.flux.drop(index=filter_name)
         
        self.N_points = len(self.flux)
        self.N_Np     = self.N_points-self.free_para
        
        # printing filters to be fitted and not_fitted
        if verbose:
            _t1, _t2            = pd.DataFrame(),pd.DataFrame()    
            _t1['to_be_fitted'] = self.flux.index.values
            _t1['wavelength']   = self.flux.wavelength.values
            _t1                 = _t1.set_index('wavelength')
            _t2['not_fitted']   = self.not_fitted.index.values
            _t2['wavelength']   = self.not_fitted.wavelength.values
            _t2                 = _t2.set_index('wavelength')
            _filter_table       = pd.concat([_t1,_t2],sort=True)    
            print('\n    RUNNING: drop_filters ')
            print(_filter_table.sort_index().fillna(''))
            print('    Filters to fit: %d' %len(self.flux))
            
    def visualise_A_component_SED(self,excess_cutoff=0.5, save_plot=True, show_plot=True, folder_path='plots/binary_SEDs/'):
        '''
        Plotting SED of single compoent as fitted by VOSA.
        - Default thresold for identifying excess = 0.5
        - Saving by default
        '''
        f, axes = plt.subplots(figsize=(6,6),nrows = 2, ncols = 1)
        [axi.set_axis_off() for axi in axes.ravel()]
        axes[0] = f.add_axes([0.15, 0.35, 0.8, 0.55])
        axes[1] = f.add_axes([0.15, 0.15, 0.8, 0.2])
        ####################### SED
        axes[0].scatter(self.not_fitted['wavelength'], self.not_fitted['flux'], color='k', marker='o',label ='No Fit', s=30, facecolors='none',zorder=2)
        axes[0].errorbar(self.flux['wavelength'], self.flux['flux'], yerr=self.flux['error'],color='k', label='Obs',fmt='none',lw=2, capsize=4,zorder=3)
        axes[0].plot(self.flux['wavelength'], self.flux['model_flux_A'], color='green', linestyle='-',label ='Model', lw=1)

        ########## Fractional residual
        # 3sigma errors
        axes[1].errorbar(self.flux['wavelength'], self.flux['flux']-self.flux['flux'], yerr=3*self.flux['error_fraction'],color='0.5', label='3-$\sigma$',fmt='none',lw=1, capsize=0) 
        # 1sigma error
        axes[1].errorbar(self.flux['wavelength'], self.flux['flux']-self.flux['flux'], yerr=self.flux['error_fraction'],color='k', label='1-$\sigma$',fmt='none',lw=2, capsize=4) 
        axes[1].errorbar(self.not_fitted['wavelength'], self.not_fitted['flux']-self.not_fitted['flux'], yerr=self.not_fitted['error_fraction'],color='0.5', label='',fmt='none',lw=2, capsize=4) 
        
        _fractional_residual_A = (self.flux['flux']-self.flux['model_flux_A'])/self.flux['flux']
        axes[1].plot(self.flux['wavelength'], _fractional_residual_A,label='',marker='',color='r',lw=1, linestyle='--')

        
        flux_A_excess = self.flux[((self.flux['flux']-self.flux['model_flux_A'])/self.flux['flux']>excess_cutoff)]
        axes[1].scatter(flux_A_excess['wavelength'], (flux_A_excess['flux']-flux_A_excess['model_flux_A'])/flux_A_excess['flux'], color='red',marker='o',label ='Excess (%d)' %len(flux_A_excess), s=30)

        ####################### Titles and labels
        label_A = 'A (' + str(self.Te_A) + ' K, logg=' + str(self.logg_A) + ')'
        axes[0].set_title(self.name + '       ' + label_A, x=0, y=1, ha='left')
        axes[0].set_ylabel('Flux (erg s$^{-1}$ cm$^{-2}$ $\AA$$^{-1}$)')
        axes[1].set_ylabel('Residual')
        axes[1].set_xlabel('Wavelength ($\AA$)')
        ####################### axes range and scales
        axes[0].set_xscale('log')
        axes[0].set_yscale('log')
        axes[1].set_xscale('log')
        wave_min = self.flux_all['wavelength'].min()
        wave_max = self.flux_all['wavelength'].max()
        axes[0].set_xlim([wave_min/1.2,wave_max*1.2])
        axes[1].set_xlim([wave_min/1.2,wave_max*1.2])

        flux_min = self.flux_all['flux'].min()
        flux_max = self.flux_all['flux'].max()
        axes[0].set_ylim([flux_min/100,flux_max*5])

        axes[1].set_ylim([-0.2,1.5])
        ####################### decoration    
        axes[1].axhline(excess_cutoff, lw=1, c='pink',zorder=0)
        plt.setp(axes[0].get_xticklabels(),visible=False)
        axes[0].grid()
        axes[1].grid()
        axes[0].tick_params(which='both', direction='out', length=4)
        axes[1].tick_params(which='both', direction='out', length=4)
        axes[0].legend(scatterpoints=1, loc='upper center', ncol=5,frameon=False,handletextpad=0.3, borderpad=0.1)
        axes[1].legend(scatterpoints=1, loc='upper right', ncol=3,frameon=False,handletextpad=0.3, borderpad=0.1)

        if save_plot: ########## Saving file   
            if not os.path.exists(folder_path):
                os.makedirs(folder_path)
            plt.savefig (folder_path + '/A_SED_' +self.name+'_'+str(self.Te_A)+'_logg'+str(self.logg_A)+'.jpg', format='jpg', dpi=300)
        if not show_plot:
                plt.close()
            
    def read_model_file(self,verbose=False):
        '''
        Reads model file from the 'DIR_MODELS' folder as a DataArray.
        
        Returns
        -------
        da_model_B : DataArray
            nD array of all available SED models. The dimentions are FilterID, Te, logg (and more if applicable).
        '''
        self.model_file_name_B = get_model_file_name(self.model_to_fit_B)

        with xr.open_dataarray(self.model_file_name_B) as da_model_B:
            # editing the "Wavelengths" attribute
            df1                             = pd.DataFrame()
            df1['FilterID']                 = da_model_B['FilterID']
            df1                             = df1.set_index('FilterID')
            df1['Wavelength']               = da_model_B.attrs['Wavelengths']
            df2                             = pd.DataFrame()
            df2['FilterID']                 = self.flux.index
            df2                             = df2.set_index('FilterID')
            df2['Wavelength']               = df1['Wavelength']

            da_model_B.attrs['Wavelengths'] = df2['Wavelength'].values
            da_model_B = da_model_B.sel(FilterID=self.flux.index)

        if verbose:
            print('    Model for B comp.: ',da_model_B.coords)
            print('\n    Provide contrains on ', da_model_B.dims[1:], 'for B compoent faster execution.')
            print('        (By default all parameter space is used)')

        self.da_model_B = da_model_B

    def constrain_fitting_parameters(self, limits, verbose=False):    
        '''
        Cropping the model DataArray according to givel limits
        
        Returns
        -------
        da_model_B : DataArray
            Cropped ``da_model_B``
        '''
        if self.model_to_fit_B[:6] == 'Kurucz':
            self.da_model_B = self.da_model_B.sel(Te=slice(limits['Te_B'][0],limits['Te_B'][1]))
            self.da_model_B = self.da_model_B.sel(logg=slice(limits['logg_B'][0],limits['logg_B'][1]))
            self.da_model_B = self.da_model_B.sel(MH=slice(limits['MH_B'][0],limits['MH_B'][1]))
            self.da_model_B = self.da_model_B.sel(MH=slice(limits['alpha_B'][0],limits['alpha_B'][1]))
        if self.model_to_fit_B == 'Koester':
            self.da_model_B = self.da_model_B.sel(Te=slice(limits['Te_B'][0],limits['Te_B'][1]))
            self.da_model_B = self.da_model_B.sel(logg=slice(limits['logg_B'][0],limits['logg_B'][1]))
            
        if verbose: print('    da_model_B:\n',self.da_model_B.coords)
        
        self.create_star_dataarrays(verbose=verbose)

    def create_star_dataarrays(self, verbose=False):
        '''
        Create observed flux and flux error DataArrays with dimentions of ``da_model_B`` for vectorised arithmatic later.
        
        Returns
        -------
        da_res_A : DataArray
            nD array in the shape of ``da_model_B``. All filters have same value: residual after subtracting A component from observed SED
        da_obs_error : DataArray
            nD array in the shape of ``da_model_B``. All filters have same value: observed error
        da_obs_error_2percent : DataArray
            nD array in the shape of ``da_model_B``. All filters have same value: observed error (min 2%)
        da_obs_error_10percent : DataArray
            nD array in the shape of ``da_model_B``. All filters have same value: observed error (min 10%)
        '''    
        da_res_A               = self.da_model_B.copy()
        da_obs_error           = self.da_model_B.copy()
        da_obs_error_2percent  = self.da_model_B.copy()
        da_obs_error_10percent = self.da_model_B.copy()

        for filter_name in self.flux.index:
            da_res_A                = da_res_A.where(da_res_A.FilterID!=filter_name,self.flux['flux'][filter_name]-self.flux['model_flux_A'][filter_name])
            da_obs_error            = da_obs_error.where(da_obs_error.FilterID!=filter_name,self.flux['error'][filter_name])
            da_obs_error_2percent   = da_obs_error_2percent.where(da_obs_error_2percent.FilterID!=filter_name,self.flux['error_2percent'][filter_name])
            da_obs_error_10percent  = da_obs_error_10percent.where(da_obs_error_10percent.FilterID!=filter_name,self.flux['error_10percent'][filter_name])

        self.da_res_A                 = da_res_A
        self.da_obs_error           = da_obs_error
        self.da_obs_error_2percent  = da_obs_error_2percent
        self.da_obs_error_10percent = da_obs_error_10percent

        self.da_res_A.attrs['long_name'] = 'Residual after A component'
        self.da_obs_error.attrs['long_name'] = 'Error'
        self.da_obs_error_2percent.attrs['long_name'] = 'Error (2\% min)'
        self.da_obs_error_10percent.attrs['long_name'] = 'Error (10\% min)'
        if verbose: 
            print('    da observed flux:\n',da_res_A.head(2))
            print('    da observed error:\n',da_obs_error.head(2))
            print('    da observed error (2\% min):\n',da_obs_error_2percent.head(2))
            print('    da observed error (10\% min):\n',da_obs_error_10percent.head(2))

    def create_sf_list(self, R_B_min=0.00001, R_B_max=1., log_sf_stepsize=0.01):
        '''
        Create list of scaling factors
        
        Parameters
        ----------
        R_B_min : float
            Minimum expected radius of the B component
        R_B_max : float
            Maximum expected radius of the B component
        log_sf_stepsize : float
            Stepsize for ``sf`` in log space
        
        Returns
        -------
        sf_list : list
            List of scaling factor (NOT in log space) within ``R_B_min`` to ``R_B_max`` with stepsize of ``log_sf_stepsize``.
        '''
        sf_B_min = calc_sf(R_B_min,self.distance)
        sf_B_max = calc_sf(R_B_max,self.distance)
        self.sf_list_B = 10**np.arange(np.log10(sf_B_min),np.log10(sf_B_max), log_sf_stepsize)

    def calculate_chi2(self, cycle, refit=True, verbose=False, trim=True, minimising_param='chi2'):
        '''
        Calculate chi2
        
        - Takes ``da_model_B`` and scales it for each sf in ``sf_list_B``.
        - Subtracts the ``da_model_B`` from ``da_res_A`` and calculates corresponding chi2_i for each datapoint.
        - Sums the chi2_i for each model and saves the chi2 values into DataArray (``da_chi2_stacked``). 
        - Simultaneously calculates vgfb2
        - Converts ``da_chi2_stacked`` into a dataframe (``df_chi2``)
        - Calculates R [Rsun] and L [Lsun] for each fit
        - Saves ``df_chi2`` into a .csv file
        
        Returns
        -------
        da_chi2 : DataFrame
            chi2 sorted DataFrame with Te [K], logg, (MH), sf, R [Rsun], L [Lsun], vgfb2.
        '''
        chi2_file_name = DIR_OUTPUTS + 'chi_files/'+self.type +'_'+self.name+'_'+self.model_to_fit_A+'_'+self.model_to_fit_B+'_'+str(cycle)+'.csv'

        # If the chi2 file exists, it will be read (depends on refit = True or False)
        if os.path.isfile(chi2_file_name): 
            if not refit:
                df_chi2      = pd.read_csv(chi2_file_name)
                if verbose:     
                    print('    Reading %s'%chi2_file_name)
                    print('    WARNING! give "refit=True" if you want to rerun the fitting process.')
                    print(df_chi2.head())
                self.df_chi2 = df_chi2
                return 
            if refit: print('    WARNING: '+chi2_file_name+' file will be overwritten.')

        da_model = self.da_model_B
        if self.model_to_fit_B[:6]=='Kurucz':
            nd_arr = np.empty((len(da_model.Te), len(da_model.logg), len(da_model.MH), len(self.sf_list_B)))
            nd_arr.fill(np.nan)
            da_chi2_stacked = xr.DataArray(nd_arr, coords={'Te'  :da_model.Te,
                                                           'logg':da_model.logg,
                                                           'MH'  :da_model.MH,
                                                           'sf'  :self.sf_list_B})
            da_vgf2_stacked  = da_chi2_stacked.copy()
            da_vgfb2_stacked = da_chi2_stacked.copy()
        if self.model_to_fit_B=='Koester':
            nd_arr = np.empty((len(da_model.Te), len(da_model.logg), len(self.sf_list_B)))
            nd_arr.fill(np.nan)
            da_chi2_stacked = xr.DataArray(nd_arr, coords={'Te'  :da_model.Te,
                                                           'logg':da_model.logg,
                                                           'sf'  :self.sf_list_B})
            da_vgf2_stacked  = da_chi2_stacked.copy()
            da_vgfb2_stacked = da_chi2_stacked.copy()

        for idx in range (len(self.sf_list_B)):
            print_progress(idx, len(self.sf_list_B),10,  message='    Calculating chi2')
            sf          = self.sf_list_B[idx]
            da_residual = self.da_res_A - (da_model * sf)

            da_chi2_i   = (da_residual/self.da_obs_error)**2
            da_chi2     = da_chi2_i.sum('FilterID',min_count=1)
            da_chi2_stacked.loc[dict(sf=sf)] = da_chi2

            da_vgf2_i   = (da_residual/self.da_obs_error_2percent)**2
            da_vgf2     = da_vgf2_i.sum('FilterID',min_count=1)
            da_vgf2_stacked.loc[dict(sf=sf)] = da_vgf2

            da_vgfb2_i   = (da_residual/self.da_obs_error_10percent)**2
            da_vgfb2     = da_vgfb2_i.sum('FilterID',min_count=1)
            da_vgfb2_stacked.loc[dict(sf=sf)] = da_vgfb2

        df_chi2          = da_chi2_stacked.to_dataframe(name='chi2').reset_index()
        df_vgf2          = da_vgf2_stacked.to_dataframe(name='vgf2').reset_index()
        df_vgfb2         = da_vgfb2_stacked.to_dataframe(name='vgfb2').reset_index()

        df_chi2['R']     = calc_radius(df_chi2['sf'], self.distance)
        df_chi2['L']     = calc_luminosity(df_chi2['R'], df_chi2['Te'])
        df_chi2['vgf2']  = df_vgf2['vgf2'].values
        df_chi2['vgfb2'] = df_vgfb2['vgfb2'].values

        if minimising_param not in ['chi2', 'vgf2', 'vgfb2']:
            raise Exception('    Error: The minimising_param should be one of the following: chi2, vgf2, vgfb2')
        df_chi2 = df_chi2.sort_values(minimising_param)
        df_chi2.reset_index(drop=True, inplace=True)

        if verbose: 
            print('\n    Calculated chi2 for %d models. \n    Best 5 models:'%len(df_chi2))
            print(df_chi2.head())
        if trim: df_chi2 = df_chi2.head(5000)
        self.df_chi2 = df_chi2

        if not os.path.exists(DIR_OUTPUTS + 'chi_files/'): os.makedirs(DIR_OUTPUTS + 'chi_files/')
        if verbose: print('    Saving %s'%chi2_file_name)
        df_chi2.to_csv(chi2_file_name, index=False, header=True, sep=',')

    def get_parameters_from_chi2_minimization(self, verbose=False):
        '''
        Estimates best fit parameters from least chi2 fit
        
        - model_flux_B, Te_B, sf_B, R_B, L_B, chi2, vgf, vgfb
        - Updates these parameter to parent class object
        - Raises warning if any datapoints have significantly (3-sigma) higher chi2_i than average
        '''
        self.Te_B                 = self.df_chi2.Te[0]
        self.logg_B               = self.df_chi2.logg[0]

        if self.model_to_fit_B[:6]   =='Kurucz':
            self.MH_B                 = self.df_chi2.MH[0]
            self.alpha_B              = self.df_chi2.alpha[0]
            best_fit_flux             = self.da_model_B.sel(Te=self.Te_B).sel(logg=self.logg_B).sel(MH=self.MH_B).sel(alpha=self.alpha_B)

        if self.model_to_fit_B       =='Koester':
            self.MH_B                 = np.nan
            self.alpha_B              = np.nan
            best_fit_flux             = self.da_model_B.sel(Te=self.Te_B).sel(logg=self.logg_B)

        self.sf_B                 = self.df_chi2.sf[0]
        self.R_B                  = self.df_chi2.R[0]
        self.L_B                  = self.df_chi2.L[0]

        self.flux['model_flux_B'] = self.sf_B * best_fit_flux
        _residual_flux            = self.flux['flux'] - self.flux['model_flux_A'] - self.flux['model_flux_B']
        self.flux['chi2_i']       = _residual_flux**2 / self.flux['error']**2
        self.flux['vgf2_i']       = _residual_flux**2 / self.flux['error_2percent']**2
        self.flux['vgfb2_i']      = _residual_flux**2 / self.flux['error_10percent']**2
        self.chi2                 = self.flux['chi2_i'].sum()
        self.chi2_r               = self.chi2/self.N_Np
        self.vgf2                 = self.flux['vgf2_i'].sum()/self.N_Np
        self.vgfb2                = self.flux['vgfb2_i'].sum()/self.N_Np

        with xr.open_dataarray(self.model_file_name_B) as da_model_all:
            da_model_all = da_model_all.sel(FilterID=self.flux_all.index)
            if self.model_to_fit_B[:6]=='Kurucz':
                da_model_all = da_model_all.sel(Te=self.Te_B).sel(logg=self.logg_B).sel(MH=self.MH_B).sel(alpha=self.alpha_B)
            if self.model_to_fit_B=='Koester':
                da_model_all = da_model_all.sel(Te=self.Te_B).sel(logg=self.logg_B)
            df_model_all = da_model_all.to_dataframe(name='model_flux_B')
        self.flux_all['model_flux_B'] = df_model_all['model_flux_B']*self.sf_B    

        if self.model_to_fit_B[:6]=='Kurucz':
            if verbose: print('\n    Fitting parameters: T=%d,logg=%.2f,MH=%.2f,alpha=%.1f,sf=%.2e,R=%f,L=%f with chi2=%.2f'%(self.Te_B,self.logg_B, self.MH_B, self.alpha_B, self.sf_B, self.R_B, self.L_B,self.chi2))
        if self.model_to_fit_B       =='Koester':
            if verbose: print('\n    Fitting parameters: T=%d,logg=%.2f,sf=%.2e,R=%f,L=%f with chi2=%.2f'%(self.Te_B,self.logg_B, self.sf_B, self.R_B, self.L_B,self.chi2))

        # Printing the filtes with too large chi2 values i.e. 3sigma away from other chi2 values
        abs_zscore = np.abs(zscore(self.flux['chi2_i']))
        outliers = self.flux[abs_zscore>=3].index.values
        if len(outliers)>0: print('\n    Based on chi2, I recommend removal of following filters: ', outliers)
  
    def plot_skeletal_SED(self):
        '''
        Builds the skeletal plots (along with labels, limits etc) for binary SED.
        
        Returns
        -------
        f : Figure
        axes : axes.Axes or array of axes
        '''
        ###################### initialising
        f, axes = plt.subplots(figsize=(12,6),nrows = 3, ncols = 3)
        [axi.set_axis_off() for axi in axes.ravel()]

        axes[0,0] = f.add_axes([0.06, 0.44, 0.49, 0.50])
        axes[1,0] = f.add_axes([0.06, 0.27, 0.49, 0.17])
        axes[2,0] = f.add_axes([0.06, 0.10, 0.49, 0.17])

        axes[0,1] = f.add_axes([0.63, 0.66, 0.30, 0.28])
        axes[1,1] = f.add_axes([0.63, 0.38, 0.30, 0.28])
        axes[2,1] = f.add_axes([0.63, 0.10, 0.30, 0.28])

        axes[0,2] = f.add_axes([0.91, 0.10, 0.02, 0.28])
        ####################### SED
        axes[0,0].scatter(self.not_fitted['wavelength'], self.not_fitted['flux'], color='k', marker='o',label ='No Fit', s=30, facecolors='none',zorder=2)
        axes[0,0].errorbar(self.flux['wavelength'], self.flux['flux'], yerr=self.flux['error'],color='k', label='Obs',fmt='none',lw=2, capsize=4,zorder=3)
        ####################### residual
        axes[1,0].errorbar(self.flux['wavelength'], self.flux['flux']-self.flux['flux'], yerr=self.flux['error_fraction'],color='k', label='Obs',fmt='none',lw=2, capsize=4) 
        axes[1,0].errorbar(self.not_fitted['wavelength'], self.not_fitted['flux']-self.not_fitted['flux'], yerr=self.not_fitted['error_fraction'],color='0.5', label='',fmt='none',lw=2, capsize=4) 
        ####################### isochrones and WD cooling curves
        iso        = pd.read_csv(DIR_MODELS + 'master_isochrone.csv')
        iso_8      = iso[iso.logAge==8]
        iso_9      = iso[iso.logAge==9]
        iso_10     = iso[iso.logAge==10]
        Bergeron_WD= pd.read_csv(DIR_MODELS + 'master_Bergeron_WD.csv')
        WD_02      = Bergeron_WD[(Bergeron_WD.mass==0.2) & (Bergeron_WD.spectral_type=='DA')]
        WD_13      = Bergeron_WD[(Bergeron_WD.mass==1.3) & (Bergeron_WD.spectral_type=='DA')]

        axes[0,1].plot(10**(iso_8.logTe),10**(iso_8.logL), label='',c='0',lw=0.5, rasterized = True, zorder=1)
        axes[0,1].plot(10**(iso_9.logTe),10**(iso_9.logL), label='',c='0',lw=0.5, rasterized = True, zorder=1)
        axes[0,1].plot(10**(iso_10.logTe),10**(iso_10.logL), label='',c='0',lw=0.5, rasterized = True, zorder=1)
        axes[0,1].plot(WD_02.Teff,10**WD_02.logL, label='',c='0', ls=(0,(5,10)),lw=0.5, rasterized = True, zorder=1)
        axes[0,1].plot(WD_13.Teff,10**WD_13.logL, label='',c='0', ls=(0,(5,10)),lw=0.5, rasterized = True, zorder=1)
        ####################### Labels
        axes[0,0].set_ylabel('Flux (erg s$^{-1}$ cm$^{-2}$ $\AA$$^{-1}$)')
        axes[1,0].set_ylabel('Fractional\nResidual')
        axes[2,0].set_ylabel('$\chi^2_i$')
        axes[2,0].set_xlabel('Wavelength ($\AA$)')
        axes[0,1].set_ylabel('L ($L_{\odot}$)')
        axes[1,1].set_ylabel('R ($R_{\odot}$)')
        axes[2,1].set_ylabel('$\chi^2_r$')
        axes[2,1].set_xlabel('Te (K)')
        ####################### axes range and scales
        axes[0,0].set_xscale('log')
        axes[0,0].set_yscale('log')
        axes[1,0].set_xscale('log')
        axes[2,0].set_xscale('log')
        axes[0,1].set_yscale('log')
        axes[1,1].set_yscale('log')
        axes[0,1].set_xscale('log')
        axes[1,1].set_xscale('log')
        axes[2,1].set_xscale('log')
        axes[2,1].set_yscale('log')

        wave_min = self.flux_all['wavelength'].min()
        wave_max = self.flux_all['wavelength'].max()

        axes[0,0].set_xlim([wave_min/1.2,wave_max*1.2])
        axes[1,0].set_xlim([wave_min/1.2,wave_max*1.2])
        axes[2,0].set_xlim([wave_min/1.2,wave_max*1.2])

        T_min = min(self.df_chi2.Te.min(), 4500)
        T_max = max(self.df_chi2.Te.max(), 40000)

        for ax in [axes[0,1],axes[1,1],axes[2,1]]:
            ax.set_xlim(T_max*1.1,T_min/1.4)

        flux_min = self.flux_all['flux'].min()
        flux_max = self.flux_all['flux'].max()
        axes[0,0].set_ylim([flux_min/100,flux_max*5])

        L_min = np.min([self.L_A, self.L_B])
        L_max = np.max([self.L_A, self.L_B])
        axes[0,1].set_ylim(L_min/5e2,L_max*5e2)

        plt.setp(axes[0,0].get_xticklabels(),visible=False)
        plt.setp(axes[1,0].get_xticklabels(),visible=False)
        plt.setp(axes[0,1].get_xticklabels(),visible=False)
        plt.setp(axes[1,1].get_xticklabels(),visible=False)

        ####################### decoration    
        for i,j in product(range(3),range(2)):
            axes[i,j].tick_params(which='both', direction='out', length=4)
            axes[i,j].grid(which='both', axis='x', zorder=0)
            axes[i,j].grid(which='major', axis='y', zorder=0)
        
        axes[0,2].tick_params(which='both', direction='out', length=4)

        return f, axes
    
    def plot_fitted_SED(self,cycle,save_plot=True,show_plot=True,folder_path='plots/binary_SEDs/'):
        '''
        Plotting binary SED
        
        - Creates SED plot after observed flux is fitted with B component.
        - Includes HR diagram with isochrone of logAge 8/9/10, WD models for mass 0.2/1.3.
        - Includes R_B vs Te_B plot
        - Includes chi2 vs Te_B plot. This plot should have explicit global minima in a good fit.
        '''
        f, axes = self.plot_skeletal_SED()
 
        axes[0,0].plot(self.flux_all['wavelength'], self.flux_all['model_flux_A'], color='r', linestyle='--',label ='A', lw=1)
        axes[0,0].plot(self.flux_all['wavelength'], self.flux_all['model_flux_B'], color='b', linestyle='--',label ='B', lw=1)
        axes[0,0].plot(self.flux_all['wavelength'], self.flux_all['model_flux_A']+self.flux_all['model_flux_B'], color='g', linestyle='-',label ='Total', lw=1)
        ########## Fractional residual
        _fractional_residual_A = (self.flux['flux']-self.flux['model_flux_A'])/self.flux['flux']
        _fractional_residual_B = (self.flux['flux']-self.flux['model_flux_A']-self.flux['model_flux_B'])/self.flux['flux']
        axes[1,0].plot(self.flux['wavelength'], _fractional_residual_A,label='',marker='',color='r',lw=1, linestyle='--')
        axes[1,0].plot(self.flux['wavelength'], _fractional_residual_B,label='',marker='',color='green',lw=1, linestyle='-')
        ########## chi2_i
        axes[2,0].plot(self.flux['wavelength'], self.flux['chi2_i'],label='',marker='o',color='green', linestyle='-',lw=1)
        ########## L vs T (HR diagram)
        if len(self.df_chi2)>10_000: 
            top_chi2        = np.nanpercentile(self.df_chi2.chi2, 1)
        else: 
            top_chi2        = self.df_chi2.chi2.max()
        top_df_chi2     = self.df_chi2[self.df_chi2.chi2<=top_chi2]
        Te_data_top     = top_df_chi2.Te
        L_data_top      = top_df_chi2.L
        vgfb2_data_top  = top_df_chi2.vgfb2
        R_data_top      = top_df_chi2.R
        chi2_r_data_top = top_df_chi2.chi2/self.N_Np

        axes[0,1].scatter(Te_data_top,L_data_top, marker='.',c=vgfb2_data_top ,label='',s=5,rasterized = True, zorder=2, cmap='summer')
        axes[0,1].scatter(self.Te_A,self.L_A, marker='s', label='A',c='r', s=40,rasterized = True, zorder=3)
        axes[0,1].scatter(self.Te_B,self.L_B, marker='s', label='B',c='b', s=40,rasterized = True, zorder=3)
        ########## R vs T    
        cs1 = axes[1,1].scatter(Te_data_top,R_data_top,c=vgfb2_data_top, cmap='summer',s=10,rasterized = True, zorder=2, norm=matplotlib.colors.LogNorm())
        axes[1,1].scatter(self.Te_B,self.R_B, marker='s', label='',c='b', s=40,rasterized = True, zorder=3)
        ########## chi2 vs T
        cs2 = axes[2,1].scatter(Te_data_top,chi2_r_data_top,c=vgfb2_data_top, cmap='summer',s=10,rasterized = True,zorder=2, norm=matplotlib.colors.LogNorm())  
        axes[2,1].scatter(self.Te_B,self.chi2_r, marker='s', label='',c='b', s=40,rasterized = True, zorder=3)
        ########## colorbar
        f.colorbar(cs2, cax=axes[0,2])
        cs1.set_clim(vgfb2_data_top.min(),vgfb2_data_top.max())
        cs2.set_clim(vgfb2_data_top.min(),vgfb2_data_top.max())
        axes[0,2].set_ylabel('$vgf_b^2$')
        axes[0,2].yaxis.set_label_position("right")
        ####################### Titles and labels
        if self.model_to_fit_A[:6]=='Kurucz':
            label_A = '%d $K$, '%(self.Te_A) +\
                    'logg=%.2f, MH=%.2f'%(self.logg_A, self.MH_A)
        if self.model_to_fit_A=='Koester':
            label_A = '%d $K$, '%(self.Te_A) +\
                    'logg=%.2f'%(self.logg_A)

        if self.model_to_fit_B[:6]=='Kurucz':
            label_B = '%d $K$, %.3f $R_{\odot}$, %.4f $L_{\odot}$, logg=%.2f, MH=%.2f, alpha=%.1f'%(self.Te_B, self.R_B, self.L_B, self.logg_B, self.MH_B, self.alpha_B)
        if self.model_to_fit_B=='Koester':
            label_B = '%d $K$, %.3f $R_{\odot}$, %.4f $L_{\odot}$, logg=%.2f'%(self.Te_B, self.R_B, self.L_B, self.logg_B)
        axes[0,0].set_title(self.name+'    A (' + self.model_to_fit_A + ': ' + label_A +
                            ')    B (' + self.model_to_fit_B + ': ' + label_B + ')', x=0, y=1, ha='left')
        axes[2,0].set_title('$\chi^2$ = %.2f\n$\chi_r^2$ = %.2f\n$vgf_b^2$ = %.2f'%(self.chi2, self.chi2_r, self.vgfb2),
                             x=0.98,y=0.9, ha='right', va='top')
        ####################### decoration and saving
        axes[1,1].set_ylim(R_data_top.min()/1.1,R_data_top.max()*1.1)

        axes[0,0].legend(scatterpoints=1, loc='upper center', ncol=5,frameon=False,handletextpad=0.3, borderpad=0.1)
        axes[0,1].legend(scatterpoints=1)
        if save_plot:   
            if not os.path.exists(folder_path):
                os.makedirs(folder_path)
            plt.savefig (folder_path +self.name+'_' + self.model_to_fit_B +'_'+str(self.Te_B)+'_'+str(cycle)+'.jpg', format='jpg', dpi=300)
            plt.savefig (folder_path +self.name+'_' + self.model_to_fit_B +'_'+str(self.Te_B)+'_'+str(cycle)+'.pdf', format='pdf', dpi=300)
        if not show_plot:
            plt.close() 

    def calculate_chi2_noisy(self, cycle, total_iterations=100, verbose=False, refit=True, percentile_threshold=5, plot=False, minimising_param='chi2'):
        '''
        Calculate chi2 after adding noise to the observed data
        
        - Takes ``da_model_B`` and scales it for each sf in ``sf_list_B_noisy``.
        - Creates ``total_iterations`` versions of ``da_res_A`` as ``da_obs_noisy``.
        - Subtracts the ``da_model_B`` from ``da_obs_noisy`` and calculates corresponding chi2_i for each datapoint.
        - Sums the chi2_i for each model and saves the chi2 values into DataArray (``da_chi2_stacked``). 
        - Simultaneously calculates vgfb2
        - Converts ``da_chi2_stacked`` into a dataframe (``df_chi2_noisy``)
        - Saves ``df_chi2_noisy`` into a .csv file
        
        Returns
        -------
        df_chi2_noisy : DataFrame
            chi2 sorted DataFrame with Te [K], logg, (MH), sf.            
        '''
        top_chi2           = np.nanpercentile(self.df_chi2.chi2, percentile_threshold)
        top_df_chi2        = self.df_chi2[self.df_chi2.chi2<top_chi2]
        sf_min             = top_df_chi2.sf.min()
        sf_max             = top_df_chi2.sf.max()
        self.sf_list_B_noisy = 10**np.arange(np.log10(sf_min),np.log10(sf_max), 0.01) 

        chi2_file_name = DIR_OUTPUTS + 'chi_files/'+self.type +'_'+self.name+'_'+self.model_to_fit_B+'_noisy_'+str(cycle)+'.csv'

        # If the chi2 file exists, it will be read (depends on refit = True or False)
        if os.path.isfile(chi2_file_name): 
            if not refit:
                df_chi2_noisy      = pd.read_csv(chi2_file_name)
                if verbose:     
                    print('    Reading %s'%chi2_file_name)
                    print('    WARNING! give "refit=True" if you want to rerun the fitting process.')
                    print(df_chi2_noisy.head())
                if plot: 
                    size_arr    = 20 * df_chi2_noisy['logg']**0
                    size_arr[0] = 100
                    sns.pairplot(df_chi2_noisy, vars=['Te','logg','sf','chi2'],corner=True, diag_kind="hist", plot_kws={"s": size_arr})
                self.df_chi2_noisy = df_chi2_noisy
                return 
            if refit: print('    WARNING: '+chi2_file_name+' file will be overwritten.')

        da_model = self.da_model_B
        if self.model_to_fit_B[:6]=='Kurucz':
            nd_arr = np.empty((len(da_model.Te), len(da_model.logg), len(da_model.MH), len(da_model.alpha), len(self.sf_list_B_noisy)))
            nd_arr.fill(np.nan)
            da_chi2_stacked = xr.DataArray(nd_arr, coords={'Te'   :da_model.Te,
                                                           'logg' :da_model.logg,
                                                           'MH'   :da_model.MH,
                                                           'alpha':da_model.alpha,
                                                           'sf'   :self.sf_list_B_noisy})
        if self.model_to_fit_B=='Koester':
            nd_arr = np.empty((len(da_model.Te), len(da_model.logg), len(self.sf_list_B_noisy)))
            nd_arr.fill(np.nan)
            da_chi2_stacked = xr.DataArray(nd_arr, coords={'Te'  :da_model.Te,
                                                           'logg':da_model.logg,
                                                           'sf'  :self.sf_list_B_noisy})

        da_obs_noisy                    = self.da_model_B.copy()
        da_obs_noisy.attrs['long_name'] = 'Noisy flux'

        df_chi2_noisy = pd.DataFrame()

        for seed in range (total_iterations):
            print_progress(seed, total_iterations=total_iterations, step=2, message='    Calculating noisy chi2')

            np.random.seed(seed)
            noise = np.random.normal(0, 1, [self.N_points]) 
            if seed == 0:    # No noise added for 0th iteration
                noisy_flux = self.flux['flux'] - self.flux['model_flux_A']
            else:            # Adding Gaussian noise to the flux 
                noisy_flux = self.flux['flux'] + self.flux['error']*noise - self.flux['model_flux_A']

            for filter_name in self.flux.index:
                da_obs_noisy = da_obs_noisy.where(da_obs_noisy.FilterID!=filter_name,noisy_flux[filter_name])

            for idx in range (len(self.sf_list_B_noisy)):
                sf = self.sf_list_B_noisy[idx]
                da_residual = da_obs_noisy - (da_model * sf)

                da_chi2_i   = (da_residual/self.da_obs_error)**2
                da_chi2     = da_chi2_i.sum('FilterID',min_count=1)
                da_chi2_stacked.loc[dict(sf=sf)] = da_chi2
                
                if minimising_param != 'chi2':
                    da_vgf2_i   = (da_residual/self.da_obs_error_2percent)**2
                    da_vgf2     = da_vgf2_i.sum('FilterID',min_count=1)
                    da_vgf2_stacked.loc[dict(sf=sf)] = da_vgf2

                    da_vgfb2_i  = (da_residual/self.da_obs_error_10percent)**2
                    da_vgfb2    = da_vgfb2_i.sum('FilterID',min_count=1)
                    da_vgfb2_stacked.loc[dict(sf=sf)] = da_vgfb2

            df_chi2         = da_chi2_stacked.to_dataframe(name='chi2').reset_index()
            if minimising_param not in ['chi2', 'vgf2', 'vgfb2']:
                raise Exception('    Error: The minimising_param should be one of the following: chi2, vgf2, vgfb2')
            df_chi2 = df_chi2.sort_values(minimising_param)
            df_chi2.reset_index(drop=True, inplace=True)

            df_chi2_noisy = pd.concat([df_chi2_noisy, df_chi2.head(1)], ignore_index = True)

        if verbose: print('\n',df_chi2_noisy.head())
        if plot: 
            size_arr    = 20 * df_chi2_noisy['logg']**0
            size_arr[0] = 100
            sns.pairplot(df_chi2_noisy, vars=['Te','logg','sf','chi2'],corner=True, diag_kind="hist",plot_kws={"s": size_arr})

        self.df_chi2_noisy = df_chi2_noisy

        if not os.path.exists(DIR_OUTPUTS + 'chi_files/'): os.makedirs(DIR_OUTPUTS + 'chi_files/')
        if verbose: print('    Saving %s'%chi2_file_name)
        df_chi2_noisy.to_csv(chi2_file_name, index=False, header=True, sep=',')

    def get_parameters_from_noisy_chi2(self, verbose=False):
        '''
        Estimates best fit parameters from noisy least chi2 fit
        
        - model_flux_B, Te_B, sf_B, R_B, L_B, chi2, vgf, vgfb
        - Errors in Te_B, sf_B, logg_B, MH_B, R_B, L_B based on their noisy distribution
        - Updates these parameter to parent class object
        '''
        self.Te_median_B, self.Te_error_lower_B, self.Te_error_upper_B       = get_realistic_errors_from_iterations(self.df_chi2_noisy['Te'],self.da_model_B['Te'].values,para_type='Te')
        self.logg_median_B, self.logg_error_lower_B, self.logg_error_upper_B = get_realistic_errors_from_iterations(self.df_chi2_noisy['logg'],self.da_model_B['logg'].values,para_type='logg')
        self.sf_median_B, self.sf_error_lower_B, self.sf_error_upper_B       = get_realistic_errors_from_iterations(self.df_chi2_noisy['sf'],self.sf_list_B_noisy,para_type='sf')

        self.R_median_B      = calc_radius(self.sf_median_B, self.distance)
        # error contribution from noisy fitting
        self.R_error_lower_B = self.R_median_B - calc_radius(self.sf_median_B-self.sf_error_lower_B, self.distance)
        self.R_error_upper_B = calc_radius(self.sf_median_B+self.sf_error_upper_B, self.distance) - self.R_median_B
        # error after including distance error,     deltaR = R * deltaD/D
        self.R_error_lower_B = addition_in_quadrature(self.R_error_lower_B, self.R_median_B * self.distance_err/self.distance)
        self.R_error_upper_B = addition_in_quadrature(self.R_error_upper_B, self.R_median_B * self.distance_err/self.distance)

        self.L_median_B      = calc_luminosity(self.R_median_B, self.Te_median_B)
        # error contribution from noisy fitting
        self.L_error_lower_B = self.L_median_B - calc_luminosity(self.R_median_B-self.R_error_lower_B, self.Te_median_B-self.Te_error_lower_B)
        self.L_error_upper_B = calc_luminosity(self.R_median_B+self.R_error_upper_B, self.Te_median_B+self.Te_error_upper_B) - self.L_median_B
        # error after including distance error,     deltaL = L * (2 * deltaD/D)
        self.L_error_lower_B = addition_in_quadrature(self.L_error_lower_B, self.L_median_B * 2. * self.distance_err/self.distance)
        self.L_error_upper_B = addition_in_quadrature(self.L_error_upper_B, self.L_median_B * 2. * self.distance_err/self.distance)

        if self.Te_B!=self.Te_median_B:
            print ('    WARNING! Median of noisy "Te" is not same as "Te".')

        if self.model_to_fit_B[:6]  =='Kurucz':
            self.MH_median_B, self.MH_error_lower_B, self.MH_error_upper_B = get_realistic_errors_from_iterations(self.df_chi2_noisy['MH'],self.da_model['MH'].values,para_type='MH')
            self.alpha_median_B, self.alpha_error_lower_B, self.alpha_error_upper_B = get_realistic_errors_from_iterations(self.df_chi2_noisy['alpha'],self.da_model['alpha'].values,para_type='alpha')
            best_fit_flux_median_B         = self.da_model_B.sel(Te=self.Te_median_B).sel(logg=self.logg_median_B).sel(MH=self.MH_median_B).sel(alpha=self.alpha_median_B)

        if self.model_to_fit_B       =='Koester':
            self.MH_median_B, self.MH_error_lower_B, self.MH_error_upper_B  = np.nan,np.nan,np.nan
            self.alpha_median_B, self.alpha_error_lower_B, self.alpha_error_upper_B  = np.nan,np.nan,np.nan
            best_fit_flux_median_B           = self.da_model_B.sel(Te=self.Te_median_B).sel(logg=self.logg_median_B)

        self.flux['model_flux_median_B'] = self.sf_median_B * best_fit_flux_median_B
        _residual_flux                 = self.flux['flux'] - self.flux['model_flux_A'] - self.flux['model_flux_median_B']
        self.flux['chi2_i_median']     = _residual_flux**2 / self.flux['error']**2
        self.flux['vgf2_i_median']     = _residual_flux**2 / self.flux['error_2percent']**2
        self.flux['vgfb2_i_median']    = _residual_flux**2 / self.flux['error_10percent']**2
        self.chi2_median               = self.flux['chi2_i_median'].sum()
        self.chi2_r_median             = self.chi2_median/self.N_Np
        self.vgf2_median               = self.flux['vgf2_i_median'].sum()/self.N_Np
        self.vgfb2_median              = self.flux['vgfb2_i_median'].sum()/self.N_Np

        with xr.open_dataarray(self.model_file_name_B) as da_model_all:
            da_model_all = da_model_all.sel(FilterID=self.flux_all.index)
            if self.model_to_fit_B[:6]=='Kurucz':
                da_model_all = da_model_all.sel(Te=self.Te_median_B).sel(logg=self.logg_median_B).sel(MH=self.MH_median_B).sel(alpha=self.alpha_median_B)
            if self.model_to_fit_B=='Koester':
                da_model_all = da_model_all.sel(Te=self.Te_median_B).sel(logg=self.logg_median_B)
            df_model_all = da_model_all.to_dataframe(name='model_flux_median_B')
        self.flux_all['model_flux_median_B'] = df_model_all['model_flux_median_B']*self.sf_median_B

        if verbose: 
            print('\n    Fitting parameters of from noisy iterations:')
            print('        Te   (A) %d (+%d-%d)  |  (B) %d (+%d-%d)'%(self.Te_A, self.Te_error_lower_A, self.Te_error_upper_A,self.Te_median_B, self.Te_error_lower_B, self.Te_error_upper_B))
            print('        logg (A) %.2f (+%.2f-%.2f)  |  (B) %.2f (+%.2f-%.2f)'%(self.logg_A, self.logg_error_lower_A, self.logg_error_upper_A,self.logg_median_B, self.logg_error_lower_B, self.logg_error_upper_B))
            print('        MH   (A) %.2f (+%.2f-%.2f)  |  (B) %.2f (+%.2f-%.2f)'%(self.MH_A, self.MH_error_lower_A, self.MH_error_upper_A,self.MH_median_B, self.MH_error_lower_B, self.MH_error_upper_B))
            print('        alpha(A) %.2f (+%.2f-%.2f)  |  (B) %.2f (+%.2f-%.2f)'%(self.alpha_A, self.alpha_error_lower_A, self.alpha_error_upper_A,self.alpha_median_B, self.alpha_error_lower_B, self.alpha_error_upper_B))
            print('        sf   (A) %.2e (+%.2e-%.2e)  |  (B) %.2e (+%.2e-%.2e)'%(self.sf_A, self.sf_error_lower_A, self.sf_error_upper_A,self.sf_median_B, self.sf_error_lower_B, self.sf_error_upper_B))
            print('        R    (A) %.4f (+%.4f-%.4f)  |  (B) %.4f (+%.4f-%.4f)'%(self.R_A, self.R_error_lower_A, self.R_error_upper_A,self.R_median_B, self.R_error_lower_B, self.R_error_upper_B))
            print('        L    (A) %.4f (+%.4f-%.4f)  |  (B) %.4f (+%.4f-%.4f)'%(self.L_A, self.L_error_lower_A, self.L_error_upper_A,self.L_median_B, self.L_error_lower_B, self.L_error_upper_B))
            print('        chi2 = %.1f, vgf2 = %.1f, vgfb2 = %.1f'%(self.chi2_median, self.vgf2_median, self.vgfb2_median))

    def plot_fitted_noisy_SED(self,cycle,plot_noisy_SEDs=True, save_plot=True,show_plot=True, folder_path='plots/binary_SEDs/'):
        '''
        Creates noisy SED

        - Includes HR diagram with isochrone of logAge 8/9/10, WD models for mass 0.2/1.3.
        - Includes radius vs Te plot
        - Includes chi2 vs Te plot. This plot should have explicit global minima in a good fit.
        '''
        f, axes = self.plot_skeletal_SED()

        axes[0,0].plot(self.flux_all['wavelength'], self.flux_all['model_flux_A'], color='r', linestyle='--',label ='A', lw=1)
        axes[0,0].plot(self.flux_all['wavelength'], self.flux_all['model_flux_B'], color='g', linestyle='--',label ='B', lw=1)
        axes[0,0].plot(self.flux_all['wavelength'], self.flux_all['model_flux_A']+self.flux_all['model_flux_B'], color='green', linestyle='-',label ='Total', lw=1)

        ########## plotting noisy fits
        if plot_noisy_SEDs:
            for idx, sf_noisy in enumerate(self.df_chi2_noisy['sf']):
                print_progress(idx, total_iterations=100, step=5, message='    Plotting noisy SEDs')
                Te_noisy   = self.df_chi2_noisy['Te'][idx]
                logg_noisy   = self.df_chi2_noisy['logg'][idx]
                if self.model_to_fit_B[:6]  =='Kurucz':
                    MH_noisy     = self.df_chi2_noisy['MH'][idx]
                    alpha_noisy  = self.df_chi2_noisy['alpha'][idx]
                    flux_B_noisy = self.da_model_B.sel(Te=Te_noisy).sel(logg=logg_noisy).sel(MH=MH_noisy).sel(alpha=alpha_noisy).values * sf_noisy 
                if self.model_to_fit_B  =='Koester':
                    flux_B_noisy = self.da_model_B.sel(Te=Te_noisy).sel(logg=logg_noisy).values * sf_noisy
                axes[0,0].plot(self.flux['wavelength'], flux_B_noisy, color='cyan', linestyle='-',label='', lw=0.2, alpha=0.2,zorder=0)
                residual_noisy = (self.flux['flux'] - flux_B_noisy - self.flux['model_flux_A'])/self.flux['flux']
                axes[1,0].plot(self.flux['wavelength'], residual_noisy, color='cyan', linestyle='-',label='', lw=0.2, alpha=0.2,zorder=0)   
        ########## Fractional residual
        _fractional_residual_A = (self.flux['flux']-self.flux['model_flux_A'])/self.flux['flux']
        _fractional_residual_B = (self.flux['flux']-self.flux['model_flux_A']-self.flux['model_flux_B'])/self.flux['flux']
        axes[1,0].plot(self.flux['wavelength'], _fractional_residual_A,label='',marker='',color='r',lw=1, linestyle='--')
        axes[1,0].plot(self.flux['wavelength'], _fractional_residual_B,label='',marker='',color='green',lw=1, linestyle='-')
        ########## chi2_i
        axes[2,0].plot(self.flux['wavelength'], self.flux['chi2_i'],label='',marker='o',color='green', linestyle='-',lw=1)
        ########## L vs T (HR diagram)
        self.df_chi2_noisy['R']    = calc_radius(self.df_chi2_noisy['sf'], self.distance)
        self.df_chi2_noisy['L']    = calc_luminosity(self.df_chi2_noisy['R'], self.df_chi2_noisy['Te'])
        Te_data     = self.df_chi2_noisy.Te
        L_data      = self.df_chi2_noisy.L
        R_data      = self.df_chi2_noisy.R
        chi2_r_data = self.df_chi2_noisy.chi2/self.N_Np
        axes[0,1].scatter(Te_data,L_data, c=chi2_r_data, cmap='summer', marker='.',label='',s=10,rasterized = True, zorder=2)
        axes[0,1].scatter(self.Te_A,self.L_A, marker='s', label='A',c='r', s=40,rasterized = True, zorder=3)
        axes[0,1].scatter(self.Te_B,self.L_B, marker='s', label='B',c='b', s=40,rasterized = True, zorder=3)
        ########## R vs T
        axes[1,1].scatter(Te_data,R_data, c=chi2_r_data, cmap='summer', s=10,rasterized = True, zorder=2)
        axes[1,1].scatter(self.Te_B,self.R_B, marker='s', label='',c='b', s=40,rasterized = True, zorder=3)
        ########## chi2 vs T
        cs2 = axes[2,1].scatter(Te_data,chi2_r_data,c=chi2_r_data, cmap='summer',s=10,rasterized = True,zorder=2)  
        axes[2,1].scatter(self.Te_B,self.chi2_r, marker='s', label='',c='b', s=40,rasterized = True, zorder=3)
        ########## colorbar
        f.colorbar(cs2, cax=axes[0,2])
        axes[0,2].set_ylabel('$\chi^2_r$')
        axes[0,2].yaxis.set_label_position("right")
        ####################### errorbars
        axes[0,1].axvspan(self.Te_A-self.Te_error_lower_A, self.Te_A+self.Te_error_upper_A, alpha=0.2, color='hotpink', zorder=1)
        axes[0,1].axhspan(self.L_A-self.L_error_lower_A, self.L_A+self.L_error_upper_A, alpha=0.2, color='hotpink', zorder=1)        

        axes[0,1].axvspan(self.Te_B-self.Te_error_lower_B, self.Te_B+self.Te_error_upper_B, alpha=0.2, color='dodgerblue', zorder=1)
        axes[0,1].axhspan(self.L_B-self.L_error_lower_B, self.L_B+self.L_error_upper_B, alpha=0.2, color='dodgerblue', zorder=1)   

        axes[1,1].axvspan(self.Te_B-self.Te_error_lower_B, self.Te_B+self.Te_error_upper_B, alpha=0.2, color='dodgerblue', zorder=1)
        axes[1,1].axhspan(self.R_B-self.R_error_lower_B, self.R_B+self.R_error_upper_B, alpha=0.2, color='dodgerblue', zorder=1)    

        axes[2,1].axvspan(self.Te_B-self.Te_error_lower_B, self.Te_B+self.Te_error_upper_B, alpha=0.2, color='dodgerblue', zorder=1)
        ####################### Titles and labels
        if self.model_to_fit_A[:6]=='Kurucz':
            label_A = '%d $K$, '%(self.Te_A) +\
                    'logg=%.2f, MH=%.2f'%(self.logg_A, self.MH_A)
        if self.model_to_fit_A=='Koester':
            label_A = '%d $K$, '%(self.Te_A) +\
                    'logg=%.2f'%(self.logg_A)

        if self.model_to_fit_B[:6]=='Kurucz':
            label_B = '%d$^{+%d}_{-%d}$ $K$, '%(self.Te_B, self.Te_error_lower_B, self.Te_error_upper_B) +\
                      '%.4f$^{+%.4f}_{-%.4f}$ $R_{\odot}$, '%(self.R_B, self.R_error_lower_B, self.R_error_upper_B) +\
                      '%.4f$^{+%.4f}_{-%.4f}$ $L_{\odot}$, '%(self.L_B, self.L_error_lower_B, self.L_error_upper_B) +\
                      'logg=%.2f, MH=%.2f, alpha=%.1f'%(self.logg_B, self.MH_B, self.alpha_B)
        if self.model_to_fit_B=='Koester':
            label_B = '%d$^{+%d}_{-%d}$ $K$, '%(self.Te_B, self.Te_error_lower_B, self.Te_error_upper_B) +\
                      '%.4f$^{+%.4f}_{-%.4f}$ $R_{\odot}$, '%(self.R_B, self.R_error_lower_B, self.R_error_upper_B) +\
                      '%.4f$^{+%.4f}_{-%.4f}$ $L_{\odot}$, '%(self.L_B, self.L_error_lower_B, self.L_error_upper_B) +\
                    '  logg=%.2f'%(self.logg_B)
        # axes[0,0].set_title(self.name+'    A (' + self.model_to_fit_A + ': ' + label_A + 
        #                     ')    B (' + self.model_to_fit_B + ': ' + label_B + ')', 
        #                     x=0, y=1, ha='left')
        axes[0,0].set_title('%s    A (%s: %s)    B (%s: %s)'%(self.name, self.model_to_fit_A, label_A, self.model_to_fit_B, label_B), 
                            x=0, y=1, ha='left')
        axes[2,0].set_title('$\chi^2$ = %.2f\n$\chi_r^2$ = %.2f\n$vgf_b^2$ = %.2f'%(self.chi2, self.chi2_r,self.vgfb2),
                            x=0.98,y=0.9, ha='right', va='top')
        ####################### decoration and saving
        axes[1,1].set_ylim(self.R_B-5*self.R_error_lower_B, self.R_B+5*self.R_error_upper_B)

        axes[1,1].set_yscale('linear')
        axes[2,1].set_yscale('linear')

        axes[0,0].legend(scatterpoints=1, loc='upper center', ncol=5,frameon=False,handletextpad=0.3, borderpad=0.1)
        axes[0,1].legend(scatterpoints=1)
        if save_plot:   
            if not os.path.exists(folder_path):
                os.makedirs(folder_path)
            plt.savefig (folder_path + self.name+'_' + self.model_to_fit_B +'_'+str(self.Te_median_B)+'_'+str(cycle)+'_noisy.jpg', format='jpg', dpi=300)#,bbox_inches='tight')
            # plt.savefig (folder_path + self.name+'_' + self.model_to_fit_B +'_'+str(self.Te_median_B)+'_'+str(cycle)+'_noisy.pdf', format='pdf', dpi=300)#,bbox_inches='tight')
        if not show_plot:
                plt.close()         
        
    def save_object(self, cycle, save_log=True):
        '''
        Saving class object and log
        
        - Saves the SingleStar object as a pickel.
        - Saves class attributes in a .csv log file         
        '''
        _self_copy = copy.deepcopy(self)
        _self_copy.time  = strftime("%Y-%m-%d %H:%M:%S", gmtime())
        _self_copy.cycle = cycle
        try:
            for attrName in ['da_model_B', 'da_res_A', 'da_obs_error', 'da_obs_error_2percent', 'da_obs_error_10percent', 'df_chi2', 'df_chi2_noisy']:
                delattr(_self_copy, attrName)
        except: pass

        file_name = DIR_OUTPUTS + 'pickels/BinaryStar_%s_%d.pkl'%(self.name,cycle)
        if not os.path.exists(DIR_OUTPUTS + 'pickels/'): os.makedirs(DIR_OUTPUTS + 'pickels/')
        with open(file_name, 'wb') as outp:
            pickle.dump(_self_copy, outp, pickle.HIGHEST_PROTOCOL)

        if save_log: 
            log_file_name              = 'log_binary_fitting.csv'

            # remove arrays and dataframes from the log file
            try:
                for attrName in ['flux', 'flux_all','not_fitted', 'sf_list_B', 'sf_list_B_noisy']:
                    delattr(_self_copy, attrName)
            except: pass
            _self_copy.filters_to_drop = ' '.join(_self_copy.__dict__['filters_to_drop'])

            df_log = pd.DataFrame(_self_copy.__dict__,index=[0])

            if not os.path.isfile(log_file_name):
                print('    Creating %s and saving log'%log_file_name)
            else:
                print('    Saving log in %s'%log_file_name)
                df_past_log = pd.read_csv(log_file_name)
                df_log = pd.concat([df_past_log,df_log])

            df_log.to_csv(log_file_name, index=False)
