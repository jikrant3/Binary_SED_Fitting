import copy
import logging
import pathlib
import os
import sys
import time

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import zscore
from astropy import constants as const
from astropy import units as u
import xarray as xr

DIR_OUTPUTS = 'outputs/'
DIR_MODELS = '/Users/vikrantjadhav/Documents/work/models_and_tools/models/'
if not os.path.exists(DIR_MODELS):
    raise FileNotFoundError('Rename the DIR_MODELS to your model folder.')

################################################################################
# set up logging to file
if 'logger' not in locals():
    if not os.path.exists('data'):
        os.makedirs('data')
    logging.basicConfig(level=logging.DEBUG,
                        format='\n%(asctime)s %(name)-12s %(levelname)-8s %(funcName)s: \n%(message)s',
                        filename='data/log_%s.txt' % (
                            time.strftime("%Y-%m-%d", time.gmtime())),
                        filemode='a')
    # # defined a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler(stream=sys.stdout)
    # console.terminator = ''
    console.setLevel(logging.WARNING)
    formatter = logging.Formatter(
        '%(asctime)s ----- %(levelname)-8s ----- %(funcName)s\n%(message)s', datefmt='%H:%M:%S')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)
    # defined logger to save all logs
    logger = logging.getLogger(__name__)
    logging.getLogger('matplotlib').setLevel(logging.WARNING)
################################################################################

matplotlib.rcParams['xtick.top'] = True
matplotlib.rcParams['ytick.right'] = True
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
matplotlib.rcParams['grid.alpha'] = 0.5

################################################################################


@u.quantity_input
def calc_sf(radius: u.m, distance: u.m):
    """
    Calculate scaling factor for given radius and distance

    Parameters
    ----------
    radius : Quantity
        Radius [Rsun]
    distance : Quantity
        distance [pc]

    Returns
    -------
    sf : Quantity
        Scaling factor
    """
    return ((radius/(distance))**2).decompose()


@u.quantity_input
def calc_radius(sf, distance: u.m):
    """
    Radius for given scaling factor and distance

    Parameters
    ----------
    sf : float or list
        Scaling factor
    distance : Quantity
        distance [pc]

    Returns
    -------
    radius : Quantity
        radius [Rsun]
    """
    return (sf**0.5 * distance).to(u.solRad)


@u.quantity_input
def calc_luminosity(radius: u.m, Te: u.K):
    """
    Calculate luminosity from radius and temperature

    Parameters
    ----------
    radius : Quantity
        radius [Rsun]
    Te : Quantity
        Temperature [K]

    Returns
    -------
    luminosity : Quantity
        luminosity [Lsun]
    """

    L = (4 * np.pi * const.sigma_sb * radius**2 * Te**4).to(u.solLum)
    return L.to(u.solLum)


def load_data(file_name, mode):
    """
    Load star data.

    Parameters
    ----------
    file_name : str
        File name containing data.
    mode : str
        Data loading mode ('csv' or 'vosa').

    Returns
    -------
    data : DataFrame
        Loaded data.
    """
    if mode == 'csv':
        data = pd.read_csv(file_name, engine='python', header=0)

    if mode == 'vosa':
        data = gather_extinction_corrected_data_from_VOSA(file_name)

    data.set_index('FilterID', inplace=True)
    return data


def gather_extinction_corrected_data_from_VOSA(file_name, load_params=False):
    """
    Gather extinction-corrected data from VOSA single-fit file.

    Parameters
    ----------
    file_name : str
        File name containing VOSA data.
    load_params : bool, optional
        Load additional parameters. Defaults to False.

    Returns
    -------
    data : DataFrame
        Extinction-corrected data.
    """
    # Reading parameters derived by VOSA
    data = pd.read_csv(file_name, engine='python', comment='#',
                       delim_whitespace=True, skipinitialspace=True, header=None)
    data.columns = ['FilterID', 'wavelength', 'Obs.Flux', 'Obs.Error',
                    'flux', 'error', 'model_flux', 'Fitted', 'Excess', 'FitExc', 'UpLim']

    # Removing filters noted as "upper limit"
    logger.debug('%s removed due to upper limit' %
                 (data[data['UpLim'] == '1']['FilterID'].values))
    data = data[data['UpLim'] == '---']
    data = data.drop(columns=['Obs.Flux', 'Obs.Error',
                     'Fitted', 'Excess', 'FitExc', 'UpLim'])
    if not load_params:
        data = data[['FilterID', 'wavelength', 'flux', 'error']]
    if load_params:
        data = data[['FilterID', 'wavelength', 'flux', 'error', 'model_flux']]
    return data


def save_fig(plot_name, show_plot, folder, format='jpg'):
    """
    Save a matplotlib figure to a specified folder and handle display options.

    Parameters
    ----------
    plot_name : str
        The name of the file to save the plot as, including the file extension (e.g., 'plot.png').
    show_plot : bool
        If True, the plot will be shown after saving. If False, the plot will be closed after saving.
    folder : str
        The directory in which to save the plot.
    format : str
        Format of the plot (jpg, pdf etc.)
    """
    if not os.path.exists(folder):
        os.makedirs(folder)
    plt.savefig(plot_name+'.'+format, dpi=300, bbox_inches='tight')
    if not show_plot:
        plt.close()


def fancy_text(text, buffer=12):
    """
    Create a fancy title with a specified buffer around the text.

    Parameters
    ----------
    text : str
        The text to be displayed as the title.
    buffer : int, optional
        The buffer of padding around the text. Default is 12.

    Returns
    -------
    str
        A formatted string with the fancy title.
    """
    thick = '=' * (len(text) + 4 * buffer)
    middle = '%s%s%s%s%s' % ('-' * buffer, ' ' * buffer,
                             text, ' ' * buffer, '-' * buffer)
    thin = '-' * (len(text) + 4 * buffer)
    return '%s\n%s\n%s\n%s\n%s' % (thick, thin, middle, thin, thick)


def estimate_runtime(current_iteration, total_iterations, start_time, message='Processing'):
    """
    Estimates time of completion.

    Parameters
    ----------
    current_iteration : int
        The current iteration.
    total_iterations : int
        Total number of iterations.
    start_time : float
        The start time of the process.
    message : str, optional
        Custom message for progress. Defaults to 'Processing'.
    """
    end_time = time.perf_counter()
    total_time = (end_time-start_time) * total_iterations/current_iteration
    logger.warning('%s: ETA ~ %d s' % (message, total_time))


class Star:
    """
    Class for working with individual stars and their SED fitting.

    Attributes
    ----------
    name : str
        Name of the source.
    distance : float
        Distance to the source in parsecs (pc).
    e_distance : float
        Error in the distance measurement in parsecs (pc).
    filters_to_drop : list, optional
        List of filters to be dropped before fitting the SED. Default is an empty list.
    wavelength_range : list, optional
        Minimum and maximum values of the filters to consider for fitting. Default is [0, 1,000,000,000].
    r_limits : str, optional
        Limits for the radius. Default is 'blackbody'.
    run_name : str, optional
        Name of the run. Default is '0'.
    component : str, optional
        Component identifier (should be either 'A', 'B' or 'C'). Default is 'A'.
    _type : str, optional
        Type of the source. Default is 'Star'.
    model : Model, optional
        Model to fit to the data.
    data_all : pd.DataFrame
        A copy of the original data with an added 'fitted' column indicating whether the filter is used. Mmodel flux, residuals, and statistical measures are added in :func:`Fitter.get_parameters_from_chi2_minimization` and :func:`Fitter.get_parameters_from_noisy_chi2_minimization`.
    data : pandas.DataFrame, optional
        Dataframe containing the observational data. It is cropped to the fitted filters based on :attr:`filters_to_drop` and :attr:`wavelength_range`.
    data_not_fitted : DataFrame
        DataFrame containing only removed filters based on :attr:`filters_to_drop` and :attr:`wavelength_range`..
    N_points : int
        Updated number of data points.
    filters_to_drop_all : array-like
        List of filters that are to be dropped based on :attr:`filters_to_drop` and :attr:`wavelength_range`.
    N : int
        The number of filters that are fitted.
    Te_blackbody : float
        Effective temperature from the blackbody fit [K].
    log_sf_blackbody : float
        Logarithm of the scaling factor from the blackbody fit.
    sf_list : list
        List of scaling factors (NOT in log space).
    _minimising_param : str
        Parameter used for minimization (Only 'chi2' is implimented).
    df_chi2 : pd.DataFrame
        DataFrame with fitting parameters for all fits.
    df_chi2_trimmed : pd.DataFrame
        DataFrame with fitting parameters for the top fits after trimming (default: top 1000).
    df_chi2_noisy : pd.DataFrame
        DataFrame with best fitting parameters for each noisy iteration.
    N_Np : float
        Number of data points minus the number of free parameters in the model.
    chi2_r : float
        Reduced chi^2 value indicating goodness of fit.
    vgf2 : float
        Variance of fractional residuals squared (vgf^2) for the fit.
    vgfb2 : float
        Variance of fractional residuals squared with broadened error (vgfb^2) for the fit.
    <dim> : float
        Updates various best fit parameters available in star.df_chi2 (e.g., chi2, sf, R, L, Te, MH, logg).
    chi2_median : float
        Median of the chi^2 values from noisy fits.
    chi2_r_median : float
        Median of the reduced chi^2 values from noisy fits.
    vgf2_median : float
        Median of the vgf^2 values from noisy fits.
    vgfb2_median : float
        Median of the vgfb^2 values from noisy fits.
    <dim>_median : float
        Updates the median best fit parameter values available in :attr:`Star.df_chi2_noisy` (sf, R, L, Te, MH, logg, etc.)
    <dim>_error_lower : float
        Updates the 16th percentile error for best fit parameters available in :attr:`Star.df_chi2_noisy`.
    <dim>_error_upper : float
        Updates the 84th percentile error for best fit parameters available in :attr:`Star.df_chi2_noisy`.

    """

    def __init__(self,
                 name,
                 distance,
                 e_distance,
                 filters_to_drop=[],
                 wavelength_range=[0, 1_000_000_000],
                 data=None,
                 model=None,
                 r_limits='blackbody',
                 run_name='0',
                 component='A',
                 _type='Star'):
        """
        Initializes a Star object using name, distance, model_to_fit, etc.

        Parameters
        ----------
        name : str
            Name of the source.
        distance : float
            Distance to the source in parsecs (pc).
        e_distance : float
            Error in the distance measurement in parsecs (pc).
        filters_to_drop : list, optional
            List of filters to be dropped before fitting the SED. Default is an empty list.
        wavelength_range : list, optional
            Minimum and maximum values of the filters to consider for fitting. Default is [0, 1,000,000,000].
        r_limits : str, optional
            Limits for the radius. Default is 'blackbody'.
        run_name : str, optional
            Name of the run. Default is '0'.
        component : str, optional
            Component identifier (should be either 'A', 'B' or 'C'). Default is 'A'.
        _type : str, optional
            Type of the source. Default is 'Star'.
        data : pandas.DataFrame, optional
            Dataframe containing the observational data.
        model : Model, optional
            Model to fit to the data.

        Raises
        ------
        NotImplementedError
            If 'component' is not 'A', 'B', or 'C'.

        Notes
        -----
        If `data` is provided, it will be sorted by wavelength and processed.
        If `model` is also provided, it will be configured with the processed data.
        """
        if component not in ['A', 'B', 'C', 'Total']:
            raise NotImplementedError(
                "'component' should be either 'A', 'B' or 'C'.")

        logger.info(fancy_text('%s %s' % (name, component)))

        self.name = name
        self.distance = distance
        self.e_distance = e_distance
        self.filters_to_drop = filters_to_drop
        self.wavelength_range = wavelength_range
        self.r_limits = r_limits
        self.run_name = run_name
        self.component = component
        self._type = _type

        if data is not None:
            data = data.sort_values(by='wavelength')
            self.data = data
            self.process_data()
            self.drop_filters()
            self.create_sf_list(r_limits)

            if model is not None:
                self.model = model
                self.model.filter_list = self.data.index.values
                self.model.crop_based_on_filter_list()
                self.model.create_star_dataarrays(self)

    def process_data(self):
        """
        Process the data by validating the data and calculating necessary parameters.

        Raises
        ------
        ValueError
            If duplicate filters are found in the data.

        Notes
        -----
        *Updates:* :attr:`Star.data`, :attr:`Star.data_not_fitted`, :attr:`Star.N_points`

        The method includes several steps to clean and preprocess the data before fitting:
        - Ensures no duplicate filters are present in the data.
        - Filters with zero errors are replaced with a value calculated as 110% of the maximum error fraction.
        - Errors are adjusted to ensure a minimum error fraction of 2% for VGF and 10% for VGFB calculations.
        - Logarithmic values of wavelength and flux are computed to facilitate fitting.
        """
        if len(self.data.index) > len(self.data.index.unique()):
            u, c = np.unique(self.data.index, return_counts=True)
            dup = u[c > 1]
            raise ValueError(
                'Some filters are repeated. Remove the duplicates %s.' % dup)

        # Replacing zeros in errors with 110% of max error
        self.data['error_fraction'] = self.data['error']/self.data['flux']
        self.data.loc[self.data.error == 0, 'error'] = self.data['flux'] * \
            (self.data['error_fraction'].max()*1.1)
        # Recalculating errors_fraction
        self.data['error_fraction'] = self.data['error']/self.data['flux']

        # error modification for calculating vgf (minimum error = 2%) and vgfb (minimum error = 10%)
        self.data['error_2percent'] = np.where(
            self.data['error_fraction'] < 0.02, 0.02, self.data['error_fraction'])*self.data['flux']
        self.data['error_10percent'] = np.where(
            self.data['error_fraction'] < 0.10, 0.10, self.data['error_fraction'])*self.data['flux']

        self.data['log_wavelength'] = np.log10(self.data['wavelength'])
        self.data['log_flux'] = np.log10(self.data['flux'])
        self.data['e_log_flux'] = self.data['error'] / \
            self.data['flux']/np.log(10)

    def drop_filters(self):
        """
        Drops unnecessary filters based on wavelength range.

        Notes
        -----
        *Updates:* :attr:`Star.data_all`, :attr:`Star.filters_to_drop_all`, :attr:`Star.data`, :attr:`Star.N`

        This function updates self.data and self.data_not_fitted DataFrames by filtering out 
        entries outside the specified wavelength range or with negative wavelengths.

        The function also logs the filters that are fitted and not fitted.
        """
        self.data_all = self.data.copy()
        _excess_filters = self.data[(self.data['wavelength'] < self.wavelength_range[0]) |
                                    (self.data['wavelength'] > self.wavelength_range[1]) |
                                    (self.data['wavelength'] < 0)].index
        self.filters_to_drop_all = np.concatenate(
            (self.filters_to_drop, _excess_filters))

        self.data_all['fitted'] = 1
        self.data_all.loc[self.filters_to_drop_all, 'fitted'] = 0
        self.data = self.data_all[self.data_all['fitted'] == 1]
        self.data_not_fitted = self.data_all[self.data_all['fitted'] == 0]

        self.N = len(self.data)

        # logging filters to be fitted and not-fitted
        _filter_table = self.data_all[['wavelength', 'fitted']].copy()
        _filter_table['Fitted'] = np.where(
            _filter_table.fitted == 1, _filter_table.index, np.nan)
        _filter_table['Not fitted'] = np.where(
            _filter_table.fitted == 0, _filter_table.index, np.nan)
        _filter_table = _filter_table[[
            'wavelength', 'Fitted', 'Not fitted']].set_index('wavelength')
        logger.info(_filter_table.sort_index().fillna(''))
        logger.info('Filters: used/all = %d/%d' %
                    (len(self.data), len(self.data_all)))

    def fit_blackbody(self, p0=[5000., -20], plot=False, show_plot=True, folder=None):
        """
        Fits a blackbody model to the SED.

        Parameters
        ----------
        p0 : list, optional
            Initial guess for the blackbody model parameters. Default is [5000., -20].
        plot : bool, optional
            Whether to plot the blackbody fit. Default is False.
        show_plot : bool, optional
            Whether to show the plot. Default is True.
        folder : str, optional
            Folder to save the plot. No plot is saved if `folder` is `None`.

        Notes
        -----
        *Updates:* :attr:`Star.Te_blackbody`, :attr:`Star.log_sf_blackbody`
        """
        Fitter.blackbody(self, p0=p0)
        if plot:
            fig, ax = plt.subplots(figsize=(8, 6), nrows=3, ncols=1)
            [axi.set_axis_off() for axi in ax.ravel()]
            ax[0] = fig.add_axes([0.1, 0.40, 0.85, 0.55])
            ax[1] = fig.add_axes([0.1, 0.25, 0.85, 0.15])
            ax[2] = fig.add_axes([0.1, 0.10, 0.85, 0.15])

            self.plot_sed_obs(ax=ax[0])
            self.plot_blackbody_model(ax=ax[0])
            # Residual
            bb_flux = 10**Fitter.get_blackbody_spectrum_loglog(
                self.data_all['log_wavelength'], self.Te_blackbody, self.log_sf_blackbody)
            fractional_residual = (
                self.data_all['flux'] - bb_flux)/self.data_all['flux']
            ax[1].plot(self.data_all['wavelength'],
                       fractional_residual, 'r--', label='')
            # chi2
            bb_flux = 10**Fitter.get_blackbody_spectrum_loglog(
                self.data['log_wavelength'], self.Te_blackbody, self.log_sf_blackbody)
            chi2_i = (self.data['flux']-bb_flux)**2 / self.data['error']**2
            ax[2].plot(self.data['wavelength'], chi2_i, 'r--',
                       label='$χ^2$ = %.2f' % (chi2_i.sum()))
            # Titles and labels
            ax[0].set_title(self.name + ' (blackbody fit)',
                            x=0, y=1, ha='left')
            ax[0].set_ylabel('log(Flux) (erg s$^{-1}$ cm$^{-2}$ $Å$$^{-1}$)')
            ax[1].set_ylabel('Fractional\nResidual')
            ax[1].set_xlabel('log(Wavelength) (Å)')
            # decoration
            plt.setp(ax[0].get_xticklabels(), visible=False)
            plt.setp(ax[1].get_xticklabels(), visible=False)
            # taking xlim and ylims from one axis to another
            ax[0].tick_params(which='both', direction='out', length=4)
            ax[1].tick_params(which='both', direction='out', length=4)
            ax[0].legend()
            ax[2].legend()
            for axes in ax.ravel():
                axes.grid()
                axes.set_xscale('log')
            ax[0].set_yscale('log')
            ax[2].set_xlim(ax[0].get_xlim())
            if folder is not None:
                plot_name = folder+'%s_blackbody' % (self.name)
                save_fig(plot_name, show_plot, folder)

    def create_sf_list(self, r_limits):
        """
        Create list of scaling factors.
        The limits on radius are either taken from the used (via r_limits parameter).
        Or calculated automatically from blackbody fitting.

        Parameters
        ----------
        r_limits : str or list
            Limits for the radius. Can be 'blackbody' to use blackbody fitting results or a list of two values for radius limits.

        Notes
        -----
        *Updates:* :attr:`Star.sf_list`
        """
        LOG_SF_STEPSIZE = 0.01
        LOG_SF_FLEXIBILITY = 2

        if r_limits == 'blackbody':
            mode = 'blackbody'
            if self.component in ['B', 'C']:
                logger.warning('It is recommended to give r_limits as [r_min, r_max] for %s %s' % (
                    self.name, self.component))
        else:
            if not isinstance(r_limits, list):
                raise ValueError(
                    'r_limits must be either "blackbody" or a list of two values, e.g. [0.001, 1]')
            if len(r_limits) != 2:
                raise ValueError(
                    'r_limits must be either "blackbody" or a list of two values, e.g. [0.001, 1]')

            mode = 'radius limits'

        if mode == 'blackbody':
            if not hasattr(self, 'Te_blackbody'):
                self.fit_blackbody()
            self.sf_list = 10**np.arange(self.log_sf_blackbody-LOG_SF_FLEXIBILITY,
                                         self.log_sf_blackbody+LOG_SF_FLEXIBILITY,
                                         LOG_SF_STEPSIZE)
        if mode == 'radius limits':
            sf_B_min = calc_sf(r_limits[0] * u.solRad, self.distance*u.pc)
            sf_B_max = calc_sf(r_limits[1] * u.solRad, self.distance*u.pc)
            self.sf_list = 10**np.arange(np.log10(sf_B_min),
                                         np.log10(sf_B_max),
                                         LOG_SF_STEPSIZE)

        logger.debug('Using mode = %s\nsf_list = [%e,...,%e] %d' % (mode,
                                                                    self.sf_list[0], self.sf_list[-1], len(self.sf_list)))

    def fit_chi2(self, refit=True, _trim=1000, _minimising_param='chi2'):
        """
        Fits the SED using a chi-squared minimization approach.

        Parameters
        ----------
        refit : bool, optional
            Whether to refit the model. Default is True.
        _trim : int, optional
            Number of data points to trim for fitting. Default is 1000.
        _minimising_param : str, optional
            Parameter used for minimizing (only 'chi2' supported as of now). Default is 'chi2'.

        Notes
        -----
        *Updates:* :attr:`Star._minimising_param`
        """
        self._minimising_param = _minimising_param
        Fitter.calculate_chi2(self, self.model, refit=refit,
                              _trim=_trim, _minimising_param=_minimising_param)
        Fitter.get_parameters_from_chi2_minimization(self, self.model)

    def fit_noisy_chi2(self, total_iterations=100, refit=True, plot=False, _percentile_threshold=10):
        """
        Fits the SED using a chi-squared minimization approach. The SED is fitted total_iterations times using bootstrapping.

        Parameters
        ----------
        total_iterations : int, optional
            Number of iterations used for bootstrapping. Default is 100.
        refit : bool, optional
            Whether to refit the model. Default is True.
        _percentile_threshold : int, optional
            Percentile threshold used for setting up the extent of the noisy parameter grid. Default is 10.
        plot : bool, optional
            Whether to plot the blackbody fit. Default is False.
        """
        Fitter.calculate_noisy_chi2(self, self.model, total_iterations=total_iterations, refit=refit,
                                    plot=plot, _percentile_threshold=_percentile_threshold)
        Fitter.get_parameters_from_noisy_chi2_minimization(self, self.model)

    def plot_sed_obs(self, ax=None):
        """
        Plot the observed SED.

        Parameters
        ----------
        ax : matplotlib.axes, optional
            Matplotlib axes to plot on. If None, a new subplot is created.
        """
        if ax is None:
            fig, ax = plt.subplots()
        ax.scatter(self.data_not_fitted['wavelength'], self.data_not_fitted['flux'],
                   color='k', marker='o', label='No Fit', s=30, facecolors='none', zorder=2)
        ax.errorbar(self.data['wavelength'], self.data['flux'], yerr=self.data['error'],
                    color='k', label='Obs', fmt='none', lw=2, zorder=3, capsize=4)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set(xlabel='Wavelength (Å)',
               ylabel='Flux (erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$)')
        ax.grid()

    def plot_sed_fitted(self, ax=None, color='C1', ls='--', median=False):
        """
        Plot the fitted SED.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Matplotlib axes to plot on. If None, a new subplot is created.
        color : str, optional
            Color of the plot. Default is 'C1'.
        ls : str, optional
            Linestyle of the plot. Default is '--'.
        median : bool, optional
            Whether to use the median values. Default is False.
        """
        if ax is None:
            fig, ax = plt.subplots()
        y = self.data_all['model_flux_median'] if median else self.data_all['model_flux']
        ax.plot(self.data_all['wavelength'], y,
                color=color, label='%s' % self.component, zorder=3, ls=ls)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set(xlabel='Wavelength (Å)',
               ylabel='Flux (erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$)')
        ax.grid()

    def plot_noisy_seds(self, ax=None, color='C1', ls='--'):
        """
        Plot the noisy SEDs.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Matplotlib axes to plot on. If None, a new subplot is created.
        color : str, optional
            Color of the plot. Default is 'C1'.
        ls : str, optional
            Linestyle of the plot. Default is '--'.
        """
        if ax is None:
            fig, ax = plt.subplots()
        for idx, sf_noisy in enumerate(self.df_chi2_noisy['sf']):
            flux_noisy = self.model.da * sf_noisy
            for dim in self.model.da.dims[1:]:
                flux_noisy = flux_noisy.sel(
                    {dim: self.df_chi2_noisy[dim][idx]})
            ax.plot(self.data['wavelength'], flux_noisy, color=color,
                    linestyle=ls, label='', lw=0.2, alpha=0.2, zorder=0)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set(xlabel='Wavelength (Å)',
               ylabel='Flux (erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$)')
        ax.grid()

    def plot_residual(self, ax=None, color='C1', ls='--', median=False):
        """
        Plot the residuals.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Matplotlib axes to plot on. If None, a new subplot is created.
        color : str, optional
            Color of the plot. Default is 'C1'.
        ls : str, optional
            Linestyle of the plot. Default is '--'.
        median : bool, optional
            Whether to use the median values. Default is False.
        """
        if ax is None:
            fig, ax = plt.subplots()
        y = self.data_all['residual_flux_median'] if median else self.data_all['residual_flux']
        ax.plot(self.data_all['wavelength'], self.data_all['residual_flux'],
                color=color, marker='.', label='', zorder=2, ls=ls)
        ax.set_xscale('log')
        ax.grid()

    def plot_fractional_residual(self, ax=None, FR_cutoff=None, color='C1', ls='--', median=False):
        """
        Plot the fractional residuals.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Matplotlib axes to plot on. If None, a new subplot is created.
        FR_cutoff : float, optional
            Cutoff value for identifying data with fractional residual excess.
        color : str, optional
            Color of the plot. Default is 'C1'.
        ls : str, optional
            Linestyle of the plot. Default is '-'.
        median : bool, optional
            Whether to use the median values. Default is False.
        """
        if ax is None:
            fig, ax = plt.subplots()
        y = self.data_all['fractional_residual_median'] if median else self.data_all['fractional_residual']
        ax.plot(self.data_all['wavelength'], y,
                label='', marker='', color=color, ls=ls)
        y = self.data_not_fitted['fractional_residual_median'] if median else self.data_not_fitted['fractional_residual']
        ax.scatter(self.data_not_fitted['wavelength'], y,
                   color=color, marker='o', facecolor='none', edgecolor=color)
        ax.set_xscale('log')
        ax.set(xlabel='Wavelength (Å)', ylabel='FR')
        ax.grid()
        if FR_cutoff is not None:
            if not median:
                _excess = self.data_all[self.data_all['fractional_residual']
                                        > FR_cutoff]
                y = _excess['fractional_residual']
            if median:
                _excess = self.data_all[self.data_all['fractional_residual_median']
                                        > FR_cutoff]
                y = _excess['fractional_residual_median']
            ax.scatter(_excess['wavelength'], y,
                       label='Excess', marker='|', color='r')
            ax.axhline(FR_cutoff, c='0.5', lw=0.5)

    def plot_ewr(self, ax=None,  deviation_cutoff=None, color='C1', ls='--', median=False):
        """
        Plot the Error Weighted Residual.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Matplotlib axes to plot on. If None, a new subplot is created.
        deviation_cutoff : float, optional
            Cutoff value for deviation in EWR.
        color : str, optional
            Color of the plot. Default is 'C1'.
        ls : str, optional
            Linestyle of the plot. Default is '-'.
        median : bool, optional
            Whether to use the median values. Default is False.
        """
        if ax is None:
            fig, ax = plt.subplots()
        y = self.data_all['ewr_median'] if median else self.data_all['ewr']
        ax.plot(self.data_all['wavelength'], y,
                label='', marker='', color=color, ls=ls)
        y = self.data_not_fitted['ewr_median'] if median else self.data_not_fitted['ewr']
        ax.scatter(self.data_not_fitted['wavelength'], y,
                   color=color, marker='o', facecolor='none', edgecolor=color)
        ax.set_xscale('log')
        ax.set(xlabel='Wavelength (Å)', ylabel='EWR')
        ax.grid()
        if deviation_cutoff is not None:
            if not median:
                _excess = self.data_all[abs(
                    self.data_all['ewr']) > deviation_cutoff]
                y = _excess['ewr']
            if median:
                _excess = self.data_all[self.data_all['ewr_median']
                                        > deviation_cutoff]
                y = _excess['ewr_median']
            ax.scatter(_excess['wavelength'], y,
                       label='Excess', marker='|', color='r')
            ax.axhline(deviation_cutoff, c='0.5', lw=0.5)
            ax.axhline(-deviation_cutoff, c='0.5', lw=0.5)

    def plot_chi2(self, ax=None, color='C1', ls='--', median=False):
        """
        Plot the chi-squared values.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Matplotlib axes to plot on. If None, a new subplot is created.
        color : str, optional
            Color of the plot. Default is 'C1'.
        ls : str, optional
            Linestyle of the plot. Default is '-'.
        median : bool, optional
            Whether to use the median values. Default is False.
        """
        if ax is None:
            fig, ax = plt.subplots()
        y = self.data['chi2_i_median'] if median else self.data['chi2_i']
        ax.plot(self.data['wavelength'], y,
                label='', marker='', color=color, ls=ls)
        # ax.scatter(self.data_not_fitted['wavelength'],self.data_not_fitted['chi2_i'],
        #            color=color, marker='o', facecolor='none', edgecolor=color)
        ax.set_xscale('log')
        ax.set(xlabel='Wavelength (Å)', ylabel='$χ^2_i$')
        ax.grid()

    def plot_hrd(self, ax=None, color='C1', median=False, cax=None):
        """
        Plot the Hertzsprung-Russell Diagram (HRD).

        Parameters
        ----------
        ax : matplotlib.axes, optional
            Matplotlib axes to plot on. If None, a new subplot is created.
        color : str, optional
            Color of the plot.
        median : bool, optional
            Whether to use the median values.
        cax : matplotlib.axes, optional
            Matplotlib axes for the colorbar.
        """
        if ax is None:
            fig, ax = plt.subplots()
        if not median:
            ax.scatter(self.Te, self.L, marker='s', label=self.component,
                       color=color, s=40, rasterized=True, zorder=3)
            p = ax.scatter(self.df_chi2_trimmed['Te'], self.df_chi2_trimmed['L'], c=self.df_chi2_trimmed['chi2'],
                           marker='.', cmap='summer', norm=matplotlib.colors.LogNorm())
        if median:
            ax.scatter(self.Te_median, self.L_median, marker='s',
                       label=self.component, color=color, s=40, rasterized=True, zorder=3)
            p = ax.scatter(self.df_chi2_noisy['Te'], self.df_chi2_noisy['L'], c=self.df_chi2_noisy['chi2'],
                           marker='.', cmap='summer', norm=matplotlib.colors.LogNorm())
        ax.set(xlabel='Te (K)', ylabel='L ($L_⊙$)')
        Model.plot_isochrone_and_wd(ax=ax)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.invert_xaxis()
        ax.grid()
        if cax is not False:
            plt.colorbar(p, ax=ax, cax=cax, pad=0, label='$χ^2$')

    def plot_temperature_radius(self, ax=None, color='C1', median=False, cax=None):
        """
        Plot the Temperature-Radius distribution.

        Parameters
        ----------
        ax : matplotlib.axes, optional
            Matplotlib axes to plot on. If None, a new subplot is created.
        color : str, optional
            Color of the plot.
        median : bool, optional
            Whether to use the median values.
        cax : matplotlib.axes, optional
            Matplotlib axes for the colorbar.
        """
        if ax is None:
            fig, ax = plt.subplots()
        if not median:
            ax.scatter(self.Te, self.R, marker='s', label='Best fit',
                       color=color, s=40, rasterized=True, zorder=3)
            p = ax.scatter(self.df_chi2_trimmed['Te'], self.df_chi2_trimmed['R'], c=self.df_chi2_trimmed['chi2'],
                           marker='.', cmap='summer', norm=matplotlib.colors.LogNorm())
        if median:
            ax.scatter(self.Te_median, self.R_median, marker='s',
                       label='Best fit', color=color, s=40, rasterized=True, zorder=3)
            p = ax.scatter(self.df_chi2_noisy['Te'], self.df_chi2_noisy['R'], c=self.df_chi2_noisy['chi2'],
                           marker='.', cmap='summer', norm=matplotlib.colors.LogNorm())
        ax.set(xlabel='Te (K)', ylabel='R ($R_⊙$)')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.invert_xaxis()
        ax.grid()
        if cax is not False:
            plt.colorbar(p, ax=ax, cax=cax, pad=0, label='$χ^2$')

    def plot_temperature_chi2(self, ax=None, color='C1', median=False, cax=None):
        """
        Plot the Temperature-Chi-squared distribution.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Matplotlib axes to plot on. If None, a new subplot is created.
        color : str, optional
            Color of the plot. Default is 'C1'.
        median : bool, optional
            Whether to use the median values. Default is False.
        cax : matplotlib.axes.Axes, optional
            Matplotlib axes for the colorbar.
        """
        if ax is None:
            fig, ax = plt.subplots()
        if not median:
            ax.scatter(self.Te, self.chi2, marker='s', label='Best fit',
                       color=color, s=40, rasterized=True, zorder=3)
            p = ax.scatter(self.df_chi2_trimmed['Te'], self.df_chi2_trimmed['chi2'], c=self.df_chi2_trimmed['chi2'],
                           marker='.', cmap='summer', norm=matplotlib.colors.LogNorm())
        if median:
            ax.scatter(self.Te_median, self.chi2_median, marker='s',
                       label='Best fit', color=color, s=40, rasterized=True, zorder=3)
            p = ax.scatter(self.df_chi2_noisy['Te'], self.df_chi2_noisy['chi2'], c=self.df_chi2_noisy['chi2'],
                           marker='.', cmap='summer', norm=matplotlib.colors.LogNorm())
        ax.set(xlabel='Te (K)', ylabel='$χ^2$')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.invert_xaxis()
        ax.grid()
        if cax is not False:
            p.set_clim(ax.get_ylim())
            plt.colorbar(p, ax=ax, cax=cax, pad=0, label='$χ^2$')

    def plot_temperature_vgfb2(self, ax=None, color='C1', median=False, cax=None):
        """
        Plot the Temperature-vgf_b^2 distribution.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Matplotlib axes to plot on. If None, a new subplot is created.
        color : str, optional
            Color of the plot. Default is 'C1'.
        median : bool, optional
            Whether to use the median values. (Not supported for vgfb2)
        cax : matplotlib.axes.Axes, optional
            Matplotlib axes for the colorbar.

        Raises
        ------
        NotImplementedError
            If median mode is True (not supported for vgfb2).
        """
        if median is not False:
            raise NotImplementedError(
                'Median mode not supported for vgfb2 yet.')

        if ax is None:
            fig, ax = plt.subplots()
        if not median:
            ax.scatter(self.Te, self.vgfb2, marker='s', label='Best fit',
                       color=color, s=40, rasterized=True, zorder=3)
            p = ax.scatter(self.df_chi2_trimmed['Te'], self.df_chi2_trimmed['vfgb2'], c=self.df_chi2_trimmed['chi2'],
                           marker='.', cmap='summer', norm=matplotlib.colors.LogNorm())
            ax.set(xlabel='Te (K)', ylabel='$vgf_b^2$')
        # if median:
        #     ax.scatter(self.Te_median,self.vgfb2_median, marker='s', label='Best fit',color=color, s=40,rasterized = True, zorder=3)
        #     p = ax.scatter(self.df_chi2_noisy['Te'], self.df_chi2_noisy['vgfb2'], c=self.df_chi2_noisy['chi2'],
        #             marker='.',cmap='summer', norm=matplotlib.colors.LogNorm())
        #     ax.set(xlabel='Te (K)',ylabel='$χ^2$')
        self._plot_helper_right(ax=ax, p=p, cax=cax)

    def plot(self, add_noisy_seds=False,
             show_plot=True, folder=None,
             FR_cutoff=0.5, _median=False):
        """
        Plot detailed properties of a single SED.

        Parameters
        ----------
        add_noisy_seds : bool, optional
            Whether to add noisy SEDs to the plot.
        show_plot : bool, optional
            Whether to display the plot.
        folder : str, optional
            Folder to save the plot. No plot is saved if `folder` is `None`.
        FR_cutoff : float, optional
            Cutoff value for identifying data with fractional residual excess.
        median : bool, optional
            Whether to use the median values for plotting.

        Raises
        ------
        NotImplementedError
            If component is other than 'A', 'B', and 'C'.

        Notes
        -----
        This method plots various detailed properties of a single SED, including observed and fitted SEDs,
        fractional residuals, error weighted residuals, chi-squared values, HR diagram, and temperature-related
        distributions.
        """
        if self.component not in ['A', 'B', 'C']:
            raise NotImplementedError(
                "Components other than 'A', 'B', and 'C' are not supported yet. Make the plot yourself by modifying the following function.")
        if self.component == 'A':
            color, ls = 'C1', '--'
        if self.component == 'B':
            color, ls = 'C0', '-.'
        if self.component == 'C':
            color, ls = 'C4', ':'
        fig, ax = self.plot_sed_skeleton_detailed()

        self.plot_sed_obs(ax[0, 0])
        self.plot_sed_fitted(ax[0, 0], color=color, ls=ls, median=_median)
        if add_noisy_seds:
            self.plot_noisy_seds(ax=ax[0, 0], color=color, ls=ls)

        self.plot_fractional_residual(
            ax[1, 0], color=color, ls=ls, FR_cutoff=FR_cutoff, median=_median)
        self.plot_ewr(ax[2, 0], color=color, ls=ls, median=_median)
        self.plot_chi2(ax[3, 0], color=color, ls=ls, median=_median)

        self.plot_hrd(ax[0, 1], color=color, median=_median, cax=ax[0, 2])
        self.plot_temperature_radius(
            ax[1, 1], color=color, median=_median, cax=ax[0, 2])
        self.plot_temperature_chi2(
            ax[2, 1], color=color, median=_median, cax=ax[0, 2])

        self.update_xylimits_detailed(ax)
        ax[0, 0].grid()
        ax[0, 0].legend()

        title = '%s (%s) [%s]\n' % (self.name, self.model.name, self.run_name)
        dims = np.concatenate([self.df_chi2.columns[:-6], ['R', 'L']])
        for dim in dims:
            title += '%s = %.2E' % (dim, getattr(self, dim))
            if hasattr(self, dim+'_error_upper'):
                ten_factor = 10**np.floor(np.log10(getattr(self, dim)))
                title += '$^{+%.2f}_{-%.2f}$' % (getattr(self, dim+'_error_upper')/ten_factor,
                                                 getattr(self, dim+'_error_lower')/ten_factor)
            title += '   '
        title += '\n $χ^2_r$ %.1f, vgfb$^2$ %.2f' % (
            self.chi2_r, self.vgfb2)
        ax[0, 0].set_title(title, loc='left')
        if folder is not None:
            plot_name = folder+'%s_single_%s' % (self.name, self.run_name)
            save_fig(plot_name, show_plot, folder)

    def plot_public(self, median=False, add_noisy_seds=False, ax=None,
                    FR_cutoff=0.5, show_plot=True, folder=None):
        """
        Plot a simplified version of various properties of a single SED for public view.

        Parameters
        ----------
        median : bool, optional
            Whether to use the median values.
        add_noisy_seds : bool, optional
            Whether to add noisy SEDs to the plot.
        ax : matplotlib.axes, optional
            Matplotlib axes to plot on. If None, a new subplot is created.
        FR_cutoff : float, optional
            Cutoff value for identifying data with fractional residual excess.
        show_plot : bool, optional
            Whether to display the plot.
        folder : str, optional
            Folder to save the plot. No plot is saved if `folder` is `None`.

        Raises
        ------
        NotImplementedError
            If component is other than 'A', 'B', and 'C'.
        """
        if self.component not in ['A', 'B', 'C']:
            raise NotImplementedError(
                "Components other than 'A', 'B', and 'C' are not supported yet. Make the plot yourself by modifying the following function.")
        if self.component == 'A':
            color, ls = 'C1', '--'
        if self.component == 'B':
            color, ls = 'C0', '-.'
        if self.component == 'C':
            color, ls = 'C4', ':'
        if ax is None:
            fig, ax = self.plot_sed_skeleton_public()

        self.plot_sed_obs(ax[0])
        self.plot_sed_fitted(ax[0], color=color, ls=ls, median=median)
        if add_noisy_seds:
            self.plot_noisy_seds(ax=ax[0], color=color, ls=ls)

        self.plot_fractional_residual(
            ax[1], color=color, ls=ls, FR_cutoff=FR_cutoff, median=median)
        self.plot_ewr(ax[2], color=color, ls=ls, median=median)
        self.update_xylimits_public(ax)
        ax[0].grid()
        ax[0].legend()
        title = self.name
        title += '  %d K, %.2e R$_⊙$, %.2e L$_⊙$' % (
            self.Te, self.R, self.L)
        ax[0].set_title(title, loc='left')

        if folder is not None:
            plot_name = folder + \
                '%s_single_%s_public.jpg' % (self.name, self.run_name)
            save_fig(plot_name, show_plot, folder)

    def plot_sed_skeleton_detailed(self):
        """
        Initialize and return the skeleton for detailed SED plots.

        Returns
        -------
        fig : matplotlib.figure.Figure
            Matplotlib figure object.
        ax : np.ndarray
            Matplotlib axes array.
        """
        # initialising
        fig, ax = plt.subplots(figsize=(12, 6), nrows=4, ncols=3)
        [axi.set_axis_off() for axi in ax.ravel()]

        ax[0, 0] = fig.add_axes([0.06, 0.43, 0.49, 0.50])
        ax[1, 0] = fig.add_axes([0.06, 0.32, 0.49, 0.11])
        ax[2, 0] = fig.add_axes([0.06, 0.21, 0.49, 0.11])
        ax[3, 0] = fig.add_axes([0.06, 0.10, 0.49, 0.11])

        ax[0, 1] = fig.add_axes([0.63, 0.66, 0.30, 0.28])
        ax[1, 1] = fig.add_axes([0.63, 0.38, 0.30, 0.28])
        ax[2, 1] = fig.add_axes([0.63, 0.10, 0.30, 0.28])

        ax[0, 2] = fig.add_axes([0.92, 0.10, 0.01, 0.28])
        # Labels
        ax[0, 0].set_ylabel('Flux (erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$)')
        ax[1, 0].set_ylabel('FR')
        ax[2, 0].set_ylabel('EWR')
        ax[3, 0].set_ylabel('$χ^2_i$')
        ax[2, 0].set_xlabel('Wavelength (Å)')
        ax[0, 1].set_ylabel('L ($L_⊙$)')
        ax[1, 1].set_ylabel('R ($R_⊙$)')
        ax[2, 1].set_ylabel('$χ^2_r$')
        ax[2, 1].set_xlabel('Te (K)')
        return fig, ax

    def plot_sed_skeleton_public(self):
        """
        Initialize and return the skeleton for public SED plots.

        Returns
        -------
        fig : matplotlib.figure.Figure
            Matplotlib figure object.
        ax : np.ndarray
            Matplotlib axes array.
        """
        # initialising
        fig, ax = plt.subplots(figsize=(6, 6), nrows=3, ncols=1)
        [axi.set_axis_off() for axi in ax.ravel()]

        ax[0] = fig.add_axes([0.10, 0.41, 0.85, 0.52])
        ax[1] = fig.add_axes([0.10, 0.26, 0.85, 0.15])
        ax[2] = fig.add_axes([0.10, 0.11, 0.85, 0.15])
        # Labels
        ax[0].set_ylabel('Flux (erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$)')
        ax[1].set_ylabel('FR')
        ax[2].set_ylabel('EWR')
        ax[2].set_xlabel('Wavelength (Å)')
        return fig, ax

    def update_xylimits_detailed(self, ax):
        """
        Update the x and y limits for public SED plots.

        Parameters
        ----------
        ax : np.ndarray
            Matplotlib axes array.
        """
        # ax range and scales
        wave_min = self.data_all['wavelength'].min()
        wave_max = self.data_all['wavelength'].max()

        for axes in [ax[1, 0], ax[2, 0], ax[3, 0]]:
            axes.set_xlim(ax[0, 0].get_xlim())
        ax[1, 0].set_ylim([-0.2, 1.1])

        T_min = min(self.df_chi2.Te.min(), 4500)
        T_max = max(self.df_chi2.Te.max(), 40000)
        for axes in [ax[0, 1], ax[1, 1], ax[2, 1]]:
            axes.set_xlim(T_max*1.1, T_min/1.4)

        ax[0, 1].set_ylim(self.L/1e3, self.L*1e3)

        for axes in [ax[0, 0], ax[1, 0], ax[2, 0], ax[0, 1], ax[1, 1]]:
            plt.setp(axes.get_xticklabels(), visible=False)
        plt.setp(ax[0, 2].get_yticklabels(which='both'), visible=False)

    def update_xylimits_public(self, ax):
        """
        Update the x and y limits for public SED plots.

        Parameters
        ----------
        ax : np.ndarray
            Matplotlib axes array.
        """
        # ax range and scales
        wave_min = self.data_all['wavelength'].min()
        wave_max = self.data_all['wavelength'].max()

        for axes in [ax[1], ax[2]]:
            axes.set_xlim(ax[0].get_xlim())
        ax[1].set_ylim([-0.2, 1.1])

        T_min = min(self.df_chi2.Te.min(), 4500)
        T_max = max(self.df_chi2.Te.max(), 40000)
        # Hide unnecessary ticks
        plt.setp(ax[0].get_xticklabels(), visible=False)
        plt.setp(ax[1].get_xticklabels(), visible=False)

    def plot_blackbody_model(self, ax=None):
        """
        Plot the blackbody model spectrum.

        Parameters
        ----------
        ax : matplotlib.axes, optional
            Matplotlib axes to plot on. If None, a new subplot is created.
        """
        if ax is None:
            fig, ax = plt.subplots()
        x = np.linspace(self.data_all['log_wavelength'].min(
        ), self.data_all['log_wavelength'].max(), 100)
        y = Fitter.get_blackbody_spectrum_loglog(
            x, self.Te_blackbody, self.log_sf_blackbody)
        ax.plot(10**x, 10**y, color='pink', ls='--',
                label='Blackbody spectrum', zorder=0)
        x = self.data_all['log_wavelength']
        y = Fitter.get_blackbody_spectrum_loglog(
            x, self.Te_blackbody, self.log_sf_blackbody)
        ax.plot(10**x, 10**y, 'r--', label='Blackbody SED (T=%d, log_sf=%.2f)' %
                (self.Te_blackbody, self.log_sf_blackbody), zorder=1)
        ax.loglog()

    def _reduce_Star(self, level=0):
        """
        Reduce the Star object for saving purposes.

        Parameters
        ----------
        level : int, optional
            Reduction level.
            0: Only remove dataframes.
            1: Remove arrays and dataframes.

        Returns
        -------
        star_copy : Star
            Reduced copy of the Star object.

        Notes
        -----
        This method creates a reduced copy of the Star object, primarily for saving to storage formats like pickle files.
        At `level` 0, dataframes (:attr:`Star.df_chi2`, :attr:`Star.df_chi2_trimmed`, :attr:`Star.df_chi2_noisy`) and 
        bigger attributes (:attr:`Star.model`, :attr:`Star.sf_list`, :attr:`Star.sf_list_noisy`) are removed.
        At `level` 1, additional attributes (:attr:`Star.data`, :attr:`Star.data_all`, :attr:`Star.data_not_fitted`) are also removed.
        Lists within the object are converted to strings for efficient saving.
        If any list-like attribute contains more than one element, a ValueError is raised.

        Examples
        --------
        >>> reduced_star = star_instance._reduce_Star(level=1)
        """
        star_copy = copy.deepcopy(self)
        star_copy.time = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
        if level == 0:
            for attrName in ['model', 'sf_list', 'sf_list_noisy', 'df_chi2', 'df_chi2_trimmed', 'df_chi2_noisy']:
                try:  # remove dataframes for the pickle file
                    delattr(star_copy, attrName)
                except:
                    logger.warning('%s attribute not in Star (%s %s)' % (
                        attrName, star_copy.name, star_copy.component))
        if level == 1:
            for attrName in ['model', 'sf_list', 'sf_list_noisy', 'df_chi2', 'df_chi2_trimmed', 'df_chi2_noisy',
                             'data', 'data_all', 'data_not_fitted']:
                try:  # remove arrays and dataframes for the summary file
                    delattr(star_copy, attrName)
                except:
                    logger.warning('%s attribute not in Star (%s %s)' % (
                        attrName, star_copy.name, star_copy.component))

            star_copy.filters_to_drop = ' '.join(
                star_copy.__dict__['filters_to_drop'])
            star_copy.filters_to_drop_all = ' '.join(
                star_copy.__dict__['filters_to_drop_all'])
            star_copy.wavelength_range = ' '.join(
                map(str, star_copy.__dict__['wavelength_range']))
            if isinstance(star_copy.r_limits, list):
                star_copy.r_limits = ' '.join(
                    map(str, star_copy.__dict__['r_limits']))

            # types = [(key, type(star_copy.__dict__[key])) for key in star_copy.__dict__.keys()]
            # print(types)
            # types = [type(star_copy.__dict__[key]) for key in star_copy.__dict__.keys()]
            # print(set(types))
            for key in star_copy.__dict__.keys():
                if type(star_copy.__dict__[key]) == list:
                    if len(star_copy.__dict__[key]) > 1:
                        raise ValueError(
                            'Some parameter (e.g. %s) in the object is list-like. Convert the list into string for saving (similar to "filter_to_drop" in above lines).' % key)
        return star_copy

    def save_summary(self):
        """
        Save a summary of the current state of the instance to a log file.

        Notes
        -----
        - Logs the current state of the instance into a CSV file ('data/log_single_fitting.csv') for future reference.
        - If the log file does not exist, it creates a new one. Otherwise, it appends to the existing log file.
        - Uses the :func:`Star._reduce_Star` method to create a reduced copy of the instance suitable for saving.
        """
        log_file_name = 'data/log_single_fitting.csv'
        _self_copy = self._reduce_Star(level=1)

        df_log = pd.DataFrame(_self_copy.__dict__, index=[0])

        if not os.path.isfile(log_file_name):
            logger.info('Creating %s and saving log' % log_file_name)
        else:
            logger.info('Saving log in %s' % log_file_name)
            df_past_log = pd.read_csv(log_file_name)
            df_log = pd.concat([df_past_log, df_log])

        df_log.to_csv(log_file_name, index=False)


class Model:
    """
    Class to work with synthetic photometry models. 

    - The models are saved as xarray DataArrays.
    - Also includes simple tools to access and plot isochrones and WD cooling curves 
    """

    def __init__(self, name, limits, filter_list=None, star=None):
        """
        Initialize a Model object.

        Parameters
        ----------
        name : str
            Name of the model ('kurucz', 'koester', 'kurucz_uvblue').
        limits : dict
            Dictionary specifying the limits for model parameters.
        filter_list : list, optional
            List of filters to use in creating the filter_list.
        star : Star, optional
            Star object to extract filter_list if not provided.

        Attributes
        ----------
        da : xarray.DataArray
            Multi-dimensional array representing spectral energy distribution (SED) models.
            Dimensions typically include FilterID, Te (effective temperature), logg (surface gravity), etc.
            The specific dimensions depend on the model file read.
        da_obs : xarray.DataArray
            DataArray of observed flux values with dimensions matching :attr:`Model.da`.
        da_obs_error : xarray.DataArray
            DataArray of observed flux error values with dimensions matching :attr:`Model.da`.
        da_obs_error_2percent : xarray.DataArray
            DataArray of observed flux error values (minimum 2%) with dimensions matching :attr:`Model.da`.
        da_obs_error_10percent : xarray.DataArray
            DataArray of observed flux error values (minimum 10%) with dimensions matching :attr:`Model.da`.

        Raises
        ------
        NotImplementedError
            If an unsupported model name is provided.

        Notes
        -----
        - Initializes a Model object with a specified model name, parameter limits, and optional filter list.
        - If a Star object is provided, extracts filter_list from its data; otherwise, uses the provided filter_list.
        - Uses :func:`Model.read_model_file` to load model data from a pre-defined file based on the model name.
        - Initializes :attr:`Model.free_params` based on the number of dimensions in the model.
        """
        if name not in ['kurucz', 'koester', 'kurucz_uvblue']:  # 'uvblue', 'levenhagen'
            raise NotImplementedError('%s model is not supported.' % name)

        self.name = name.lower()
        self.limits = limits
        self.model_file_name = '%s%s_synthetic_photometry.nc' % (
            DIR_MODELS, self.name)

        self.filter_list = filter_list

        if filter_list is not None:
            if len(filter_list) <= 1:
                raise ValueError('The number of filters must be more than 1.')

        if star is not None:
            self.filter_list = star.data.index.values

        self.read_model_file()
        self.constrain_fitting_parameters()
        if star is not None:
            self.create_star_dataarrays(star)
        """    
        free_params = model_dimentions     + 1 
                    = (len(sel.dims) - 1)  + scaling_factor
        Here, for simplicity, the filter dimention is treated as a proxy to the scaling_factor
        """
        _shape = np.array(self.da.shape)
        self.free_params = len(_shape[_shape > 1])

    def read_model_file(self):
        """
        Reads model file from the 'DIR_MODELS' folder as a DataArray.

        Notes
        -----
        *Updates:* :attr:`Model.da`

        - Updates :attr:`Model.da` with the contents of the model file specified by :attr:`Model.model_file_name`.
        - Calls :func:`Model.crop_based_on_filter_list` to adjust :attr:`Model.da` based on the provided filter list.
        - If DIR_MODELS is different from the current path of :attr:`Model.model_file_name`, updates it accordingly.
        """
        # Updating path of the model file if DIR_MODELS has been changed
        if pathlib.Path(self.model_file_name).parent != pathlib.Path(DIR_MODELS):
            self.model_file_name = DIR_MODELS + \
                self.model_file_name.split('/')[-1]

        with xr.open_dataarray(self.model_file_name) as da:
            self.da = da
            self.crop_based_on_filter_list()
        logger.debug(self.da.coords)

    def crop_based_on_filter_list(self):
        """
        Crop the model to only include specified filters.

        Notes
        -----
        *Updates:* :attr:`Model.da`

        Updates self.da.attrs['Wavelengths'] to reflect the wavelengths corresponding
        to the selected filters in self.filter_list.
        """
        if self.filter_list is None:
            return

        # Creating a DataFrame to map filter IDs to wavelengths
        df1 = pd.DataFrame()
        df1['FilterID'] = self.da['FilterID']
        df1 = df1.set_index('FilterID')
        df1['Wavelength'] = self.da.attrs['Wavelengths']

        # Creating a DataFrame for the selected filters
        df2 = pd.DataFrame()
        df2['FilterID'] = self.filter_list
        df2 = df2.set_index('FilterID')
        df2['Wavelength'] = df1['Wavelength']

        # Updating the Wavelengths attribute and selecting the data
        self.da.attrs['Wavelengths'] = df2['Wavelength'].values
        self.da = self.da.sel(FilterID=self.filter_list)

    def constrain_fitting_parameters(self):
        """
        Crop the model DataArray according to specified limits for each parameter dimension.

        Raises
        ------
        ValueError
            If the provided limits do not match all and only the model parameter dimensions,
            or if any of the provided limits exceed the range of the corresponding model parameter.

        Notes
        -----
        *Updates:* :attr:`Model.da`
        """
        limits = self.limits
        dims = np.array(self.da.dims[1:])
        if not np.array_equal(np.sort(list(limits)), np.sort(dims)):
            raise ValueError(
                'Provide limits for all (and only) model parameters %s.' % dims)
        for dim in dims:
            dim_range = [self.da[dim].min().to_numpy().item(),
                         self.da[dim].max().to_numpy().item()]
            if (limits[dim][0] < dim_range[0]) | (limits[dim][1] > dim_range[1]):
                raise ValueError(
                    '%s limits are outside model range %s.' % (dim, dim_range))
        for dim in dims:
            self.da = self.da.sel({dim: slice(limits[dim][0], limits[dim][1])})
        logger.debug(self.da.coords)

    def create_star_dataarrays(self, star):
        """
        Create observed flux and flux error DataArrays matching the dimensions of `self.da`.

        Parameters
        ----------
        star : Star
            Star object containing observed flux and error data.

        Notes
        -----
        *Updates:* :attr:`Model.da_obs`, :attr:`Model.da_obs_error`, :attr:`Model.da_obs_error_2percent`, :attr:`Model.da_obs_error_10percent`

        """

        da_obs = self.da.copy()
        da_obs_error = self.da.copy()
        da_obs_error_2percent = self.da.copy()
        da_obs_error_10percent = self.da.copy()

        for filter_name in star.data.index:
            da_obs = da_obs.where(da_obs.FilterID != filter_name,
                                  star.data['flux'][filter_name])
            da_obs_error = da_obs_error.where(da_obs_error.FilterID != filter_name,
                                              star.data['error'][filter_name])
            da_obs_error_2percent = da_obs_error_2percent.where(da_obs_error_2percent.FilterID != filter_name,
                                                                star.data['error_2percent'][filter_name])
            da_obs_error_10percent = da_obs_error_10percent.where(da_obs_error_10percent.FilterID != filter_name,
                                                                  star.data['error_10percent'][filter_name])

        self.da_obs = da_obs
        self.da_obs_error = da_obs_error
        self.da_obs_error_2percent = da_obs_error_2percent
        self.da_obs_error_10percent = da_obs_error_10percent

        self.da_obs_error.attrs['long_name'] = 'Error'
        self.da_obs_error_2percent.attrs['long_name'] = 'Error (2 precent minimum)'
        self.da_obs_error_10percent.attrs['long_name'] = 'Error (10 precent minimum)'

        logger.debug(da_obs.head(2).to_numpy().squeeze())
        logger.debug(da_obs_error.head(2).to_numpy().squeeze())
        logger.debug(da_obs_error_2percent.head(2).to_numpy().squeeze())
        logger.debug(da_obs_error_10percent.head(2).to_numpy().squeeze())

    @staticmethod
    def load_isochrone():
        """
        Load Parsec ([M/H] = 0) isochrone data from a CSV file and extract subsets based on logAge values.

        Returns
        -------
        iso : pandas.DataFrame
            Complete isochrone data.
        iso_8 : pandas.DataFrame
            Subset of isochrone data for logAge = 8.
        iso_9 : pandas.DataFrame
            Subset of isochrone data for logAge = 9.
        iso_10 : pandas.DataFrame
            Subset of isochrone data for logAge = 10.
        """
        iso = pd.read_csv(DIR_MODELS + 'master_isochrone.csv')
        iso_8 = iso[iso.logAge == 8]
        iso_9 = iso[iso.logAge == 9]
        iso_10 = iso[iso.logAge == 10]
        return iso, iso_8, iso_9, iso_10

    @staticmethod
    def load_wd_cooling_curves():
        """
        Load Bergeron white dwarf cooling curves from a CSV file.

        Returns
        -------
        WD_02 : pandas.DataFrame
            White dwarf cooling curve for mass = 0.2 and spectral type = 'DA'.
        WD_05 : pandas.DataFrame
            White dwarf cooling curve for mass = 0.5 and spectral type = 'DA'.
        WD_13 : pandas.DataFrame
            White dwarf cooling curve for mass = 1.3 and spectral type = 'DA'.
        """
        Bergeron_WD = pd.read_csv(DIR_MODELS + 'master_Bergeron_WD.csv')
        WD_02 = Bergeron_WD[(Bergeron_WD.mass == 0.2) &
                            (Bergeron_WD.spectral_type == 'DA')]
        WD_05 = Bergeron_WD[(Bergeron_WD.mass == 0.5) &
                            (Bergeron_WD.spectral_type == 'DA')]
        WD_13 = Bergeron_WD[(Bergeron_WD.mass == 1.3) &
                            (Bergeron_WD.spectral_type == 'DA')]
        return WD_02, WD_05, WD_13

    @staticmethod
    def plot_isochrone_and_wd(ax=None):
        """
        Plot isochrone and white dwarf cooling curves.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Matplotlib axes to plot on. If None, a new subplot is created.

        Notes
        -----
        - Parsec Isochrones for log(age) of 8, 9, and 10 and [M/H] = 0.
        - Bergeron model white dwarf cooling curves for mass of 0.2, 0.5, 1.3 and spectral type DA.
        """
        if ax is None:
            fig, ax = plt.subplots()
        _, iso_8, iso_9, iso_10 = Model.load_isochrone()
        WD_02, WD_05, WD_13 = Model.load_wd_cooling_curves()
        ax.plot(10**(iso_8.logTe), 10**(iso_8.logL), label='',
                c='0', lw=0.5, rasterized=True, zorder=1)
        ax.plot(10**(iso_9.logTe), 10**(iso_9.logL), label='',
                c='0', lw=0.5, rasterized=True, zorder=1)
        ax.plot(10**(iso_10.logTe), 10**(iso_10.logL), label='',
                c='0', lw=0.5, rasterized=True, zorder=1)
        ax.plot(WD_02.Teff, 10**WD_02.logL, label='', c='0',
                ls=(0, (5, 10)), lw=0.5, rasterized=True, zorder=1)
        ax.plot(WD_05.Teff, 10**WD_05.logL, label='', c='0',
                ls=(0, (5, 10)), lw=0.5, rasterized=True, zorder=1)
        ax.plot(WD_13.Teff, 10**WD_13.logL, label='', c='0',
                ls=(0, (5, 10)), lw=0.5, rasterized=True, zorder=1)


class Fitter:
    """
    Collections of functions which deal with SED fitting. 

    - Fitting blackbody spectrum
    - Functions to calculate chi2 values for a parameter grid.
    - Functions to find best fit values and their errors
    - At present, only chi2 minimization is implimented.
    """
    @staticmethod
    def get_blackbody_spectrum_loglog(log_wavelength, Te, log_scaling_factor):
        """
        Compute the log-log blackbody spectrum flux.

        Adapted from the IUE RDAF 1989.

        Parameters
        ----------
        log_wavelength : float or list
            Logarithm of the wavelength in Angstroms.
        Te : float
            Temperature of the source in Kelvin.
        log_scaling_factor : float
            Logarithmic scaling factor to normalize flux (for a star: log10(radius/distance^2)).

        Returns
        -------
        log_bbflux : float or list
            Logarithm of the blackbody flux in ergs/cm2/s/A for the given temperature and wavelength range.
        """
        c1 = 3.7417749e-5                    # =2*!dpi*h*c*c
        c2 = 1.4387687                       # =h*c/k

        wave = (10**log_wavelength)/1.e8     # angstroms to cm

        bbflux = c1 / (wave**5 * (np.exp(c2/wave/Te)-1.))
        bbflux = bbflux*1.e-8                # ergs/cm2/s/a

        log_bbflux = np.log10(bbflux)+log_scaling_factor

        return log_bbflux

    @staticmethod
    def blackbody(star, p0=[5000., -20]):
        """
        Fit blackbody parameters to observational data using `scipy.optimize.curve_fit`.

        Parameters
        ----------
        star : Star object
            Star data with observed flux.
        p0 : list, optional
            Initial guess for the fit parameters [Temperature, Log Scaling Factor].

        Notes
        -----
        *Updates:* :attr:`Star.Te_blackbody`, :attr:`Star.log_sf_blackbody`
        """
        x = star.data['log_wavelength']
        y = star.data['log_flux']
        err = star.data['e_log_flux']
        nan_filter = np.isfinite(x+y+err)
        popt, pcov = curve_fit(Fitter.get_blackbody_spectrum_loglog,
                               x[nan_filter], y[nan_filter], p0=p0, sigma=err[nan_filter])
        star.Te_blackbody = popt[0]
        star.log_sf_blackbody = popt[1]
        logger.info('Fit parameters: T=%d K, log_sf=%.2f' % tuple(popt))

    @staticmethod
    def plot_chi2_matrix(df, model, alpha=0.3):
        """
        Plot scatter matrix for chi2 DataFrame.

        Parameters
        ----------
        df : DataFrame
            DataFrame containing chi2 values.
        model : Model object
            Model information.
        alpha : float, optional
            Transparency of points in the scatter matrix plot (default is 0.3).
        """
        df['log_chi2'] = np.log10(df['chi2'])
        df['log_R'] = np.log10(df['R'])
        df = df[np.concatenate(
            (np.array(model.da.dims[1:]), ['log_chi2', 'log_R']))]
        pd.plotting.scatter_matrix(df, alpha=alpha, figsize=(8, 8))

    @staticmethod
    def calculate_chi2(star, model, refit=True, plot=False, _trim=None, _minimising_param='chi2'):
        """
        Calculate chi2 values for different model parameters.

        Parameters
        ----------
        star : Star object
            Star data.
        model : Model object
            Model information.
        refit : bool, optional
            Whether to refit the data (default is True).
        plot : bool, optional
            Whether to plot the scatter matrix (default is False).
        _trim : int, optional
            Number of top results to trim the DataFrame (default is None).
        _minimising_param : str, optional
            Parameter over which to determine the best fit ('chi2', 'vgf2', 'vgfb2'). Currently, only 'chi2' is implemented.

        Notes
        -----
        *Updates:* :attr:`Star.df_chi2`, :attr:`Star.df_chi2_trimmed`

        - Calculates chi2 values based on model parameters and observed data.
        - Saves the results into DataFrame `star.df_chi2`.
        - Optionally plots the scatter matrix using `plot_chi2_matrix` if `plot=True`.
        - Trims `star.df_chi2` to top `_trim` results if `_trim` is provided.
        - Saves `star.df_chi2` into a CSV file.
        - Uses `calc_radius` and `calc_luminosity` functions to calculate radius (R) and luminosity (L) based on model parameters.

        Raises
        ------
        NotImplementedError
            If `_minimising_param` is not 'chi2'. ('vgf2' and 'vgfb2' are to be implemented).
        """
        if not os.path.exists(DIR_OUTPUTS+'chi_files/'):
            os.makedirs(DIR_OUTPUTS+'chi_files/')
        chi2_file_name = '%s%s_%s_%s_%s_%s.csv' % (DIR_OUTPUTS+'chi_files/', star.name, star.component,
                                                   model.name, star.run_name, star._type)

        # If the chi2 file exists, it will be read (depends on refit = True or False)
        if os.path.isfile(chi2_file_name):
            if not refit:
                df_chi2 = pd.read_csv(chi2_file_name)
                logger.debug('Loading previous fit (%s)' % chi2_file_name)
                logger.warning(
                    'Give "refit=True" if you want to rerun the fitting process.')
                logger.info(df_chi2.iloc[:5, :-2])
                star.df_chi2 = df_chi2
                star.df_chi2_trimmed = df_chi2.head(5000)
                if plot:
                    Fitter.plot_chi2_matrix(star.df_chi2_trimmed, model)
                return
            if refit:
                logger.info('%s file will be overwritten.' % chi2_file_name)

        shape = np.concatenate((model.da.shape[1:], [len(star.sf_list)]))
        nd_arr = np.full(shape, np.nan)
        _coords = {}
        for dim in model.da.dims[1:]:
            _coords[dim] = model.da[dim]
        _coords['sf'] = star.sf_list
        da_chi2_stacked = xr.DataArray(nd_arr, coords=_coords)
        da_vgf2_stacked = da_chi2_stacked.copy()
        da_vgfb2_stacked = da_chi2_stacked.copy()

        start_time = time.perf_counter()
        for idx in range(len(star.sf_list)):
            if idx == 50:
                estimate_runtime(idx, len(star.sf_list),
                                 start_time, message='Calculating chi2')
            sf = star.sf_list[idx]
            da_residual = model.da_obs - (model.da * sf)

            da_chi2_i = (da_residual/model.da_obs_error)**2
            da_chi2 = da_chi2_i.sum('FilterID', min_count=1)
            da_chi2_stacked.loc[dict(sf=sf)] = da_chi2

            da_vgf2_i = (da_residual/model.da_obs_error_2percent)**2
            da_vgf2 = da_vgf2_i.sum('FilterID', min_count=1)
            da_vgf2_stacked.loc[dict(sf=sf)] = da_vgf2

            da_vgfb2_i = (da_residual/model.da_obs_error_10percent)**2
            da_vgfb2 = da_vgfb2_i.sum('FilterID', min_count=1)
            da_vgfb2_stacked.loc[dict(sf=sf)] = da_vgfb2

        df_chi2 = da_chi2_stacked.to_dataframe(name='chi2').reset_index()
        df_vgf2 = da_vgf2_stacked.to_dataframe(name='vgf2').reset_index()
        df_vgfb2 = da_vgfb2_stacked.to_dataframe(name='vgfb2').reset_index()

        df_chi2['R'] = calc_radius(
            df_chi2['sf'].values, star.distance*u.pc).value
        df_chi2['L'] = calc_luminosity(
            df_chi2['R'].values*u.solRad, df_chi2['Te'].values*u.K).value
        df_chi2['vgf2'] = df_vgf2['vgf2'].values
        df_chi2['vgfb2'] = df_vgfb2['vgfb2'].values
        if _minimising_param != 'chi2':
            raise NotImplementedError(
                'The _minimising_param other than chi2 has not been implemented yet.')
        df_chi2 = df_chi2.sort_values(_minimising_param)

        df_chi2.reset_index(drop=True, inplace=True)
        logger.info(df_chi2.iloc[:5, :-2])

        star.df_chi2 = df_chi2
        if _trim is not None:
            star.df_chi2_trimmed = df_chi2.head(_trim)
            star.df_chi2 = star.df_chi2_trimmed
        else:
            star.df_chi2_trimmed = df_chi2.head(1000)

        logger.debug('Saving %s' % chi2_file_name)
        star.df_chi2.to_csv(chi2_file_name, index=False, header=True, sep=',')

    @staticmethod
    def calculate_noisy_chi2(star, model, total_iterations=100, refit=True, plot=False,
                             _percentile_threshold=10):
        """
        Calculate chi2 values for noise added data.

        Parameters
        ----------
        star : Star object
            Star data.
        model : Model object
            Model information.
        total_iterations : int, optional
            Number of iterations to perform noisy chi2 calculations (default is 100).
        refit : bool, optional
            Whether to refit the data if a previous result exists (default is True).
        plot : bool, optional
            Whether to plot the scatter matrix of noisy chi2 results (default is False).
        _percentile_threshold : int, optional
            Percentile threshold for filtering chi2 values from the original chi2 results (default is 10).

        Notes
        -----
        *Updates:* :attr:`Star.df_chi2_noisy`

        - Uses `star.df_chi2` to filter out top performing parameters based on `_percentile_threshold`.
        - Generates `total_iterations` versions of noisy observed data (`da_obs_noisy`) by adding Gaussian noise to the original flux.
        - Computes chi2 values for each noisy iteration and stores the best fitting parameters in `star.df_chi2_noisy`.
        - Saves `star.df_chi2_noisy` into a CSV file named based on star and model information.

        Raises
        ------
        AttributeError
            If `star.df_chi2` does not exist or `_percentile_threshold` is not within the valid range [0, 100].
        """
        top_chi2 = np.nanpercentile(star.df_chi2.chi2, _percentile_threshold)
        top_df_chi2 = star.df_chi2[star.df_chi2.chi2 < top_chi2]
        sf_min = top_df_chi2.sf.min()
        sf_max = top_df_chi2.sf.max()
        star.sf_list_noisy = 10**np.arange(np.log10(sf_min),
                                           np.log10(sf_max), 0.01)

        chi2_file_name = '%s%s_%s_%s_noisy_%s_%s.csv' % (DIR_OUTPUTS+'chi_files/', star.name, star.component,
                                                         model.name, star.run_name, star._type)
        # If the chi2 file exists, it will be read (depends on 'refit' = True or False)
        if os.path.isfile(chi2_file_name):
            if not refit:
                df_chi2_noisy = pd.read_csv(chi2_file_name)
                logger.debug('Loading previous fit (%s)' % chi2_file_name)
                logger.warning(
                    'Give "refit=True" if you want to rerun the fitting process.')
                logger.info(df_chi2_noisy.iloc[:5, :])
                star.df_chi2_noisy = df_chi2_noisy
                if plot:
                    Fitter.plot_chi2_matrix(df_chi2_noisy, model, alpha=1)
                return
            if refit:
                logger.info('%s file will be overwritten.' % chi2_file_name)

        shape = np.concatenate((model.da.shape[1:], [len(star.sf_list_noisy)]))
        nd_arr = np.full(shape, np.nan)
        _coords = {}
        for dim in model.da.dims[1:]:
            _coords[dim] = model.da[dim]
        _coords['sf'] = star.sf_list_noisy
        da_chi2_stacked = xr.DataArray(nd_arr, coords=_coords)

        da_obs_noisy = model.da.copy()
        da_obs_noisy.attrs['long_name'] = 'Noisy flux'
        df_chi2_noisy = pd.DataFrame()

        start_time = time.perf_counter()
        for seed in range(total_iterations):
            if seed == 2:
                estimate_runtime(seed, total_iterations, start_time,
                                 message='Calculating noisy chi2')
            rng = np.random.default_rng(seed)
            noise = rng.normal(0, 1, [star.N])
            # adding gaussian noise to the flux (except seed 0)
            if seed == 0:
                noisy_flux = star.data['flux']
            else:
                noisy_flux = star.data['flux'] + star.data['error']*noise

            for filter_name in star.data.index:
                da_obs_noisy = da_obs_noisy.where(
                    da_obs_noisy.FilterID != filter_name, noisy_flux[filter_name])

            for idx in range(len(star.sf_list_noisy)):
                sf = star.sf_list_noisy[idx]
                da_residual = da_obs_noisy - (model.da * sf)

                da_chi2_i = (da_residual/model.da_obs_error)**2
                da_chi2 = da_chi2_i.sum('FilterID', min_count=1)
                da_chi2_stacked.loc[dict(sf=sf)] = da_chi2
            df_chi2 = da_chi2_stacked.to_dataframe(name='chi2').reset_index()

            df_chi2 = df_chi2.sort_values('chi2').reset_index(drop=True)
            df_chi2_noisy = pd.concat(
                [df_chi2_noisy, df_chi2.head(1)], ignore_index=True)

        df_chi2_noisy['R'] = calc_radius(
            df_chi2_noisy['sf'].values, star.distance*u.pc).value
        df_chi2_noisy['L'] = calc_luminosity(
            df_chi2_noisy['R'].values*u.solRad, df_chi2_noisy['Te'].values*u.K).value
        logger.info(df_chi2_noisy.iloc[:5, :])
        if plot:
            Fitter.plot_chi2_matrix(df_chi2_noisy, model, alpha=1)
        star.df_chi2_noisy = df_chi2_noisy

        logger.debug('Saving %s' % chi2_file_name)
        df_chi2_noisy.to_csv(chi2_file_name, index=False, header=True, sep=',')

    @staticmethod
    def get_parameters_from_chi2_minimization(star, model):
        """
        Estimate best fit parameters from the least chi2 fit.

        Parameters
        ----------
        star : Star object
            Star data containing chi2 results and observed data.
        model : Model object
            Model information used for fitting.

        Notes
        -----
        *Updates:* :attr:`Star.data_all`, :attr:`Star.data`, :attr:`Star.data_not_fitted`, :attr:`Star.N_Np`, :attr:`Star.chi2_r`, :attr:`Star.vgf2`, :attr:`Star.vgfb2`, :attr:`Star.<dim>`

        This method retrieves and updates parameters from the best chi2 fit found in star.df_chi2.
        It calculates additional statistical metrics for the fitted data.
        It raises a warning if any data points have chi2 values significantly higher than the average (3-sigma threshold).
        """
        fitting_results = {}
        for dim in star.df_chi2.columns:
            value = star.df_chi2[dim][0]
            setattr(star, dim, value)
            fitting_results[dim] = value
        fitting_results = "\n".join("{}\t{}".format(k, v)
                                    for k, v in fitting_results.items())
        logger.info(fitting_results)
        with xr.open_dataarray(model.model_file_name) as da_all:
            da_all = da_all.sel(FilterID=star.data_all.index)
            for dim in model.da.dims[1:]:
                da_all = da_all.sel({dim: getattr(star, dim)})
            df_all = da_all.to_dataframe(name='model_flux')
        star.data_all['model_flux'] = df_all['model_flux']*star.sf
        star.data_all['residual_flux'] = star.data_all['flux'] - \
            star.data_all['model_flux']
        star.data_all['fractional_residual'] = star.data_all['residual_flux'] / \
            star.data_all['flux']
        star.data_all['ewr'] = star.data_all['residual_flux'] / \
            star.data_all['error']
        star.data_all['chi2_i'] = star.data_all['residual_flux']**2 / \
            star.data_all['error']**2
        star.data_all['vgf2_i'] = star.data_all['residual_flux']**2 / \
            star.data_all['error_2percent']**2
        star.data_all['vgfb2_i'] = star.data_all['residual_flux']**2 / \
            star.data_all['error_10percent']**2

        star.data = star.data_all.loc[star.data_all['fitted'] == 1]
        star.data_not_fitted = star.data_all.loc[star.data_all['fitted'] == 0]

        star.N_Np = star.N - model.free_params
        star.chi2_r = star.chi2/star.N_Np
        star.vgf2 = star.data['vgf2_i'].sum()/star.N_Np
        star.vgfb2 = star.data['vgfb2_i'].sum()/star.N_Np

        # Printing the filters with too large chi2 values i.e. 3sigma away from other chi2 values
        abs_zscore = np.abs(zscore(star.data['chi2_i']))
        outliers = star.data[abs_zscore >= 3].index.values
        if len(outliers) > 0:
            _chi_i_values = star.data[abs_zscore >= 3].chi2_i.values
            logger.warning(
                'Based on chi2, I recommend removal of following filters: %s; chi2=%s' % (outliers, _chi_i_values))

    @staticmethod
    def get_realistic_errors_from_iterations(parameter_values, grid_steps, component,
                                             para_type='parameter'):
        """
        Estimates errors using spread in the noisy fits and edge cases.

        - If the best fit parameter is near a boundary:
            - Exaggerates errors and issues a warning.
        - If spread in iterations is less than step size:
            - Keeps errors similar to the step size.
        - Otherwise:
            - Calculates the 16nd, 50th, and 84th percentiles of the parameter_values.
            - Adjusts percentiles to grid_steps boundaries if needed.
            - Computes lower and upper bound errors in the parameter.

        Parameters
        ----------
        parameter_values : list
            List of parameter values obtained from noisy fits.
        grid_steps : list
            List of grid steps available for the parameter.
        component : str
            Component information related to the parameter.
        para_type : str, optional
            Type of the parameter (default is 'parameter').

        Returns
        -------
        para_50 : float
            Median of the parameter_values.
        error_lower : float
            Lower bound error in the parameter.
        error_upper : float
            Upper bound error in the parameter.

        Notes
        -----
        This method estimates realistic errors based on noisy fits and grid steps:
        - It calculates percentiles and adjusts them to grid_steps boundaries if necessary.
        - Issues warnings if the best fit parameter is close to model boundaries.
        - Handles cases where spread in iterations is minimal or zero.

        """
        def find_nearest_with_index(array, value):
            """
            Finds the nearest value available in an array
            """
            array = np.asarray(array)
            idx = (np.abs(array - value)).argmin()
            return idx, array[idx]

        if len(grid_steps) == 1:
            para_50, error_lower, error_upper = np.percentile(
                parameter_values, 50, interpolation='nearest'), np.nan, np.nan
            return para_50, error_lower, error_upper

        para_16 = np.percentile(parameter_values, 15.9, interpolation='lower')
        para_50 = np.percentile(parameter_values, 50.0,
                                interpolation='nearest')
        para_84 = np.percentile(parameter_values, 84.1, interpolation='higher')

        error_lower = para_50 - para_16
        error_upper = para_84 - para_50

        median_index, _ = find_nearest_with_index(grid_steps, para_50)

        if median_index == 0:
            logger.warning(
                '%s_%s : The best fit value is at lower limit of the model.' % (para_type, component))
            if error_upper > 0:
                # Lower error is kept as 3x of upper limit errors
                error_lower = 3 * error_upper
            else:
                error_upper = grid_steps[median_index+1] - \
                    grid_steps[median_index]
                error_lower = 3 * error_upper

        if median_index == len(grid_steps)-1:
            logger.warning(
                '%s_%s : The best fit value is at upper limit of the model.' % (para_type, component))
            if error_lower > 0:
                # upper erro is kept as 3x of lower limit errors
                error_upper = 3 * error_lower
            else:
                error_lower = grid_steps[median_index] - \
                    grid_steps[median_index-1]
                error_upper = 3 * error_lower

        if (median_index > 0) & (median_index < len(grid_steps)-1):
            if error_lower == 0:
                error_lower = grid_steps[median_index] - \
                    grid_steps[median_index-1]
            if error_upper == 0:
                error_upper = grid_steps[median_index+1] - \
                    grid_steps[median_index]

        return para_50, error_lower, error_upper

    @staticmethod
    def get_parameters_from_noisy_chi2_minimization(star, model):
        """
        Estimate best fit parameters from noisy least chi2 fit.

        Parameters
        ----------
        star : Star object
            Star data containing noisy chi2 results.
        model : Model object
            Model information used for fitting.

        Notes
        -----
        *Updates:* :attr:`Star.data_all`, :attr:`Star.data`, :attr:`Star.data_not_fitted`, :attr:`Star.chi2_median`, :attr:`Star.chi2_r_median`, :attr:`Star.vgf2_median`, :attr:`Star.vgfb2_median`, :attr:`Star.<dim>_median`, :attr:`Star.<dim>_error_lower`, :attr:`Star.<dim>_error_upper`

        This method estimates the best fit parameters by calculating the median and errors from the noisy chi2 fits:
        - It uses Fitter.get_realistic_errors_from_iterations() to estimate realistic errors based on noisy fits and grid steps.
        - Updates various attributes in star related to best fit parameters, errors, and fit quality metrics.
        - Logs warnings if the best fit parameters differ from their median value.
        """
        fitting_results = {}
        for dim in star.df_chi2_noisy.columns:
            if dim in model.da.dims[1:]:
                median, error_lower, error_upper = Fitter.get_realistic_errors_from_iterations(
                    star.df_chi2_noisy[dim], model.da[dim].values, component=star.component, para_type=dim)
            elif dim == 'chi2':
                continue
            elif dim == 'sf':
                median, error_lower, error_upper = Fitter.get_realistic_errors_from_iterations(
                    star.df_chi2_noisy[dim], star.sf_list_noisy, component=star.component, para_type=dim)
            elif dim == 'R':
                R_list = calc_radius(star.sf_list_noisy,
                                     star.distance * u.pc).value
                median, error_lower, error_upper = Fitter.get_realistic_errors_from_iterations(
                    star.df_chi2_noisy[dim], R_list, component=star.component, para_type=dim)
                # error after including distance error,     deltaR = R * deltaD/D
                error_lower = np.hypot(
                    error_lower, median * star.e_distance/star.distance)
                error_upper = np.hypot(
                    error_upper, median * star.e_distance/star.distance)
            elif dim == 'L':
                median, _lower, _upper = np.percentile(
                    star.df_chi2_noisy[dim], [50, 31.7, 68.3])
                error_lower, error_upper = median-_lower, _upper-median
                # error after including distance error,     deltaL = L * (2 * deltaD/D)
                error_lower = np.hypot(
                    error_lower, median * 2. * star.e_distance/star.distance)
                error_upper = np.hypot(
                    error_upper, median * 2. * star.e_distance/star.distance)
            value = '%s(-%s,+%s)' % (median, error_lower, error_upper)
            setattr(star, dim+'_median', median)
            setattr(star, dim+'_error_lower', error_lower)
            setattr(star, dim+'_error_upper', error_upper)
            fitting_results[dim] = value

        fitting_results = "\n".join("{}\t{}".format(k, v)
                                    for k, v in fitting_results.items())
        logger.info(fitting_results)

        if star.Te != star.Te_median:
            logger.warning(
                'Te_%s (%d) != Te_median_%s (%d) : Proceed with caution!' % (star.component, star.Te, star.component, star.Te_median))

        with xr.open_dataarray(model.model_file_name) as da_all:
            da_all = da_all.sel(FilterID=star.data_all.index)
            for dim in model.da.dims[1:]:
                da_all = da_all.sel({dim: getattr(star, dim+'_median')})
            df_all = da_all.to_dataframe(name='model_flux_median')

        star.data_all['model_flux_median'] = df_all['model_flux_median']*star.sf_median
        star.data_all['residual_flux_median'] = star.data_all['flux'] - \
            star.data_all['model_flux_median']
        star.data_all['fractional_residual_median'] = star.data_all['residual_flux_median'] / \
            star.data_all['flux']
        star.data_all['ewr_median'] = star.data_all['residual_flux_median'] / \
            star.data_all['error']
        star.data_all['chi2_i_median'] = star.data_all['residual_flux_median']**2 / \
            star.data_all['error']**2
        star.data_all['chi2_i_median'] = star.data_all['residual_flux_median']**2 / \
            star.data_all['error']**2
        star.data_all['vgf2_i_median'] = star.data_all['residual_flux_median']**2 / \
            star.data_all['error_2percent']**2
        star.data_all['vgfb2_i_median'] = star.data_all['residual_flux_median']**2 / \
            star.data_all['error_10percent']**2

        star.data = star.data_all.loc[star.data_all['fitted'] == 1]
        star.data_not_fitted = star.data_all.loc[star.data_all['fitted'] == 0]

        star.chi2_median = star.data['chi2_i_median'].sum()
        star.chi2_r_median = star.chi2_median/star.N_Np
        star.vgf2_median = star.data['vgf2_i_median'].sum()/star.N_Np
        star.vgfb2_median = star.data['vgfb2_i_median'].sum()/star.N_Np


class System:
    """
    Class for working with multi-component SEDs. At present, maximum 3 components are supported.

    Attributes
    ----------
    name : str
        Name of the system.
    distance : float
        Distance to the system [pc].
    e_distance : float
        Error in the distance measurement [pc].
    run_name : str
        Name of the run.
    filters_to_drop : list
        List of filters to drop from the data.
    Total : Star
        Star object representing the total fitted flux.
    components : list
        List of components within the system.
    <component> : Star
        Star object for a given component (e.g. 'A', 'B', 'C' or 'Total').
    Total.data_all : pd.DataFrame
        Data for all filters.
    Total.data : pd.DataFrame
        Data for fitted filters.
    Total.data_not_fitted : pd.DataFrame
        Data for non-fitted filters.
    Total.N_Np : float
        Number of data points minus number of free parameters.
    Total.chi2 : float
        Total chi^2 of the fit.
    Total.chi2_r : float
        Reduced chi^2 of the fit.
    Total.vgf2 : float
        Total vgf^2 of the fit.
    Total.vgfb2 : float
        Total vgfb^2 of the fit.
    """

    def __init__(self, name, distance, e_distance, data, run_name='0', filters_to_drop=[]):
        """
        Initialize a System object.

        Parameters
        ----------
        name : str
            The name of the system.
        distance : float
            The distance to the system [pc].
        e_distance : float
            The error in the distance measurement [pc].
        data : DataFrame
            Input data for the system.
        run_name : str, optional
            Name of the run (default is '0').
        filters_to_drop : list, optional
            List of filters to drop from the data (default is []).

        Notes
        -----
        Initializes a Star object (:attr:`System.Total`) representing the total fitted flux for the system based on the input data.
        """
        self.name = name
        self.distance = distance
        self.e_distance = e_distance
        self.run_name = run_name
        self.filters_to_drop = filters_to_drop

        logger.debug(data)
        self.Total = Star(name=name,
                          distance=distance,
                          e_distance=e_distance,
                          filters_to_drop=filters_to_drop,
                          run_name=run_name,
                          _type='System',
                          data=data.copy(),
                          component='Total')
        self.components = []

    def setup_A_component(self, model, wavelength_range=[0, 1_000_000_000],
                          r_limits='blackbody'):
        """
        Set up the A component of the system.

        Parameters
        ----------
        model : Model object
            Model information for the A component.
        wavelength_range : list, optional
            Range of wavelengths to consider for fitting. Default is [0, 1,000,000,000].
        r_limits : str, optional
            R limits for setting the component. Default is 'blackbody'.

        Raises
        ------
        AttributeError
            If 'A' component is already defined.

        Notes
        -----
        *Updates*: :attr:`System.A`, :attr:`System.components`

        Initializes the A component of the system as a Star object (`self.A`), using the provided model
        information (`model`), wavelength range (`wavelength_range`), and R limits (`r_limits`). Updates
        `self.A` with the appropriate data and attributes, and adds 'A' to the list of components 
        (`self.components`) within the system.
        """
        if 'A' in self.components:
            raise AttributeError("'A' component is already defined.")

        self.A = Star(name=self.name,
                      distance=self.distance,
                      e_distance=self.e_distance,
                      filters_to_drop=self.filters_to_drop,
                      data=self.Total.data_all.copy(),
                      run_name=self.run_name,
                      _type='System',

                      model=model,
                      wavelength_range=wavelength_range,
                      r_limits=r_limits,
                      component='A')
        self.components.append('A')

    def create_residual_star(self, component, model,
                             wavelength_range=[0, 1_000_000_000], r_limits='blackbody'):
        """
        Create a residual star object for a specific component based on residual flux.

        Parameters
        ----------
        component : str
            Identifier for the component for which residual star is to be created.
        model : Model object
            Model information for the residual star.
        wavelength_range : list, optional
            Range of wavelengths to consider for fitting. Default is [0, 1,000,000,000].
        r_limits : str, optional
            R limits for setting the component. Default is 'blackbody'.

        Raises
        ------
        AttributeError
            If the `component` is already defined in the System.

        Notes
        -----
        *Updates:* :attr:`System.Total`, :attr:`System.<component>`, :attr:`System.components`

        This method computes residual flux from `self.Total.data_all` by subtracting flux from all fitted components present in the System. 
        """
        if component in self.components:
            raise AttributeError(
                "'%s' component is already defined." % component)
        self.update_Total_flux()
        data_res = self.Total.data_all[[
            'wavelength', 'residual_flux', 'error']]
        data_res.columns = ['wavelength', 'flux', 'error']
        star_res = Star(name=self.name,
                        distance=self.distance,
                        e_distance=self.e_distance,
                        filters_to_drop=self.filters_to_drop,
                        wavelength_range=wavelength_range,
                        run_name=self.run_name,
                        _type='System',

                        data=data_res,
                        model=model,
                        r_limits=r_limits,
                        component=component)

        setattr(self, component, star_res)
        self.components.append(component)

    def update_Total_flux(self, median=False):
        """
        Update the self.Total flux based on all fitted component fluxes.

        Parameters
        ----------
        median : bool, optional
            Use median parameters for calculations.

        Notes
        -----
        *Updates:* :attr:`System.Total.data_all`, :attr:`System.Total.data`, :attr:`System.Total.data_not_fitted`, :attr:`System.Total.N_Np`, :attr:`System.Total.chi2`, :attr:`System.Total.chi2_r`, :attr:`System.Total.vgf2`, :attr:`System.Total.vgfb2`
        """
        for idx, component in enumerate(self.components):
            star = getattr(self, component)
            if idx == 0:
                if median:
                    self.Total.data_all['model_flux'] = star.data_all['model_flux_median']
                else:
                    self.Total.data_all['model_flux'] = star.data_all['model_flux']
                free_params = star.N - star.N_Np
            else:
                if median:
                    self.Total.data_all['model_flux'] += star.data_all['model_flux_median']
                else:
                    self.Total.data_all['model_flux'] += star.data_all['model_flux']
                free_params += star.N - star.N_Np

        self.Total.data_all['residual_flux'] = self.Total.data_all['flux'] - \
            self.Total.data_all['model_flux']
        self.Total.data_all['fractional_residual'] = self.Total.data_all['residual_flux'] / \
            self.Total.data_all['flux']
        self.Total.data_all['ewr'] = self.Total.data_all['residual_flux'] / \
            self.Total.data_all['error']
        self.Total.data_all['chi2_i'] = self.Total.data_all['residual_flux']**2 / \
            self.Total.data_all['error']**2
        self.Total.data_all['vgf2_i'] = self.Total.data_all['residual_flux']**2 / \
            self.Total.data_all['error_2percent']**2
        self.Total.data_all['vgfb2_i'] = self.Total.data_all['residual_flux']**2 / \
            self.Total.data_all['error_10percent']**2

        self.Total.data_all['fitted'] = star.data_all['fitted'] == 1
        self.Total.data = self.Total.data_all.loc[self.Total.data_all['fitted'] == 1]
        self.Total.data_not_fitted = self.Total.data_all.loc[self.Total.data_all['fitted'] == 0]

        self.Total.N_Np = self.Total.N - free_params
        self.Total.chi2 = self.Total.data['chi2_i'].sum()
        self.Total.chi2_r = self.Total.chi2/self.Total.N_Np
        self.Total.vgf2 = self.Total.data['vgf2_i'].sum()/self.Total.N_Np
        self.Total.vgfb2 = self.Total.data['vgfb2_i'].sum()/self.Total.N_Np

    def plot(self, add_noisy_seds=False, show_plot=True, folder=None, FR_cutoff=0.5, _median=False):
        """
        Plot detailed information about the system.

        Parameters
        ----------
        add_noisy_seds : bool, optional
            Include noisy SEDs in the plot. Default is False.
        show_plot : bool, optional
            Whether to display the plot. Default is True.
        folder : str, optional
            Folder to save the plot. No plot is saved if `folder` is `None`.
        FR_cutoff : float, optional
            Fractional residual cutoff for plotting. Default is 0.5.
        _median : bool, optional
            Use median for calculations. Default is False.
        """
        fig, ax = self.Total.plot_sed_skeleton_detailed()

        title = '%s [%s]' % (self.name, self.run_name)
        for idx, component in enumerate(self.components):
            star = getattr(self, component)
            if star.component == 'A':
                color, ls = 'C1', '--'
            if star.component == 'B':
                color, ls = 'C0', '-.'
            if star.component == 'C':
                color, ls = 'C4', ':'
            star.plot_sed_fitted(ax[0, 0], color=color, ls=ls, median=_median)
            if add_noisy_seds:
                star.plot_noisy_seds(ax=ax[0, 0], color=color, ls=ls)

            star.plot_fractional_residual(
                ax[1, 0], color=color, ls=ls, FR_cutoff=FR_cutoff, median=_median)
            star.plot_ewr(ax[2, 0], color=color, ls=ls, median=_median)
            star.plot_chi2(ax[3, 0], color=color, ls=ls, median=_median)

            star.plot_hrd(ax[0, 1], color=color, median=_median, cax=ax[0, 2])
            star.plot_temperature_radius(
                ax[1, 1], color=color, median=_median, cax=ax[0, 2])
            star.plot_temperature_chi2(
                ax[2, 1], color=color, median=_median, cax=ax[0, 2])

            star.update_xylimits_detailed(ax)
            title += '\n%s (%s):' % (component, star.model.name)
            dims = np.concatenate([star.df_chi2.columns[:-6], ['R', 'L']])
            for dim in dims:
                title += '  %s = %.2E' % (dim, getattr(star, dim))
                if hasattr(star, dim+'_error_upper'):
                    ten_factor = 10**np.floor(np.log10(getattr(star, dim)))
                    title += '$^{+%.2f}_{-%.2f}$' % (getattr(star, dim+'_error_upper')/ten_factor,
                                                     getattr(star, dim+'_error_lower')/ten_factor)
        title += '\n $χ^2_r$ %.1f, vgfb$^2$ %.2f' % (
            star.chi2_r, star.vgfb2)
        self.update_Total_flux()
        # Don't know why the next line changes xlim
        ax[0, 0].set_xlim(ax[0, 0].get_xlim())
        self.Total.plot_sed_obs(ax[0, 0])
        self.Total.plot_sed_fitted(ax[0, 0], color='g', ls='-')
        self.Total.plot_fractional_residual(ax[1, 0], color='g', ls='-')
        self.Total.plot_ewr(ax[2, 0], color='g', ls='-')
        self.Total.plot_chi2(ax[3, 0], color='g', ls='-')
        ax[0, 0].grid()
        ax[0, 0].legend()
        ax[0, 1].legend()
        ax[0, 0].set_title(title, loc='left')
        if folder is not None:
            plot_name = folder+'%s_system_%s' % (self.name, self.run_name)
            save_fig(plot_name, show_plot, folder)

    def plot_public(self, add_noisy_seds=False, median=False, show_plot=True, folder=None):
        """
        Plot summarised information about the system.

        Parameters
        ----------
        add_noisy_seds : bool, optional
            Include noisy SEDs in the plot. Default is False.
        median : bool, optional
            Use median for calculations. Default is False.
        show_plot : bool, optional
            Whether to display the plot. Default is True.
        folder : str, optional
            Folder to save the plot. No plot is saved if `folder` is `None`.
        """
        fig, ax = self.Total.plot_sed_skeleton_public()
        for idx, component in enumerate(self.components):
            star = getattr(self, component)
            if star.component == 'A':
                color, ls = 'C1', '--'
            if star.component == 'B':
                color, ls = 'C0', '-.'
            if star.component == 'C':
                color, ls = 'C4', ':'
            star.plot_sed_fitted(ax[0], color=color, ls=ls, median=median)

            if add_noisy_seds:
                star.plot_noisy_seds(ax=ax[0], color=color, ls=ls)

            if idx == 0:
                star.plot_fractional_residual(
                    ax[1], color=color, ls=ls, median=median)
                star.plot_ewr(ax[2], color=color, ls=ls, median=median)
                label = '   %d K' % star.Te
            else:
                label += ' + %d K' % star.Te
            star.update_xylimits_public(ax)

        self.update_Total_flux()
        ax[0].set_xlim(ax[0].get_xlim())
        # Don't know why the next line changes xlim
        self.Total.plot_sed_obs(ax[0])
        self.Total.plot_sed_fitted(ax[0], color='g', ls='-')
        self.Total.plot_fractional_residual(ax[1], color='g', ls='-')
        self.Total.plot_ewr(ax[2], color='g', ls='-')
        ax[0].grid()
        ax[0].legend()
        ax[0].set_title(self.name + label)

        if folder is not None:
            plot_name = folder + \
                '%s_system_%s_public' % (self.name, self.run_name)
            save_fig(plot_name, show_plot, folder)

    def save_summary(self):
        """
        Save a summary of the current state of the instance to a CSV file.

        Notes
        -----
        If the log file (`data/log_starsystem_fitting.csv`) does not exist, creates a new one. Otherwise, appends to the existing log file.
        """
        log_file_name = 'data/log_starsystem_fitting.csv'
        dict = {'run_name': self.run_name}
        for idx, component in enumerate(self.components):
            star = getattr(self, component)
            _dict = star._reduce_Star(level=1).__dict__
            _dict = {key+'_%s' %
                     (component): val for key, val in _dict.items()}
            dict = dict | _dict
        df_log = pd.DataFrame(dict, index=[0])

        if not os.path.isfile(log_file_name):
            logger.info('Creating %s and saving log' % log_file_name)
        else:
            logger.info('Saving log in %s' % log_file_name)
            df_past_log = pd.read_csv(log_file_name)
            df_log = pd.concat([df_past_log, df_log])

        df_log.to_csv(log_file_name, index=False)
