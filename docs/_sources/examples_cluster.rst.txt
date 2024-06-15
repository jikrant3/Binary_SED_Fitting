Example cluster fit
===================

⚠️🔴 ONLY FOR ADVANCE USERS 🔴⚠️ 

Example of how multiple stars in a cluster can be fitted together 

- Create a star_param dictionary which stores the necessary parameters (e.g. filters_to_drop, wavelength_range_A, model_name, r_limit…)

- Plots are not shown in the notebook (``show_plot=False``) to reduce clutter

Cluster and stellar parameters
------------------------------

.. code:: ipython3

    import src.binary_sed_fitting as bsf
    import importlib
    importlib.reload(bsf)
    import warnings
    warnings.filterwarnings("ignore")
    ################################################################################
    distance    = 831. # pc
    e_distance  = 11.  # pc

.. code:: ipython3

    optical_IR_wave_range = [3000,1_000_000_000]
    all_wave_range = [1,1_000_000_000]
    star_params = {'WOCS2002':{'filters_to_drop':['KPNO/Mosaic.I'],
                               'wavelength_range_A':optical_IR_wave_range,
                               'wavelength_range_B':all_wave_range,
                               },
                   'Y1168':{'filters_to_drop':['2MASS/2MASS.J', '2MASS/2MASS.H', '2MASS/2MASS.Ks',
                                                  'WISE/WISE.W1','WISE/WISE.W2'],
                               'wavelength_range_A':[2000,10_000_000],
                               'wavelength_range_B':all_wave_range,
                               },
                   'WOCS1007'   :{'filters_to_drop':[],
                               'wavelength_range_A':optical_IR_wave_range,
                               'wavelength_range_B':all_wave_range,
                               }}

Fitting single SEDs
-------------------

.. code:: ipython3

    bsf.console.setLevel(bsf.logging.WARNING)
    
    refit = False
    run_name = '0'
    
    for name in star_params.keys():
        file_name   = 'data/extinction_corrected_flux_files/%s.csv'%name
        data        = bsf.load_data(file_name, mode='csv')
        ################################################################################
        if name == 'Y1168':
            model_name = 'koester'
            limits = {'Te'   : [5000, 80000],
                      'logg' : [ 6.5,   9.5]}
        else:
            model_name = 'kurucz'
            limits = {'Te'   : [3500, 9000],
                      'logg' : [   3,    5],
                      'MH'   : [ 0.0,  0.0],
                      'alpha': [ 0.0,  0.0]}
        model = bsf.Model(model_name, limits=limits)
        ################################################################################
        star = bsf.Star(name=name, 
                        distance=distance, 
                        e_distance=e_distance,
                        filters_to_drop=star_params[name]['filters_to_drop'], 
                        wavelength_range=star_params[name]['wavelength_range_A'],
                        data=data, 
                        model=model, 
                        r_limits='blackbody',     
                        run_name=run_name)
        ################################################################################
        star.fit_chi2(refit=refit)
        star.fit_noisy_chi2(refit=refit)
        ################################################################################
        star.plot_public(add_noisy_seds=False,
                         show_plot=False, 
                         folder='plots/M67/',
                         FR_cutoff=0.5)
        star.save_summary()


.. parsed-literal::

    22:37:55 ----- WARNING  ----- calculate_chi2
    Give "refit=True" if you want to rerun the fitting process.
    22:37:55 ----- WARNING  ----- get_parameters_from_chi2_minimization
    Based on chi2, I recommend removal of following filters: ['GAIA/GAIA3.Grvs']; chi2=[19.2204403]
    22:37:55 ----- WARNING  ----- calculate_noisy_chi2
    Give "refit=True" if you want to rerun the fitting process.
    22:37:55 ----- WARNING  ----- calculate_chi2
    Give "refit=True" if you want to rerun the fitting process.
    22:37:55 ----- WARNING  ----- calculate_noisy_chi2
    Give "refit=True" if you want to rerun the fitting process.
    22:37:55 ----- WARNING  ----- get_realistic_errors_from_iterations
    logg_A : The best fit value is at upper limit of the model.
    22:37:56 ----- WARNING  ----- calculate_chi2
    Give "refit=True" if you want to rerun the fitting process.
    22:37:56 ----- WARNING  ----- calculate_noisy_chi2
    Give "refit=True" if you want to rerun the fitting process.
    

Identifying stars with UV excess
--------------------------------

-  Manually

.. code:: ipython3

    stars_with_uv_excess = ['WOCS2002', 'WOCS1007']

Fitting double SEDs
-------------------

.. code:: ipython3

    bsf.console.setLevel(bsf.logging.WARNING)
    
    refit = False
    run_name = '0'
    
    for name in stars_with_uv_excess:
        file_name   = 'data/extinction_corrected_flux_files/%s.csv'%name
        data        = bsf.load_data(file_name, mode='csv')
        ################################################################################
        model_name = 'kurucz_uvblue'
        limits = {'Te'   : [3500, 9000],
                  'logg' : [   3,    5],
                  'MH'   : [ 0.0,  0.0]}
        model_A = bsf.Model(model_name, limits=limits)
    
        model_name = 'koester'
        limits = {'Te'   : [5000, 80000],
                  'logg' : [ 6.5,  9.5]}
        model_B = bsf.Model(model_name, limits=limits)
        ################################################################################
        star_system = bsf.System(name=name,
                                distance=distance,
                                e_distance=e_distance,
                                data=data,
                                run_name=run_name,
                                filters_to_drop=star_params[name]['filters_to_drop'])
        ################################################################################
        star_system.setup_A_component(model=model_A, 
                                    wavelength_range=star_params[name]['wavelength_range_A'],
                                    r_limits='blackbody')
    
        star_system.A.fit_chi2(refit=refit)
        star_system.A.fit_noisy_chi2(refit=refit)
        ################################################################################
        star_system.create_residual_star(component='B', 
                                        model=model_B,
                                        wavelength_range=star_params[name]['wavelength_range_B'], 
                                        r_limits=[0.001,1.0])
        star_system.B.fit_chi2(refit=refit)
        star_system.B.fit_noisy_chi2(refit=refit)
        ################################################################################
        star_system.plot(add_noisy_seds=False, 
                         FR_cutoff=0.5,
                         folder='plots/',
                         show_plot=False)
        star_system.save_summary()


.. parsed-literal::

    22:37:56 ----- WARNING  ----- calculate_chi2
    Give "refit=True" if you want to rerun the fitting process.
    22:37:56 ----- WARNING  ----- calculate_noisy_chi2
    Give "refit=True" if you want to rerun the fitting process.
    22:37:56 ----- WARNING  ----- get_realistic_errors_from_iterations
    logg_A : The best fit value is at upper limit of the model.
    22:37:56 ----- WARNING  ----- calculate_chi2
    Give "refit=True" if you want to rerun the fitting process.
    22:37:56 ----- WARNING  ----- get_parameters_from_chi2_minimization
    Based on chi2, I recommend removal of following filters: ['GALEX/GALEX.NUV']; chi2=[2582.21594585]
    22:37:56 ----- WARNING  ----- calculate_noisy_chi2
    Give "refit=True" if you want to rerun the fitting process.
    22:37:56 ----- WARNING  ----- get_realistic_errors_from_iterations
    logg_B : The best fit value is at upper limit of the model.
    22:37:56 ----- WARNING  ----- get_parameters_from_noisy_chi2_minimization
    Te_B (14750) != Te_median_B (14500) : Proceed with caution!
    22:37:57 ----- WARNING  ----- calculate_chi2
    Give "refit=True" if you want to rerun the fitting process.
    22:37:57 ----- WARNING  ----- calculate_noisy_chi2
    Give "refit=True" if you want to rerun the fitting process.
    22:37:57 ----- WARNING  ----- get_realistic_errors_from_iterations
    logg_A : The best fit value is at lower limit of the model.
    22:37:58 ----- WARNING  ----- calculate_chi2
    Give "refit=True" if you want to rerun the fitting process.
    22:37:58 ----- WARNING  ----- get_parameters_from_chi2_minimization
    Based on chi2, I recommend removal of following filters: ['GALEX/GALEX.NUV']; chi2=[754.58954018]
    22:37:58 ----- WARNING  ----- calculate_noisy_chi2
    Give "refit=True" if you want to rerun the fitting process.
    22:37:58 ----- WARNING  ----- get_realistic_errors_from_iterations
    Te_B : The best fit value is at upper limit of the model.
    
