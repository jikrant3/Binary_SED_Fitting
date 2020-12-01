#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
from itertools import product
import linecache
import scipy.stats as st
import os
import sys

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import ticker

###############################################################################
def get_number_from_string(line): 
    # Read a float number from sentence
    for t in line.split():
        try:
            l = float(t)
        except ValueError:
            pass
    return l


def read_A_comp_from_VOSA(STAR_NAME):
    # Reading single component fitted data from VOSA
    file_A = DIR_OBS + STAR_NAME +'/bestfitp/'+ STAR_NAME +'.bfit.phot.dat'
    flux_A = pd.read_csv(file_A, engine='python', comment='#',delim_whitespace= True, skipinitialspace=True, header=None)
    # Reading the column names and fit parameters automatically
    flux_A_col = pd.read_csv(file_A, engine='python',delim_whitespace= True,skipinitialspace=True,skiprows=39, nrows=1,escapechar='#', header=None)
    flux_A_col=flux_A_col.drop([0], axis=1)
    flux_A.columns = flux_A_col.values[0]
    flux_A = flux_A.set_index('FilterID')
    
    Teff_A = get_number_from_string(linecache.getline(file_A, 9))
    logg_A = get_number_from_string(linecache.getline(file_A, 10))
    sf_A = get_number_from_string(linecache.getline(file_A, 13))
    L_A = get_number_from_string(linecache.getline(file_A, 37))
    print('A-Component fit parameters: T=%d, logg=%f, SF=%.2e, L=%f' %(Teff_A, logg_A, sf_A, L_A))  
    print ('Total filters: %d' %len(flux_A.index.values))           # Check if the filters are same as required
    return flux_A, Teff_A,logg_A,logg_A,sf_A,L_A
    '''
    flux_A_col.values[0]     # The column names in the file are as follows: 
            'Filter', 'Wave', 'ObsFlux', 'ObsErr', 'CorFlux','CorErr', 'ModFlux', 'Fitted', 'Excess', 'FitExc','UpLim'
            We are interested in 
                      CorFlux, CorErr (extinction corrected observed flux)
                      ModFlux (Model Flux of the cooler component fitted by VOSA)
            Residual will be 'CorFlux-ModFlux'
            This residual will be fitted with B-component
    '''

def drop_filters(STAR_NAME,flux_A):
    '''
    Some filters will have to be removed while fitting.
        In case of IR excess, you can remove Wise filters,
        or you want to specifically remove something, you can do it here.
    e.g. Removing I,W3,W4 filter in WOCS2002 and none in BSS2
    '''    
    # Select a set of filters which were Fitted or given excess
    # Removing points with upper limits or FitExcess
    flux_A['UpLim'] = flux_A['UpLim'].replace(['---'],0)
    flux_A['UpLim'] = flux_A['UpLim'].replace(['1'],1)
    flux_A['FitExc'] = flux_A['FitExc'].replace(['---'],0)
    flux_A['FitExc'] = flux_A['FitExc'].replace(['1'],1)
    flux_A['to_be_fitted'] = flux_A.Fitted+flux_A.Excess+flux_A.FitExc+flux_A.UpLim

    not_fitted_A = flux_A[(flux_A['to_be_fitted']!=1)]
    flux_A = flux_A[(flux_A['to_be_fitted']==1)]

    if (STAR_NAME == 'WOCS2002'):
        for filter_name in ['GALEX/GALEX.FUV','WISE/WISE.W3']:
            not_fitted_A = pd.concat([not_fitted_A, flux_A[(flux_A.index==filter_name)]])
            flux_A=flux_A.drop(index=filter_name)
    if (STAR_NAME=='BSS2'):
        pass
    
    # printing filters to be fitted and not_fitted
    _t1, _t2 = pd.DataFrame(),pd.DataFrame()    
    _t1['to_be_fitted'] = flux_A.index.values
    _t1['Wavelength'] = flux_A.Wavelength.values
    _t1 = _t1.set_index('Wavelength')
    _t2['not_fitted_A'] = not_fitted_A.index.values
    _t2['Wavelength'] = not_fitted_A.Wavelength.values
    _t2 = _t2.set_index('Wavelength')
    _filter_table = pd.concat([_t1,_t2],sort=True)    
    print (_filter_table.sort_index().fillna(''))

    N_points = len(flux_A)
    N_Np = N_points-FREE_PARA
    return flux_A,not_fitted_A, N_points, N_Np
    
def read_model_file(model, logg_B, Z_B):
    # loading the synthetic flux file
    if model =='Koe':
        flux_model = pd.read_csv('models/Koe_logg'+logg_B+'.csv', engine='python')
    if model == 'Kr':
        flux_model = pd.read_csv('models/Kr_logg'+logg_B+'_Z'+Z_B+'.csv', engine='python')
    if model == 'Tlusty':
        flux_model = pd.read_csv('models/Tlusty_logg'+logg_B+'.csv', engine='python')

    flux_model = flux_model.sort_values('wave')
    flux_model = flux_model.set_index('filter')
    return flux_model


def add_A_comp_to_model(flux_A,flux_model):
    # combining observed and model fluxes for easier intercomparison
    # Adding the observed corrected flux and errors to the combined file
    flux_model['CorFlux']=flux_A['Flux']
    flux_model['CorErr']=flux_A['Error']
    flux_model['flux_A']=flux_A['FluxMod']

    #Replacing zeros in errors with 110% of max error
    flux_model['CorErr_frac']=flux_model['CorErr']/flux_model['CorFlux']
    flux_model['CorErr'] = flux_model['CorErr'].replace(0, flux_model['CorFlux']*(flux_model['CorErr_frac'].max()+0.1)) 

    # Removing rows without datapoints in filters
    flux_model = flux_model[flux_model['CorFlux'].notna()]
    return flux_model


def initializing_sf_T(mode, sf_min, sf_max):
    # Initialising the scaling factor grid and temperature grid for different modes
    if mode=='test_Kr_Koe': # Testing of working of Kurucz (already fitted)+Koster(to be fitted) model binary.
        scale_list = np.geomspace(sf_min, sf_max, num=3)
        temp_B_list = ['10000','22000','30000']
    if mode=='rough_Kr_Koe': # rough fitting of Kurucz (already fitted)+Koster(to be fitted) model binary. 
        scale_list = np.geomspace(sf_min, sf_max, num=20)
        temp_B_list = ['05000','06000','07000','08000','09000','10000','11000','12000','13000','14000','15000','16000','17000','18000','19000','20000','21000','22000','23000','24000','25000','26000','27000','28000','29000','30000','35000','40000','50000','60000','70000','80000']
    if mode=='finer_Kr_Koe': # finer fitting of Kurucz (already fitted)+Koster(to be fitted) model binary.
        scale_list = np.geomspace(sf_min, sf_max, num=100)
        temp_B_list = ['05000','05250','05500','05750','06000','06250','06500','06750','07000','07250','07500','07750','08000','08250','08500','08750','09000','09250','09500','09750','10000','10250','10500','10750','11000','11250','11500','11750','12000','12250','12500','12750','13000','13250','13500','13750','14000','14250','14500','14750','15000','15250','15500','15750','16000','16250','16500','16750','17000','17250','17500','17750','18000','18250','18500','18750','19000','19250','19500','19750','20000','21000','22000','23000','24000','25000','26000','27000','28000','29000','30000','32000','34000','35000','36000','38000','40000','45000','50000','60000','70000','80000']
    if mode=='rough_Kr_Kr': # rough fitting of Kurucz (already fitted)+Kurucz(to be fitted) model binary.
        scale_list = np.geomspace(sf_min, sf_max, num=20)
        temp_B_list = ['5000','6000','7000','8000','9000','10000','11000','12000','13000','14000','15000','16000','17000','18000','19000','20000','21000','22000','23000','24000','25000','26000','27000','28000','29000','30000','32000','34000','36000','38000']
    if mode=='finer_Kr_Kr': # finer fitting of Kurucz (already fitted)+Kurucz(to be fitted) model binary.
        scale_list = np.geomspace(sf_min, sf_max, num=100)
        temp_B_list = ['5000','5250','5500','5750','6000','6250','6500','6750','7000','7250','7500','7750','8000','8250','8500','8750','9000','9250','9500','9750','10000','10250','10500','10750','11000','11250','11500','11750','12000','12250','12500','12750','13000','14000','15000','16000','17000','18000','19000','20000','21000','22000','23000','24000','25000','26000','27000','28000','29000','30000','31000','32000','33000','34000','35000','36000','37000','38000','39000']
    if mode=='test_Kr_Kr': # Testing of working of Kurucz (already fitted)+Kurucz(to be fitted) model binary.
        scale_list = np.geomspace(sf_min, sf_max, num=3)
        temp_B_list = ['5000','5250','5500']
    if mode=='rough_Kr_Tlusty': # rough fitting of Kurucz (already fitted)+Tlusty(to be fitted) model binary.
        scale_list = np.geomspace(sf_min, sf_max, num=20)
        temp_B_list = ['15000','17000','19000','21000','23000','25000','27000','29000','30000','32500' ,'35000','37500','40000','42500','45000','47500','50000','52500','55000']
    if mode=='finer_Kr_Tlusty': # finer fitting of Kurucz (already fitted)+Tlusty(to be fitted) model binary.
        scale_list = np.geomspace(sf_min, sf_max, num=100)
        temp_B_list = ['15000','16000','17000','18000','19000','20000','21000','22000','23000','24000','25000','26000','27000','28000','29000','30000','32500' ,'35000','37500','40000','42500','45000','47500','50000','52500','55000']
    if mode=='test_Kr_Tlusty': # Testing of working of Kurucz (already fitted)+Tlusty(to be fitted) model binary.
        scale_list = np.geomspace(sf_min, sf_max, num=3)
        temp_B_list = ['15000','16000','17000']
#     print  ('len(SF),  len(T),  len(T)*len(SF) = ','',len(scale_list),len(temp_B_list) , len(temp_B_list)*len(scale_list))
    # Return scale factor list and temperature list
    return scale_list,temp_B_list    


def calculating_chi2(scale_list, temp_B_list, fitting_required, mode):
    ChiSqr_File_B = 'outputs/chi_files/'+STAR_NAME+'_ChiSqur_logg_B' + str(logg_B) + '_' +Z_B+'_'+ mode +'_' +model+'_'+str(cycle)  + '.csv'
    
    if fitting_required == 0:    # Read chi2 from previously generated file
        data_chi = pd.read_csv(ChiSqr_File_B, engine='python', dtype={'Temp':'object','SF':'float','ChiSqr':'float'})
        
    elif fitting_required == 1:
        data_chi = pd.DataFrame()
        idx=0
        count = len(scale_list)*len(temp_B_list)
        for counter, (sf, temp) in enumerate(product(scale_list, temp_B_list),1):
            flux_model['Total'] = flux_model['flux_A'] + sf*flux_model[temp]
            flux_model['ResidualFlux']=flux_model['CorFlux'] - flux_model['Total']
            flux_model['chisqu']= flux_model['ResidualFlux']**2 / flux_model['CorErr']**2
            Sumchi = flux_model['chisqu'].sum()

            data_chi.loc[counter, 'Temp'] = temp
            data_chi.loc[counter, 'SF'] = sf
            data_chi.loc[counter, 'ChiSqr'] = Sumchi

            if (counter%100==0):  # printing progress
                print ('\r Calculating chi2:  %d/%d  (%d%%)' %(counter, count, 100*counter/count), end='')
                
        data_chi=data_chi.sort_values('ChiSqr')
        data_chi = data_chi.reset_index(drop=True)
        if not os.path.exists('outputs/chi_files/'):
            os.makedirs('outputs/chi_files/')
        data_chi.to_csv(ChiSqr_File_B, index=False, header=True, sep=',')
    print ('\r',data_chi.head())
    return data_chi
   
    
def minimizing_chi2(flux_model,data_chi):  
    # Getting best fit parameters from least chi2 fit
    Teff_B = data_chi.Temp[0]#.astype(np.float16)[0]
    sf_B = data_chi.SF[0]
    R_B = radius(sf_B)
    L_B = lumi(R_B,Teff_B)
    flux_model['flux_B']=sf_B * flux_model[str(Teff_B)]
    flux_model['Total'] = flux_model['flux_B']+flux_model['flux_A']
    flux_model['ResidualFlux']=(flux_model['CorFlux'] - flux_model['Total'])
    flux_model['chisqu']= flux_model['ResidualFlux']**2 / flux_model['CorErr']**2
    return Teff_B,sf_B, R_B, L_B

    
def create_plots(mode, save=1):
    iso = pd.read_csv('data/example_isochrone.txt',engine='python',delimiter= ',', header=0)
    ###################### initialising
    label_A = 'A (' + str(Teff_A) + ' K, logg=' + str(logg_A) + ')'
    label_B = 'B (' + str(Teff_B) + ' K, logg=' + str(logg_B)+', R='+str(round(R_B,3)) + ', L='+str(round(L_B,3)) + ')'
    f, axes = plt.subplots(figsize=(12,6),nrows = 3, ncols = 3)
    [axi.set_axis_off() for axi in axes.ravel()]
    axes[0][0] = f.add_axes([0.06, 0.44, 0.49, 0.50])
    axes[1][0] = f.add_axes([0.06, 0.27, 0.49, 0.17])
    axes[2][0] = f.add_axes([0.06, 0.10, 0.49, 0.17])

    axes[0][1] = f.add_axes([0.63, 0.66, 0.30, 0.28])
    axes[1][1] = f.add_axes([0.63, 0.38, 0.30, 0.28])
    axes[2][1] = f.add_axes([0.63, 0.10, 0.30, 0.28])
    
    axes[0][2] = f.add_axes([0.91, 0.10, 0.02, 0.56])
    ####################### SED
    axes[0][0].plot(flux_model['wave'], flux_model['flux_A'], color='orange', linestyle='-.',label ='A', lw=1)
    axes[0][0].plot(flux_model['wave'], flux_model['flux_B'], color='dodgerblue', linestyle=(0, (5, 5)),label ='B', lw=1)
    axes[0][0].plot(flux_model['wave'], flux_model['Total'], color='green', linestyle='-',label ='Model', lw=1)
    axes[0][0].scatter(not_fitted_A['Wavelength'], not_fitted_A['Flux'], color='orange', marker='o',label ='No Fit', s=30)

    
    matplotlib.rcParams.update({'errorbar.capsize': 4})
    axes[0][0].errorbar(flux_model['wave'], flux_model['CorFlux'], yerr=flux_model['CorErr'],color='k', label='Obs',fmt='none',lw=2)
    ########## Fractional residual
    axes[1][0].plot(flux_model['wave'], (flux_model['CorFlux']-flux_model['flux_A'])/flux_model['CorFlux'],label='',marker='',color='orange', linestyle='-.',lw=1)
    axes[1][0].plot(flux_model['wave'], (flux_model['CorFlux']-flux_model['flux_A']-flux_model['flux_B'])/flux_model['CorFlux'],label='',marker='',color='green',lw=1, linestyle='-')
    axes[1][0].errorbar(flux_model['wave'], flux_model['CorFlux']-flux_model['CorFlux'], yerr=flux_model['CorErr_frac'],color='k', label='Obs',fmt='none',lw=2)
    ########## hi2_i
    axes[2][0].plot(flux_model['wave'], flux_model['chisqu'],label='',marker='o',color='green', linestyle='-',lw=1)
    ####################### L vs T (HR diagram)
    x_data = data_chi.Temp.astype(np.float).head(100)
    y_data = lumi2(radius(data_chi['SF'].head(100)),x_data) 
    axes[0][1].scatter(x_data,y_data, marker='.', label='',c='dodgerblue', s=1,rasterized = True, zorder=1)
    axes[0][1].scatter(Teff_A,L_A, marker='s', label='A-comp',c='r', s=40,rasterized = True, zorder=2)
    axes[0][1].scatter(int(Teff_B),L_B, marker='o', label='B-comp',c='b', s=40,rasterized = True, zorder=3)
    axes[0][1].scatter(10**(iso.logTe),10**(iso.logL), marker='.', label='',c='0.5', s=1,rasterized = True, zorder=0)
    ########## R vs T    
    x_data = data_chi.Temp.astype(np.float).head(100)
    y_data = radius(data_chi['SF'].head(100))
    c_data = data_chi['ChiSqr'].head(100)
    cs = axes[1][1].scatter(x_data,y_data,c=c_data, cmap='hot',rasterized = True, zorder=2)
    axes[1][1].scatter(data_chi.Temp.astype(np.int),radius(data_chi['SF']), marker='.', c='0.5', s=1,rasterized = True, zorder=0)
    ########## chi2 vs T
    y_data = data_chi['ChiSqr'].head(100)
    axes[2][1].scatter(x_data,y_data,c=c_data, cmap='hot',rasterized = True,zorder=2)  
    ####################### colorbar
    f.colorbar(cs, cax=axes[0][2])
    cs.set_clim(c_data.min(),c_data.max())
    axes[0][2].set_ylabel('$\chi^2$   (for Best 100 fits)')
    axes[0][2].yaxis.set_label_position("right")
    ####################### best fit lines
    axes[1][1].axvline(int(Teff_B), ls=(0, (5, 10)), lw=2, c='g',zorder=1)
    axes[2][1].axvline(int(Teff_B), ls=(0, (5, 10)), lw=2, c='g',zorder=1)
    axes[1][1].axhline(radius(data_chi.SF[0]), ls=(0, (5, 10)), lw=2, c='g',zorder=1)
    axes[2][1].axhline(data_chi.ChiSqr[0], ls=(0, (5, 10)), lw=2, c='g',zorder=1)
    ####################### Titles and labels
    axes[0][0].set_title(STAR_NAME+'       ' +mode + '       ' + label_A +'       '+ label_B, x=0, y=1, ha='left')
    axes[2][0].set_title('$\chi^2$ = '+str(round(flux_model['chisqu'].sum(),1))+'\n$\chi_r^2$ = '+str(round(flux_model['chisqu'].sum()/N_Np,2)),x=0.98,y=0.9, ha='right', va='top')

    axes[0][0].set_ylabel('Flux (erg s$^{-1}$ cm$^{-2}$ $\AA$$^{-1}$)')
    axes[1][0].set_ylabel('Residual')
    axes[2][0].set_ylabel('$\chi^2_i$')
    axes[2][0].set_xlabel('Wavelength ($\AA$)')
    axes[0][1].set_ylabel('L ($L_{\odot}$)')
    axes[1][1].set_ylabel('R ($R_{\odot}$)')
    axes[2][1].set_ylabel('$\chi^2$')
    axes[2][1].set_xlabel('Temp (K)')
    ####################### axes range and scales
    axes[0][0].set_xscale('log')
    axes[0][0].set_yscale('log')
    axes[1][0].set_xscale('log')
    axes[2][0].set_xscale('log')
    axes[0][1].set_yscale('log')
    axes[1][1].set_yscale('log')
    axes[0][1].set_xscale('log')
    axes[1][1].set_xscale('log')
    axes[2][1].set_xscale('log')

    wave_min = min(not_fitted_A['Wavelength'].min(),flux_model['wave'].min())
    wave_max = max(not_fitted_A['Wavelength'].max(),flux_model['wave'].max())
    
    axes[0][0].set_xlim([wave_min/1.2,wave_max*1.2])
    axes[1][0].set_xlim([wave_min/1.2,wave_max*1.2])
    axes[2][0].set_xlim([wave_min/1.2,wave_max*1.2])

    axes[0][1].set_xlim([data_chi.Temp.astype(np.int).max()*1.1,data_chi.Temp.astype(np.int).min()/1.5])
    axes[1][1].set_xlim([data_chi.Temp.astype(np.int).max()*1.1,data_chi.Temp.astype(np.int).min()/1.5])
    axes[2][1].set_xlim([data_chi.Temp.astype(np.int).max()*1.1,data_chi.Temp.astype(np.int).min()/1.5])

    flux_min = min(flux_model['CorFlux'].min(),  not_fitted_A['Flux'].min())
    flux_max = max(flux_model['CorFlux'].max(),  not_fitted_A['Flux'].max())
    axes[0][0].set_ylim([flux_min/5,flux_max*5])
    axes[0][1].set_ylim(min(L_B,L_A)/10,max(L_A,L_B)*10)
    axes[1][1].set_ylim(radius(sf_min),radius(sf_max))

    plt.setp(axes[0][0].get_xticklabels(),visible=False)
    plt.setp(axes[1][0].get_xticklabels(),visible=False)
    plt.setp(axes[0][1].get_xticklabels(),visible=False)
    plt.setp(axes[1][1].get_xticklabels(),visible=False)
    axes[0][1].set_xticks([5000, 10000, 20000, 40000,80000])
    axes[0][1].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    axes[1][1].set_xticks([5000, 10000, 20000, 40000,80000])
    axes[1][1].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    axes[2][1].set_xticks([5000, 10000, 20000, 40000,80000])
    axes[2][1].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ####################### decoration    
    for i,j in product(range(3),range(2)):
        axes[i][j].tick_params(which='both', direction='out', length=4)
        axes[i][j].grid()
    axes[0][2].tick_params(which='both', direction='out', length=4)
     
    axes[0][0].legend(scatterpoints=1, loc='upper center', ncol=5,frameon=False,handletextpad=0.3, borderpad=0.1)
    axes[0][1].legend(scatterpoints=1)
    if save==1: ########## Saving file   
        if not os.path.exists('outputs/'+mode):
            os.makedirs('outputs/'+mode)
        plt.savefig ('outputs/'+mode+'/'+STAR_NAME+'_'+str(Teff_B)+'_logg'+str(logg_B)+'_Z'+Z_B+'_'+model+'_'+str(cycle)+'.png', format='png', dpi=300)#,bbox_inches='tight')
        # plt.savefig ('outputs/'+mode+'/'+STAR_NAME+'_'+str(Teff_B)+'_logg'+str(logg_B)+'_Z'+Z_B+'_'+model+'_'+str(cycle)+'.pdf', format='pdf', dpi=300)#,bbox_inches='tight')
    plt.show()


def radius(sf):
    # Radius for given scaling factor and distance (in pc)
    return sf**0.5 * (DISTANCE*44353566.0)
    # returns R in Rsun

    
def scaling_factor(radius):
    # scaling factor for given radius (in Rsun) and distance (in pc)
    return (radius/(DISTANCE*44353566.0))**2

    
def lumi(r,t):
    sigma = 5.67e-8  #W m−2 K−4
    # r in Rsun, T in K
    return (sigma * 4 * 3.141592 * (r*6.957e+8)**2 * float(t)**4)/3.828e+26 
    # returns L in Lsun

def lumi2(r,t):  # for list of temperatures
    sigma = 5.67e-8  #W m−2 K−4
    # r in Rsun, T in K
    return (sigma * 4 * 3.141592 * (r*6.957e+8)**2 * t**4)/3.828e+26 
    # returns L in Lsun
    
def save_log(): 
    # Save the log in a file
    if not os.path.isfile('log_file.csv'):
        print ('Creating log_file.csv and saving log')
        file_object = open('log_file.csv', 'a')
        header = 'cycle,name,Teff_A_K,sf_A,logg_A,R_A_rsun,L_A_lsun,Teff_B_K,sf_B,logg_B,R_B_rsun,L_B_lsun,Z_B,e_R_B_rsun,e_L_B_lsun,chi2,chi2_r,model,N_points,N_Np\n'
        file_object.write(header)    
    else:
        file_object = open('log_file.csv', 'a')
        print ('Saving log in log_file.csv')
    details = (str(cycle) +',' + STAR_NAME+ ', '+
               str(Teff_A)+  ', '+str(sf_A)+  ', '+str(logg_A)+  ', '+ str(radius(sf_A))+', '+ str(lumi(radius(sf_A),Teff_A)) + ', '+
               Teff_B+ ', '+str(data_chi.SF[0])+ ', '+str(logg_B) + ', '+str(R_B)+  ', '+str(L_B) +', ' +str(Z_B) +', '+
               str(R_B*DISTANCE_ERR/DISTANCE)+',---,'+             # Add e_Teff_B and e_L_B_lsun manually
               str(data_chi.ChiSqr[0])+ ', '+str(data_chi['ChiSqr'][0]/N_Np)+','+
               model+ ', ' +str(N_points)+','+str(N_Np)+ '\n' )    # Append at the end of file
    file_object.write(details)
    file_object.close()


# In[2]:


STAR_NAME,logg_B,Z_B,model,cycle,fitting_required = 'WOCS2002', '7.0', '00', 'Koe', 4, 1
DIR_OBS, DISTANCE, DISTANCE_ERR, FREE_PARA = 'data/vosa_results_38873/objects/', 831.76, 11, 2+1
double_fitting = 1
'''
    STAR_NAME          Name of the star
    logg_B             logg of the star-to-be-fitted (B-component)
    Z_B                metallicity of B-component
    model              Spectral model for B-component 
                            Here Kurucz (Kr; Castelli et al 1997, AA 318, 841) and Koester (Koe; Trembley & Bergeron 2010, ApJ696, 1755) are used
                            One can add new models in "models" directory and use them
    cycle              You can use different cycle to try different things
    fitting_required   '1' if you are fitting for the first time 
                            '0' if you already have a chi2 file and just want to plot SEDs/analyse the fits
    DIR_OBS            Directory with VOSA fit results for A-component
    DISTANCE           Distance in pc
    DISTANCE_ERR       Error in distance in pc
    FREE_PARA          '2' for A-component (Kr_Temp and Kr_logg) and '1' for B-component (_Temp)
                            One has to keep it same as the number of free parameters all sed fittings
    double_fitting     '1' for double component fitting, '0' for single component fitting (---planned---)
'''

# Read data from A component fitted with VOSA
flux_A, Teff_A,logg_A,logg_A,sf_A,L_A = read_A_comp_from_VOSA(STAR_NAME)
# Remove some filters from the above file
'''
You can manually specify filters not to be fitted in drop_filters()
'''
flux_A, not_fitted_A,N_points, N_Np = drop_filters(STAR_NAME,flux_A)
# read model flux file
flux_model = read_model_file(model, logg_B, Z_B)

if (double_fitting == 0): # Planning to add a way to fit single components, basically keep "flux_A=0"
    # Upload the flux and filter names
    pass
if (double_fitting == 1):
    flux_model = add_A_comp_to_model(flux_A,flux_model)


# In[3]:


###############################################################################
mode ='rough_Kr_'+model

# Expected upper and lower bounds for B-component
#       R_WD ~ R_earth = 0.009 R_sun,    R_MS ~ 1 Rsun
R_min,  R_max  = 0.001,                 2           # in Rsun
sf_min, sf_max = scaling_factor(R_min), scaling_factor(R_max)   

scale_list, temp_B_list = initializing_sf_T(mode,sf_min, sf_max)
data_chi                = calculating_chi2(scale_list, temp_B_list,fitting_required, mode)
Teff_B,sf_B, R_B, L_B   = minimizing_chi2(flux_model,data_chi)
print ('B-component fitting parameters: T=%d,sf=%.2e,R=%f,L=%f'%(int(Teff_B),sf_B, R_B, L_B))
create_plots(mode, save=1)


# In[4]:


###############################################################################
mode='finer_Kr_'+model     # limit the fits to smaller sf range. Also use all available temperatures
sf_min = data_chi['SF'].head(100).min()
sf_max = data_chi['SF'].head(100).max()

scale_list, temp_B_list = initializing_sf_T(mode,sf_min, sf_max)
data_chi = calculating_chi2(scale_list, temp_B_list,fitting_required, mode)
Teff_B,sf_B, R_B, L_B   = minimizing_chi2(flux_model,data_chi)
print ('B-component fitting parameters: T=%d,sf=%.2e,R=%f,L=%f'%(int(Teff_B),sf_B, R_B, L_B))
create_plots(mode)
################################################################################
save_log()

