#!/usr/bin/env python
# coding: utf-8

# This file has 2 parts, run the one you want.
# 
# VOSA requires following format for uploading the photometric information files
#  ----------------------------------------------------------------------------
# | object  | RA  | DEC | dis | Av | filter | flux | error | pntopts | objopts |
# | ---     | --- | --- | --- | ---| ---    | ---  | ---   | ---     | ---     |
# | ---     | --- | --- | --- | ---| ---    | ---  | ---   | ---     | ---     |
# 
# Identify the filters in the files you are uploading (create the "filter_list" accordingly)
# 
# Create a file with "name, ra, dec, magnitudes and magnitude errors".
# This "photomety_file" be converted to VOSA format
#     Note: This code is for a cluster, hence distance and extinction is kept constant
#     
# The VOSA_input.txt file has magnitudes in the "flux" column. So while uploading to VOSA, keep "file type" as "magnitudes"

# In[1]:


# ---------------------------------------------------------------------- #
# -------------------- Creating VOSA upload files ---------------------- #
# ---------------------------------------------------------------------- #
import pandas as pd

DISTANCE = '831.76+-11'  # in pc
A_V = 0.1736    # = 3.1*E(B-V)
e_A_V = 0

photomety_file = 'data/example_photomety_file.csv'    # Has name, ra, dec, magnitude, errors
photometry = pd.read_csv(photomety_file, engine='python')
print (list(photometry.columns.values))


# In[3]:


# Make the three lists. 1) VOSA names for filters, 2)names of the magnitude colums, 3)names of the error column
filter_list = ['KPNO/Mosaic.B', 'KPNO/Mosaic.V','KPNO/Mosaic.I', 'KPNO/Mosaic.U','KPNO/Mosaic.R','GALEX/GALEX.FUV','GALEX/GALEX.NUV','Astrosat/UVIT.F148W','Astrosat/UVIT.F154W','Astrosat/UVIT.F169M']
mag_list = ['B', 'V', 'I','U','I','R','GALEX_NUV','GALEX_FUV','F148W','F154W','F169M']
err_list = ['e_B', 'e_V', 'e_I','e_U','e_I','e_R','e_GALEX_NUV','e_GALEX_FUV','e_F148W','e_F154W','e_F169M']

''' 
If one wants to calculate flux in jy instead of magnitudes:

# create a zero_point_list with "AB/VEGA magnitude zeropoints" (http://svo2.cab.inta-csic.es/theory/fps/)
zero_point_list = [3954.5,3632,2384.1,1681.2,2945.8,3631,3631,3631,3631,3631]

flux_jy = zero_point_list[j]* 10**(-0.4*photometry[mag_list[j]][i])
e_flux_jy = zero_point_list[j]* 10**(-0.4*(photometry[mag_list[j]][i]-photometry[err_list[j]][i])) - flux_jy
calculate flux_jy and e_flux_jy in each loop.
'''

# combining data from all stars to make the VOSA upload file 
VOSA = pd.DataFrame(columns = ['object', 'RA', 'DEC','dis','Av','filter','flux','error','pntopts','objopts']) 
for i in range (0,len(photometry)):
    for j in range (0,len(filter_list)):
        VOSA = VOSA.append({'object': photometry['name'][i], 
                            'RA':photometry['ra'][i], 
                            'DEC':photometry['dec'][i],
                            'dis':DISTANCE,
                            'Av':str(A_V)+'+-'+str(e_A_V),
                            'filter':filter_list[j],
                            'flux':photometry[mag_list[j]][i],
                            'error':photometry[err_list[j]][i],
                            'pntopts':'---',
                            'objopts': '---'},ignore_index = True)
VOSA.fillna('---', inplace=True)

print (VOSA)
# VOSA.to_csv('data/example_VOSA_input_file.txt', header=None, index=None, sep=' ')
'''
Now upload the file at http://svo2.cab.inta-csic.es/theory/vosa/index.php?action=myfiles&otype=star&seeall= 
Make sure to change the File type: To magnitude or Flux (jy) 
Keep SED Type: Flux vs Lambda 
Select the file and search through VO for all possible detections Look at the SEDs, 
Possibly remove some telescopes (SDSS creates problems most of the times)
'''


# In[18]:


# ---------------------------------------------------------------------- #
# --------------- Creating Synthetic Photometry files ------------------ #
#         Download synthetic photometry for all filters from             #  
#        http://svo2.cab.inta-csic.es/theory/newov2/syph.php             #
#                 unzip the files in "models" folder                     #
# ---------------------------------------------------------------------- #
import os
import pandas as pd
# For Kr models:
#          logg is like: 3.0, 4.5, 5.0 etc.
#          Z is m05 (-0.5), p00 (0), p05(0.5) etc.
Kr_logg_list=['0.0','0.5','1.0','1.5','2.0','2.5','3.0','3.5','4.0','4.5','5.0']
Kr_Z_list = ['m25','m20','m15','m10','m05','p00','p05']
# I dont need all models, so I have used shorter lists
Kr_logg_list=['3.0','3.5','4.0','4.5','5.0']
Kr_Z_list = ['m05','p00','p05']
# For Koe models:
#          logg is like: 700, 750, 850 etc. (meaning 7,8,9 etc.)
#          Z does not matter 
Koe_logg_list = ['650','675','700','725','750','775','800','825','850','875','900','925','950']
Koe_logg_list = ['700','800','900']

def create_synthetic_photometry(model,logg,Z):
    if model=='Koe':     # List of all available temperatures for Koester model
        temp_list=['05000','05250','05500','05750','06000','06250','06500','06750','07000','07250','07500','07750','08000','08250','08500','08750','09000','09250','09500','09750','10000','10250','10500','10750','11000','11250','11500','11750','12000','12250','12500','12750','13000','13250','13500','13750','14000','14250','14500','14750','15000','15250','15500','15750','16000','16250','16500','16750','17000','17250','17500','17750','18000','18250','18500','18750','19000','19250','19500','19750','20000','21000','22000','23000','24000','25000','26000','27000','28000','29000','30000','32000','34000','35000','36000','38000','40000','45000','50000','60000','70000','80000']
    if model == 'Kr':    # List of all available temperatures for Kurucz model
        temp_list=['3500','3750','4000','4250','4500','4750','5000','5250','5500','5750','6000','6250','6500','6750','7000','7250','7500','7750','8000','8250','8500','8750','9000','9250','9500','9750','10000','10250','10500','10750','11000','11250','11500','11750','12000','12250','12500','12750','13000','14000','15000','16000','17000','18000','19000','20000','21000','22000','23000','24000','25000','26000','27000','28000','29000','30000','31000','32000','33000','34000','35000','36000','37000','38000','39000','40000','41000','42000','43000','44000','45000','46000','47000','48000','49000','50000']

    flux_model =  pd.DataFrame()

    for i in range (0,len(temp_list)):
        if model=='Koe':
            _model_files = 'models/koester2_phot_1605778413.1972/koester2_da'+str(temp_list[i])+'_'+logg+'.dk.phot.dat'   # Check the names for logg, Z, model etc. 
        if model == 'Kr':
            _model_files= 'models/Kurucz_phot_1605787543.3471/Kurucz_f'+Z+'k2odfnew.pck.teff='+str(temp_list[i])+'..logg='+logg+'0000.phot.dat'
        if not os.path.isfile(_model_files):
            #print (temp_list[i]),
            pass
        else:
            _data_flux = pd.read_csv(_model_files,comment='#',engine='python', header=None, delim_whitespace= True,skipinitialspace=True)
            _data_flux.columns = ['filter', 'wave', 'flux']
            flux_model[temp_list[i]]=_data_flux['flux']    # Saving the flux from above file in a new dataframe as a single column

    # adding the filter and wavelength as new columns    
    flux_model['filter'] = _data_flux['filter']
    flux_model['wave'] = _data_flux['wave']
    flux_model = flux_model.sort_values(by=['wave'])
    if model == 'Koe':
        output_name = 'models/'+model+'_logg'+str(int(logg)/100.0)+'.csv'
    if model =='Kr':
        output_name = 'models/'+model+'_logg'+logg+'_Z'+Z+'.csv'
    flux_model.to_csv(output_name, index=False)
    print (output_name)

# for logg in Kr_logg_list:
#     for Z in Kr_Z_list:
#         print (logg,Z)
#         create_synthetic_photometry('Kr',logg,Z)
for logg in Koe_logg_list:
        print (logg)
        Z='--'
        create_synthetic_photometry('Koe',logg,Z)

