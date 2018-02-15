'''This will crawl the preprocessor files and make plots 
and distributions of threshold and El. '''

import numpy as np
import os
import pandas
import allensdk.core.json_utilities as ju
import matplotlib.pyplot as plt
import sys
relative_path=os.path.dirname(os.getcwd())
sys.path.append(os.path.join(relative_path, 'libraries'))
from data_library import check_and_organize_data, get_file_path_endswith
from pub_plot_library import distribution_plot, draw_unity, color_dict

#---------------------------------------------------------------
#------------NO SPECIFICATIONS NEEDED---------------------------
#----------this file only works with mouse data-----------------
#-----------but could be altered to look at human---------------
#---------------------------------------------------------------


def extract_values(folder, the_end_file_match):
    '''extract threshold values from a config file with a specified ending
    inputs:
        folder: string 
            Path to second tier data folder with the data files of a specific neuron inside
        the_end_file_match: string
            Specifies the end of desired file name. Used for grabbing different model level config files.
    outputs:
        Dictionary of values
    '''
    dictionary={}
    dictionary['th_NOT_opt']={'from_zero': np.nan, 'absolute': np.nan}
    dictionary['th_opt']={'from_zero': np.nan, 'absolute': np.nan}
    dictionary['th_coeff']=np.nan
    dictionary['El_reference']=np.nan

    if np.any([f.endswith(the_end_file_match) for f in os.listdir(folder)]):
        file=get_file_path_endswith(folder, the_end_file_match)
        config_dict=ju.read(file)
        
        dictionary['th_NOT_opt']['absolute']=config_dict['th_inf']+config_dict['El_reference']
        dictionary['th_opt']['absolute']=config_dict['th_inf']*config_dict['coeffs']['th_inf']+config_dict['El_reference']
        
        dictionary['th_NOT_opt']['from_zero']=config_dict['th_inf']
        dictionary['th_opt']['from_zero']=config_dict['th_inf']*config_dict['coeffs']['th_inf']
        
        dictionary['th_coeff']=config_dict['coeffs']['th_inf']
        dictionary['El_reference']=config_dict['El_reference']
    else:
        if 'GLIF1' in the_end_file_match:
            print 'THERE IS NO LIF MODEL FOR NEURON', specimen_ID
            raise Exception('there should be a GLIF1 in every level')
        pass
    
    return dictionary

# load data out of configuration files
data_path=os.path.join(relative_path,'mouse_struc_data_dir')
folders=[os.path.join(data_path, f) for f in  os.listdir(data_path)]

all_neurons=[]   
for folder in folders:
    specimen_ID=os.path.basename(folder)[:9]
    cre=os.path.basename(folder)[10:]

    #-read threshold values from config files.
    #LIF 
    LIF_dict=extract_values(folder, '_GLIF1_neuron_config.json')
    LIFR_dict=extract_values(folder, '_GLIF2_neuron_config.json')
    LIFASC_dict=extract_values(folder, '_GLIF3_neuron_config.json')
    LIFRASC_dict=extract_values(folder, '_GLIF4_neuron_config.json')
    LIFRASCAT_dict=extract_values(folder, '_GLIF5_neuron_config.json')

    # simple check to make sure the non optimized threshold values are the same for every file within a neuron
    if np.any([(LIF_dict['th_NOT_opt']['absolute'] !=th and ~np.isnan(th)) for th in 
               [LIFR_dict['th_NOT_opt']['absolute'], LIFASC_dict['th_NOT_opt']['absolute'], LIFRASC_dict['th_NOT_opt']['absolute'], LIFRASCAT_dict['th_NOT_opt']['absolute']]]): 
        raise Exception('threshold values should all be the same')

    all_neurons.append([specimen_ID,                                            #0
                      cre,                                                      #1
                      LIF_dict['th_NOT_opt']['absolute']*1e3,                   #2
                      LIF_dict['El_reference']*1e3,                             #3
                      LIF_dict['th_NOT_opt']['from_zero']*1e3,                  #4
                      LIF_dict['th_opt']['absolute']*1e3,                       #5
                      LIF_dict['th_opt']['from_zero']*1e3,                      #6
                      LIFR_dict['th_opt']['absolute']*1e3,                      #7
                      LIFR_dict['th_opt']['from_zero']*1e3,                     #8                        
                      LIFASC_dict['th_opt']['absolute']*1e3,                    #9                           
                      LIFASC_dict['th_opt']['from_zero']*1e3,                   #10         
                      LIFRASC_dict['th_opt']['absolute']*1e3,                   #11                            
                      LIFRASC_dict['th_opt']['from_zero']*1e3,                  #12 
                      LIFRASCAT_dict['th_opt']['absolute']*1e3,                 #13                              
                      LIFRASCAT_dict['th_opt']['from_zero']*1e3])               #14


(cre_dict)=check_and_organize_data(all_neurons)  #organizes data into format used for distribution plotting

percentile_dict=distribution_plot(cre_dict, 3, 2, 
                                  xlabel='Resting potential (mV)', ylabel=r'Measured $\theta_{\infty}$ (mV)')
percentile_dict=distribution_plot(cre_dict, 3, 4, 
                                  xlabel='Resting potential (mV)', ylabel=r'Measured $\Delta V$ $\theta_{\infty}$ (mV)') 
plt.figure()
for cre in cre_dict:
    for neuron in cre_dict[cre]:
        plt.plot(neuron[2], neuron[5], '.', ms=12, color=color_dict[cre])
    plt.xlabel(r'Measured $\theta_{\infty}$ (mV)')
    plt.ylabel(r'GLIF$_1$ $\theta_{\infty}$ (mV)')
#    plt.title('absolute threshold')
draw_unity()

plt.figure()
for cre in cre_dict:
    for neuron in cre_dict[cre]:
        plt.plot(neuron[4], neuron[6], '.', ms=12, color=color_dict[cre])
    plt.xlabel(r'Measured $\Delta V$ $\theta_{\infty}$ (mV)')
    plt.ylabel(r'GLIF$_1 \Delta V$ $\theta_{\infty}$ (mV)')
#    plt.title(r'$\Delta V$  Threshold (mV)')
draw_unity()

plt.figure()
for cre in cre_dict:
    for neuron in cre_dict[cre]:
        plt.plot(neuron[4], neuron[8], '.', ms=12, color=color_dict[cre])
    plt.xlabel(r'Measured $\Delta V$ $\theta_{\infty}$ (mV)')
    plt.ylabel(r'GLIF$_2 \Delta V$ $\theta_{\infty}$ (mV)')
#    plt.title(r'$\Delta V$  Threshold (mV)')
draw_unity()

plt.figure()
for cre in cre_dict:
    for neuron in cre_dict[cre]:
        plt.plot(neuron[4], neuron[10], '.', ms=12, color=color_dict[cre])
    plt.xlabel(r'Measured $\Delta V$ $\theta_{\infty}$ (mV)')
    plt.ylabel(r'GLIF$_3 \Delta V$ $\theta_{\infty}$ (mV)')
#    plt.title(r'$\Delta V$  Threshold (mV)')
draw_unity()

plt.figure()
for cre in cre_dict:
    for neuron in cre_dict[cre]:
        plt.plot(neuron[4], neuron[12], '.', ms=12, color=color_dict[cre])
    plt.xlabel(r'Measured $\Delta V$ $\theta_{\infty}$ (mV)')
    plt.ylabel(r'GLIF$_4 \Delta V$ $\theta_{\infty}$ (mV)')
#    plt.title(r'$\Delta V$  Threshold (mV)')
draw_unity()

plt.figure()
for cre in cre_dict:
    for neuron in cre_dict[cre]:
        plt.plot(neuron[4], neuron[14], '.', ms=12, color=color_dict[cre])
    plt.xlabel(r'Measured $\Delta V$ $\theta_{\infty}$ (mV)')
    plt.ylabel(r'GLIF$_5 \Delta V$ $\theta_{\infty}$ (mV)')
#    plt.title(r'$\Delta V$  Threshold (mV)')
draw_unity()

plt.show()
    