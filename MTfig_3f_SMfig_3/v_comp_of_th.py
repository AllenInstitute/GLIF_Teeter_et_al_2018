'''Plots the parameters for the voltage component of the threshold.
'''

import os
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams.update({'font.size': 16})
import allensdk.core.json_utilities as ju
import sys
relative_path=os.path.dirname(os.getcwd())
sys.path.append(os.path.join(relative_path, 'libraries'))
from data_library import check_and_organize_data, get_file_path_endswith, get_pp_path
from pub_plot_library import distribution_plot

data_path=os.path.join(relative_path,'mouse_struc_data_dir')
folders=[os.path.join(data_path, f) for f in  os.listdir(data_path)]

all_neurons=[]   
for folder in folders:
    specimen_ID=os.path.basename(folder)[:9]
    cre=os.path.basename(folder)[10:]
    try:
        get_file_path_endswith(folder, '_GLIF5_neuron_config.json') #checks if the file is there, if not the values should not be used
    except: 
        continue
    pp_file=get_pp_path(folder)
    pp_dict=ju.read(pp_file)
    if pp_dict['threshold_adaptation']['a_voltage_comp_of_thr_from_fitab'] is not None and pp_dict['threshold_adaptation']['b_voltage_comp_of_thr_from_fitab'] is not None:
        all_neurons.append([specimen_ID,
                        cre,
                        pp_dict['threshold_adaptation']['a_voltage_comp_of_thr_from_fitab']/pp_dict['threshold_adaptation']['b_voltage_comp_of_thr_from_fitab'],
                        np.log10(pp_dict['threshold_adaptation']['b_voltage_comp_of_thr_from_fitab'])])
    else:
        raise Exception('this should not have made it though the exclusions')

        
cre_dict=check_and_organize_data(all_neurons)

percentile_dict=distribution_plot(cre_dict, 2, 3, xlabel=r'$a_v/b_v$', ylabel=r'log$_{10}(b_v)$')

#plt.annotate('Voltage Component of Threshold', xy=(.5, .93),
#                     xycoords='figure fraction',
#                     horizontalalignment='center', verticalalignment='bottom',
#                     fontsize=22)

plt.show()