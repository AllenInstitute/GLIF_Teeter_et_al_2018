'''Creates plots for spike component of voltage.
'''

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 16})
import os
from scipy import signal
import allensdk.core.json_utilities as ju
import sys
relative_path=os.path.dirname(os.getcwd())
sys.path.append(os.path.join(relative_path, 'libraries'))
from glif_sdk.find_spikes import find_spikes_list, find_spikes_ssq_list
from glif_sdk.threshold_adaptation import calc_spike_component_of_threshold_from_multiblip
from glif_sdk.preprocess_neuron import estimate_dv_cutoff
import allensdk.ephys.ephys_extractor as efex
from data_access import load_sweeps 
from pub_plot_library import distribution_plot
from data_library import check_and_organize_data, get_sweep_num_by_name, get_file_path_endswith
from allensdk.core.cell_types_cache import CellTypesCache
ctc = CellTypesCache(manifest_file=os.path.join(relative_path,'cell_types_manifest.json'))

#---------------------------------------------------------------
#------------NO SPECIFICATIONS NEEDED---------------------------
#---------------------------------------------------------------


RESTING_POTENTIAL = 'slow_vm_mv'
def make_plots(specimen_id):
    dir_name=os.path.join(relative_path, 'mouse_nwb/specimen_'+ str(specimen_id))
    the_sweeps=ctc.get_ephys_sweeps(int(specimen_id),  os.path.join(dir_name, 'ephys_sweeps.json'))
    ssq_triple_sweeps=get_sweep_num_by_name(the_sweeps, 'Short Square - Triple')
    if specimen_id == 512322162: #the last sweep on this neuron is strange
        ssq_triple_sweeps=ssq_triple_sweeps[:-1]
    
    # put data in the format required for functions below
    multi_ssq_data={'current':[], 'voltage':[]}
    voltage=[]
    for s in ssq_triple_sweeps[:-1]:
        data=ctc.get_ephys_data(specimen_id, os.path.join(dir_name, 'ephys.json')).get_sweep(s)
        voltage.append(data['response'])
        multi_ssq_data['current'].append(data['stimulus'])
        sr=data['sampling_rate']
    dt=1./sr
    
    #--bessel the traces
    bessel={ 'N': 4, 'freq': 10000 }
    for trace in voltage:
        filt_coeff = (bessel["freq"]) / (sr / 2.) # filter fraction of Nyquist frequency
        b, a = signal.bessel(bessel["N"], filt_coeff, "low")
        multi_ssq_data['voltage'].append(signal.filtfilt(b, a, trace, axis=0))
    
    multi_ssq_dv_cutoff, multi_ssq_thresh_frac = estimate_dv_cutoff(multi_ssq_data['voltage'], dt, 
                                                                    efex.SHORT_SQUARE_TRIPLE_WINDOW_START,
                                                                    efex.SHORT_SQUARE_TRIPLE_WINDOW_END)
    (a_spike_component_of_threshold, b_spike_component_of_threshold, \
        mean_voltage_first_spike_of_blip) = calc_spike_component_of_threshold_from_multiblip(multi_ssq_data, 
                                                                                             dt, 
                                                                                             multi_ssq_dv_cutoff,
                                                                                             multi_ssq_thresh_frac,
                                                                                             MAKE_PLOT=True, 
                                                                                             SHOW_PLOT=True, 
                                                                                             BLOCK=True,
                                                                                             PUBLICATION_PLOT=True)





# load data out of configuration files
data_path=os.path.join(relative_path,'mouse_struc_data_dir')
folders=[os.path.join(data_path, f) for f in  os.listdir(data_path)]

all_neurons=[]   
for folder in folders:
    specimen_ID=os.path.basename(folder)[:9]
    cre=os.path.basename(folder)[10:]
    try:
        file=get_file_path_endswith(folder, '_GLIF2_neuron_config.json')
    except: 
        continue
    neuron_dict=ju.read(file)
    all_neurons.append([specimen_ID,
                        cre,
                        neuron_dict['threshold_reset_method']['params']['a_spike']*1.e3,
                        1./neuron_dict['threshold_reset_method']['params']['b_spike']*1.e3])

(cre_dict)=check_and_organize_data(all_neurons)

percentile_dict=distribution_plot(cre_dict, 2, 3, xlabel=r'$\delta \Theta_s (mV)$', ylabel=r'$1/b_s (ms)$')

##--------------plot examples-----------------------------

make_plots(474637203) #htr3 

make_plots(512322162) #ctgf

