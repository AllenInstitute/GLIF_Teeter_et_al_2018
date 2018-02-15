'''
creates plots in manuscript
'''

from scipy import signal
import os
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 16})
import allensdk.core.json_utilities as ju
import sys
relative_path=os.path.dirname(os.getcwd())
sys.path.append(os.path.join(relative_path, 'libraries'))
from glif_sdk.spike_cutting import calc_spike_cut_and_v_reset_via_expvar_residuals
from pub_plot_library import distribution_plot
from data_library import check_and_organize_data, get_sweep_num_by_name, get_pp_path
from allensdk.core.cell_types_cache import CellTypesCache
ctc = CellTypesCache(manifest_file=os.path.join(relative_path,'cell_types_manifest.json'))

#---------------------------------------------------------------
#------------NO SPECIFICATIONS NEEDED---------------------------
#---------------------------------------------------------------

def make_plots(specimen_id):    
    dir_name=os.path.join(relative_path, 'mouse_nwb/specimen_'+ str(specimen_id))
    the_sweeps=ctc.get_ephys_sweeps(specimen_id,  os.path.join(dir_name, 'ephys_sweeps.json'))
    noise1_sweeps=get_sweep_num_by_name(the_sweeps, 'Noise 1')
    
    # get data
    noise1={'current':[], 'voltage':[], 'bio_spike_times':[]}
    voltage=[]
    for s in noise1_sweeps:
        #data=ctc.get_ephys_data(specimen_id).get_sweep(s)
        data=ctc.get_ephys_data(specimen_id, os.path.join(dir_name, 'ephys.nwb')).get_sweep(s)

        spike_times=ctc.get_ephys_data(specimen_id, os.path.join(dir_name, 'ephys.nwb')).get_spike_times(s)
        voltage.append(data['response'])
        noise1['current'].append(data['stimulus'])
        noise1['bio_spike_times'].append(spike_times)
        sr=data['sampling_rate']
    dt=1./sr
    
    #--bessel the traces
    bessel={ 'N': 4, 'freq': 10000 }
    for trace in voltage:
        filt_coeff = (bessel["freq"]) / (sr / 2.) # filter fraction of Nyquist frequency
        b, a = signal.bessel(bessel["N"], filt_coeff, "low")
        noise1['voltage'].append(signal.filtfilt(b, a, trace, axis=0))
    
    calc_spike_cut_and_v_reset_via_expvar_residuals(noise1['current'], noise1['voltage'], 
                                                    dt, 0, 0, 
                                                    MAKE_PLOT=True, 
                                                    PUBLICATION_PLOT=True, 
                                                    SHOW_PLOT=True, 
                                                    BLOCK=True)


# load data out of configuration files
data_path=os.path.join(relative_path,'mouse_struc_data_dir')
folders=[os.path.join(data_path, f) for f in  os.listdir(data_path)]

all_neurons=[]   
for folder in folders:
    specimen_ID=os.path.basename(folder)[:9]
    pp_file=get_pp_path(folder)
    pp_dict=ju.read(pp_file)
    cre=os.path.basename(folder)[10:]
    all_neurons.append([specimen_ID,
              cre,
              pp_dict['spike_cut_length']['no deltaV shift']['slope'],
              pp_dict['spike_cut_length']['no deltaV shift']['intercept']*1000.,
              pp_dict['spike_cut_length']['no deltaV shift']['length']*1000.*pp_dict['dt_used_for_preprocessor_calculations']])

(cre_dict)=check_and_organize_data(all_neurons)
              
distribution_plot(cre_dict, 3 , 2, ylabel='Slope (mV/mS)', xlabel='Intercept (mV)')
distribution_plot(cre_dict, 4 , 2, ylabel='Slope (mV/mS)', xlabel='Spike cut length (ms)')

specimen_id=474637203 #htr3
make_plots(specimen_id)

specimen_id=512322162 #ctgf
make_plots(specimen_id)

