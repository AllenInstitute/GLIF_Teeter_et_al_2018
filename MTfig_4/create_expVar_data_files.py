'''Written by Corinne Teeter.  Calculates explained variance at several temporal
resolutions for Figure 4b of the manuscript.
'''

from allensdk.core.nwb_data_set import NwbDataSet
#from allensdk.api.queries.glif_api import GlifApi
#glif_api = GlifApi()
import allensdk.core.json_utilities as ju
import numpy as np
import os
import sys
relative_path=os.path.dirname(os.getcwd())
sys.path.append(os.path.join(relative_path, 'libraries'))
import expVarOfSpikeTrains as expVar
from data_library import get_sweep_num_by_name, get_model_nwb_path_from_folder, get_ev_from_folder
from allensdk.core.cell_types_cache import CellTypesCache
ctc = CellTypesCache(manifest_file=os.path.join(relative_path,'cell_types_manifest.json'))

sigma=np.array([.0001, .001, .004, .01, .020, .1, 1.0])


def exVar(data_spike_ind_list, model_spike_ind, sigma, dt, data_length):
    '''Calculates explained variance
    data_spike_ind_list: list of numpy arrays containing spike indicies for each data trace 
    model_spike_ind: list containing a singly array of data
    sigma: vector of time resolution (seconds) of which one would calculate explained variance 
    dt: time step size in seconds
    data_length: integer specifying the length of the data traces
    returns list of pairwise explained variance
    '''
    pwExVarData = []
    pwExVarDataWModel = []
    for sig in sigma:
        pwExVarData.append(expVar.fromSpikesToPWExpVar_ofDataSet(data_spike_ind_list, dt, sig, data_length))
        pwExVarDataWModel.append(expVar.fromSpikesToPWExpVar_ofDataWModel(data_spike_ind_list, model_spike_ind, dt, sig, data_length))               
    
    pwRatio=(np.array(pwExVarDataWModel) / np.array(pwExVarData)).tolist() 
              
    return pwExVarData, pwExVarDataWModel, pwRatio
    
def calc_ev(ew, folder, s, sweeps, stim_len, data_spike_times, dt):
    '''
    '''
    print ew, folder
    #convert data times to indicies
    data_spike_ind=[]
    for d in data_spike_times:
        data_spike_ind.append((d/dt).astype(int))
    
    #get model data    
    path=get_model_nwb_path_from_folder(ew, folder, s)  #get nwb file path
    if isinstance(path, basestring):
        model=NwbDataSet(path)
        model_spike_ind=[]
        for sw in sweeps:
            spikes=(model.get_spike_times(sw)/dt).astype(int)
            model_spike_ind.append(spikes)               
        #check to make sure all spike time arrays are the same for the model
        for ii in range(1,len(model_spike_ind)):
            if not np.array_equal(model_spike_ind[ii], model_spike_ind[ii-1]):
                print 'MODEL SPIKE TIMES SHOULD BE THE SAME AND THEY ARE NOT!', os.path.basename(folder)[:9]
                print len(model_spike_ind), model_spike_ind
#                raise Exception('model spike times should be the same and they are not')
        return exVar(data_spike_ind, [model_spike_ind[0]], sigma, dt, stim_len)
    else:
        return np.nan


def main(specimen_id):
    #find the folder name for specimen
    for dir in folders:
        sp_id=int(os.path.basename(dir)[:9])
        if sp_id == specimen_id:
            folder=dir
    
    #---get spike times from the data
    all_sweeps=ctc.get_ephys_sweeps(specimen_id)
    sweeps=get_sweep_num_by_name(all_sweeps, 'Noise 2')
    nwb=ctc.get_ephys_data(specimen_id)
    data=[]
    spike_times=[]
    for s in sweeps:
        spike_times.append(nwb.get_spike_times(s))
        data.append(nwb.get_sweep(s))    
    dt=1./data[0]['sampling_rate']
    stim_len=len(data[0]['stimulus'])
    
    # calculate explained variance
    ev_LIF=calc_ev('_GLIF1_neuron_config.json', folder, '(LIF)', sweeps, stim_len, spike_times, dt)
    ev_LIFR=calc_ev('_GLIF2_neuron_config.json', folder, '(LIF-R)', sweeps, stim_len, spike_times, dt)
    ev_LIFASC=calc_ev('_GLIF3_neuron_config.json', folder, '(LIF-ASC)', sweeps, stim_len, spike_times, dt)  
    ev_LIFRASC=calc_ev('_GLIF4_neuron_config.json', folder, '(LIF-R-ASC)', sweeps, stim_len, spike_times, dt)  
    ev_LIFRASCAT=calc_ev('_GLIF5_neuron_config.json', folder, '(LIF-R-ASC-A)', sweeps, stim_len, spike_times, dt)

    return [ev_LIF, ev_LIFR, ev_LIFASC, ev_LIFRASC, ev_LIFRASCAT]


data_path=os.path.join(relative_path,'mouse_struc_data_dir')
folders=np.sort([os.path.join(data_path, f) for f in  os.listdir(data_path)])

try: os.mkdir('json_data')
except: pass

ev_list=main(474637203) #use this if you want to recalculate (it take a few minutes)
ju.write('json_data/474637203htr3_ev_data2.json', ev_list)

ev_list=main(512322162)  #use this if you want to recalculate (it take a few minutes)
ju.write('json_data/512322162ctgf_ev_data2.json', ev_list)