'''Written by Corinne Teeter 6-18-17. Will calculate the explained variance before optimization.
This code produces a spread sheet which can then be used.'''

How is this code used?

import allensdk.internal.model.AIC as AIC
from allensdk.api.queries.glif_api import GlifApi
from allensdk.core.cell_types_cache import CellTypesCache
import allensdk.core.json_utilities as ju
from allensdk.model.glif.glif_neuron import GlifNeuron
from allensdk.core.nwb_data_set import NwbDataSet
import matplotlib.pyplot as plt
import numpy as np
import os
from for_release.data_library import get_sweep_num_by_name, get_model_nwb_path_from_folder, get_file_path_endswith
import for_release.expVarOfSpikeTrains as expVar
import pandas as pd
ctc = CellTypesCache(manifest_file=os.path.join(relative_path,'cell_types_manifest.json'))
glif_api = GlifApi()
# Depricate this and use what is in data library
#def get_model_spike_ind(ew, folder, model_string, sweeps):
#    path=get_model_nwb_path_from_folder(ew, folder, model_string)  #get nwb file path
#    if isinstance(path, basestring):
#        model=NwbDataSet(path)
#        model_spike_ind=[]
#        for sw in sweeps:
#            model_spike_ind.append((model.get_spike_times(sw)/dt).astype(int))
#        #check to make sure all spike time arrays are the same for the model
#        for ii in range(1,len(model_spike_ind)):
#            if not np.array_equal(model_spike_ind[ii], model_spike_ind[ii-1]):
#                print 'MODEL SPIKE TIMES SHOULD BE THE SAME AND THEY ARE NOT!', os.path.basename(folder)[:9]
##                raise Exception('model spike times should be the same and they are not')
#        return model_spike_ind[0]
#    else:
#        return np.nan

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
              
    return pwRatio
# load data out of configuration files
folder_path='/home/corinnet/workspace/GLIF_stripped/analysis/publication_plots/NOT_for_release/all_data'
folders=np.sort([os.path.join(folder_path, f) for f in  os.listdir(folder_path)])


sp_ids=[]
cres=[]
all_neurons=[]   
for folder in folders[:2]:
    specimen_id=int(os.path.basename(folder)[:9])
    sp_ids.append(specimen_id)
    print specimen_id
    cre=os.path.basename(folder)[10:]
    cres.append(cre)
    
    # get the file and open it
    
    #---get spike times from the data
    the_sweeps=ctc.get_ephys_sweeps(specimen_id)
    noise_sweeps=get_sweep_num_by_name(the_sweeps, 'Noise 2')
    data=ctc.get_ephys_data(specimen_id)
    noise_data=[]
    noise_spike_times=[]
    for s in noise_sweeps:
        noise_spike_times.append(data.get_spike_times(s))
        noise_data.append(data.get_sweep(s))    
    dt=1./noise_data[0]['sampling_rate']
    stim=noise_data[0]['stimulus']
    stim_len=len(stim)
    
    #convert data times to indicies
    noise1_spike_ind=[]
    for d in noise_spike_times:
        noise1_spike_ind.append((d/dt).astype(int))

    # First get the models for the neuron
    try:
        file=get_file_path_endswith(folder, '_GLIF1_neuron_config.json')
        neuron_config=ju.read(file)
        neuron_config['dt']=dt
        neuron = GlifNeuron.from_dict(neuron_config)
        print 'running', specimen_id, 'after'
        model = neuron.run(stim)
        GLIF1_after=exVar(noise1_spike_ind, [model['grid_spike_times']/dt], [.01], dt, stim_len)
    
        neuron_config['coeff'][th_inf]=1.0
        neuron = GlifNeuron.from_dict(neuron_config)
        print 'running', specimen_id, 'before'
        model = neuron.run(stim)    
        GLIF1_before=exVar(noise1_spike_ind, [model['grid_spike_times']/dt], [.01], dt, stim_len)
    except:
        GLIF1_after=np.NAN
        GLIF2_before=np.NAN
    #----get spike times of models
    
df=pd.DataFrame({'specimen_id':sp_ids,
                 'cre':cres,
                 'GLIF1_before':GLIF1_before,
                 'GLIF1_after':GLIF1_after})    
                    
#df.to_csv('expVar_beforeOpt.csv')    