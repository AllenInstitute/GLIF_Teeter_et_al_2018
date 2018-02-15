'''Written by Corinne Teeter. Will calculate the AIC for the fit models with the training data
This code produces a spread sheet of the AIC for the different model levels which can then be used to
make difference plots. Note that the AIC calculation we are using leaves out a constant (this constant is 
left out in many stats packages).  The constant is irrelevant because the differences in AIC are considered.'''

import allensdk.core.json_utilities as json_utilities
from allensdk.core.nwb_data_set import NwbDataSet
import numpy as np
import os
import sys
relative_path=os.path.dirname(os.getcwd())
sys.path.append(os.path.join(relative_path, 'libraries'))
from data_library import get_sweep_num_by_name, get_model_nwb_path_from_folder
import glif_sdk.aic as AIC
import expVarOfSpikeTrains as expVar
import pandas as pd
from allensdk.api.queries.glif_api import GlifApi
from allensdk.core.cell_types_cache import CellTypesCache
ctc = CellTypesCache(manifest_file=os.path.join(relative_path,'cell_types_manifest.json'))
glif_api = GlifApi()

# Depricate this and use what is in data_library
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

def calc_aic(model_num, sweeps, avgConvolvedData, dt, sigma, stim_len):
    '''
    Calculates AIC where an agreement in spike times (as opposed to subthreshold voltage) is the desired 
    result of the model. The SS is calculated by the difference between the convolved GLIF spike times and 
    the average of the data convolved spike times. 
    Inputs:
        model_num: integer
            specifies GLIF model 1 though 5
        sweeps: integer
            specifies sweeps of the desired data
        avgConvolvedData: array
            average of data spike times convolved with a Gaussian
        dt: float
            size of time step
        sigma: float
            standard deviation of the Gaussian used in convolution
        stim_len: integer 
            number of time steps in the stimulus
    Returns the AIC or a np.nan if there are no spikes in the model
        
    '''
    pairs=[['_GLIF1_neuron_config.json', '(LIF)'],
           ['_GLIF2_neuron_config.json', '(LIF-R)'],
           ['_GLIF3_neuron_config.json', '(LIF-ASC)'],  
           ['_GLIF4_neuron_config.json', '(LIF-R-ASC)'],  
           ['_GLIF5_neuron_config.json', '(LIF-R-ASC-A)']]
    pair_index=model_num-1
    glif_spike_ind=get_model_spike_ind(pairs[pair_index][0], folder, pairs[pair_index][1], sweeps)
    if isinstance(glif_spike_ind, np.ndarray):
        glif_convolved=expVar.makeConvolvedSpikeTrains([glif_spike_ind], dt, sigma, stim_len, convolveType='same', plot=False)[0]
        SS=np.sum((glif_convolved-avgConvolvedData)**2)
        return AIC.AIC(SS, model_num, stim_len*dt/sigma)
    elif np.isnan(glif_spike_ind):
        return np.nan

# get list of data located in data folder
folders=np.sort([os.path.join(relative_path, 'mouse_struc_data_dir', dir) for dir in  os.listdir(os.path.join(relative_path, 'mouse_struc_data_dir'))])

aic_LIF=[]
aic_LIFR=[]
aic_LIFASC=[]
aic_LIFRASC=[]
aic_LIFRASCAT=[]

sp_ids=[]
cres=[]
all_neurons=[]   
for folder in folders:
    specimen_id=int(os.path.basename(folder)[:9])
    sp_ids.append(specimen_id)
    print specimen_id
    cre=os.path.basename(folder)[10:]
    cres.append(cre)

    #---get spike times from the data
    the_sweeps=ctc.get_ephys_sweeps(specimen_id)
    noise1_sweeps=get_sweep_num_by_name(the_sweeps, 'Noise 1')
    data=ctc.get_ephys_data(specimen_id)
    noise1_data=[]
    noise1_spike_times=[]
    for s in noise1_sweeps:
        noise1_spike_times.append(data.get_spike_times(s))
        noise1_data.append(data.get_sweep(s))    
    dt=1./noise1_data[0]['sampling_rate']
    stim_len=len(noise1_data[0]['stimulus'])
    
    # convert data times to indicies
    noise1_spike_ind=[]
    for d in noise1_spike_times:
        noise1_spike_ind.append((d/dt).astype(int))
    
    # calculating convolved trains because they will be needed to calculate a SS between the traces
    sigma=.01
    convolved_data=expVar.makeConvolvedSpikeTrains(noise1_spike_ind, dt, sigma, stim_len, convolveType='same', plot=False)
    avgConvolvedData = np.mean(np.array(convolved_data), axis=0)
    
    # calculate the AIC
    aic_LIF.append(calc_aic(1, noise1_sweeps, avgConvolvedData, dt, sigma, stim_len))
    aic_LIFR.append(calc_aic(2, noise1_sweeps, avgConvolvedData, dt, sigma, stim_len))
    aic_LIFASC.append(calc_aic(3, noise1_sweeps, avgConvolvedData, dt, sigma, stim_len))
    aic_LIFRASC.append(calc_aic(4, noise1_sweeps, avgConvolvedData, dt, sigma, stim_len))
    aic_LIFRASCAT.append(calc_aic(5, noise1_sweeps, avgConvolvedData, dt, sigma, stim_len))

df=pd.DataFrame({'specimen_id':sp_ids,
                 'cre':cres,
                 'aic_LIF':aic_LIF, 
                 'aic_LIFR':aic_LIFR, 
                 'aic_LIFASC':aic_LIFASC, 
                 'aic_LIFRASC':aic_LIFRASC, 
                 'aic_LIFRASCAT':aic_LIFRASCAT})    
                    
df.to_csv('new_st_file')    