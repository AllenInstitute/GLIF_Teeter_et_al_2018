'''Written by Corinne Teeter. Will calculate the AIC for the subthreshold voltage of the fit models with 
the training data.  Note that some parameters were fit to subthreshold data but the data was not 
optimized to fit subthreshold data. 
This code produces a spread sheet of the AIC for the different model levels which can then be used to
make difference plots. Note that the AIC constant we are using leaves out a constant which is widely 
done in stats packages.  It is the differences in AIC that matters.'''

import os
import sys
relative_path=os.path.dirname(os.getcwd())
sys.path.append(os.path.join(relative_path, 'libraries'))
import numpy as np
import glif_sdk.aic as AIC
import pandas as pd
from data_library import get_sweep_num_by_name, get_file_path_endswith
from allensdk.api.queries.glif_api import GlifApi
from allensdk.core.cell_types_cache import CellTypesCache
import allensdk.core.json_utilities as json_utilities
ctc = CellTypesCache(manifest_file=os.path.join(relative_path,'cell_types_manifest.json'))
glif_api = GlifApi()


def grab_diff_v_from_folder(ew, folder):
    '''get the voltage difference if the file exists
    inputs:
        ew: string that is to be matched with a file in the folder in order to return a 
            value. i.e. if GLIF2 is being requested but there is not GLIF2 file in the
            folder a nan will be returned regardless of whether there is a value of explain 
            variance in the database. For example this would happen if the model was excluded from 
            analysis because of an aberrant parameter.  
        folder: path to the structured folder used in the rest of analysis
    returns:
        either np.nan or the explained variance ratio for the requested model.
        
    '''
    try:
        file=get_file_path_endswith(folder, ew)
        contents=json_utilities.read(file)
        RSS_of_voltage_diff=contents['noise1']['RSS_of_voltage_diff']
        num_data_points_wo_spike_shape=contents['noise1']['num_data_points_wo_spike_shape']
    except:
        RSS_of_voltage_diff=np.nan
        num_data_points_wo_spike_shape=np.nan
        
    return RSS_of_voltage_diff, num_data_points_wo_spike_shape

def calc_aic(model_num, SS, dt, sigma, stim_len):
    '''
    Calculates AIC where an agreement subthreshold voltage (as opposed to spike times) is the desired 
    result of the model. The SS is calculated by the difference between voltage of the data and the
    model for each subthreshold data point in voltage traces during the forced spike paradigm (so that 
    spikes are aligned).
    Inputs:
        model_num: integer
            specifies GLIF model 1 though 5
        SS: float
            sum of squares between v_data and v_model
        dt: float
            size of time step
        sigma: float
            standard deviation of the Gaussian used in convolution
        stim_len: integer 
            number of time steps in the stimulus
    Returns the AIC or a np.nan if there are no spikes in the model
        
    '''
    return AIC.AIC(SS, model_num, stim_len*dt/sigma)

def main():    
    # get list of data located in data folder
    folders=np.sort([os.path.join(relative_path, 'data', dir) for dir in  os.listdir(os.path.join(relative_path, 'data'))])
    
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
        
        #---get dt from the data
        #TODO: depricate after one working run though
#        cell_type_directory='/local2/workspace/cell_types/' #pointing to directory not the default ctc directory which will download to working directory
#        sweeps_file=os.path.join(cell_type_directory,'specimen_'+ str(specimen_id), 'ephys_sweeps.json')
#        data_file=os.path.join(cell_type_directory,'specimen_'+ str(specimen_id), 'ephys_sweeps.nwb')
#
#        the_sweeps=ctc.get_ephys_sweeps(specimen_id, sweeps_file)
#        data=ctc.get_ephys_data(specimen_id, data_file)
#        
#        noise1_sweeps=get_sweep_num_by_name(the_sweeps, 'Noise 1')
#        noise1_data=[]
        the_sweeps=ctc.get_ephys_sweeps(specimen_id)
        noise1_sweeps=get_sweep_num_by_name(the_sweeps, 'Noise 1')
        data=ctc.get_ephys_data(specimen_id)
        noise1_data=[]
        for s in noise1_sweeps:
            noise1_data.append(data.get_sweep(s))    
        dt=1./noise1_data[0]['sampling_rate']
        stim_len=len(noise1_data[0]['stimulus'])
        
        # using same sigma as was used in the spike time AIC calculation
        sigma=.01
        
        SS_of_model_data_v_diff_LIF, n_LIF=grab_diff_v_from_folder('_GLIF1_subthr_v.json', folder)
        SS_of_model_data_v_diff_LIFR, n_LIFR=grab_diff_v_from_folder('_GLIF2_subthr_v.json', folder)
        SS_of_model_data_v_diff_LIFASC, n_LIFASC=grab_diff_v_from_folder('_GLIF3_subthr_v.json', folder)
        SS_of_model_data_v_diff_LIFRASC, n_LIFRASC=grab_diff_v_from_folder('_GLIF4_subthr_v.json', folder)
        SS_of_model_data_v_diff_LIFRASCAT, n_LIFRASCAT=grab_diff_v_from_folder('_GLIF5_subthr_v.json', folder)
        
        # calculate the AIC
        aic_LIF.append(calc_aic(1, SS_of_model_data_v_diff_LIF, dt, sigma, n_LIF))
        aic_LIFR.append(calc_aic(2, SS_of_model_data_v_diff_LIFR, dt, sigma, n_LIFR))
        aic_LIFASC.append(calc_aic(3, SS_of_model_data_v_diff_LIFASC, dt, sigma, n_LIFASC))
        aic_LIFRASC.append(calc_aic(4, SS_of_model_data_v_diff_LIFRASC, dt, sigma, n_LIFRASC))
        aic_LIFRASCAT.append(calc_aic(5, SS_of_model_data_v_diff_LIFRASCAT, dt, sigma, n_LIFRASCAT))
        
    df=pd.DataFrame({'specimen_id':sp_ids,
                 'cre':cres,
                 'aic_LIF':aic_LIF, 
                 'aic_LIFR':aic_LIFR, 
                 'aic_LIFASC':aic_LIFASC, 
                 'aic_LIFRASC':aic_LIFRASC, 
                 'aic_LIFRASCAT':aic_LIFRASCAT})    
                    
#    df.to_csv('aic_subthresh_v_noise1.csv')    
df.to_csv('new_subv_file.csv')
    
if __name__ =='__main__':
    main()