'''Calculates subthreshold differences between the model and data in the forced spike paradigm
to quantify how well models reproduce subthreshold voltage. Also calculates the point by point 
variance in subthreshold voltage of the data itself.
'''

import numpy as np
import os
import sys
relative_path=os.path.dirname(os.getcwd())
sys.path.append(os.path.join(relative_path, 'libraries'))
import time
import matplotlib.pyplot as plt
from data_library import get_sweep_num_by_name, get_model_spike_ind_from_nwb, get_file_path_endswith, convert_spike_times_to_ind
from allensdk.model.glif.glif_neuron import GlifNeuron
from allensdk.internal.model.glif.glif_optimizer_neuron import GlifOptimizerNeuron
from allensdk.core.cell_types_cache import CellTypesCache
ctc = CellTypesCache(manifest_file=os.path.join(relative_path,'cell_types_manifest.json'))
import allensdk.core.json_utilities as json_utilities

#---------------------------------------------------------------------------
#----------------------SPECIFY SPECIES--------------------------------------
#---------------------------------------------------------------------------
species="mouse"
#species="human"

'''Note that this code can be run by specifying the system argument input or 
via the command line. This is done so that one or a few neurons could be run at
once or a motivated individual could set this up to run on their cluster.  If 
no inputs are specified it will run the first 1000 neurons.  Please see  
if __name__ == '__main__' below.

Example: python calc_all_exp_var_RUN.py 0 20
'''
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------

def extract_data(sweep_name, the_sweeps, nwb):
    '''loads data for desired sweeps from the specified nwb file.
    inputs:
        sweep name: string
            string specifying sweeps to look for in the_sweeps dictionary
        the_sweeps: dictionary
            describes sweeps in experiment returned by the allensdk ctc.get_ephys_sweeps
        nwb: object
            experiment data returned by allensdk ctc.get_ephys_data
    returns:
        sweeps: list of integers
            sweep numbers corresponding to specified sweep_name
        data: list of dictionaries
            data extracted from nwb
        spike_ind: list of arrays
            each array has the indices of the spikes
        spike_times:  list of arrays
            each array has the times of the spikes
        dt: float
            time step of first sweep
    '''
        
    sweeps=get_sweep_num_by_name(the_sweeps, sweep_name)
    data=[]
    spike_times=[]
    for s in sweeps:
        spike_times.append(nwb.get_spike_times(s))
        data.append(nwb.get_sweep(s))    
    dt=1./data[0]['sampling_rate']
    stim=data[0]['stimulus']
    spike_ind=convert_spike_times_to_ind(spike_times, dt)
    
    return sweeps, data, spike_ind, spike_times, dt 

def running(stimulus_type, file, data_dict):
    '''runs the forced spike model run
    input:
        stimulus_type: string
            specifies whether to run noise 1 or noise two.
            can be 'noise1' or 'noise2'
        file: 
        data_dict:
    '''
    if stimulus_type=='noise1':
        sweep_key='noise1_sweeps'
        dt_key='n1_dt'
        data_key='noise1_data'
        spike_ind_key='noise1_spike_ind'
    elif stimulus_type=='noise2':
        sweep_key='noise2_sweeps'
        dt_key='n2_dt'
        data_key='noise2_data'
        spike_ind_key='noise2_spike_ind'
    else:
        raise Exception('Stimulus type has not been properly specified')
   
#        model_GLIF1_n2_after=get_model_spike_ind(end_with, folder, model_string, data_dict[sweep_key], data_dict[dt_key])[0]
    
    neuron_config=json_utilities.read(file)
    neuron_config['dt']=data_dict[dt_key] #set dt to be loaded into model run to be the dt the data was collected at
    forced_spike_neuron=GlifOptimizerNeuron.from_dict(neuron_config)  #class initialization to run model in forced spike paradigm

    # run the model in the forced spike paradigm for all noise sweeps
    forced_spike_output=[]
    for data, spike_ind in zip(data_dict[data_key], data_dict[spike_ind_key]):
        forced_spike_output.append(forced_spike_neuron.run_with_biological_spikes(data['stimulus'], data['response'], spike_ind))    
    
    # adjust the model voltage traces by El
    for ii in range(len(forced_spike_output)):
        forced_spike_output[ii]['voltage']=forced_spike_output[ii]['voltage']+neuron_config['El_reference']
    
    # calculated the difference in voltage between the model and data (eliminate nans due to spike removal for the calculation of variance)
    all_diff_v=np.array([]) #for combining difference in voltage (without nans) for all sweeps 
    all_data_v=np.array([]) #for combining voltage (without nans) for all sweeps   
    all_squares_diff_v=np.array([])
    for model,data, spike_ind in zip(forced_spike_output, data_dict[data_key], data_dict[spike_ind_key]):
        #diff_subsection=model['voltage'][spike_ind[0]-2:spike_ind[0]+50]-data['response'][spike_ind[0]-2:spike_ind[0]+50]
        diff_whole=model['voltage']-data['response']
        non_nans=np.where(~np.isnan(model['voltage']))[0]
        diff_no_nan=diff_whole[non_nans]
        data_no_nan=data['response'][non_nans]
        all_diff_v=np.append(all_diff_v, diff_no_nan)
        all_data_v=np.append(all_data_v, data_no_nan)
        all_squares_diff_v=np.append(all_squares_diff_v, diff_no_nan**2)
            
    out={stimulus_type: {'var_of_voltage_difference': np.var(all_diff_v), 
                         'var_of_voltage_data': np.var(all_data_v),
                         'RSS_of_voltage_diff': np.sum(all_squares_diff_v),
                         'num_data_points_wo_spike_shape': len(all_data_v)}}
    
    return out


def cycle(folder, end_with, model_string, data_dict):
    '''Calculates the squared error of the subthreshold voltage specified folder and model.
    This function is here for the repetitive nature of the code.  
    inputs:
        folder: string
            path to folder containing model configurations.
        end_with: string
            end of file searching for:  options '_GLIF1_neuron_config.json',_GLIF2_neuron_config.json' etc.
        model_string: string
            string searching for in model name: options '(LIF)', '(LIF-R)', '(LIF-ASC)', '(LIF-R_ASC)', '(LIF-R_ASC_A')
        data_dict: dictionary
            contains data returned by extract_data
        stimulus_type: string
            can be 'noise1' or 'noise2'
    output:
        writes RSS, and variance of the voltage difference between the model and the data, the variance
        of the data itself and the number of data points in subthreshold data considering into file ending
        with '*_subthr_v.json' in folder
    '''

    try:   
        file=get_file_path_endswith(folder, end_with)
    except:
        return

    specimen_id=int(os.path.basename(folder)[:9])
    cre=os.path.basename(folder)[10:]
    
    # run forced spike protocol and save subthreshold voltage difference to file
    out=running('noise1', file, data_dict)
    out_to_update=running('noise2', file, data_dict)
    out.update(out_to_update)
    output_file_name=os.path.join(folder, str(specimen_id)+'_'+cre+end_with[:7]+'subthr_v.json')
    json_utilities.write(output_file_name, out) 

def main(folder):
    '''Calculates the squared error between the subthreshold voltages of the model and data of the 
    noise 2 stimulus in the forced spike paradigm.
    input:
        folder: string
            path to specimen id folder in the data folder
    output:
        Writes voltage squared error values to a .json file in data folder
    '''

    specimen_id=int(os.path.basename(folder)[:9])
    print specimen_id
    cre=os.path.basename(folder)[10:]
    
    if species=='mouse':
        nwb_dir=os.path.join(relative_path,'mouse_nwb', 'specimen_'+ str(specimen_id))
    elif species=='human':
        nwb_dir=os.path.join(relative_path,'human_nwb', 'specimen_'+ str(specimen_id))    
    the_sweeps=ctc.get_ephys_sweeps(specimen_id,  os.path.join(nwb_dir, 'ephys_sweeps.json'))
    nwb=ctc.get_ephys_data(specimen_id, os.path.join(nwb_dir, 'ephys.nwb'))
    
    # get data
    data_dict={}
    data_dict['noise1_sweeps'], data_dict['noise1_data'], data_dict['noise1_spike_ind'], \
        data_dict['noise1_spike_times'], data_dict['n1_dt']= extract_data('Noise 1', the_sweeps, nwb)
    
    data_dict['noise2_sweeps'], data_dict['noise2_data'], data_dict['noise2_spike_ind'], \
        data_dict['noise2_spike_times'], data_dict['n2_dt']= extract_data('Noise 2', the_sweeps, nwb)

    
    start_time=time.time()
#    cycle(folder, '_GLIF1_neuron_config.json', '(LIF)', data_dict, 'noise2')
    cycle(folder, '_GLIF1_neuron_config.json', '(LIF)', data_dict)
    print 'GLIF1 done at',(time.time()-start_time)/60., 'min'
    cycle(folder, '_GLIF2_neuron_config.json', '(LIF-R)', data_dict)
    print 'GLIF2 done at',(time.time()-start_time)/60., 'min'
    cycle(folder, '_GLIF3_neuron_config.json','(LIF-ASC)', data_dict)
    print 'GLIF3 done at',(time.time()-start_time)/60., 'min'
    cycle(folder, '_GLIF4_neuron_config.json', '(LIF-R-ASC)', data_dict)
    print 'GLIF4 done at',(time.time()-start_time)/60., 'min'
    cycle(folder, '_GLIF5_neuron_config.json', '(LIF-R-ASC-A)',data_dict)    
    print 'GLIF5 done at',(time.time()-start_time)/60., 'min'


if __name__ == '__main__':
    
    '''this code can be run by specifying the system argument inputs or 
    via the command line. This is done so that one or a few neurons could be run at
    once or a motivated individual could set this up to run on their cluster.  If 
    no inputs are specified it will run the first 1000 neurons.  
    
    Example: python calc_all_exp_var_RUN.py 0 20'''
   
   # look for system arguments 
    try:
        start=int(sys.argv[1]) #index of first file to run 
        end=int(sys.argv[2]) #index +1 of last file to run
    # if no system arguments use defaults
    except:
        start=0
        end=1000
    
    # specify structured data directory by species    
    if species=='mouse':
        folders=np.sort([os.path.join(relative_path, 'mouse_struc_data_dir', dir) for dir in  os.listdir(os.path.join(relative_path, 'mouse_struc_data_dir'))])
    elif species=='human':
        folders=np.sort([os.path.join(relative_path, 'human_struc_data_dir', dir) for dir in  os.listdir(os.path.join(relative_path, 'human_struc_data_dir'))])

    # run code
    for ii in range(start, end):
        print ii, "out of", end
        print folders[ii]
        main(folders[ii])
        

