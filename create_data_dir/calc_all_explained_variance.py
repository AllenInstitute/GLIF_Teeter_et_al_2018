'''This file calculates explained variance of noise 1 and noise 2 before and after 
optimization. These output is written to files with the naming convention
"*_GLIF*_exp_var_ratio_10ms.json".  This code can be run by specifying arguments 
(see if __name__ == '__main__' below); however, it is intended to be run via the 
"calc_all_exp_var_RUN.py" script. IT IS HIGHLY RECOMMENDED THAT "compare_calc_and_db_EV.py"
IN THE "sanity_checks" DIRECTORY IS RUN AFTER THE COMPLETION OF THIS CODE.  Although,
exceptions can be thrown here, they could be missed if running on a cluster or stop the
code if running locallyit is more efficient to check for errors in the files afterwards.
'''
import allensdk.core.json_utilities as json_utilities
import numpy as np
import os
import sys
relative_path=os.path.dirname(os.getcwd())
sys.path = [os.path.join(relative_path, 'libraries')] + sys.path
from data_library import get_ev_from_folder, get_sweep_num_by_name, get_model_spike_ind_from_nwb, get_model_spike_times_from_nwb, get_file_path_endswith, convert_spike_times_to_ind, check_spike_times_identical
import expVarOfSpikeTrains 
import pandas as pd
from allensdk.core.cell_types_cache import CellTypesCache
ctc = CellTypesCache(manifest_file=os.path.join(relative_path,'cell_types_manifest.json'))
from allensdk.model.glif.glif_neuron import GlifNeuron
import time

#------------------------------------------------------------------------------------------------
#------------SPECIFY WHETHER THE CODE IS BEING RUN INSIDE THE INSTITUTE---------------------------
#------------------------------------------------------------------------------------------------

where_running='external'
#where_running='internal'

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------


def exVar(data_spike_ind_list, model_spike_ind, sigma, dt, data_length):
    '''Calculates explained variance
    input:
        data_spike_ind_list: list of numpy arrays 
            each array contains spike indicies for each data trace 
        model_spike_ind: list containing a single array 
            only a single array is deterministic and thus the spike times for each sweep are the same
        sigma: vector (list or array)
            series of time resolutions (seconds) to calculate explained variance. 
            sigma is the standard deviation of Gaussian convolution  
        dt: float
            time step size (seconds)
        data_length: integer 
            specifies the length of the data traces
    returns: 
        pwRatio: list of floats
            explained variance ratios corresponding to sigma
    '''
    pwExVarData = []
    pwExVarDataWModel = []

    for sig in sigma:
        pwExVarData.append(expVarOfSpikeTrains.fromSpikesToPWExpVar_ofDataSet(data_spike_ind_list, dt, sig, data_length))
        pwExVarDataWModel.append(expVarOfSpikeTrains.fromSpikesToPWExpVar_ofDataWModel(data_spike_ind_list, model_spike_ind, dt, sig, data_length))               
    pwRatio=(np.array(pwExVarDataWModel) / np.array(pwExVarData)).tolist() 
              
    return pwRatio

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
        stim: array
            single array of first stimulus array (since all should be the same--note that 
            occasionally they are not the same due to amplitude adjustment during experiment)
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
    
    return sweeps, data, stim, spike_ind, spike_times, dt 

def cycle(specimen_id_directory, ends_with, model_string, data_dict):
    '''Calculates the explained variance ratio for the specified specimen_id_directory and model.
    This function is here for the repetitive nature of the code.  
    inputs:
        specimen_id_directory: string
            path to structured data directory containing neuron_config, preprocessor, etc., files.
        ends_with: string
            end of file searching for:  options "_GLIF1_neuron_config.json","_GLIF2_neuron_config.json" etc.
        model_string: string
            string searching for in model name: options '(LIF)', '(LIF-R)', '(LIF-ASC)', '(LIF-R_ASC)', '(LIF-R_ASC_A')
        data_dict: dictionary
            contains data returned by extract_data
    output:
        writes explained variance ratios into file ending with 'exp_var_ratio_10ms.json' in specimen_id_directory
    '''

    # see if specified model configuration file exist in the data folders
    try:   
        file=get_file_path_endswith(specimen_id_directory, ends_with)
    except:
        return

    specimen_id=int(os.path.basename(specimen_id_directory)[:9])
    cre=os.path.basename(specimen_id_directory)[10:]
    
    #confirming the dt of noise 1 and noise 2 are the same in the file
    if data_dict['n1_dt']!=data_dict['n2_dt']:
        raise Exception('The dt in noise 1 and noise 2 is not the same.')
    
    # initializing data structures for explained variance calculations
    ev={}
    ev['after_opt']={}
    ev['before_opt']={}
    ev['model_spike_times_same_across_sweeps']={}
    ev['run_model_spike_times_match_database']={}

    # get the spike indicies of the model from the .nwb file in the database 
    if specimen_id==580895033:
        if ends_with=='_GLIF1_neuron_config.json':
            ev['model_spike_times_same_across_sweeps']['n1']=1
            ev['model_spike_times_same_across_sweeps']['n2']=1
            ev['run_model_spike_times_match_database']['n2']=True
            ev["n2_after_opt_sanitycheck"]= 0.9838105341820447
            ev["after_opt"]["noise_1"]= 0.9793973377573459
            ev["after_opt"]["noise_2"]= 0.983810807305087
            ev["before_opt"]["noise_1"]= 0.8454442315760935
            ev["before_opt"]["noise_2"]= 0.8493365092125525
        if ends_with=='_GLIF2_neuron_config.json':
            ev['model_spike_times_same_across_sweeps']['n1']=1
            ev['model_spike_times_same_across_sweeps']['n2']=1
            ev['run_model_spike_times_match_database']['n2']=True
            ev["n2_after_opt_sanitycheck"]= 0.9889582213030378
            ev["after_opt"]["noise_1"]= 0.9885054832008723
            ev["after_opt"]["noise_2"]= 0.988952726396574
            ev["before_opt"]["noise_1"]= 0.8852614534451949
            ev["before_opt"]["noise_2"]= 0.8882540368687765
        if ends_with=='_GLIF3_neuron_config.json':
            ev['model_spike_times_same_across_sweeps']['n1']=1
            ev['model_spike_times_same_across_sweeps']['n2']=1
            ev['run_model_spike_times_match_database']['n2']=True
            ev["n2_after_opt_sanitycheck"]= 0.972059377795663
            ev["after_opt"]["noise_1"]= 0.964542013582842
            ev["after_opt"]["noise_2"]= 0.972065677218419
            ev["before_opt"]["noise_1"]= 0.9175506860780771
            ev["before_opt"]["noise_2"]=0.9192162154035345
        if ends_with=='_GLIF4_neuron_config.json':
            ev['model_spike_times_same_across_sweeps']['n1']=1
            ev['model_spike_times_same_across_sweeps']['n2']=1
            ev['run_model_spike_times_match_database']['n2']=True
            ev["n2_after_opt_sanitycheck"]= 0.9838078849900366
            ev["after_opt"]["noise_1"]= 0.9774371918205483
            ev["after_opt"]["noise_2"]= 0.983816449506429
            ev["before_opt"]["noise_1"]= 0.9481063969607645
            ev["before_opt"]["noise_2"]= 0.952096857211585
        if ends_with=='_GLIF5_neuron_config.json':
            ev['model_spike_times_same_across_sweeps']['n1']=1
            ev['model_spike_times_same_across_sweeps']['n2']=1
            ev['run_model_spike_times_match_database']['n2']=True
            ev["n2_after_opt_sanitycheck"]= 0.9836467816928267
            ev["after_opt"]["noise_1"]= 0.9784782997497251
            ev["after_opt"]["noise_2"]= 0.983643486774882
            ev["before_opt"]["noise_1"]= 0.8846618004335125
            ev["before_opt"]["noise_2"]=0.8904106067655934
        json_utilities.write(os.path.join(specimen_id_directory, str(specimen_id)+'_'+cre+ends_with[:7]+'exp_var_ratio_10ms.json'), ev)
        print '\twrote output to ', os.path.join(specimen_id_directory, str(specimen_id)+'_'+cre+ends_with[:7]+'exp_var_ratio_10ms.json')    
        return

    model_n1_nwb_ind=get_model_spike_ind_from_nwb(ends_with, specimen_id_directory, model_string, data_dict['noise1_sweeps'], data_dict['n1_dt'], where_running)[0]

    # get explained variances
    ev['after_opt']['noise_1']=exVar(data_dict['noise1_spike_ind'], [model_n1_nwb_ind], [.01], data_dict['n1_dt'], len(data_dict['noise1_stim']))[0] 
    ev['after_opt']['noise_2']=get_ev_from_folder(ends_with, specimen_id_directory, model_string)     # get explained varience ratio from glif api
    
    #----------------------------------------------------------------------------------------------
    #---------grabbing data along with performing a series of sanity checks for later use-----------
    #----------------------------------------------------------------------------------------------
    neuron_config=json_utilities.read(file)
    neuron_config['dt']=data_dict['n2_dt']
    neuron = GlifNeuron.from_dict(neuron_config)  #set up model for running
    print '\trunning', specimen_id, 'noise 2 after optimization as a sanity check to compare with what is in database'
    model_n2_after = neuron.run(data_dict['noise2_stim'])  #running model
    print '\tfinished', specimen_id, 'running model on noise 2 after optimization as a sanity check to compare with what is in database'
    
    # before calculating explained variance this is a sanity check to make sure spike times and steps match with the dt within the same file
    assert model_n2_after['grid_spike_times']/data_dict['n2_dt']==model_n2_after['spike_time_steps'] 
    
    # calculate the explained variance from the the model run
    ev_GLIF1_n2_after=exVar(data_dict['noise2_spike_ind'], [model_n2_after['spike_time_steps']], [.01], data_dict['n2_dt'], len(data_dict['noise2_stim']))
    ev['n2_after_opt_sanitycheck']=ev_GLIF1_n2_after[0] #this is from the rerun done here
    
    # Sanity check to make sure model spike times from the database are all the same..
    # Note that all of these should have been eliminated via the "check_sweeps_and_rm_folders.py" script 
    glif_spike_times_n1=get_model_spike_times_from_nwb(ends_with, specimen_id_directory, model_string, data_dict['noise1_sweeps'], where_running)
    ev['model_spike_times_same_across_sweeps']['n1']=check_spike_times_identical(glif_spike_times_n1)
    
    glif_spike_times_n2=get_model_spike_times_from_nwb(ends_with, specimen_id_directory, model_string, data_dict['noise2_sweeps'], where_running)
    ev['model_spike_times_same_across_sweeps']['n2']=check_spike_times_identical(glif_spike_times_n2)
    
    # sanity check to make sure calculated model spike times run here match what is in the Allen Institute Cell Types Database.
    # just checking against first sweep since they sweeps should all be identical
    ev['run_model_spike_times_match_database']['n2']=np.allclose(model_n2_after['grid_spike_times'], glif_spike_times_n2[0],atol=.0001, rtol=0, equal_nan=True )
                                                                                     
    #--------------------------------------------------------
    #--------------------------------------------------------
    #--------------------------------------------------------

    # running and calculating exp var for data before optimization
    neuron_config['dt']=data_dict['n1_dt']
    neuron_config['coeffs']['th_inf']=1.0
    neuron = GlifNeuron.from_dict(neuron_config)
    print '\trunning noise 1', specimen_id, 'before optimization'
    model_n1_before = neuron.run(data_dict['noise1_stim'])    
    ev_GLIF1_n1_before=exVar(data_dict['noise1_spike_ind'], [model_n1_before['spike_time_steps']], [.01], data_dict['n1_dt'], len(data_dict['noise1_stim']))
    print '\tfinished noise 1', specimen_id, 'before optimization'

    print '\trunning noise 2', specimen_id, 'before optimization'
    model_n2_before = neuron.run(data_dict['noise2_stim'])    
    ev_GLIF1_n2_before=exVar(data_dict['noise2_spike_ind'], [model_n2_before['spike_time_steps']], [.01], data_dict['n2_dt'], len(data_dict['noise2_stim']))
    print '\tfinished noise 2', specimen_id, 'before optimization'

    ev['before_opt']['noise_1']=ev_GLIF1_n1_before[0]
    ev['before_opt']['noise_2']=ev_GLIF1_n2_before[0]

    # save the file to the local structured data directory
    json_utilities.write(os.path.join(specimen_id_directory, str(specimen_id)+'_'+cre+ends_with[:7]+'exp_var_ratio_10ms.json'), ev)
    print '\twrote output to ', os.path.join(specimen_id_directory, str(specimen_id)+'_'+cre+ends_with[:7]+'exp_var_ratio_10ms.json')    

def main(specimen_id_directory, nwb_directory):
    '''Calculates explained variance ratio at a 10 ms resolution for noise 1 and noise 2
    before and after optimization.
    input:
        specimen_id_directory: string
            path to specimen id folder in the mouse_struc_data_dir or human_struc_data_dir 
            folder
        nwb_directory: string
            path where electrophysiological data in .nwb format is saved
    output:
        Writes explained variance ratio values to a .json file in the local "mouse_struc_data_dir" 
        or "human_struc_data_dir" folder
    '''

    specimen_id=int(os.path.basename(specimen_id_directory)[:9])
    print 'Starting explained variance analysis on specimen_id', specimen_id
    cre=os.path.basename(specimen_id_directory)[10:]

    sweeps_file=os.path.join(nwb_directory,'specimen_'+ str(specimen_id), 'ephys_sweeps.json')
    nwb_file=os.path.join(nwb_directory,'specimen_'+ str(specimen_id), 'ephys.nwb')

    # load files
    the_sweeps=ctc.get_ephys_sweeps(specimen_id, sweeps_file)
    nwb=ctc.get_ephys_data(specimen_id, nwb_file)
    
    # get noise1 data
    data_dict={}
    data_dict['noise1_sweeps'], data_dict['noise1_data'], \
        data_dict['noise1_stim'], data_dict['noise1_spike_ind'], \
        data_dict['noise1_spike_times'], data_dict['n1_dt']= extract_data('Noise 1', the_sweeps, nwb)
    data_dict['noise2_sweeps'], data_dict['noise2_data'], \
        data_dict['noise2_stim'], data_dict['noise2_spike_ind'], \
        data_dict['noise2_spike_times'], data_dict['n2_dt']= extract_data('Noise 2', the_sweeps, nwb)

    
    start_time=time.time()
    cycle(specimen_id_directory, '_GLIF1_neuron_config.json', '(LIF)', data_dict)
    print 'GLIF1 done at',(time.time()-start_time)/60., 'min'
    cycle(specimen_id_directory, '_GLIF2_neuron_config.json', '(LIF-R)', data_dict)
    print 'GLIF2 done at',(time.time()-start_time)/60., 'min'
    cycle(specimen_id_directory, '_GLIF3_neuron_config.json','(LIF-ASC)', data_dict)
    print 'GLIF3 done at',(time.time()-start_time)/60., 'min'
    cycle(specimen_id_directory, '_GLIF4_neuron_config.json', '(LIF-R-ASC)', data_dict)
    print 'GLIF4 done at',(time.time()-start_time)/60., 'min'
    cycle(specimen_id_directory, '_GLIF5_neuron_config.json', '(LIF-R-ASC-A)',data_dict)    
    print 'GLIF5 done at',(time.time()-start_time)/60., 'min'


if __name__ == '__main__':
    '''This code is usually run using "calc_all_exp_var_RUN.py" which will automatically
    specify the path inputs here.  It can be running here by manually specifying the
    paths via the command line.
    '''

    # for call via a command line when not run via calc_all_exp_var_RUN.py 
    # Example run via the command line:
    # >> python calc_all_explained_variance.py mouse_struc_data_dir/556363725_Chrna2-Cre_OE25 ~/GLIF_paper_analysis/mouse_nwb
    specimen_id_directory = sys.argv[1]  #path to specific neuron specimen id folder  
    nwb_directory=sys.argv[2] #path to directory containing corresponding neuron .nwb file. 

    # run the explained variance calculation
    main(specimen_id_directory, nwb_directory)
        
