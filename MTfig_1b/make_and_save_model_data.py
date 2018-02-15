'''This code runs and pickles data for use in mech_zoom_ex_fig_1.py.
This code takes ~30 min to run. 
'''
import numpy as np
import pickle
import os
import allensdk.core.json_utilities as ju
from allensdk.model.glif.glif_neuron import GlifNeuron
import sys
import time
relative_path=os.path.dirname(os.getcwd())
sys.path.append(os.path.join(relative_path, 'libraries'))
from data_library import get_file_path_endswith, get_sweep_num_by_name
from allensdk.core.cell_types_cache import CellTypesCache
ctc = CellTypesCache(manifest_file=os.path.join(relative_path,'cell_types_manifest.json'))

#---------------------------------------------------------------
#------------NO SPECIFICATIONS NEEDED---------------------------
#---------------------------------------------------------------

def get_model(path, EW):
    '''Runs the model for a specified neuron and model
    inputs:
        path: string
            folder path with files for the neuron
        EW: string
            end of file searching for:  options '_GLIF1_neuron_config.json',_GLIF2_neuron_config.json' etc.
    returns:
        run_data: dictionary
            contains data from the model run
            
   '''

    specimen_id=int(os.path.basename(path)[:9])
    file=get_file_path_endswith(path, EW)

    # load data
    dir_name=os.path.join(relative_path, 'mouse_nwb/specimen_'+ str(specimen_id))
    all_sweeps=ctc.get_ephys_sweeps(specimen_id,  os.path.join(dir_name, 'ephys_sweeps.json'))
    #all_sweeps=ctc.get_ephys_sweeps(specimen_id)
    sweeps=get_sweep_num_by_name(all_sweeps, 'Noise 2')
    
    noise2_sweeps = get_sweep_num_by_name(all_sweeps, 'Noise 2')
#    noise2_data=ctc.get_ephys_data(specimen_id).get_sweep(noise2_sweeps[0])
    noise2_data=ctc.get_ephys_data(specimen_id, os.path.join(dir_name, 'ephys.nwb')).get_sweep(noise2_sweeps[0])

    # run model with current
    stimulus2=noise2_data['stimulus']
    neuron_config=ju.read(file)
    neuron_config['dt']=1./noise2_data['sampling_rate'] #reset dt to the stimulus dt not the optimization dt
    neuron = GlifNeuron.from_dict(neuron_config)
    1/noise2_data['sampling_rate']
    run_data = neuron.run(stimulus2)
    run_data['time']=np.arange(0, len(run_data['voltage']))*neuron_config['dt']
    run_data['El_reference']=neuron_config['El_reference']    
    run_data['stimulus']=noise2_data['stimulus']

    return run_data

def make_and_save_models(specimen_id):
    '''Runs models and creates resulting voltage waveforms and saves them to a pickle file
    inputs:
        specimen_id: integer
            specifies neuron to be run
    outputs:
        pickle files
    '''
    
    global start_time #grab start_time from outside this module
    
    # finding the folder associated with the desired specimen_id 
    for dir in folders:
        sp_id=int(os.path.basename(dir)[:9])
        if sp_id == specimen_id:
            folder=dir
    cre=os.path.basename(folder)[10:]
    
    try:
        os.makedirs('pkl_data')
    except: pass
    
    print 'running LIF'
    LIF_model=get_model(folder, '_GLIF1_neuron_config.json')
    pickle.dump(LIF_model, open("pkl_data/"+str(specimen_id)+cre+"_LIF_model.pkl", "wb" ))
    print 'GLIF1 done at',(time.time()-start_time)/60., 'min'

    print 'running LIFR'
    LIFR_model=get_model(folder, '_GLIF2_neuron_config.json')
    pickle.dump(LIFR_model, open("pkl_data/"+str(specimen_id)+cre+"_LIFR_model.pkl", "wb" ))
    print 'GLIF2 done at',(time.time()-start_time)/60., 'min'

    print 'running LIFASC'
    LIFASC_model=get_model(folder, '_GLIF3_neuron_config.json')
    pickle.dump(LIFASC_model, open("pkl_data/"+str(specimen_id)+cre+"_LIFASC_model.pkl", "wb" ))
    print 'GLIF3 done at',(time.time()-start_time)/60., 'min'

    print 'running LIFRASC'
    LIFRASC_model=get_model(folder, '_GLIF4_neuron_config.json')
    pickle.dump(LIFRASC_model, open("pkl_data/"+str(specimen_id)+cre+"_LIFRASC_model.pkl", "wb" ))
    print 'GLIF4 done at',(time.time()-start_time)/60., 'min'
    
    print 'running LIFRASCAT'
    LIFRASCAT_model=get_model(folder, '_GLIF5_neuron_config.json')
    pickle.dump(LIFRASCAT_model, open("pkl_data/"+str(specimen_id)+cre+"_LIFRASCAT_model.pkl", "wb" ))
    print 'GLIF5 done at',(time.time()-start_time)/60., 'min'


path=os.path.join(relative_path,'mouse_struc_data_dir')
folders=[os.path.join(path, f) for f in  os.listdir(path)]

start_time=time.time()

#--make and save models
specimen_id=474637203
make_and_save_models(specimen_id)

specimen_id=512322162
make_and_save_models(specimen_id)


