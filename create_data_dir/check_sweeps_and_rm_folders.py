'''This code will crawl the structured data directory and remove neuron (specimen id) directories
that:
1. do not have 2 or more noise 1 and noise 2 sweeps
2. do not have a GLIF1 model file
3. the spike times within the model file for the sweeps of a noise type (noise 1 or noise 2) are not the same

It will replace the "mouse_rm_bad_sweep.log" or a "human_rm_bad_sweep.log" file that provides explanations
of the removed neurons.
'''

import logging
import allensdk.core.json_utilities as json_utilities
import numpy as np
import datetime
import os
import shutil
import sys
relative_path=os.path.dirname(os.getcwd())
sys.path = [os.path.join(relative_path, 'libraries')] + sys.path
from data_library import get_ev_from_folder, get_sweep_num_by_name, get_model_spike_times_from_nwb, get_file_path_endswith, convert_spike_times_to_ind, check_spike_times_identical
import expVarOfSpikeTrains 
import pandas as pd
from allensdk.core.cell_types_cache import CellTypesCache
ctc = CellTypesCache(manifest_file=os.path.join(relative_path,'cell_types_manifest.json'))

#---------------------------------------------------------------------------
#----------------------SPECIFY SPECIES--------------------------------------
#---------------------------------------------------------------------------
species="mouse"
#species="human"

#------------------------------------------------------------------------------------------------
#------------SPECIFY WHETHER THE CODE IS BEING RUN INSIDE THE INSTITUTE---------------------------
#------------------------------------------------------------------------------------------------

where_running='external'
#where_running='internal'

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------

if species=="mouse":
    structured_data_directory='mouse_struc_data_dir'
    nwb_directory=os.path.join(relative_path,'mouse_nwb')
    log_filename="mouse_rm_bad_sweep.log"
elif species=="human":
    structured_data_directory='human_struc_data_dir'  #note that this is pointing at the local directory not the universal directory that is located one directory upstream
    nwb_directory=os.path.join(relative_path,'human_nwb')
    log_filename="human_rm_bad_sweep.log"
else:
    raise Exception("Species unrecognized")

try: os.remove(log_filename)
except: pass
logging.basicConfig(level=logging.WARNING, filename=log_filename,filemode='a+')
now=datetime.datetime.now()
logging.warning(now.strftime("%Y-%m-%d %H:%M"))

def check_more_than_two_sweeps(sweeps, nwb):
    '''Checks if ther are more than two sweeps in the electrophysiological data
    (.nwb file). This function may be irrelevant because one could just test the 
    length of the "sweeps" input. However, this function remains here because it 
    will raise an Exception if the number of sweeps listed in the in the 
    "ephys_sweeps.json" file does not match the number of sweeps in the .nwb 
    data file that contains the electrophysiological data.  This should not happen
     but there is no harm in performing this check.
    Inputs:
        sweeps: list of integers
            corresponds to sweep numbers
        nwb: NwbDataSet object
            contains electrophysiological data loaded from the nwb file using ctc.get_ephys_data()
    Returns:
        Boolian where True specifies that there are at least two data sweeps
        
    '''
    data=[]
    spike_times=[]
    for s in sweeps:
        spike_times.append(nwb.get_spike_times(s))
        data.append(nwb.get_sweep(s)) 
    if len(data)!=len(sweeps):
        raise exception
    if len(data)<2:
        return False
    else:
        return True
    
# sort the data so that specifying start and end integers works
folders=np.sort([os.path.join(structured_data_directory, f) for f in  os.listdir(structured_data_directory)])

sp_id_to_remove=[]
for ii, specimen_id_directory in enumerate(folders):
    print 'looking at', ii, 'of', len(folders)  #prints the sequential file number being run
    specimen_id=int(os.path.basename(specimen_id_directory)[:9])

    sweeps_file=os.path.join(nwb_directory,'specimen_'+ str(specimen_id), 'ephys_sweeps.json')
    nwb_file=os.path.join(nwb_directory,'specimen_'+ str(specimen_id), 'ephys.nwb')

    # load files
    the_sweeps=ctc.get_ephys_sweeps(specimen_id, sweeps_file)
    nwb=ctc.get_ephys_data(specimen_id, nwb_file) 
    n1_sweeps=get_sweep_num_by_name(the_sweeps, 'Noise 1')
    n2_sweeps=get_sweep_num_by_name(the_sweeps, 'Noise 2')
    
    # check to see if their are at least two noise 1 and noise 2 sweeps in the data nwb file
    if not check_more_than_two_sweeps(n1_sweeps, nwb):
        print specimen_id ,"has less than two noise_1 sweeps"
        logging.warning(str(specimen_id) +" has less than two noise_1 sweeps")
        sp_id_to_remove.append(specimen_id)
        continue
    if not check_more_than_two_sweeps(n2_sweeps, nwb):
        print specimen_id ,"has less than two noise_2 sweeps"
        logging.warning(str(specimen_id) +" has less than two noise_2 sweeps")
        sp_id_to_remove.append(specimen_id)
        continue
    
    # check if there is at least a level 1 GLIF model in the structured data directory    
    if not get_file_path_endswith(specimen_id_directory, '_GLIF1_neuron_config.json'):
        print specimen_id, "has no model file!!!!!!!"
        logging.warning(str(specimen_id) +" has no model file")
        sp_id_to_remove.append(specimen_id)
        continue

    # check to see that the spike times in the model files are the same (if not, it is likely that the was more than one stimulus amplitude played to the neuron) 
    glif_spike_times_n1=get_model_spike_times_from_nwb('_GLIF1_neuron_config.json', specimen_id_directory, '(LIF)', n1_sweeps, where_running)
    if check_spike_times_identical(glif_spike_times_n1):
        pass
    else: 
        print specimen_id , "noise 1 has inconsistent model sweep times"
        logging.warning(str(specimen_id)+"noise 1 has inconsistent model sweep times")
        sp_id_to_remove.append(specimen_id)
    
    glif_spike_times_n2=get_model_spike_times_from_nwb('_GLIF1_neuron_config.json', specimen_id_directory, '(LIF)', n2_sweeps, where_running)
    if check_spike_times_identical(glif_spike_times_n2):
        pass
    else: 
        print specimen_id , "noise 2 has inconsistent model sweep times"
        logging.warning(str(specimen_id)+"noise 2 has inconsistent model sweep times")
        sp_id_to_remove.append(specimen_id)
        
# remove the directories that do not pass the data requirements
print "removing", len(np.unique(sp_id_to_remove)), "out of", len(folders), 'folders' 
logging.warning("removing "+str(len(np.unique(sp_id_to_remove)))+" out of "+str(len(folders))+" folders") 
for specimen_id_directory in folders:
    specimen_id=int(os.path.basename(specimen_id_directory)[:9])
    if specimen_id in np.unique(sp_id_to_remove):
        cre=os.path.basename(specimen_id_directory)[10:]
        shutil.rmtree(os.path.join(structured_data_directory,str(specimen_id)+'_'+cre))



    