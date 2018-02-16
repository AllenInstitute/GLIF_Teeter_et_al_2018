'''Note this code will go through and check to make sure all models 
of a given stimulus have the same spike times since the model is 
deterministic.  It is unlikely that it will be needed as these 
data are eliminated via the "check_sweeps_and_rm_folders.py" 
script in the "create_data_dir" directory.
'''

import os
import sys
import numpy as np
relative_path=os.path.dirname(os.getcwd())
sys.path.append(os.path.join(relative_path, 'libraries'))
from data_library import get_model_spike_times_from_nwb, get_sweep_num_by_name, check_spike_times_identical
from allensdk.core.cell_types_cache import CellTypesCache
ctc = CellTypesCache(manifest_file=os.path.join(relative_path,'cell_types_manifest.json'))

#------------------------------------------------------------------------------------------------
#------------SPECIFY WHETHER THE CODE IS BEING RUN INSIDE THE INSTITUTE---------------------------
#------------------------------------------------------------------------------------------------

where_running='external'
#where_running='internal'


# load data out of configuration files
data_path=os.path.join(relative_path,'create_data_dir/human_data')
folders=np.sort([os.path.join(data_path, f) for f in  os.listdir(data_path)])

all_neurons=[]   
for ii, folder in enumerate(folders):
# sort the data so that specifying start and end integers works
    specimen_id=int(os.path.basename(folder)[:9])
    cre=os.path.basename(folder)[10:]
    the_sweeps=ctc.get_ephys_sweeps(specimen_id)
    noise1_sweeps=get_sweep_num_by_name(the_sweeps, 'Noise 1')
    noise2_sweeps=get_sweep_num_by_name(the_sweeps, 'Noise 2')
    pairs=[['_GLIF1_neuron_config.json', '(LIF)'],
           ['_GLIF2_neuron_config.json', '(LIF-R)'],
           ['_GLIF3_neuron_config.json', '(LIF-ASC)'],  
           ['_GLIF4_neuron_config.json', '(LIF-R-ASC)'],  
           ['_GLIF5_neuron_config.json', '(LIF-R-ASC-A)']]
    for pair in pairs: #TODO: update this to use new function or check files
        glif_spike_ind_n1=get_model_spike_times_from_nwb(pair[0], folder, pair[1], noise1_sweeps, where_running)
        if check_spike_times_identical(glif_spike_ind_n1):
            pass
        else:
            print ii, specimen_id, 'has unmatching noise 1 spike times'
            print glif_spike_ind_n1
        glif_spike_ind_n2=get_model_spike_times_from_nwb(pair[0], folder, pair[1], noise1_sweeps, where_running)
        if check_spike_times_identical(glif_spike_ind_n2):
            pass
        else:
            print ii, specimen_id, 'has unmatching noise 2 spike times'
            print glif_spike_ind_n2
