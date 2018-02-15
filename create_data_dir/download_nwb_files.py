'''Grabs the ephys.nwb data files and the ephys_sweeps.nwb files 
for the data within the structured data directory created by 
'put_neuron_config_in_folder.py.  The files will be save in 
a directory named either "mouse_nwb" or "human_nwb" in the root
directory (for use with other scripts). Note that the nwb files 
this will consume ~39G of memory or more for mouse and ~4G or 
more for human.
'''
import os
import numpy as np
import sys
relative_path=os.path.dirname(os.getcwd())
from allensdk.core.cell_types_cache import CellTypesCache
ctc = CellTypesCache(manifest_file=os.path.join(relative_path,'cell_types_manifest.json'))

#-------------------------------------------------
# !!!!!! specify species 'mouse or human' !!!!!!!!
#-------------------------------------------------
species='mouse'
#species='human'
#-------------------------------------------------
#-------------------------------------------------

if species =='mouse':
    structured_data_directory='mouse_struc_data_dir'  #Note that this is accessing data local to this folder not the global 'mouse_data' folder saved one directory up for convenience
elif species=='human':
    structured_data_directory='human_struc_data_dir'
else:
    raise Exception('species not recognized')

# sorting folders into an order (not necessary)
folders=np.sort([os.path.join(structured_data_directory, f) for f in  os.listdir(structured_data_directory)])

def make_the_directory(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)

for folder in folders:
    specimen_id=int(os.path.basename(folder)[:9])
    
    if species =='mouse':
        dir_name=os.path.join(relative_path, 'mouse_nwb/specimen_'+ str(specimen_id))
    elif species=='human':
        dir_name=os.path.join(relative_path, 'human_nwb/specimen_'+ str(specimen_id))
    else:
        raise Exception('species not recognized')    
    
    make_the_directory(dir_name)
    ctc.get_ephys_sweeps(specimen_id,  os.path.join(dir_name, 'ephys_sweeps.json'))
    ctc.get_ephys_data(specimen_id, os.path.join(dir_name, 'ephys.nwb'))

