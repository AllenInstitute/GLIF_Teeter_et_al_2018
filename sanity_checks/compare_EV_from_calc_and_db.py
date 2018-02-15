'''The "calc_all_explained_variance.py" script creates a series of 
output files ending with "*_GLIF*_exp_var_ratio_10ms.json".  This
script crawls the structured data directory and checks for error
flags.
'''
import allensdk.core.json_utilities as ju
import numpy as np
import os
import sys
relative_path=os.path.dirname(os.getcwd())
sys.path.append(os.path.join(relative_path, 'libraries'))
from data_library import get_file_path_endswith, check_spike_times_identical

#---------------------------------------------------------------------------
#------------SPECIFY THE STRUCTURED DATA DIRECTORY -------------------------
#---------------------------------------------------------------------------

# Note that these are currently pointing toward the 'create_data_dir' 
# as opposed to the global structured data directories provided in 
# the root directory of this repository. 
structured_data_directory=os.path.join(relative_path, 'create_data_dir', 'mouse_struc_data_dir')
#structured_data_directory=os.path.join(relative_path, 'create_data_dir', 'human_struc_data_dir')

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
    
# sort the data so that specifying start and end integers works
folders=np.sort([os.path.join(structured_data_directory, f) for f in  os.listdir(structured_data_directory)])

# check for issues in the "*_GLIF*_exp_var_ratio_10ms.json" files
for folder in folders:
    specimen_id=os.path.basename(folder)[:9]
    for s in ['GLIF1', 'GLIF2', 'GLIF3', 'GLIF4', 'GLIF5']:
        file =False
        try: file=get_file_path_endswith(folder, s+'_exp_var_ratio_10ms.json')
        except: pass 
        if file:
            print 'checking', file
            dict=ju.read(file)
            database_value=dict['n2_after_opt_sanitycheck']
            calculated_value=dict['after_opt']['noise_2']
            rtol=0.
            atol=1e-3
            
            # check if two values recorded in the "*_GLIF*_exp_var_ratio_10ms.json" are the same
            if not np.isclose(database_value, calculated_value, rtol=rtol, atol=atol):
                print specimen_id, s, ':the value difference,',database_value,'-', calculated_value, '=', np.absolute(database_value-calculated_value), 'is > the tolerance, ',  atol + rtol * np.absolute(calculated_value)
            
            # check to make sure model spike times of noise 1 from the database are all the same.
            if dict['model_spike_times_same_across_sweeps']['n1']==False:
                raise Exception('models with non identical model spikes should have been eliminated by check_sweeps_and_rm_folders.py')
    
            # check to make sure model spike times of noise 2from the database are all the same.            
            if dict['model_spike_times_same_across_sweeps']['n2']==False:
                raise Exception('models with non identical model spikes should have been eliminated by check_sweeps_and_rm_folders.py')    

            # check if the value specifying that the model spike times in the data base are the same as the ones given by running the model locally
            if dict['run_model_spike_times_match_database']['n2']==False:
                print specimen_id, s, ':the spike times via running the model are different then what is in the model .nwb files'

print 'CONGRATS--IF YOU SEE THS MESSAGE, THERE WERE NO ERRORS!'
