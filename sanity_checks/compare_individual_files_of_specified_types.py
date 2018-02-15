'''Compare two structured data directories.  This is useful 
characterize the difference in directories if running the code 
from scratch.  
'''

import numpy as np
import os
import sys
relative_path=os.path.dirname(os.getcwd())
sys.path = [os.path.join(relative_path, 'libraries')] + sys.path
from data_library import get_file_path_endswith
import dictdiffer
import pprint
from test_library import get_directory_structure
import pprint
import allensdk.core.json_utilities as ju


#------------------------------------------------------------------------------------
#------------SPECIFY THE TWO DIRECTORIES YOU WOULD LIKE TO COMPARE-------------------
#------------------------------------------------------------------------------------

struct_data_dir1=os.path.join(relative_path, 'create_data_dir', 'mouse_struc_data_dir')
struct_data_dir2=os.path.join(relative_path, 'mouse_struc_data_dir')
#end='_subthr_v.json'
end='_exp_var_ratio_10ms.json'
#end='_neuron_config.json'
#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------     
    
# sorting folders into an order (not necessary)
folders=np.sort([os.path.join(struct_data_dir1, f) for f in  os.listdir(struct_data_dir1)])
found1_flag=0
for specimen_id_directory in folders:
    specimen_id=os.path.basename(specimen_id_directory)[:9]
    for ends_with in ['_GLIF1'+end, 
                     '_GLIF2'+end,
                     '_GLIF3'+end,
                     '_GLIF4'+end,
                     '_GLIF5'+end]:
        # see if specified model configuration file exist in the data folders
        try:   
            json_file1=get_file_path_endswith(specimen_id_directory, ends_with)
        except:
            continue
        try:   
            json_file2=get_file_path_endswith(os.path.join(struct_data_dir2,os.path.basename(specimen_id_directory)), ends_with)
        except:
            continue
        
        # if both files exist get their contents
        d1=ju.read(json_file1)
        d2=ju.read(json_file2)
        
        # look for differences in contents with a tolerance
        result=list(dictdiffer.diff(d1,d2, tolerance=.001))
        pp=pprint.PrettyPrinter(indent=4)
        if result!=[]:
            found1_flag==found1_flag+1
            pp.pprint(result)    
if found1_flag<1:
    print "there were no significant differences found between the specified file types"
