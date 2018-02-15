'''Compare the skeletons of two structured data directories to see 
if the same neuron directories are present.  Note it does't look to 
if the same files within the neuron (specimen id) directories are 
the same.  Note that compare_inner_dir_structure.py will compare if
the same files within the directories are the same.
'''

import numpy as np
import os
import sys
relative_path=os.path.dirname(os.getcwd())
sys.path = [os.path.join(relative_path, 'libraries')] + sys.path
import dictdiffer
import pprint
from test_library import get_directory_structure


#------------------------------------------------------------------------------------
#------------SPECIFY THE TWO DIRECTORIES YOU WOULD LIKE TO COMPARE-------------------
#------------------------------------------------------------------------------------

struct_data_dir1=os.path.join(relative_path, 'create_data_dir', 'mouse_struc_data_dir')
struct_data_dir2=os.path.join(relative_path, 'mouse_struc_data_dir')

#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------

def remove_down_stream_dict_data(d):
    for neuron_dir in d:
        for file in d[neuron_dir].keys():
            d[neuron_dir].pop(file)
    return d
            
    
# puts the directory stuctures into two dictionaries
d1=get_directory_structure(struct_data_dir1)
d2=get_directory_structure(struct_data_dir2)

# git rid of upstream key for comparison
d1=d1.pop(d1.keys()[0])
d2=d2.pop(d2.keys()[0])

# get rid of data files with in the structured data directories
# so that the skeleton directory structure can be compared 
# without confusing it with inner files that may have not already been 
# such as the processed such as the "*_GLIF*_subthr_v.json" files.  
# This will show you if more neurons are being processed.
d1=remove_down_stream_dict_data(d1)
d2=remove_down_stream_dict_data(d2)

# assess how directories are different
pp=pprint.PrettyPrinter(indent=4)
result=list(dictdiffer.diff(d1,d2))
if result==[]:
    print 'the directory skeletons are the same!'
else:
    pp.pprint(result)
