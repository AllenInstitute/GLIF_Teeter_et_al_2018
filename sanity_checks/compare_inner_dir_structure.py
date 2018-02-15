'''Compare two structured data directories.  This is useful 
characterize the difference in directories if running the code 
from scratch.  
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
    
# puts the directory stuctures into two dictionaries
d1=get_directory_structure(struct_data_dir1)
d2=get_directory_structure(struct_data_dir2)

# git rid of upstream key for comparison
d1=d1.pop(d1.keys()[0])
d2=d2.pop(d2.keys()[0])

# assess how directories are different
pp=pprint.PrettyPrinter(indent=4)
result=list(dictdiffer.diff(d1,d2))
if list(result)==[]:
    print 'the structured data directories contains the same files!'
else:
    pp.pprint(result)
