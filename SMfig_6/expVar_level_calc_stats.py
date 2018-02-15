'''
Evaluate statistical significance of different model levels.  A Wilcoxon is done with a 
Benjamini_Hochberg correction for alpha error inflation.
Can visualize output with expVar_level_box_plots.py

'''

import os
import numpy as np
import allensdk.core.json_utilities as ju
import scipy.stats as stats
import pickle
import sys
relative_path=os.path.dirname(os.getcwd())
sys.path.append(os.path.join(relative_path, 'libraries'))
from data_library import get_file_path_endswith
import pandas as pd
import EV_level_calc_stats

data_path=os.path.join(relative_path,'mouse_struc_data_dir') 
neuron_data, cre_data_dict, stats_out=EV_level_calc_stats.main(data_path, 'before_opt')

try:
    os.makedirs('saved_data')
except: pass
pickle.dump((neuron_data,cre_data_dict, stats_out), open("saved_data/stats_out2.pkl", "wb" ))
ju.write('saved_data/stats_out2.json',[neuron_data,cre_data_dict, stats_out]) #to enable viewing if desired