'''In order to assess if the distributions of the differences between levels are different
'''
import os
import numpy as np
import allensdk.core.json_utilities as ju
import scipy.stats as stats
import sys
relative_path=os.path.dirname(os.getcwd())
sys.path.append(os.path.join(relative_path, 'libraries'))
from data_library import get_ev_from_folder, pos_cre_lines, excitatory_cre_lines, inhibitory_cre_lines
import pandas as pd
import allensdk.core.json_utilities as ju

def set_up_data(folders):
    neuron_data=[]
    for folder in folders:
        name=os.path.basename(folder)
        specimen_id=name[:9]
        cre=os.path.basename(folder)[10:]
        if cre in excitatory_cre_lines:
            excit='excitatory'
        elif cre in inhibitory_cre_lines:
            excit='inhibitory'
        else:
            raise Exception('this cre-line is not classified')
        ev_LIF=np.nan
        ev_LIFR=np.nan
        ev_LIFASC=np.nan
        ev_LIFRASC=np.nan
        ev_LIFRASCAT=np.nan

        
        ev_LIF=get_ev_from_folder('_GLIF1_neuron_config.json', folder, '(LIF)')*100.    
        ev_LIFR=get_ev_from_folder('_GLIF2_neuron_config.json', folder, '(LIF-R)')*100.    
        ev_LIFASC=get_ev_from_folder('_GLIF3_neuron_config.json', folder, '(LIF-ASC)')*100.    
        ev_LIFRASC=get_ev_from_folder('_GLIF4_neuron_config.json', folder, '(LIF-R-ASC)')*100.    
        ev_LIFRASCAT=get_ev_from_folder('_GLIF5_neuron_config.json', folder, '(LIF-R-ASC-A)')*100.    
        
        neuron_data.append([specimen_id, cre, excit, ev_LIF, ev_LIFR, ev_LIFASC, ev_LIFRASC, ev_LIFRASCAT]) 

    df=pd.DataFrame(neuron_data)
    df.columns=['specimen_id', 'cre', 'excitation','ev_LIF', 'ev_LIFR', 'ev_LIFASC', 'ev_LIFRASC', 'ev_LIFRASCAT']
    
    #look at distributions of differences
    df['1-2']=df['ev_LIF']-df['ev_LIFR']
    df['1-3']=df['ev_LIF']-df['ev_LIFASC']
    df['1-4']=df['ev_LIF']-df['ev_LIFRASC']
    df['1-5']=df['ev_LIF']-df['ev_LIFRASCAT']
    df['2-3']=df['ev_LIFR']-df['ev_LIFASC']
    df['2-4']=df['ev_LIFR']-df['ev_LIFRASC']
    df['2-5']=df['ev_LIFR']-df['ev_LIFRASCAT']
    df['3-4']=df['ev_LIFASC']-df['ev_LIFRASC']
    df['3-5']=df['ev_LIFASC']-df['ev_LIFRASCAT']
    df['4-5']=df['ev_LIFRASC']-df['ev_LIFRASCAT']
    
    return df 

#---comment these in to not use saved file
#data_path=os.path.join(relative_path,'data')
#folders=[os.path.join(data_path, f) for f in  os.listdir(data_path)]
#df=set_up_data(folders)

df=pd.read_csv('saved_data/diff_df.csv')
#--------------------------------------

def get_differences(key):
    excite[key]=df[key][df['excitation']=='excitatory'].values
    excite[key] = excite[key][~np.isnan(excite[key])]
    inhibit[key]=df[key][df['excitation']=='inhibitory'].values
    inhibit[key] = inhibit[key][~np.isnan(inhibit[key])]

excite={}
inhibit={}
keys=['1-2', '1-3','1-4','1-5','2-3','2-4','2-5','3-4','3-5', '4-5']
for key in keys:
    get_differences(key)


p_list=[]
for key in keys:
    p_list.append([key, stats.mannwhitneyu(excite[key], inhibit[key])[1]])

p_df=pd.DataFrame(p_list)
p_df.columns=['level', 'p_value']
p_df['bonferroni_p']=p_df['p_value']*10.
p_df=p_df.sort('p_value').reset_index(drop=True)
p_df['Benjamini_Hochberg_p']=p_df['bonferroni_p']/(p_df.index+1)
print p_df
