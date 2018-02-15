'''Written by Corinne Teeter creates Supplementary Material Figure 4 and prints
out numbers reported in the article.
'''

import os 
import numpy as np
import allensdk.core.json_utilities as json_utilities
import matplotlib.pyplot as plt
import sys
relative_path=os.path.dirname(os.getcwd())
sys.path.append(os.path.join(relative_path, 'libraries'))
from data_library import get_pp_path, get_ev_percent_from_calculated_file
from allensdk.core.cell_types_cache import CellTypesCache
ctc = CellTypesCache(manifest_file=os.path.join(relative_path,'cell_types_manifest.json'))
from data_library import get_ev_from_folder, check_and_organize_data, get_sweep_num_by_name
from pub_plot_library import distribution_plot
import pandas as pd
import statsmodels.api as sm
from scipy.stats import mannwhitneyu as manu

df=pd.read_csv('spikecut_standard_err.csv')#, delimiter='\t')
data_path=os.path.join(relative_path,'mouse_struc_data_dir')
nwb_directory=os.path.join(relative_path,'mouse_nwb')
folders=[os.path.join(data_path, f) for f in  os.listdir(data_path)]

all_neurons=[]
std_error_list=[]
spike_length_list=[]
reciprocal_num_sp_list=[]
ev_LIFASC_list=[]
for folder in folders:
    specimen_id=os.path.basename(folder)[:9]
    cre=os.path.basename(folder)[10:]
    
    #get standard error from fitting spike reset rules
    std_err=df[df['specimen_id']==int(specimen_id)]['standard_err'].values[0]    
    std_error_list.append(std_err)
    
    #get spike cut length
    pp_file=get_pp_path(folder)
    pp_dict=json_utilities.read(pp_file)
    length=pp_dict['spike_cut_length']['no deltaV shift']['length']*1000.*pp_dict['dt_used_for_preprocessor_calculations'] #note this is converted to ms as opposed to seconds
    spike_length_list.append(length)

    #get number of spikes in noise_1
    if specimen_id=='580895033':
        # note that one could copy the data from the ephys_sweeps.json file from the archive
        # at http://download.alleninstitute.org/informatics-archive/september-2017/mouse_cell_types/glif/
        # to a directory named 'cell_data' (where the other data is automatically downloaded if one 
        # is reprocessing data from the Allen Institute Cell Types Database) and then comment 
        # in the relevant line below to get the values. However, here I save you from that necessity.
#        all_sweeps=ctc.get_ephys_sweeps(580895033, file_name=os.path.join(relative_path,'mouse_nwb/specimen_580895033/ephys_sweeps.json'))
#        num_of_spikes=np.mean([s['num_spikes'] for s in all_sweeps if s['stimulus_name'] == 'Noise 1' ])
#        recip=1./num_of_spikes
        num_of_spikes=295.333333333
        recip=0.0033860045146726866
    else:
        sweeps_file=os.path.join(nwb_directory,'specimen_'+ str(specimen_id), 'ephys_sweeps.json')
        all_sweeps=ctc.get_ephys_sweeps(int(specimen_id), sweeps_file)
        num_of_spikes=np.mean([s['num_spikes'] for s in all_sweeps if s['stimulus_name'] == 'Noise 1' ])
        recip=1./num_of_spikes
    reciprocal_num_sp_list.append(recip)
    
    ev_LIFASC=get_ev_percent_from_calculated_file('GLIF3_exp_var_ratio_10ms.json', folder, 'after_opt', 'noise_2')
    ev_LIFASC_list.append(ev_LIFASC)
    all_neurons.append([specimen_id, cre, std_err, length, num_of_spikes, recip, ev_LIFASC]) 
    
# create input for multiple linear regression
matrix=[std_error_list, spike_length_list, reciprocal_num_sp_list]

def multiple_regression(y, x):
    '''Performs regression of multi-linear regression.
    inputs:
        y: list of floats
            dependent variable. Here, these values correspond to the explained variance.
        x: list of lists of floats
            each list contains the values of the independent variables. Note that the order of
            this input matrix determines the interpretation of the output.  Here, we
            expect the variable lists to be in the following order: 
            [std_error_list, spike_length_list, reciprocal_num_sp_list]. Note that
            this function prepends the variables to the matrix , X, which is input
            to the statistics algorithm.  Therefore the output values are in the 
            reverse order. 
    output:
        prints output
    '''
    ones = np.ones(len(x[0]))
    # instantiating X by taking the 0th list of x and a row of 1's to account for variability we have not identified
    X = sm.add_constant(np.column_stack((x[0], ones)))
    # prepends additional lists of x to X
    for ele in x[1:]:
        X = sm.add_constant(np.column_stack((ele, X)))
    #resulting columns in the following order: reciprocal_num_sp, spike_length, std_error of spike cut length, ones
    the_fit = sm.OLS(y, X).fit()
    print 'pvalues of multilinear regression:'
    print '\treciprocal_num_sp, spike_length, std_error of spike cut length, ones'
    print '\t', the_fit.pvalues
    # Can print summary out put by 
    #print the_fit.summary()
    
# Perform multiple linear regression on the variable in question
multiple_regression(ev_LIFASC_list, matrix)

#--- calculate statistics of individual variables----
  
cre_dict=check_and_organize_data(all_neurons)

# create dataframes for easy median calculations
inhibit_df=pd.DataFrame(cre_dict['inhibitory'], columns=['specimen_id', 'cre', 'std_err', 'length','n',  '1/n', 'ev_LIFASC'])
excite_df=pd.DataFrame(cre_dict['excitatory'], columns=['specimen_id', 'cre', 'std_err', 'length', 'n', '1/n', 'ev_LIFASC'])

inhibit={}
excite={}
inhibit['std_err']={}
excite['std_err']={}
inhibit['length']={}
excite['length']={}
inhibit['1/n']={}
excite['1/n']={}
inhibit['n']={}
excite['n']={}

inhibit['std_err']['median']=np.median(inhibit_df['std_err'].values)
excite['std_err']['median']=np.median(excite_df['std_err'].values)

inhibit['length']['median']=np.median(inhibit_df['length'].values)
excite['length']['median']=np.median(excite_df['length'].values)

inhibit['1/n']['median']=np.median(inhibit_df['1/n'].values)
excite['1/n']['median']=np.median(excite_df['1/n'].values)

inhibit['n']['median']=np.median(inhibit_df['n'].values)
excite['n']['median']=np.median(excite_df['n'].values)

print '\nMedians for individual variables'
print '\tinhibitory', inhibit
print '\texcitatory', excite

# test whether variables are statistically different between inhibitory and excitatory neurons
_, std_err_p=manu(inhibit_df['std_err'].values, excite_df['std_err'].values)
_, length_p=manu(inhibit_df['length'].values, excite_df['length'].values)
_, recip_n_p=manu(inhibit_df['n'].values, excite_df['n'].values)
print '\np-values for significance between individual inhibitory and excitatory variables'
print '\tstd_err_p', std_err_p
print '\tlength_p', length_p
print '\tn_p', recip_n_p

# make plots  
distribution_plot(cre_dict, 2 , 6, xlabel='spike cut length standard error', ylabel='expVar GLIF3') 
distribution_plot(cre_dict, 3 , 6, xlabel='spike cut length', ylabel='expVar GLIF3') 
distribution_plot(cre_dict, 4 , 6, xlabel='ave number of noise 1 spikes', ylabel='expVar GLIF3') 
distribution_plot(cre_dict, 5 , 6, xlabel='1/ave number of noise 1 spikes', ylabel='expVar GLIF3') 
    
plt.show()
    