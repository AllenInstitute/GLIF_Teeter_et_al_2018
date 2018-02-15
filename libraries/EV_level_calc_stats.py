'''
Evaluate statistical significance of different model levels.  A Wilcoxon is done with a 
Benjamini_Hochberg correction for alpha error inflation.
Can visualize output with expVar_level_box_plots.py
Note that this is for data BEFORE optimization.
'''
import os
import numpy as np
import allensdk.core.json_utilities as ju
import scipy.stats as stats
import pickle
import sys
relative_path=os.path.dirname(os.getcwd())
sys.path.append(os.path.join(relative_path, 'libraries'))
from data_library import get_ev_percent_from_calculated_file, get_file_path_endswith, pos_cre_lines, excitatory_cre_lines, inhibitory_cre_lines, get_ev_from_folder
import pandas as pd
import logging

def do_friedman(neuron_data):
    '''Returns number of neurons in computation and the statistics returned from friedman test
    NOTE: this function requires neurons data to have explained variance in the same order.
    input:
        neuron_data: list of lists
            each list corresponds to data from one neuron          
    output:
        n: integer
            number of neurons used in the Wilcoxon
        out: list of floats
            output of Wilcoxon (second value is pvalue)'''
    a=[]
    b=[]
    c=[]
    d=[]
    e=[]
    n=0
    for neuron in neuron_data:
        if not np.any(np.isnan(neuron[2:6])):
            a.append(neuron[2])
            b.append(neuron[3])
            c.append(neuron[4])
            d.append(neuron[5])
            e.append(neuron[6])            
            n=n+1
    if len(a)>0:
        _, p_value=stats.friedmanchisquare(np.array(a),np.array(b),np.array(c),np.array(d),np.array(e))
        return n, p_value
    else:
        return np.nan, np.nan   

def do_wilcoxon(neuron_data, index_x, index_y):
    '''Returns number of neurons in computation and the statistics returned from Wilcoxon
    input:
        neuron_data: list of lists
            each list corresponds to data from one neuron
        index_x: integer
            index of neruon_data list for input to Wilcoxon
        index_y: integer
            index of neruon_data list for input to Wilcoxon            
    output:
        n: integer
            number of neurons used in the Wilcoxon
        out: list of floats
            output of Wilcoxon (second value is pvalue)
    '''
    x=[]
    y=[]
    n=0
    for neuron in neuron_data:
        if not np.isnan(neuron[index_x]) and not np.isnan(neuron[index_y]):
            x.append(neuron[index_x])
            y.append(neuron[index_y])
            n=n+1
    if len(x)>0:
        _, p_value=stats.wilcoxon(x, y)
        return n, p_value
    else:
        return np.nan, np.nan 
 

def add_to_stats_dictionary(stats_out, catagory, neuron_data):
    '''NOTE: that this module requires the neuron data to be in the correct order:
    [specimen_id, cre, ev_LIF, ev_LIFR, ev_LIFASC, ev_LIFRASC, ev_LIFRASCAT]
    '''
    
    def try_quartiles(neuron_data, index):
        try: 
            return [np.percentile(non_nan_data(neuron_data, index), 5), np.percentile(non_nan_data(neuron_data, index), 95)]
        except:
            return [np.nan, np.nan]
        
        
    
    stats_out[catagory]={'LIF':{'n':len(non_nan_data(neuron_data, 2)), 
                                'median':np.median(non_nan_data(neuron_data, 2)),
                                'mean':np.mean(non_nan_data(neuron_data, 2)), 
                                'quartiles': try_quartiles(neuron_data, 2) },
                         'LIFR':{'n':len(non_nan_data(neuron_data, 3)), 
                                'median':np.median(non_nan_data(neuron_data, 3)),
                                'mean':np.mean(non_nan_data(neuron_data, 3)), 
                                'quartiles': try_quartiles(neuron_data, 3) },
                         'LIFASC':{'n':len(non_nan_data(neuron_data, 4)), 
                                'median':np.median(non_nan_data(neuron_data, 4)),
                                'mean':np.mean(non_nan_data(neuron_data, 4)), 
                                'quartiles': try_quartiles(neuron_data, 4)},
                         'LIFRASC':{'n':len(non_nan_data(neuron_data, 5)), 
                                'median':np.median(non_nan_data(neuron_data, 5)),
                                'mean':np.mean(non_nan_data(neuron_data, 5)), 
                                'quartiles': try_quartiles(neuron_data, 5) },
                         'LIFRASCAT':{'n':len(non_nan_data(neuron_data, 6)), 
                                'median':np.median(non_nan_data(neuron_data, 6)),
                                'mean':np.mean(non_nan_data(neuron_data, 6)), 
                                'quartiles': try_quartiles(neuron_data, 6)},
                         'LIF_LIFR':{'n':np.nan, 'p_value': np.nan},
                         'LIF_LIFASC':{'n':np.nan, 'p_value': np.nan},
                         'LIF_LIFRASC':{'n':np.nan, 'p_value': np.nan},
                         'LIF_LIFRASCAT':{'n':np.nan, 'p_value': np.nan},
                              
                         'LIFR_LIFASC':{'n':np.nan, 'p_value': np.nan},
                         'LIFR_LIFRASC':{'n':np.nan, 'p_value': np.nan},
                         'LIFR_LIFRASCAT':{'n':np.nan, 'p_value': np.nan},

                         'LIFASC_LIFRASC':{'n':np.nan, 'p_value': np.nan},
                         'LIFASC_LIFRASCAT':{'n':np.nan, 'p_value': np.nan},
                              
                         'LIFRASC_LIFRASCAT':{'n':np.nan, 'p_value': np.nan}}
                                                            
    
    stats_out[catagory]['LIF_LIFR']['n'], stats_out[catagory]['LIF_LIFR']['p_value'] =do_wilcoxon(neuron_data, 2, 3)
    stats_out[catagory]['LIF_LIFASC']['n'], stats_out[catagory]['LIF_LIFASC']['p_value']=do_wilcoxon(neuron_data, 2, 4)
    stats_out[catagory]['LIF_LIFRASC']['n'], stats_out[catagory]['LIF_LIFRASC']['p_value']=do_wilcoxon(neuron_data, 2, 5)
    stats_out[catagory]['LIF_LIFRASCAT']['n'], stats_out[catagory]['LIF_LIFRASCAT']['p_value']=do_wilcoxon(neuron_data, 2, 6)

    stats_out[catagory]['LIFR_LIFASC']['n'], stats_out[catagory]['LIFR_LIFASC']['p_value']=do_wilcoxon(neuron_data, 3, 4)
    stats_out[catagory]['LIFR_LIFRASC']['n'], stats_out[catagory]['LIFR_LIFRASC']['p_value']=do_wilcoxon(neuron_data, 3, 5)
    stats_out[catagory]['LIFR_LIFRASCAT']['n'], stats_out[catagory]['LIFR_LIFRASCAT']['p_value']=do_wilcoxon(neuron_data, 3, 6)    

    stats_out[catagory]['LIFASC_LIFRASC']['n'], stats_out[catagory]['LIFASC_LIFRASC']['p_value']=do_wilcoxon(neuron_data, 4, 5)
    stats_out[catagory]['LIFASC_LIFRASCAT']['n'], stats_out[catagory]['LIFASC_LIFRASCAT']['p_value']=do_wilcoxon(neuron_data, 4, 6)

    stats_out[catagory]['LIFRASC_LIFRASCAT']['n'], stats_out[catagory]['LIFRASC_LIFRASCAT']['p_value']=do_wilcoxon(neuron_data, 5, 6)
    
    print 'Friedman test for',  catagory, do_friedman(neuron_data) #friedman just checks if there is one significant value
    
    
def non_nan_data(neuron_data, index):
    x=[]
    for neuron in neuron_data:
        if not np.isnan(neuron[index]):
            x.append(neuron[index])

    return x

def set_up_data(folders, neuron_data, opt_state):
    '''
    '''
    for folder in folders:
        name=os.path.basename(folder)
        specimen_id=name[:9]
        cre=os.path.basename(folder)[10:]

        # grab explained variance values from the structured data directory
        ev_LIF=get_ev_percent_from_calculated_file('GLIF1_exp_var_ratio_10ms.json', folder, opt_state, 'noise_2')
        ev_LIFR=get_ev_percent_from_calculated_file('GLIF2_exp_var_ratio_10ms.json', folder, opt_state, 'noise_2')          
        ev_LIFASC=get_ev_percent_from_calculated_file('GLIF3_exp_var_ratio_10ms.json', folder, opt_state, 'noise_2')
        ev_LIFRASC=get_ev_percent_from_calculated_file('GLIF4_exp_var_ratio_10ms.json', folder, opt_state, 'noise_2')
        ev_LIFRASCAT=get_ev_percent_from_calculated_file('GLIF5_exp_var_ratio_10ms.json', folder, opt_state, 'noise_2')
        
        neuron_data.append([specimen_id, cre, ev_LIF, ev_LIFR, ev_LIFASC, ev_LIFRASC, ev_LIFRASCAT]) 
    return neuron_data 


#data_path=os.path.join(relative_path,'data') #TODO: deprecate this line
def main(data_path, opt_state):
    folders=[os.path.join(data_path, f) for f in  os.listdir(data_path)]
    
    neuron_data=[]
    neuron_data=set_up_data(folders, neuron_data, opt_state)
        
    stats_out={}
    # get stats for all neurons
    add_to_stats_dictionary(stats_out, 'all_neurons', neuron_data)
    
    cre_data_dict={}
    #--get stats for cre lines
    for cre in pos_cre_lines:
        cre_data_dict[cre]=[]
        for neuron in neuron_data:
            if neuron[1]==cre:
                cre_data_dict[cre].append(neuron)
        add_to_stats_dictionary(stats_out, cre, cre_data_dict[cre])
    
    #--get stats for excitatory cre lines
    for cre in pos_cre_lines:
        cre_data_dict['excitatory']=[]
        for neuron in neuron_data:
            if neuron[1] in excitatory_cre_lines:
                cre_data_dict['excitatory'].append(neuron)
        add_to_stats_dictionary(stats_out, 'excitatory', cre_data_dict['excitatory'])
    
    #--get stats for inhibitory cre lines
    for cre in pos_cre_lines:
        cre_data_dict['inhibitory']=[]
        for neuron in neuron_data:
            if neuron[1] in inhibitory_cre_lines:
                cre_data_dict['inhibitory'].append(neuron)
        add_to_stats_dictionary(stats_out, 'inhibitory', cre_data_dict['inhibitory'])
    
    
    # count the number statistical tests and make a dictionary of just the stats tests
    num_stats_tests=0
    stats_list=[]
    for key, value in stats_out.iteritems():
        for k2, v2 in value.iteritems():
            for k3,v3 in v2.iteritems():
                if k3=='p_value' and key !='human':
                    num_stats_tests=num_stats_tests+1
                    stats_list.append({'cre': key, 'models':k2, 'p_value':stats_out[key][k2]['p_value'], 'n':stats_out[key][k2]['n']})
    if num_stats_tests != 190:
        raise Exception('number of statistical tests has changed: check data to make sure it is correct and update this Exception')
    
    # make a data frame and include alpha error corrected p values
    df=pd.DataFrame.from_dict(stats_list)
    df['bonferroni_p']=df['p_value']*num_stats_tests
    df=df.sort('p_value').reset_index(drop=True)
    df['Benjamini_Hochberg_p']=df['bonferroni_p']/(df.index+1)
    
    #now put these back into dictionary frame work
    for key, value in stats_out.iteritems():
        if key !='human':
            for comp in ['LIFASC_LIFRASC', 'LIFRASC_LIFRASCAT', 'LIF_LIFRASCAT', 'LIFR_LIFRASC', 'LIF_LIFASC', 'LIFR_LIFASC', 'LIFASC_LIFRASCAT', 'LIFR_LIFRASCAT', 'LIF_LIFRASC', 'LIF_LIFR']:
                stats_out[key][comp]['bonferroni_p']=df[(df['cre']==key) & (df['models']==comp)]['bonferroni_p'].values[0]
                stats_out[key][comp]['Benjamini_Hochberg_p']=df[(df['cre']==key) & (df['models']==comp)]['Benjamini_Hochberg_p'].values[0]

    return neuron_data, cre_data_dict, stats_out