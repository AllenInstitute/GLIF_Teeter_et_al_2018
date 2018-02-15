'''
Created by Corinne Teeter, based in summary statics code in analysis
'''

import os
import numpy as np
import allensdk.core.json_utilities as ju
import pandas as pd
import sys
relative_path=os.path.dirname(os.getcwd())
sys.path.append(os.path.join(relative_path, 'libraries'))
from data_library import get_file_path_endswith, excitatory_cre_lines, inhibitory_cre_lines, get_ev_percent_from_calculated_file, get_pp_path
from pub_plot_library import pos_cre_lines, distribution_analysis


folder_path=os.path.join(relative_path,'mouse_struc_data_dir')
folders=[os.path.join(folder_path, f) for f in  os.listdir(folder_path)]

def max_discount_nan(array):
    '''Find the maximum value in an array while ignoring nans
    '''
    temp=[]
    for value in array:
        if not np.isnan(value):
            temp.append(value)
    return np.max(np.array(temp))

def get_explained_var(folder):
    '''Finds which model files are in the folder (which means that they are not excluded) 
    and uses the glif api to get the explained variance ratio for these values. 
    '''

    out={}
    # get the explained variance for models in folders
    out['ev_LIF']=get_ev_percent_from_calculated_file('_GLIF1_exp_var_ratio_10ms.json', folder, 'after_opt', 'noise_2')    
    out['ev_LIFR']=get_ev_percent_from_calculated_file('_GLIF2_exp_var_ratio_10ms.json', folder, 'after_opt', 'noise_2')
    out['ev_LIFASC']=get_ev_percent_from_calculated_file('_GLIF3_exp_var_ratio_10ms.json', folder, 'after_opt', 'noise_2')
    out['ev_LIFRASC']=get_ev_percent_from_calculated_file('_GLIF4_exp_var_ratio_10ms.json', folder, 'after_opt', 'noise_2')
    out['ev_LIFRASCAT']=get_ev_percent_from_calculated_file('_GLIF5_exp_var_ratio_10ms.json', folder, 'after_opt', 'noise_2')

    
    return out

    
def get_optimized_th_inf(folder):

    def get_th_inf(folder, endswith):
        if np.any([f.endswith(endswith) for f in os.listdir(folder)]):
            file=get_file_path_endswith(folder, endswith)
            optimized=ju.read(file)
            return optimized['th_inf']*optimized['coeffs']['th_inf']
        else:
            return np.nan
    
    out={}
    out['th_inf_LIF']=np.nan
    out['th_inf_LIFR']=np.nan
    out['th_inf_LIFASC']=np.nan
    out['th_inf_LIFRASC']=np.nan
    out['th_inf_LIFRASCAT']=np.nan
    out['th_inf_LIF']=get_th_inf(folder, '_GLIF1_neuron_config.json')
    out['th_inf_LIFR']=get_th_inf(folder, '_GLIF2_neuron_config.json')
    out['th_inf_LIFASC']=get_th_inf(folder, '_GLIF3_neuron_config.json')
    out['th_inf_LIFRASC']=get_th_inf(folder, '_GLIF4_neuron_config.json')    
    out['th_inf_LIFRASCAT']=get_th_inf(folder, '_GLIF5_neuron_config.json')           
    return out

def num_of_no_nan(data_list):
    n=0
    for neuron in data_list:
        if np.any(np.array(neuron)=='nan'):
            pass
        else:
            n=n+1
    return n

# load data out of configuration files
neuron_data=[]   
for folder in folders:
    specimen_ID=os.path.basename(folder)[:9]
    pp_file=get_pp_path(folder)
    pp_dict=ju.read(pp_file)
    cre=os.path.basename(folder)[10:]
    #--pre assign everything as a nan incase one of the values doesnt exist
    slope=np.nan
    intercept=np.nan
    th_inf=np.nan
    El=np.nan
    R=np.nan
    R_asc=np.nan
    C=np.nan
    TC1=np.nan
    TC2=np.nan
    a_spike=np.nan
    one_over_b_spike=np.nan
    a_voltage_over_b_voltage=np.nan
    log10_b_voltage=np.nan
    spike_cut_length=np.nan
    tau=np.nan
    delta_th_inf=np.nan
    a_spike=np.nan
    one_over_b_spike=np.nan
    a_voltage_over_b_voltage=np.nan
    log10_b_voltage=np.nan
    
    #--fill in values that exist
    cre=cre 
    slope=pp_dict['spike_cutting']['NOdeltaV']['slope']
    intercept=pp_dict['spike_cutting']['NOdeltaV']['intercept']*1000.
    th_inf=pp_dict['th_inf']['via_Vmeasure']['value']*1e3
    El=pp_dict['El']['El_noise']['measured']['mean']*1e3
    R=pp_dict['resistance']['R_test_list']['mean']/1.e6
    R_asc=pp_dict['resistance']['R_fit_ASC_and_R']['mean']/1.e6
    C=pp_dict['capacitance']['C_test_list']['mean']*1.e12                      
    TC1=(1./pp_dict['asc']['k'][0])*pp_dict['asc']['amp'][0]*1.e12
    TC2=(1./pp_dict['asc']['k'][1])*pp_dict['asc']['amp'][1]*1.e12    
    spike_cut_length=pp_dict['spike_cutting']['NOdeltaV']['cut_length']*1000.*pp_dict['dt_used_for_preprocessor_calculations']
    tau=pp_dict['resistance']['R_test_list']['mean']*pp_dict['capacitance']['C_test_list']['mean']*1.e3
    delta_th_inf=(pp_dict['th_inf']['via_Vmeasure']['value']-pp_dict['El']['El_noise']['measured']['mean'])*1e3
    if np.any([f.endswith('_GLIF2_neuron_config.json') for f in os.listdir(folder)]):  #easiest way to check for spike comp of threshold adaptation 
        a_spike=pp_dict['threshold_adaptation']['a_spike_component_of_threshold']
        one_over_b_spike=1./pp_dict['threshold_adaptation']['b_spike_component_of_threshold']
    if np.any([f.endswith('_GLIF2_neuron_config.json') for f in os.listdir(folder)]):  #easiest way to check for voltage comp of threshold adaptation
        a_voltage_over_b_voltage=pp_dict['threshold_adaptation']['a_voltage_comp_of_thr_from_fitab']/pp_dict['threshold_adaptation']['b_voltage_comp_of_thr_from_fitab']
        log10_b_voltage=np.log10(pp_dict['threshold_adaptation']['b_voltage_comp_of_thr_from_fitab'])
    
    evaluated_dict=get_explained_var(folder)
    opt_th_dict=get_optimized_th_inf(folder)
    
    #--put data in a list
    neuron_data.append([
                    evaluated_dict['ev_LIF'], 
                    evaluated_dict['ev_LIFR'], 
                    evaluated_dict['ev_LIFASC'], 
                    evaluated_dict['ev_LIFRASC'], 
                    evaluated_dict['ev_LIFRASCAT'],     #0-4
                    specimen_ID,                        #5
                    cre,                                
                    slope,                              #7
                    intercept,                      
                    th_inf,                             #9
                    El,
                    R,
                    R_asc,                              #12
                    C,
                    TC1,
                    TC2,                                #15
                    a_spike,
                    one_over_b_spike,                   #17
                    a_voltage_over_b_voltage, 
                    log10_b_voltage,                   #19
                    spike_cut_length,
                    tau,                               #21
                    delta_th_inf,
                    opt_th_dict['th_inf_LIF'],          #23
                    opt_th_dict['th_inf_LIFR'],
                    opt_th_dict['th_inf_LIFASC'],       #25
                    opt_th_dict['th_inf_LIFRASC'],
                    opt_th_dict['th_inf_LIFRASCAT']])

#seperate data into cre_lines
cre_data_dict={}
#--get stats for cre lines
for cre in pos_cre_lines:
    cre_data_dict[cre]={'data':[]}
    for neuron in neuron_data:
        if neuron[6]==cre:
            cre_data_dict[cre]['data'].append(neuron)

#--get stats for excitatory cre lines
for cre in pos_cre_lines:
    cre_data_dict['excitatory']={'data':[]}
    for neuron in neuron_data:
        if neuron[6] in excitatory_cre_lines:
            cre_data_dict['excitatory']['data'].append(neuron)

#--get stats for inhibitory cre lines
for cre in pos_cre_lines:
    cre_data_dict['inhibitory']={'data':[]}
    for neuron in neuron_data:
        if neuron[6] in inhibitory_cre_lines:
            cre_data_dict['inhibitory']['data'].append(neuron)

#Find medians of all cre lines
param_key_idex_pairs=[['slope', 7],
                      ['intercept', 8],
                      ['th_inf', 9],
                      ['El', 10],
                      ['R', 11],
                      ['R_asc', 12],
                      ['C', 13],
                      ['TC1', 14],
                      ['TC2', 15],
                      ['a_spike', 16],
                      ['one_over_b_spike', 17],
                      ['a_voltage_over_b_voltage', 18],
                      ['log10_b_voltage', 19],
                      ['spike_cut_length', 20],
                      ['tau', 21],
                      ['deltath_inf', 22],
                      ['th_inf_LIF', 23],
                      ['th_inf_LIFR', 24],
                      ['th_inf_LIFASC', 25],
                      ['th_inf_LIFRASC', 26],
                      ['th_inf_LIFRASCAT', 27]]
                                        
for cre in cre_data_dict.keys():
    for pair in param_key_idex_pairs:
        cre_data_dict[cre][pair[0]]={'near_median_IDs':np.nan,'median':np.nan, 'quartiles':np.nan, 'data_column':np.nan, 'dist_tightness':np.nan}
        (cre_data_dict[cre][pair[0]]['near_median_IDs'], \
            cre_data_dict[cre][pair[0]]['median'], \
            cre_data_dict[cre][pair[0]]['quartiles'], \
            cre_data_dict[cre][pair[0]]['data_column'], \
            cre_data_dict[cre][pair[0]]['dist_tightness'])=distribution_analysis(cre_data_dict[cre]['data'], pair[1], percentile_boundries=[5, 95])
            
rows=[['inhibitory', 'inhibitory'],
          ['excitatory', 'excitatory'],
          ['Scnn1a-Tg2', 'Scnn1a-Tg2-Cre'],
          ['Nr5a1','Nr5a1-Cre'], 
          ['Scnn1a-Tg3', 'Scnn1a-Tg3-Cre'],  
          ['Rorb', 'Rorb-IRES2-Cre-D'],
          ['Cux2', 'Cux2-CreERT2'],    
          ['Ntsr1','Ntsr1-Cre_GN220'],   
          ['Ctgf','Ctgf-2A-dgCre'],
          ['Rbp4','Rbp4-Cre_KL100'],  
          ['Sst','Sst-IRES-Cre' ],
          ['Pvalb','Pvalb-IRES-Cre'],               
          ['Htr3a', 'Htr3a-Cre_NO152'],
          ['Ndnf', 'Ndnf-IRES2-dgCre'],
          ['Chat', 'Chat-IRES-Cre-neo' ], 
          ['Vip', 'Vip-IRES-Cre'],
          ['Chrna2', 'Chrna2-Cre_OE25'],
          ['Nkx2','Nkx2-1-CreERT2'] ]

columns=[ ['R (M$\\Omega$)', 'R'], 
         ['$\tau$ (ms)', 'tau'],  
         ['C (pF)', 'C'], 
         ['El (mV)', 'El'],
         ['th$_\\infty$ (mV)', 'th_inf'],
         ['$\Delta$ th$_\\infty$ (mV)', 'deltath_inf'],
         ['$R_{ASC} (M\\Omega$)', 'R_asc'], 
         ['$Q_1$ (pC)', 'TC1'],  
         ['$Q_2$ (pC)', 'TC2'], 
         ['$\delta t$ (ms)', 'spike_cut_length'],
         ['slope', 'slope'], 
         ['intercept (mV)', 'intercept'],  
         ['$a_{spike}$ (V)', 'a_spike'],  
         ['1/($b_{spike}$) (s)', 'one_over_b_spike'],
         ['$(a_v)/(b_v)$', 'a_voltage_over_b_voltage'], 
         ['$log_{10}(b_v)$', 'log10_b_voltage']]

#calculate tightness
tight_dict={}
for col in columns[0:10]:
    n=0
    tightness=0
    for index in rows:
        if not np.isnan(cre_data_dict[index[1]][col[1]]['dist_tightness']):
            n=n+1
            tightness=tightness+cre_data_dict[index[1]][col[1]]['dist_tightness']
    tight_dict[col[1]]=tightness/n 

#--create data frames for tables
df1=pd.DataFrame(index=[ind[0] for ind in rows], columns=[col[0] for col in columns[0:8]])
for index in rows:
    for ii, col in enumerate(columns[0:8]):
        if ii==0 and len(cre_data_dict[index[1]][col[1]]['data_column'])<7:
            df1=df1.drop(index[0])
        elif len(cre_data_dict[index[1]][col[1]]['data_column'])<7:
            pass
        else:
#            print index, col, len(cre_data_dict[index[1]][col[1]]['data_column']), cre_data_dict[index[1]][col[1]]['median'], cre_data_dict[index[1]][col[1]]['quartiles'][0], cre_data_dict[index[1]][col[1]]['quartiles'][1]

            value ="\shortstack{%0.3g" % cre_data_dict[index[1]][col[1]]['median']+"\\\\{[}"+"%0.3g" %cre_data_dict[index[1]][col[1]]['quartiles'][0]+", %0.3g" %cre_data_dict[index[1]][col[1]]['quartiles'][1] + '{]}}'
            df1.set_value(index[0], col[0], value)          
print df1.to_latex(escape=False)

df2=pd.DataFrame(index=[ind[0] for ind in rows], columns=[col[0] for col in columns[8:]])
for index in rows:
    for ii, col in enumerate(columns[8:]):
        if ii==0 and len(cre_data_dict[index[1]][col[1]]['data_column'])<7:
            df2=df2.drop(index[0])
        elif len(cre_data_dict[index[1]][col[1]]['data_column'])<7:
            pass
        else:
            #print index, col, len(cre_data_dict[index[1]][col[1]]['data_column']), cre_data_dict[index[1]][col[1]]['median'], cre_data_dict[index[1]][col[1]]['quartiles'][0], cre_data_dict[index[1]][col[1]]['quartiles'][1]
            value ="\shortstack{%0.3g" % cre_data_dict[index[1]][col[1]]['median']+"\\\\{[}"+"%0.3g" %cre_data_dict[index[1]][col[1]]['quartiles'][0]+", %0.3g" %cre_data_dict[index[1]][col[1]]['quartiles'][1] + '{]}}'
            df2.set_value(index[0], col[0], value)          
print df2.to_latex(escape=False)

#--find neuron that are within the percentiles for all parameters (cononical neurons)
for cre in cre_data_dict.keys():
    cre_data_dict[cre]['possible_cononical_IDs']=[]
    possible_con_neuronIDs=set.intersection(
                             set(cre_data_dict[cre]['slope']['near_median_IDs']),
                             set(cre_data_dict[cre]['intercept']['near_median_IDs']),
#                             set(cre_data_dict[cre]['th_inf']['near_median_IDs']),
                             set(cre_data_dict[cre]['El']['near_median_IDs']),
                             set(cre_data_dict[cre]['R']['near_median_IDs']),
                             set(cre_data_dict[cre]['R_asc']['near_median_IDs']),
                             set(cre_data_dict[cre]['C']['near_median_IDs']),
                             set(cre_data_dict[cre]['TC1']['near_median_IDs']),
                             set(cre_data_dict[cre]['TC2']['near_median_IDs']),
                             set(cre_data_dict[cre]['a_spike']['near_median_IDs']),
                             set(cre_data_dict[cre]['one_over_b_spike']['near_median_IDs']),
                             set(cre_data_dict[cre]['a_voltage_over_b_voltage']['near_median_IDs']),
                             set(cre_data_dict[cre]['log10_b_voltage']['near_median_IDs']),
                             set(cre_data_dict[cre]['spike_cut_length']['near_median_IDs']),
                             set(cre_data_dict[cre]['th_inf_LIF']['near_median_IDs']),
                             set(cre_data_dict[cre]['th_inf_LIFR']['near_median_IDs']),
                             set(cre_data_dict[cre]['th_inf_LIFASC']['near_median_IDs']),
                             set(cre_data_dict[cre]['th_inf_LIFRASC']['near_median_IDs']),
                             set(cre_data_dict[cre]['th_inf_LIFRASCAT']['near_median_IDs']))
    
# organize the overlaping parameter so can look for the best explained variance (cononical)
    for id in possible_con_neuronIDs:
        for neuron in neuron_data:
            if neuron[5]==id:
                cre_data_dict[cre]['possible_cononical_IDs'].append([neuron[5], neuron[0], neuron[1], neuron[2], neuron[3], neuron[4]]) 
                'there should be no neurons here that have any nans'

#find neurons that have all level models and all parameters are in the percentage intervals
np.set_printoptions(precision=3)
for cre in cre_data_dict.keys():
    n=num_of_no_nan(cre_data_dict[cre]['data'])
    if n>7:
        print cre, 'number of neurons with all entries:', n, '. However only the following neurons are within 5% and 95% of all values' 
        best_model_win_neuron=[]
        maybe_cononical=[]
        if len(cre_data_dict[cre]['possible_cononical_IDs'])>0:
            for neuron in cre_data_dict[cre]['possible_cononical_IDs']:
                if np.all(np.array(neuron)!='nan'):
                    maybe_cononical.append(neuron)
                    best_model_win_neuron.append([neuron[0], np.max(np.array(neuron[1:]))])
                    print ' ', neuron[0], '& %.3f' % neuron[1], '& %.3f' % neuron[2], '& %.3f' % neuron[3], '& %.3f' % neuron[4], '& %.3f' % neuron[5], '\\\\'
            best_row=np.where(np.max(np.array([value[1] for value in best_model_win_neuron]))==np.array([value[1] for value in best_model_win_neuron]))[0][0]
            print '\tBEST MODEL within parameter percentiles and has all level models ', best_model_win_neuron[best_row], 'which is GLIF model', 1+np.where(np.array(maybe_cononical[best_row][1:])==np.max(np.array(maybe_cononical[best_row][1:])))[0][0]
        else:
            print 'There are no models for', cre, 'that fit into the percentiles'
print '\n'

#Find highest scoring model for each cre line that has all models (can be outside of parameters)
for cre in cre_data_dict.keys():
    best_model_win_neuron=[]
    has_model=[]
    if len(cre_data_dict[cre]['data'])>=7:
        for neuron in cre_data_dict[cre]['data']:
            if np.all(~np.isnan(neuron[0:5])):  #if all five models exist
                has_model.append([neuron[5], neuron[0:5]])
                best_model_win_neuron.append([neuron[5], max_discount_nan(neuron[0:5])])    
        best_row=np.where(np.max(np.array([value[1] for value in best_model_win_neuron]))==np.array([value[1] for value in best_model_win_neuron]))[0][0]
        if has_model[best_row][0] in [neuron[0] for neuron in cre_data_dict[cre]['possible_cononical_IDs']]:
            print cre, ': BEST ExpVar of neurons with ALL MODEL LEVELS', has_model[best_row], 'in GLIF level', 1+np.where(np.array(has_model[best_row][1])==max_discount_nan(np.array(has_model[best_row][1])))[0][0], 'which IS within the parameter percentiles'
#            print cre,':', has_model[best_row][0], 'GLIF', 1+np.where(np.array(has_model[best_row][1])==max_discount_nan(np.array(has_model[best_row][1])))[0][0]
        else:
            print cre, ': BEST ExpVar of neurons with ALL MODEL LEVELS', has_model[best_row], 'in GLIF level', 1+np.where(np.array(has_model[best_row][1])==max_discount_nan(np.array(has_model[best_row][1])))[0][0], 'which IS NOT within the parameter percentiles'
#            print cre,':', has_model[best_row][0], 'GLIF', 1+np.where(np.array(has_model[best_row][1])==max_discount_nan(np.array(has_model[best_row][1])))[0][0]

print '\n'
            
#Find highest scoring model for each cre line (can be outside of parameters and not have all models)
for cre in cre_data_dict.keys():
    best_model_win_neuron=[]
    has_model=[]
    if len(cre_data_dict[cre]['data'])>=7:
        for neuron in cre_data_dict[cre]['data']:
            if np.any(~np.isnan(neuron[0:5])):  #if there are any models
                has_model.append([neuron[5], neuron[0:5]])
                best_model_win_neuron.append([neuron[5], max_discount_nan(neuron[0:5])])    
        best_row=np.where(np.max(np.array([value[1] for value in best_model_win_neuron]))==np.array([value[1] for value in best_model_win_neuron]))[0][0]
        print cre, ': BEST ExpVar OVERALL', has_model[best_row], 'in GLIF level', 1+np.where(np.array(has_model[best_row][1])==max_discount_nan(np.array(has_model[best_row][1])))[0][0]
 
print 'finished'
