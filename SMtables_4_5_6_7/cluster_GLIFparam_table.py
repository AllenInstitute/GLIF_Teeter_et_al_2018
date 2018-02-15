'''
Loads cluster csv file (cluster_ids_log_20170727_glif4.csv) 
to specify cluster identity and reports GLIF param distributions 
for the clusters in LaTeX table format.. 
'''
import os
import numpy as np
import allensdk.core.json_utilities as ju
import sys
relative_path=os.path.dirname(os.getcwd())
sys.path.append(os.path.join(relative_path, 'libraries'))
from pub_plot_library import distribution_analysis
import allensdk.core.json_utilities as ju
import pandas as pd
import csv
from data_library import get_pp_path


#load the cluster id's 
csv_data=[]
with open('cluster_ids_log_20170727_glif4.csv', 'rb') as csvfile:
    stuff = csv.reader(csvfile)
    for ii, row in enumerate(stuff):
        if ii !=0: #skip the header
            csv_data.append(row)

#find specimen ids in the different clusters
clusters=set([row[1] for row in csv_data])
spID_clust_dict={}
for clust in clusters:            
    spID_clust_dict[clust]=[]
    for row in csv_data:
        if clust == row[1]:
            spID_clust_dict[clust].append(row[0])

# load data out of configuration files
folder_path=os.path.join(relative_path,'mouse_struc_data_dir')
folders=[os.path.join(folder_path, f) for f in  os.listdir(folder_path)]

neuron_data=[]   
for folder in folders:
    specimen_id=os.path.basename(folder)[:9]
    cre=os.path.basename(folder)[10:]        
    pp_file=get_pp_path(folder)
    neuron_dict=ju.read(pp_file)
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
    
    #--fill in values that are there
    if cre == 'Ntsr1-CrePOS':
        cre ='Ntsr1-Cre_GN220POS' 
    slope=neuron_dict['spike_cutting']['NOdeltaV']['slope']
    intercept=neuron_dict['spike_cutting']['NOdeltaV']['intercept']*1000.
    th_inf=neuron_dict['th_inf']['via_Vmeasure']['value']*1e3
    El=neuron_dict['El']['El_noise']['measured']['mean']*1e3
    R=neuron_dict['resistance']['R_test_list']['mean']/1.e6
    R_asc=neuron_dict['resistance']['R_fit_ASC_and_R']['mean']/1.e6
    C=neuron_dict['capacitance']['C_test_list']['mean']*1.e12                      
    TC1=(1./neuron_dict['asc']['k'][0])*neuron_dict['asc']['amp'][0]*1.e12
    TC2=(1./neuron_dict['asc']['k'][1])*neuron_dict['asc']['amp'][1]*1.e12    
    spike_cut_length=neuron_dict['spike_cutting']['NOdeltaV']['cut_length']*1000.*neuron_dict['dt_used_for_preprocessor_calculations']
    tau=neuron_dict['resistance']['R_test_list']['mean']*neuron_dict['capacitance']['C_test_list']['mean']*1.e3
    delta_th_inf=(neuron_dict['th_inf']['via_Vmeasure']['value']-neuron_dict['El']['El_noise']['measured']['mean'])*1e3
    if (neuron_dict['threshold_adaptation']['a_spike_component_of_threshold'] is not None) and \
        (neuron_dict['threshold_adaptation']['b_spike_component_of_threshold'] is not None):
        a_spike=neuron_dict['threshold_adaptation']['a_spike_component_of_threshold']
        one_over_b_spike=1./neuron_dict['threshold_adaptation']['b_spike_component_of_threshold']
    if neuron_dict['threshold_adaptation']['a_voltage_comp_of_thr_from_fitab'] is not None and \
        neuron_dict['threshold_adaptation']['b_voltage_comp_of_thr_from_fitab'] is not None:
        a_voltage_over_b_voltage=neuron_dict['threshold_adaptation']['a_voltage_comp_of_thr_from_fitab']/neuron_dict['threshold_adaptation']['b_voltage_comp_of_thr_from_fitab']
        log10_b_voltage=np.log10(neuron_dict['threshold_adaptation']['b_voltage_comp_of_thr_from_fitab'])
    #--put data in a list
    neuron_data.append([    specimen_id,                        #0
                            cre,                                
                            slope,                              #2
                            intercept,                      
                            th_inf,                             #4
                            El,
                            R,
                            R_asc,                              #7
                            C,
                            TC1,
                            TC2,                                #10
                            a_spike,
                            one_over_b_spike,                   #12
                            a_voltage_over_b_voltage, 
                            log10_b_voltage,                   #14
                            spike_cut_length,
                            tau,                               #16
                            delta_th_inf])
        
#separate data into clusters
cluster_data_dict={}
for clust in spID_clust_dict.keys():
    cluster_data_dict[clust]={'data':[]}
    for neuron in neuron_data:
        if neuron[0] in spID_clust_dict[clust]:
            cluster_data_dict[clust]['data'].append(neuron)


#Find medians of all cre lines
param_key_idex_pairs=[['slope', 2],
                      ['intercept', 3],
                      ['th_inf', 4],
                      ['El', 5],
                      ['R', 6],
                      ['R_asc', 7],
                      ['C', 8],
                      ['TC1', 9],
                      ['TC2', 10],
                      ['a_spike', 11],
                      ['one_over_b_spike', 12],
                      ['a_voltage_over_b_voltage', 13],
                      ['log10_b_voltage', 14],
                      ['spike_cut_length', 15],
                      ['tau', 16],
                      ['deltath_inf', 17]]
                      
for cluster in cluster_data_dict.keys():
    for pair in param_key_idex_pairs:
        cluster_data_dict[cluster][pair[0]]={'near_median_IDs':np.nan,'median':np.nan, 'quartiles':np.nan, 'data_column':np.nan, 'dist_tightness':np.nan}
        (cluster_data_dict[cluster][pair[0]]['near_median_IDs'], \
            cluster_data_dict[cluster][pair[0]]['median'], \
            cluster_data_dict[cluster][pair[0]]['quartiles'], \
            cluster_data_dict[cluster][pair[0]]['data_column'], \
            cluster_data_dict[cluster][pair[0]]['dist_tightness'])=distribution_analysis(cluster_data_dict[cluster]['data'], pair[1], percentile_boundries=[5, 95])

# Get cluster IDs for each specimen id
# NOTE THESE VALUES ARE SPECIFIC FOR cluster_ids_log_20170727_glif4.csv              
rows=[['1','1_2_2_2_2'],
    ['2','1_2_2_2_1'],
    ['3','1_2_2_1_2'],
    ['4','1_2_2_1_1'],
    ['5','1_2_1'],
    ['6','1_1_2_2_2'],
    ['7','1_1_2_2_1'],
    ['8','1_1_2_1_2'],
    ['9','1_1_2_1_1'],
    ['10','1_1_1_2_2_2_2'],
    ['11','1_1_1_2_2_2_1'],
    ['12','1_1_1_2_2_1'],
    ['13','1_1_1_2_1_2_2'],
    ['14','1_1_1_2_1_2_1'],
    ['15','1_1_1_2_1_1'],
    ['16','1_1_1_1']]


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

# Output data is separated into two tables to fit onto pages
#--------------------------------1st data frame-----------------------------------------
df1=pd.DataFrame(index=[ind[0] for ind in rows], columns=[col[0] for col in columns[0:8]])
for index in rows:
    for ii, col in enumerate(columns[0:8]):
        if ii==0 and len(cluster_data_dict[index[1]][col[1]]['data_column'])<5:
            df1=df1.drop(index[0])
        elif len(cluster_data_dict[index[1]][col[1]]['data_column'])<5:
            pass
        else:
#            print index, col, len(cluster_data_dict[index[1]][col[1]]['data_column']), cluster_data_dict[index[1]][col[1]]['median'], cluster_data_dict[index[1]][col[1]]['quartiles'][0], cluster_data_dict[index[1]][col[1]]['quartiles'][1]
            value ="\shortstack{%0.3g" % cluster_data_dict[index[1]][col[1]]['median']+"\\\\{[}"+"%0.3g" %cluster_data_dict[index[1]][col[1]]['quartiles'][0]+", %0.3g" %cluster_data_dict[index[1]][col[1]]['quartiles'][1] + '{]}}'
            df1.set_value(index[0], col[0], value)          
print df1.to_latex(escape=False)

#--------------------------------2nd data frame-----------------------------------------
df2=pd.DataFrame(index=[ind[0] for ind in rows], columns=[col[0] for col in columns[8:]])
for index in rows:
    for ii, col in enumerate(columns[8:]):
        if ii==0 and len(cluster_data_dict[index[1]][col[1]]['data_column'])<5:
            df2=df2.drop(index[0])
        elif len(cluster_data_dict[index[1]][col[1]]['data_column'])<5:
            pass
        else:
#            print index, col, len(cluster_data_dict[index[1]][col[1]]['data_column']), cluster_data_dict[index[1]][col[1]]['median'], cluster_data_dict[index[1]][col[1]]['quartiles'][0], cluster_data_dict[index[1]][col[1]]['quartiles'][1]
            value ="\shortstack{%0.3g" % cluster_data_dict[index[1]][col[1]]['median']+"\\\\{[}"+"%0.3g" %cluster_data_dict[index[1]][col[1]]['quartiles'][0]+", %0.3g" %cluster_data_dict[index[1]][col[1]]['quartiles'][1] + '{]}}'
            df2.set_value(index[0], col[0], value) 
print df2.to_latex(escape=False)