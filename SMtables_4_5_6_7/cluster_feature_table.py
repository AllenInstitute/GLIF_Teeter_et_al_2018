'''
Loads cluster csv file (cluster_ids_log_20170727_glif4.csv) 
to specify cluster identity and then loads feature data csv 
(features_7_27_17.csv) and reports feature distributions for 
clusters in LaTeX table format.
'''

import os
import numpy as np
import allensdk.core.json_utilities as ju
import sys
relative_path=os.path.dirname(os.getcwd())
sys.path.append(os.path.join(relative_path, 'libraries'))
from pub_plot_library import distribution_analysis
import scipy.stats as stats
import pandas as pd
import csv

def transform_csv_read_to_data(csv_data, ind_list):
    '''If there are empty entries in the csv make them nan.  
    Change specified columns to floats (as opposed to strings).
    Inputs: 
        csv_data: list of list of strings 
            contains data in the csv file (header should be excluded here)
        ind_list: list of integers 
            specifies which columns to change to floats.
    Returns:
        csv_data: list of list of strings
            data with desired columns converted to nans and floats.
            '''
    
    for jj, row in enumerate(csv_data):
        for ii, value in enumerate(row):
            if value=='':
                row[ii]=np.NAN
            elif ii in ind_list:
                try:
                    row[ii]=float(value)
                except:
                    pass
    return csv_data

#load the cluster id's 
csv_data=[]
with open('cluster_ids_log_20170727_glif4.csv', 'rb') as csvfile:
    stuff = csv.reader(csvfile)
    for ii, row in enumerate(stuff):
        if ii !=0: #skip the header
            csv_data.append(row)

#find specimen ids in the individual GLIF clusters
clusters=set([row[1] for row in csv_data])
spID_clust_dict={}
for clust in clusters:            
    spID_clust_dict[clust]=[]
    for row in csv_data:
        if clust == row[1]:
            spID_clust_dict[clust].append(row[0])
        
# get feature data from from previously spread sheet
csv_data=[]
with open('features_7_27_17.csv', 'rU') as csvfile:
    stuff = csv.reader(csvfile)
    for ii, row in enumerate(stuff):
        if ii==0:
            header=row
        if ii !=0:
            csv_data.append(row)


#the csv file is reading as strings so the necessary rows must be changed to floats
feature_data=transform_csv_read_to_data(csv_data, range(5,21))

cluster_data_dict={}
for clust in spID_clust_dict.keys():
    cluster_data_dict[clust]={'data':[]}
    for neuron in feature_data:
        if neuron[1] in spID_clust_dict[clust]:
            cluster_data_dict[clust]['data'].append(neuron)


features=['tau', 'ri', 'vrest', 'threshold_i_long_square','threshold_v_long_square', 
          'peak_v_long_square', 'fast_trough_v_long_square', 
          'trough_v_long_square', 'upstroke_downstroke_ratio_long_square',
          'upstroke_downstroke_ratio_short_square', 'sag', 'avg_isi', 'f_i_curve_slope',
          'latency', 'adaptation', 'max_burstiness_across_sweeps']
                      
for cluster in cluster_data_dict.keys():
    for feature in features:
        cluster_data_dict[cluster][feature]={'near_median_IDs':np.nan,'median':np.nan, 'quartiles':np.nan, 'data_column':np.nan, 'dist_tightness':np.nan}
        index_of_data=header.index(feature)
        if index_of_data==20:
            pass
        (cluster_data_dict[cluster][feature]['near_median_IDs'], \
            cluster_data_dict[cluster][feature]['median'], \
            cluster_data_dict[cluster][feature]['quartiles'], \
            cluster_data_dict[cluster][feature]['data_column'], \
            cluster_data_dict[cluster][feature]['dist_tightness'])=distribution_analysis(cluster_data_dict[cluster]['data'], index_of_data, percentile_boundries=[5, 95])
            
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

# specify column titles for output table
columns=[['$\tau_\mathrm{m} (ms)$', 'tau'],
         ['R$_i$ (M$\Omega$)','ri'],
         ['V$_{rest}$ (mV)', 'vrest'],
         ['I$_{threshold}\dagger$ (pA)','threshold_i_long_square'],
         ['V$_{threshold}\dagger$ (mV)', 'threshold_v_long_square'],
         ['V$_{peak}\dagger$ (mV)','peak_v_long_square'],
         ['V$_{fast trough}\dagger$ (mV)','fast_trough_v_long_square'],
         ['V$_{trough}\dagger$','trough_v_long_square'],
         ['up:downstroke$\dagger$', 'upstroke_downstroke_ratio_long_square'],
         ['up:downstroke*', 'upstroke_downstroke_ratio_short_square'],
         ['sag (mV)', 'sag'],
         ['Avg. ISI (ms)', 'avg_isi'],
         ['\shortstack{$f$-$I$ curve slope \\\\(spikes $\textrm{s}^{-1} \textrm{pA}^{-1}$)}', 'f_i_curve_slope'],
         ['latency (s)', 'latency'],
         ['Adaptation', 'adaptation'],
         ['\shortstack{Max. burst\\\\index}', 'max_burstiness_across_sweeps']]
                                 

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
            #print index, col, len(cluster_data_dict[index[1]][col[1]]['data_column']), cluster_data_dict[index[1]][col[1]]['median'], cluster_data_dict[index[1]][col[1]]['quartiles'][0], cluster_data_dict[index[1]][col[1]]['quartiles'][1]
            value ="\shortstack{%0.3g" % cluster_data_dict[index[1]][col[1]]['median']+"\\\\{[}"+"%0.3g" %cluster_data_dict[index[1]][col[1]]['quartiles'][0]+", %0.3g" %cluster_data_dict[index[1]][col[1]]['quartiles'][1] + '{]}}'
            df1.set_value(index[0], col[0], value)

print 'First data frame'           
print df1.to_latex(escape=False)

#--------------------------------2nd data frame-----------------------------------------

df2=pd.DataFrame(index=[ind[0] for ind in rows], columns=[col[0] for col in columns[8:16]])
for index in rows:
    for ii, col in enumerate(columns[8:16]):
        if ii==0 and len(cluster_data_dict[index[1]][col[1]]['data_column'])<5:
            df2=df2.drop(index[0])
        elif len(cluster_data_dict[index[1]][col[1]]['data_column'])<5:
            pass
        else:
            #print index, col, len(cluster_data_dict[index[1]][col[1]]['data_column']), cluster_data_dict[index[1]][col[1]]['median'], cluster_data_dict[index[1]][col[1]]['quartiles'][0], cluster_data_dict[index[1]][col[1]]['quartiles'][1]
            value ="\shortstack{%0.3g" % cluster_data_dict[index[1]][col[1]]['median']+"\\\\{[}"+"%0.3g" %cluster_data_dict[index[1]][col[1]]['quartiles'][0]+", %0.3g" %cluster_data_dict[index[1]][col[1]]['quartiles'][1] + '{]}}'
            df2.set_value(index[0], col[0], value)
  
print 'Second data frame'           
print df2.to_latex(escape=False)



