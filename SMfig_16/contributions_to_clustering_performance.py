'''Takes data created via expVar_level_calc_stats.py and calculated cluster 
similarity and plots their relationship .
'''

import socket
import copy
import os
import numpy as np
import allensdk.core.json_utilities as ju
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 16})
import scipy.stats as stats
import statsmodels.api as sm
import pandas as pd
import sys
relative_path=os.path.dirname(os.getcwd())
sys.path.append(os.path.join(relative_path, 'libraries'))
from data_library import get_ev_from_folder, pos_cre_lines, excitatory_cre_lines, inhibitory_cre_lines
from pub_plot_library import color_dict
import pickle

def non_nan_data(neuron_data, index):
    x=[]
    for neuron in neuron_data:
        if not np.isnan(neuron[index]):
            x.append(neuron[index])

    return x

def color_y_axis(ax, color):
    """Color your axes."""
    for t in ax.get_yticklabels():
        t.set_color(color)
    ax.yaxis.label.set_color('red')
    return None

cluster_similarity={'binary_split':{'GLIF':{'ARI': [0.074543, 0.075066, 0.118164, 0.179298],
                                            'AVOI': [0.726094, 0.83255, 1.244296, 1.331085]},
                                    'all_features':{'ARI': 0.211155,
                                                    'AVOI':1.743444},
                                    'no_spike_features':{'ARI': 0.132005,
                                                         'AVOI':1.182007},
                                    'GLIF_w_spike':{'ARI':[0.148347, 0.204573, 0.153283, 0.181506],
                                                    'AVOI':[1.400633, 1.483436, 1.407738, 1.580448]}},

                    'affinity_prop':{'GLIF': {'ARI': [0.074988, 0.081342, 0.121039, 0.14567],
                                              'AVOI': [0.629492, 0.662893, 1.132159, 1.148379]},
                                    'all_features':{'ARI': 0.128679,
                                                    'AVOI':1.164425},
                                    'no_spike_features':{'ARI': 0.111107,
                                                         'AVOI': 1.117475},
                                    'GLIF_w_spike':{'ARI':[0.154543, 0.160928, 0.166824, 0.168998],
                                                    'AVOI':[1.168143, 1.21287, 1.297733, 1.313232]}}}
                    

v_diff_medians=np.array([0.1740242398143918, 
                         0.1205088716374338, 
                         0.23375818446206129,
                         0.089098920847300381]) # GLIF 1 though 4 calculated via sub_v_plots.py

num_param_fit=np.array([5, 7, 7, 9])

data_path=os.path.join(relative_path, 'MTfig_5_MTtable_3_SMFigs_8_9_10', 'saved_data', 'stats_out.pkl')
neuron_data, cre_data_dict, stats_out=pickle.load(open(data_path, "rb"))  

#list of mean values of explained variance for the different values 
st_expVar=np.array([np.mean(non_nan_data(neuron_data, 2)),
                    np.mean(non_nan_data(neuron_data, 3)),        
                    np.mean(non_nan_data(neuron_data, 4)),
                    np.mean(non_nan_data(neuron_data, 5))])

#----------------------------------------------------
#------------spike time performance------------------
#----------------------------------------------------
def make_plots(key, title):
    fig, ax1 = plt.subplots()
    ax1.plot(st_expVar, cluster_similarity[key]['GLIF']['ARI'], '-k.', LW=2, MS=16)
    ax2=ax1.twinx()
    ax2.plot(st_expVar, cluster_similarity[key]['GLIF']['AVOI'], '-r.', LW=2, MS=16)
    plt.title(title)
    ax1.set_xlabel('Median % Explained Variance')
    ax1.set_ylabel('Adjusted Rand Index')
    ax2.set_ylabel('Adjusted VOI score')
    color_y_axis(ax2, 'r')
    

    print 'Stats of spike time performance versus clustering'
    print stats.linregress(st_expVar, cluster_similarity[key]['GLIF']['ARI'])
    print stats.linregress(st_expVar, cluster_similarity[key]['GLIF']['AVOI'])
    
    #---------------------------------------------------
    # -----------subthreshold voltage--------------------
    #----------------------------------------------------
    
    fig, ax3 = plt.subplots()
    ax3.plot(v_diff_medians, cluster_similarity[key]['GLIF']['ARI'], '-k.', LW=2, MS=16)
    ax4=ax3.twinx()
    ax4.plot(v_diff_medians, cluster_similarity[key]['GLIF']['AVOI'], '-r.', LW=2, MS=16)
    ax3.set_xlabel('Subthreshold Voltage Difference Measure')
    ax3.set_ylabel('Adjusted Rand Index')
    ax4.set_ylabel('Adjusted VOI score')
    plt.title(title)
    color_y_axis(ax4, 'r')
    
    #---------------------------------------------------
    # -----------number of parameters fit---------------
    #----------------------------------------------------
    fig, ax5 = plt.subplots()
    ax5.plot(num_param_fit, cluster_similarity[key]['GLIF']['ARI'], '-k.', LW=2, MS=16)
    ax6=ax5.twinx()
    ax6.plot(num_param_fit, cluster_similarity[key]['GLIF']['AVOI'], '-r.', LW=2, MS=16)
    #add features
    ax5.plot(14, cluster_similarity[key]['all_features']['ARI'], 'k.', LW=2, MS=16)
    ax6.plot(14, cluster_similarity[key]['all_features']['AVOI'], 'r.', LW=2, MS=16)
    ax5.plot(10, cluster_similarity[key]['no_spike_features']['ARI'], 'k.', LW=2, MS=16)
    ax6.plot(10, cluster_similarity[key]['no_spike_features']['AVOI'], 'r.', LW=2, MS=16)

    ax5.plot(np.array(num_param_fit)+4, cluster_similarity[key]['GLIF_w_spike']['ARI'], '-k.', LW=2, MS=16)
    ax6.plot(np.array(num_param_fit)+4, cluster_similarity[key]['GLIF_w_spike']['AVOI'], '-r.', LW=2, MS=16)
    
    
    plt.title(title)
    ax5.set_xlabel('Number of Cluster Parameters')
    ax5.set_ylabel('Adjusted Rand Index')
    ax6.set_ylabel('Adjusted VOI score')
    color_y_axis(ax6, 'r')

    
    
    print 'Stats of num param fit versus clustering'
    print stats.linregress(np.append(num_param_fit, np.array([10, 14])), 
                           cluster_similarity[key]['GLIF']['ARI']+
                            [cluster_similarity[key]['no_spike_features']['ARI'],
                            cluster_similarity[key]['all_features']['ARI']])
    print stats.linregress(np.append(num_param_fit, np.array([10, 14])), 
                           cluster_similarity[key]['GLIF']['AVOI']+
                            [cluster_similarity[key]['no_spike_features']['AVOI'],
                            cluster_similarity[key]['all_features']['AVOI']])
    
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax3D = fig.add_subplot(111, projection='3d')
    ax3D.plot(st_expVar, v_diff_medians, cluster_similarity[key]['GLIF']['ARI'])
    #ax3D.plot(st_expVar, v_diff_medians, cluster_similarity[key]['GLIF']['VOI'])
    
    return ax1, ax2, ax3, ax4, ax5, ax6


###--------------Binary spitting------------------
ax1, ax2, ax3, ax4, ax5, ax6=make_plots('binary_split', 'Binary Splitting')

#annotate plots
ax2.annotate(r'GLIF$_1$', xy=(69.7, .735), horizontalalignment='left',fontsize=18)
ax2.annotate(r'GLIF$_2$', xy=(67.3, .74), horizontalalignment='left',fontsize=18)
ax2.annotate(r'GLIF$_3$', xy=(70.2, 1.12), horizontalalignment='left',fontsize=18)
ax2.annotate(r'GLIF$_4$', xy=(73.9, 1.25), horizontalalignment='left',fontsize=18)

ax4.annotate(r'GLIF$_1$', xy=(.18, .73), horizontalalignment='left',fontsize=18)
ax4.annotate(r'GLIF$_2$', xy=(.095, .8), horizontalalignment='left',fontsize=18)
ax4.annotate(r'GLIF$_3$', xy=(.22, 1.14), horizontalalignment='left',fontsize=18)
ax4.annotate(r'GLIF$_4$', xy=(.082, 1.26), horizontalalignment='left',fontsize=18)

ax6.annotate(r'GLIF$_1$', xy=(5.1, .8), horizontalalignment='left',fontsize=16)
ax6.annotate(r'GLIF$_2$', xy=(7.2, .75), horizontalalignment='left',fontsize=16)
ax6.annotate(r'GLIF$_3$', xy=(6, 1.12), horizontalalignment='left',fontsize=16)
ax6.annotate(r'GLIF$_4$', xy=(7.5, 1.4), horizontalalignment='left',fontsize=16)
ax6.annotate('features - SS', xy=(9.5, 1.08), horizontalalignment='left',fontsize=16)
ax6.annotate('all features', xy=(11.8, 1.725), horizontalalignment='left',fontsize=16) 

ax6.annotate(r'GLIF$_1$'+ '\n+ SS', xy=(9.5, 1.25), horizontalalignment='left',fontsize=16)
ax6.annotate(r'GLIF$_2$'+ '\n+ SS', xy=(9.7, 1.63), horizontalalignment='left',fontsize=16)
ax6.annotate(r'GLIF$_3$'+ '\n+ SS', xy=(11.3, 1.2), horizontalalignment='left',fontsize=16)
ax6.annotate(r'GLIF$_4$'+ '\n+ SS', xy=(12.7, 1.35), horizontalalignment='left',fontsize=16)


#--------------Affinity propagation------------------
ax1, ax2, ax3, ax4, ax5, ax6=make_plots('affinity_prop', 'Affinity Propagation')

#annotate plots
ax2.annotate(r'GLIF$_1$', xy=(69.7, .62), horizontalalignment='left',fontsize=18)
ax2.annotate(r'GLIF$_2$', xy=(67.2, .73), horizontalalignment='left',fontsize=18)
ax2.annotate(r'GLIF$_3$', xy=(70.2, 1.03), horizontalalignment='left',fontsize=18)
ax2.annotate(r'GLIF$_4$', xy=(73.9, 1.1), horizontalalignment='left',fontsize=18)

ax4.annotate(r'GLIF$_1$', xy=(.18, .62), horizontalalignment='left',fontsize=18)
ax4.annotate(r'GLIF$_2$', xy=(.095, .67), horizontalalignment='left',fontsize=18)
ax4.annotate(r'GLIF$_3$', xy=(.22, 1.05), horizontalalignment='left',fontsize=18)
ax4.annotate(r'GLIF$_4$', xy=(.082, 1.1), horizontalalignment='left',fontsize=18)

ax6.annotate(r'GLIF$_1$', xy=(5.1, .75), horizontalalignment='left',fontsize=16)
ax6.annotate(r'GLIF$_2$', xy=(7.2, .7), horizontalalignment='left',fontsize=16)
ax6.annotate(r'GLIF$_3$', xy=(6, 1.05), horizontalalignment='left',fontsize=16)
ax6.annotate(r'GLIF$_4$', xy=(8.5, 1.09), horizontalalignment='left',fontsize=16)
ax6.annotate('features - SS', xy=(9.5, 1.01), horizontalalignment='left',fontsize=16)
ax6.annotate('all features', xy=(12, 1.1), horizontalalignment='left',fontsize=16) 

ax6.annotate(r'GLIF$_1$'+ '\n+ SS', xy=(8, 1.2), horizontalalignment='left',fontsize=16)
ax6.annotate(r'GLIF$_2$'+ '\n+ SS', xy=(11.2, 1.18), horizontalalignment='left',fontsize=16)
ax6.annotate(r'GLIF$_3$'+ '\n+ SS', xy=(10, 1.3), horizontalalignment='left',fontsize=16)
ax6.annotate(r'GLIF$_4$'+ '\n+ SS', xy=(12.9, 1.21), horizontalalignment='left',fontsize=16)

plt.show()