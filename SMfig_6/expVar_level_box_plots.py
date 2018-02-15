'''Takes data created via expVar_level_calc_stats.py and makes corresponding plots.
'''

import socket
import copy
import os
import numpy as np
import allensdk.core.json_utilities as ju
import matplotlib
if socket.gethostname()=="ibs-corinnet-ux1.corp.alleninstitute.org":
    matplotlib.use('tkagg')
elif socket.gethostname()=="corinne-teeters-air.local":
    matplotlib.use('Qt4Agg')
else:
    matplotlib.use('agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 16})
import matplotlib.image as mpimg
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

#--make a figure
def make_box_plot(neuron_data, the_color, p=[1,2,3,4,5]):
    bp=plt.boxplot([non_nan_data(neuron_data, 2),
                    non_nan_data(neuron_data, 3),
                    non_nan_data(neuron_data, 4),
                    non_nan_data(neuron_data, 5),
                    non_nan_data(neuron_data, 6)], notch=0, whis=[5, 95], sym='.', positions=p, widths=.08)
    plt.setp(bp['boxes'], color=the_color, linewidth=3)
    plt.setp(bp['whiskers'], color=the_color, linewidth=3)
    plt.setp(bp['fliers'], mfc=the_color, markersize=12)
    plt.setp(bp['medians'], color=the_color, linewidth=3)
    plt.setp(bp['caps'], color=the_color, linewidth=3)

# Custom function to draw the diff bars
def label_diff(x1, x2, y, text, c):
    connectionstyle='bar, fraction='+str(.05/(x2-x1))
    props = {'connectionstyle':connectionstyle,'arrowstyle':'-','lw':2, 'color':c}
    plt.annotate(text, xy=((x2-x1)/2.+x1,y-.5),  color=c, xycoords='data', horizontalalignment='center', verticalalignment='bottom', fontsize=12)
    plt.annotate('', xy=(x1,y), xycoords='data', xytext=(x2,y),textcoords='data', arrowprops=props)

def plot_level_of_sig(stats_out, cre, level_compare, x1, x2, n, metric):
    if stats_out[cre][level_compare]['n']>5:
        if stats_out[cre][level_compare][metric]>0.05:            
            label_diff(x1, x2, n, '', color_dict[cre])
        elif stats_out[cre][level_compare][metric]<0.05 and stats_out[cre][level_compare][metric]>0.01 :            
            label_diff(x1, x2, n, '* p=%0.2g' %stats_out[cre][level_compare][metric], color_dict[cre])
        elif stats_out[cre][level_compare][metric]<0.01 :            
            label_diff(x1, x2, n, '** p=%0.2g' %stats_out[cre][level_compare][metric], color_dict[cre])
        else:
            raise Exception('p value doesnt make sense')

def make_sig_plots(key_list, metric):
    plt.xlim(.25, 5)
     
    #make plots
    n=-5    
    for cre in key_list:
        if stats_out[cre]['LIF_LIFASC']['n']>=5:
            n=n+5        

            plot_level_of_sig(stats_out, cre, 'LIF_LIFR', .5, 1.45, n+0, metric=metric)
            plot_level_of_sig(stats_out, cre, 'LIFR_LIFASC', 1.55, 2.45, n+0, metric=metric)     
            plot_level_of_sig(stats_out, cre, 'LIFASC_LIFRASC', 2.55, 3.45, n+0, metric=metric)
            plot_level_of_sig(stats_out, cre, 'LIFRASC_LIFRASCAT', 3.55, 4.45, n+0, metric=metric) 
                       
            plot_level_of_sig(stats_out, cre, 'LIF_LIFASC', .5, 2.45, n+1, metric=metric)
            plot_level_of_sig(stats_out, cre, 'LIFASC_LIFRASCAT', 2.55, 4.5, n+1, metric=metric)  
    
            plot_level_of_sig(stats_out, cre, 'LIF_LIFRASC', .5, 3.5, n+2, metric=metric)         
    
            plot_level_of_sig(stats_out, cre, 'LIFR_LIFRASCAT', 1.5, 4.5, n+3, metric=metric)                
    
            plot_level_of_sig(stats_out, cre, 'LIF_LIFRASCAT', .5, 4.5, n+4, metric=metric) 
            

    plt.ylim(-0.7, n+4)
    plt.axis('off')

neuron_data, cre_data_dict, stats_out=pickle.load(open("saved_data/stats_out.pkl", "rb"))   

#-------box plots for total, excitatory and inhibitory
metric='Benjamini_Hochberg_p'
plt.figure()
plt.subplot(2,1,2)
make_box_plot(neuron_data, 'k')
make_box_plot(cre_data_dict['inhibitory'], color_dict['inhibitory'], p=np.array([1.,2.,3.,4.,5.]) + .1)
make_box_plot(cre_data_dict['excitatory'], color_dict['excitatory'], p=np.array([1.,2.,3.,4.,5.]) - .1)
plt.xlim(.7, 5.5)
plt.xticks([1, 2,3,4,5], [r'$GLIF_1$', r'$GLIF_2$', '$GLIF_3$', r'$GLIF_4$', r'$GLIF_5$'])
plt.tick_params(axis='x', length=0)
plt.ylabel('% Explained Variance')
plt.subplot(2,1,1)
make_sig_plots(['all_neurons', 'inhibitory', 'excitatory'], metric=metric)
plt.tight_layout(h_pad=-1)
plt.annotate('All neurons', xy=(.15, .6), fontsize=20, color='k', xycoords='figure fraction',
                         horizontalalignment='right', verticalalignment='bottom')
plt.annotate('Inhibitory', xy=(.15, .75), fontsize=20, color='r', xycoords='figure fraction',
                         horizontalalignment='right', verticalalignment='bottom')
plt.annotate('Excitatory', xy=(.15, .9), fontsize=20, color='b', xycoords='figure fraction',
                         horizontalalignment='right', verticalalignment='bottom')


#-----box plots for individual excitatory cre_lines
plt.figure()
plt.subplot(2,1,2)
plt.subplot2grid((5,1),(3, 0), rowspan=2)
ii=0
for cre in excitatory_cre_lines:
    if len(cre_data_dict[cre])>=5:
        ii=ii+1
        p=np.array([1.,2.,3.,4.,5.]) + .1*ii
        make_box_plot(cre_data_dict[cre], color_dict[cre], p=p)
plt.xlim(1,6)
plt.xticks([1.4, 2.4,3.4,4.4,5.4], [r'$GLIF_1$', r'$GLIF_2$', '$GLIF_3$', r'$GLIF_4$', r'$GLIF_5$'])
plt.tick_params(axis='x', length=0)
plt.ylabel('% Explained Variance')
plt.subplot2grid((5,1),(0, 0), rowspan=3)
make_sig_plots(excitatory_cre_lines, metric=metric)
plt.tight_layout(h_pad=-2)
#get total number of sigline sections that will be plotted
tp=0
for cre in excitatory_cre_lines:
    if stats_out[cre]['LIF_LIFASC']['n']>5:
        tp=tp+1 
ii=0            
for cre in excitatory_cre_lines:
    if stats_out[cre]['LIF_LIFASC']['n']>5:
        ii=ii+1
        print cre, cre.replace('POS', '')
        plt.annotate(cre.replace('POS', ''), xy=(.15, .4+ii*.53/tp), fontsize=20, color=color_dict[cre], xycoords='figure fraction',
                         horizontalalignment='right', verticalalignment='bottom')


#---------box plots for inhibitory cre_lines
plt.figure()
plt.subplot(2,1,2)
ii=0
for cre in inhibitory_cre_lines:
    if len(cre_data_dict[cre])>=5:
        ii=ii+1
        p=np.array([1.,2.,3.,4.,5.]) + .1*ii
        make_box_plot(cre_data_dict[cre], color_dict[cre], p=p)
plt.xticks([1.3, 2.3,3.3,4.3,5.3], [r'$GLIF_1$', r'$GLIF_2$', '$GLIF_3$', r'$GLIF_4$', r'$GLIF_5$'])
plt.xlim(1,6)
plt.tick_params(axis='x', length=0)
plt.ylabel('% Explained Variance')
plt.subplot(2,1,1)
make_sig_plots(inhibitory_cre_lines, metric=metric)
plt.tight_layout(h_pad=-1)
#get total number of sigline sections that will be plotted
tp=0
for cre in inhibitory_cre_lines:
    if stats_out[cre]['LIF_LIFASC']['n']>5:
        tp=tp+1 
ii=0            
for cre in inhibitory_cre_lines:
    if stats_out[cre]['LIF_LIFASC']['n']>5:
        ii=ii+1
        plt.annotate(cre.replace('POS', ''), xy=(.16, .47+ii*.47/tp), fontsize=20, color=color_dict[cre], xycoords='figure fraction',
                         horizontalalignment='right', verticalalignment='bottom')

#---create a plot of the neurons progression.

plt.figure()
for cre_key in excitatory_cre_lines+inhibitory_cre_lines:
    median=[stats_out[cre_key]['LIF']['median'],
            stats_out[cre_key]['LIFR']['median'],
            stats_out[cre_key]['LIFASC']['median'],
            stats_out[cre_key]['LIFRASC']['median'],
            stats_out[cre_key]['LIFRASCAT']['median']]
    plt.plot([1,2,3,4,5], median, color=color_dict[cre_key], lw=2, label=cre_key)
for cre_key in ['all_neurons', 'excitatory', 'inhibitory']:
    median=[stats_out[cre_key]['LIF']['median'],
            stats_out[cre_key]['LIFR']['median'],
            stats_out[cre_key]['LIFASC']['median'],
            stats_out[cre_key]['LIFRASC']['median'],
            stats_out[cre_key]['LIFRASCAT']['median']]
    plt.plot([1,2,3,4,5], median, color=color_dict[cre_key], lw=6, label=cre_key)

plt.xlim(1,5)
plt.xticks([1,2,3,4,5])
plt.ylabel('% Median Explained Variance')
plt.xlabel('GLIF model level')


#--create table for manuscript-----
rows=[['all', 'all_neurons'],
          ['inhibitory', 'inhibitory'],
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

columns=[['n', 'n'], 
         ['$GLIF_1$', 'LIF'],
         ['$GLIF_2$', 'LIFR'],
         ['$GLIF_3$', 'LIFASC'],
         ['$GLIF_4$', 'LIFRASC'], 
         ['$GLIF_5$','LIFRASCAT']]

df1=pd.DataFrame(index=[ind[0] for ind in rows], columns=[ind[0] for ind in columns])
for index in rows:
    if stats_out[index[1]]['LIF']['n']<5:
        df1=df1.drop(index[0])
        print [index[1]], ' n:', stats_out[index[1]]['LIF']['n']
        pass
    else:
        for ii, col in enumerate(columns):
            if ii ==0:
                value=stats_out[index[1]]['LIF']['n']
            else:
                if stats_out[index[1]][col[1]]['n']<5:
                    value =np.nan
                else:
                    value ="\shortstack{%0.3g" % stats_out[index[1]][col[1]]['median']+"\\\\{[}"+"%0.3g" %stats_out[index[1]][col[1]]['quartiles'][0]+", %0.3g" %stats_out[index[1]][col[1]]['quartiles'][1] + '{]}}'
            df1.set_value(index[0], col[0], value)
print df1.to_latex(escape=False)

plt.show()