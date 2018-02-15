'''This code will take the data from the spreadsheet that is created by aic_spike_times_noise1.py 
or aic_subthresh_v_noise1.py and plots the differences between them for the different models'''

import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 16})
import numpy as np

def non_nan_data(neuron_data):
    x=[]
    for neuron in neuron_data:
        if not np.isnan(neuron):
            x.append(neuron)

    return x

def make_box_plot(df, the_color, p=[1,2,3,4,5]):
    bp=plt.boxplot([non_nan_data(df['aic_LIF'].values),
                    non_nan_data(df['aic_LIFR'].values),
                    non_nan_data(df['aic_LIFASC'].values),
                    non_nan_data(df['aic_LIFRASC'].values),
                    non_nan_data(df['aic_LIFRASCAT'].values)], notch=0, whis=[5, 95], sym='.', positions=p, widths=.08)
    plt.setp(bp['boxes'], color=the_color, linewidth=3)
    plt.setp(bp['whiskers'], color=the_color, linewidth=3)
    plt.setp(bp['fliers'], mfc=the_color, markersize=12)
    plt.setp(bp['medians'], color=the_color, linewidth=3)
    plt.setp(bp['caps'], color=the_color, linewidth=3)

df_alldata=pd.read_csv('aic_subthresh_v_noise1.csv', delimiter='\t')
make_box_plot(df_alldata, 'k')
plt.title('all data')
plt.xlabel('GLIF level')
plt.ylabel('AIC')
plt.show()

df_data_w_all_levels=df_alldata.dropna(axis=0)
make_box_plot(df_data_w_all_levels, 'k')
plt.title('data with all levels')
plt.xlabel('GLIF level')
plt.ylabel('AIC')
plt.show()

df_diff_w_lif=pd.DataFrame()
df_diff_w_lif['specimen_id']=df_data_w_all_levels['specimen_id']
df_diff_w_lif['cre']=df_data_w_all_levels['cre']
df_diff_w_lif['LIFR']=df_data_w_all_levels['aic_LIF']-df_data_w_all_levels['aic_LIFR']
df_diff_w_lif['LIFASC']=df_data_w_all_levels['aic_LIF']-df_data_w_all_levels['aic_LIFASC']
df_diff_w_lif['LIFRASC']=df_data_w_all_levels['aic_LIF']-df_data_w_all_levels['aic_LIFRASC']
df_diff_w_lif['LIFRASCAT']=df_data_w_all_levels['aic_LIF']-df_data_w_all_levels['aic_LIFRASCAT']

the_color='k'
bp=plt.boxplot([df_diff_w_lif['LIFR'].values,
                df_diff_w_lif['LIFASC'].values,
                df_diff_w_lif['LIFRASC'].values,
                df_diff_w_lif['LIFRASCAT'].values], notch=0, whis=[5, 95], sym='.', positions=[2,3,4,5], labels=[2,3,4,5], widths=.08)
plt.setp(bp['boxes'], color=the_color, linewidth=3)
plt.setp(bp['whiskers'], color=the_color, linewidth=3)
plt.setp(bp['fliers'], mfc=the_color, markersize=12)
plt.setp(bp['medians'], color=the_color, linewidth=3)
plt.setp(bp['caps'], color=the_color, linewidth=3)
plt.ylabel('AIC GLIF1-GLIF*')
plt.xlabel('GLIF level')
plt.show()

