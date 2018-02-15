'''Written by Corinne Teeter. Grabs subthreshold voltage values and explained variance from
the structured data folder and makes plots.
'''

import pickle
import os
import numpy as np
import allensdk.core.json_utilities as ju
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 16})
import sys
relative_path=os.path.dirname(os.getcwd())
sys.path.append(os.path.join(relative_path, 'libraries'))
from data_library import get_file_path_endswith, convert_list_to_numpy, get_ev_percent_from_calculated_file
from scipy import stats
from decimal import Decimal

def plotLineRegress(stats_out,xlim, color, level_string):
    '''stats_out: out put object from scipy.stats.linregress
    xlim: must be numpy array
    '''
    slope=stats_out[0]
    intercept=stats_out[1]
    r=stats_out[2]
    p=stats_out[3]
    y=slope*xlim+intercept
    print stats_out
    plt.plot(xlim, y, '-k', color=color, lw=4, label=level_string + ': slope=%.2E, int=%.2E, r=%.2f, p=%.2E' %(Decimal(slope), Decimal(intercept), r, Decimal(p)))

def grab_diff_v_from_folder(ew, folder):
    '''get the voltage difference if the file exists
    inputs:
        ew: string that is to be matched with a file in the folder in order to return a 
            value. i.e. if GLIF2 is being requested but there is not GLIF2 file in the
            folder a nan will be returned regardless of whether there is a value of explain 
            variance in the database. For example this would happen if the model was excluded from 
            analysis because of an aberrant parameter.  
        folder: path to the structured folder used in the rest of analysis
    returns:
        either or nan or the explained variance ratio for the requested model
        
    '''
    try:
        file=get_file_path_endswith(folder, ew)
        contents=ju.read(file)
        RSS_of_voltage_diff=contents['noise2']['RSS_of_voltage_diff']
        var_of_voltage_data=contents['noise2']['var_of_voltage_data']
        num_data_points_wo_spike_shape=contents['noise2']['num_data_points_wo_spike_shape']
    except:
        RSS_of_voltage_diff=np.nan
        var_of_voltage_data=np.nan
        num_data_points_wo_spike_shape=np.nan
        
    return RSS_of_voltage_diff, var_of_voltage_data, num_data_points_wo_spike_shape

#def grab_ev_from_folder(ew, folder):
#    try:
#        file=get_file_path_endswith(folder, ew)
#        contents=ju.read(file)
#        ev=contents['after_opt']['noise_2']
#    except:
#        ev=np.nan
#
#    return ev

def test_and_remove_nans(ev_list, diff_v_list):
    '''If a index is nan in ev than it should be nan in diff for a given level 
    '''
    short_ev=[]
    short_diff_v=[]
    for ev, diff in zip(ev_list, diff_v_list):
        if np.isnan(ev) and np.isnan(diff):
            continue
        if np.isnan(ev) ^ np.isnan(diff):
            raise Exception('both ev and v_diff should be nan and they are not')
        else:
            short_ev.append(ev)
            short_diff_v.append(diff)
    
    return np.array(short_ev), np.array(short_diff_v)   

def set_up_data(folders):
    
    GLIF1_dict={'ev':[], 'RSS_of_voltage_diff':[], 'var_of_voltage_data':[], 'num_data_points_wo_spike_shape':[]}
    GLIF2_dict={'ev':[], 'RSS_of_voltage_diff':[], 'var_of_voltage_data':[], 'num_data_points_wo_spike_shape':[]}
    GLIF3_dict={'ev':[], 'RSS_of_voltage_diff':[], 'var_of_voltage_data':[], 'num_data_points_wo_spike_shape':[]}
    GLIF4_dict={'ev':[], 'RSS_of_voltage_diff':[], 'var_of_voltage_data':[], 'num_data_points_wo_spike_shape':[]}
    GLIF5_dict={'ev':[], 'RSS_of_voltage_diff':[], 'var_of_voltage_data':[], 'num_data_points_wo_spike_shape':[]}

    for folder in folders:
        name=os.path.basename(folder)
        specimen_id=name[:9]
        cre=os.path.basename(folder)[10:]
        
        #grabbing explained variance
#        GLIF1_dict['ev'].append(grab_ev_from_folder('_GLIF1_exp_var_ratio_10ms.json', folder)*100.)   
#        GLIF2_dict['ev'].append(grab_ev_from_folder('_GLIF2_exp_var_ratio_10ms.json', folder)*100.)
#        GLIF3_dict['ev'].append(grab_ev_from_folder('_GLIF3_exp_var_ratio_10ms.json', folder)*100.)
#        GLIF4_dict['ev'].append(grab_ev_from_folder('_GLIF4_exp_var_ratio_10ms.json', folder)*100.)
#        GLIF5_dict['ev'].append(grab_ev_from_folder('_GLIF5_exp_var_ratio_10ms.json', folder)*100.)
        GLIF1_dict['ev'].append(get_ev_percent_from_calculated_file('_GLIF1_exp_var_ratio_10ms.json', folder, 'after_opt', 'noise_2'))
        GLIF2_dict['ev'].append(get_ev_percent_from_calculated_file('_GLIF2_exp_var_ratio_10ms.json', folder, 'after_opt', 'noise_2'))
        GLIF3_dict['ev'].append(get_ev_percent_from_calculated_file('_GLIF3_exp_var_ratio_10ms.json', folder, 'after_opt', 'noise_2'))
        GLIF4_dict['ev'].append(get_ev_percent_from_calculated_file('_GLIF4_exp_var_ratio_10ms.json', folder, 'after_opt', 'noise_2'))
        GLIF5_dict['ev'].append(get_ev_percent_from_calculated_file('_GLIF5_exp_var_ratio_10ms.json', folder, 'after_opt', 'noise_2'))


        #grabbing the subthreshold voltage values 
        RSS_of_model_data_v_diff_LIF, var_of_data_v_LIF, n_LIF=grab_diff_v_from_folder('_GLIF1_subthr_v.json', folder)
        RSS_of_model_data_v_diff_LIFR, var_of_data_v_LIFR, n_LIFR=grab_diff_v_from_folder('_GLIF2_subthr_v.json', folder)
        RSS_of_model_data_v_diff_LIFASC, var_of_data_v_LIFASC, n_LIFASC=grab_diff_v_from_folder('_GLIF3_subthr_v.json', folder)
        RSS_of_model_data_v_diff_LIFRASC, var_of_data_v_LIFRASC, n_LIFRASC=grab_diff_v_from_folder('_GLIF4_subthr_v.json', folder)
        RSS_of_model_data_v_diff_LIFRASCAT, var_of_data_v_LIFRASCAT, n_LIFRASCAT=grab_diff_v_from_folder('_GLIF5_subthr_v.json', folder)

        GLIF1_dict['RSS_of_voltage_diff'].append(RSS_of_model_data_v_diff_LIF)
        GLIF2_dict['RSS_of_voltage_diff'].append(RSS_of_model_data_v_diff_LIFR)
        GLIF3_dict['RSS_of_voltage_diff'].append(RSS_of_model_data_v_diff_LIFASC)
        GLIF4_dict['RSS_of_voltage_diff'].append(RSS_of_model_data_v_diff_LIFRASC)
        GLIF5_dict['RSS_of_voltage_diff'].append(RSS_of_model_data_v_diff_LIFRASCAT)
        
        
        GLIF1_dict['var_of_voltage_data'].append(var_of_data_v_LIF)
        GLIF2_dict['var_of_voltage_data'].append(var_of_data_v_LIFR)
        GLIF3_dict['var_of_voltage_data'].append(var_of_data_v_LIFASC)
        GLIF4_dict['var_of_voltage_data'].append(var_of_data_v_LIFRASC)
        GLIF5_dict['var_of_voltage_data'].append(var_of_data_v_LIFRASCAT)
        
        GLIF1_dict['num_data_points_wo_spike_shape'].append(n_LIF)
        GLIF2_dict['num_data_points_wo_spike_shape'].append(n_LIFR)
        GLIF3_dict['num_data_points_wo_spike_shape'].append(n_LIFASC)
        GLIF4_dict['num_data_points_wo_spike_shape'].append(n_LIFRASC)
        GLIF5_dict['num_data_points_wo_spike_shape'].append(n_LIFRASCAT)
        
        
    #convert lists to numpy arrays for potential ease in future computations    
    GLIF1_dict=convert_list_to_numpy(GLIF1_dict)
    GLIF2_dict=convert_list_to_numpy(GLIF2_dict)
    GLIF3_dict=convert_list_to_numpy(GLIF3_dict)
    GLIF4_dict=convert_list_to_numpy(GLIF4_dict)
    GLIF5_dict=convert_list_to_numpy(GLIF5_dict)
    
        
    return GLIF1_dict, GLIF2_dict, GLIF3_dict, GLIF4_dict, GLIF5_dict
    


data_path=os.path.join(relative_path,'mouse_struc_data_dir')
folders=[os.path.join(data_path, f) for f in  os.listdir(data_path)]

GLIF1_dict, GLIF2_dict, GLIF3_dict, GLIF4_dict, GLIF5_dict=set_up_data(folders)  #creates several GLIF model level dictionaries with data in lists

#plt.plot(GLIF1_dict['ev'], GLIF1_dict['var_btw_model_data_v_diff'], '.', MS=16, label='GLIF1')
#plt.plot(GLIF2_dict['ev'], GLIF2_dict['var_btw_model_data_v_diff'], '.', MS=16, label='GLIF2')
#plt.plot(GLIF3_dict['ev'], GLIF3_dict['var_btw_model_data_v_diff'], '.', MS=16, label='GLIF3')
#plt.plot(GLIF4_dict['ev'], GLIF4_dict['var_btw_model_data_v_diff'], '.', MS=16, label='GLIF4')
#plt.plot(GLIF5_dict['ev'], GLIF5_dict['var_btw_model_data_v_diff'], '.', MS=16, label='GLIF5')
#plt.legend(loc=2)
#plt.xlabel('spike time performance (exp var)')
#plt.ylabel('((vdata-vmodel)-mean(vdata-vmodel))^2/n')

# normalize data
#--before redo--delete when confirmed
#GLIF1_dict['vdiff_normalized']=GLIF1_dict['var_btw_model_data_v_diff']/GLIF1_dict['var_of_v_data']
#GLIF2_dict['vdiff_normalized']=GLIF2_dict['var_btw_model_data_v_diff']/GLIF2_dict['var_of_v_data']
#GLIF3_dict['vdiff_normalized']=GLIF3_dict['var_btw_model_data_v_diff']/GLIF3_dict['var_of_v_data']
#GLIF4_dict['vdiff_normalized']=GLIF4_dict['var_btw_model_data_v_diff']/GLIF4_dict['var_of_v_data']
#GLIF5_dict['vdiff_normalized']=GLIF5_dict['var_btw_model_data_v_diff']/GLIF5_dict['var_of_v_data']
#--------

GLIF1_dict['vdiff_normalized']=GLIF1_dict['RSS_of_voltage_diff']/(GLIF1_dict['num_data_points_wo_spike_shape']*GLIF1_dict['var_of_voltage_data'])
GLIF2_dict['vdiff_normalized']=GLIF2_dict['RSS_of_voltage_diff']/(GLIF2_dict['num_data_points_wo_spike_shape']*GLIF2_dict['var_of_voltage_data'])
GLIF3_dict['vdiff_normalized']=GLIF3_dict['RSS_of_voltage_diff']/(GLIF3_dict['num_data_points_wo_spike_shape']*GLIF3_dict['var_of_voltage_data'])
GLIF4_dict['vdiff_normalized']=GLIF4_dict['RSS_of_voltage_diff']/(GLIF4_dict['num_data_points_wo_spike_shape']*GLIF4_dict['var_of_voltage_data'])
GLIF5_dict['vdiff_normalized']=GLIF5_dict['RSS_of_voltage_diff']/(GLIF5_dict['num_data_points_wo_spike_shape']*GLIF5_dict['var_of_voltage_data'])

GLIF2_dict['ev'], GLIF2_dict['vdiff_normalized']=test_and_remove_nans(GLIF2_dict['ev'], GLIF2_dict['vdiff_normalized'])
GLIF4_dict['ev'], GLIF4_dict['vdiff_normalized']=test_and_remove_nans(GLIF4_dict['ev'], GLIF4_dict['vdiff_normalized'])
GLIF5_dict['ev'], GLIF5_dict['vdiff_normalized']=test_and_remove_nans(GLIF5_dict['ev'], GLIF5_dict['vdiff_normalized'])

## save the data
#dictionary={'GLIF1':GLIF1_dict, 'GLIF2':GLIF2_dict, 'GLIF3':GLIF3_dict, 'GLIF4':GLIF4_dict, 'GLIF5':GLIF5_dict}
#output = open('subthreshold_data.pkl', 'w')
#pickle.dump(dictionary, output)
#output.close()

#plot data
plt.figure()
plt.plot(GLIF1_dict['ev'], GLIF1_dict['vdiff_normalized'], 'b.', MS=16)
plt.plot(GLIF2_dict['ev'], GLIF2_dict['vdiff_normalized'], 'g.', MS=16)
plt.plot(GLIF3_dict['ev'], GLIF3_dict['vdiff_normalized'], 'r.', MS=16)
plt.plot(GLIF4_dict['ev'], GLIF4_dict['vdiff_normalized'], 'c.', MS=16)
plt.plot(GLIF5_dict['ev'], GLIF5_dict['vdiff_normalized'], 'm.', MS=16)
plt.legend(loc=2)

#plot medians
plt.plot(np.median(GLIF1_dict['ev']), np.median(GLIF1_dict['vdiff_normalized']), 'bs', MS=20, MEC='k',MEW=5)
plt.plot(np.median(GLIF2_dict['ev']), np.median(GLIF2_dict['vdiff_normalized']), 'gs', MS=20, MEC='k',MEW=5)
plt.plot(np.median(GLIF3_dict['ev']), np.median(GLIF3_dict['vdiff_normalized']), 'rs', MS=20, MEC='k',MEW=5)
plt.plot(np.median(GLIF4_dict['ev']), np.median(GLIF4_dict['vdiff_normalized']), 'cs', MS=20, MEC='k',MEW=5)
plt.plot(np.median(GLIF5_dict['ev']), np.median(GLIF5_dict['vdiff_normalized']), 'ms', MS=20, MEC='k',MEW=5)

x_lim=np.array(plt.gca().get_xlim())

#calc linear regression
GLIF1_regress=stats.linregress(GLIF1_dict['ev'], GLIF1_dict['vdiff_normalized'])
GLIF2_regress=stats.linregress(GLIF2_dict['ev'], GLIF2_dict['vdiff_normalized'])
GLIF3_regress=stats.linregress(GLIF3_dict['ev'], GLIF3_dict['vdiff_normalized'])
GLIF4_regress=stats.linregress(GLIF4_dict['ev'], GLIF4_dict['vdiff_normalized'])
GLIF5_regress=stats.linregress(GLIF5_dict['ev'], GLIF5_dict['vdiff_normalized'])
all_regress=stats.linregress(np.concatenate((GLIF1_dict['ev'], GLIF2_dict['ev'], GLIF3_dict['ev'], GLIF4_dict['ev'], GLIF5_dict['ev'])), 
                             np.concatenate((GLIF1_dict['vdiff_normalized'], GLIF2_dict['vdiff_normalized'], GLIF3_dict['vdiff_normalized'],GLIF4_dict['vdiff_normalized'], GLIF5_dict['vdiff_normalized'])))

plotLineRegress(GLIF1_regress, x_lim, 'b', 'GLIF1')
plotLineRegress(GLIF2_regress, x_lim, 'g', 'GLIF2')
plotLineRegress(GLIF3_regress, x_lim, 'r', 'GLIF3')
plotLineRegress(GLIF4_regress, x_lim, 'c', 'GLIF4')
plotLineRegress(GLIF5_regress, x_lim, 'm', 'GLIF5')
plotLineRegress(all_regress, x_lim, 'k', 'All')
plt.legend()
plt.xlim(x_lim)
plt.xlabel('% Explained Variance')
plt.ylabel('Suthreshold Voltage Difference Measure')

#plot means on new plot
ev_medians=[np.median(GLIF1_dict['ev']), 
          np.median(GLIF2_dict['ev']), 
          np.median(GLIF3_dict['ev']), 
          np.median(GLIF4_dict['ev']), 
          np.median(GLIF5_dict['ev'])]
v_diff_medians=[np.median(GLIF1_dict['vdiff_normalized']),
             np.median(GLIF2_dict['vdiff_normalized']),
             np.median(GLIF3_dict['vdiff_normalized']),
             np.median(GLIF4_dict['vdiff_normalized']),
             np.median(GLIF5_dict['vdiff_normalized'])]
print 'v_diff_medians', v_diff_medians
plt.figure()
plt.plot(ev_medians, v_diff_medians, 'k', linestyle='-', marker='.', LW=2, MS=16)
plt.xlabel('Median % Explained Variance')
plt.ylabel('Suthreshold Voltage Difference Measure')
plt.annotate(r'GLIF$_1$', xy=(70.7, .17), horizontalalignment='left',fontsize=18)
plt.annotate(r'GLIF$_2$', xy=(67.2, .105), horizontalalignment='left',fontsize=18)
plt.annotate(r'GLIF$_3$', xy=(73, .227), horizontalalignment='left',fontsize=18)
plt.annotate(r'GLIF$_4$', xy=(74, .085), horizontalalignment='left',fontsize=18)
plt.annotate(r'GLIF$_5$', xy=(76.6, .095), horizontalalignment='left',fontsize=18)

plt.show()