'''Creates figures related to the calculation of resistance, 
capacitance and total charge deposited by the after-spike currents.
'''

import numpy as np
import matplotlib.pyplot as plt
import allensdk.core.json_utilities as ju
import sys
import os
relative_path=os.path.dirname(os.getcwd())
sys.path.append(os.path.join(relative_path, 'libraries'))
from data_access import load_sweeps 
from data_library import check_and_organize_data, get_file_path_endswith, get_sweep_num_by_name, get_pp_path
from pub_plot_library import distribution_plot, pos_cre_lines, color_dict, draw_unity
from allensdk.core.cell_types_cache import CellTypesCache
ctc = CellTypesCache(manifest_file=os.path.join(relative_path,'cell_types_manifest.json'))

#---------------------------------------------------------------
#------------NO SPECIFICATIONS NEEDED---------------------------
#---------------------------------------------------------------

# load data out of configuration files
data_path=os.path.join(relative_path,'mouse_struc_data_dir')
folders=[os.path.join(data_path, f) for f in  os.listdir(data_path)]

all_neurons=[]   
for folder in folders:
    specimen_ID=os.path.basename(folder)[:9]
    cre=os.path.basename(folder)[10:]
    pp_file=get_pp_path(folder)
    pp_dict=ju.read(pp_file)
    all_neurons.append([specimen_ID,                                                 #0
              cre,                                                                   #1
              pp_dict['resistance']['R_test_list']['mean']/1.e6,
              pp_dict['resistance']['R_fit_ASC_and_R']['mean']/1.e6,
              pp_dict['capacitance']['C_test_list']['mean']*1.e12,              #4
              (1./pp_dict['asc']['k'][0])*pp_dict['asc']['amp'][0]*1.e12,   #total charge
              (1./pp_dict['asc']['k'][1])*pp_dict['asc']['amp'][1]*1.e12,   #total charge
              pp_dict['resistance']['R_from_lims']['value'],                    #7
              pp_dict['capacitance']['C_from_lims']['value'],                   #8           
              pp_dict['capacitance']['C_test_list']['mean']*pp_dict['resistance']['R_test_list']['mean']*1.e3])                   

(cre_dict)=check_and_organize_data(all_neurons)
            
#---------------------------------------------------
#---plotting simple cap versus resistance no asc----
#---------------------------------------------------

percentile_dict=distribution_plot(cre_dict, 2, 4, xlabel='R (MOhms)', ylabel='C (pF)')
plt.annotate('Resistance Fit Without ASC', xy=(.5, .93),
                     xycoords='figure fraction',
                     horizontalalignment='center', verticalalignment='bottom',
                     fontsize=22)

#------------------------------------------------------
#---plotting simple cap versus tau with no ASC---------
#------------------------------------------------------
percentile_dict=distribution_plot(cre_dict, 9, 4, xlabel=r'$\tau$ (ms)', ylabel='C (pF)')
plt.annotate('Tau Fit Without ASC', xy=(.5, .93),
                     xycoords='figure fraction',
                     horizontalalignment='center', verticalalignment='bottom',
                     fontsize=22)

#------------------------------------------------------
#---plotting simple cap versus resistance WITH ASC-----
#------------------------------------------------------
percentile_dict=distribution_plot(cre_dict, 3, 4, xlabel=r'$R_{ASC}$ (MOhms)', ylabel='C (pF)')
plt.annotate('Resistance Fit With ASC', xy=(.5, .93),
                     xycoords='figure fraction',
                     horizontalalignment='center', verticalalignment='bottom',
                     fontsize=22)
#---------------------------------------------------
#---plotting R with asc versus no asc---------------
#---------------------------------------------------
#--setting up figure
plt.figure()
#--plotting main scatter plot figure
for cre in pos_cre_lines: 
    for neuron in cre_dict[cre]:
        plt.plot(neuron[2], neuron[3], '.',color=color_dict[cre], ms=12)
plt.ylabel(r'$R_{ASC}$ (MOhms)')
plt.xlabel('R (MOhms)')

draw_unity()

#---------------------------------------------------
#---plotting total charge---------------------------
#---------------------------------------------------

percentile_dict=distribution_plot(cre_dict, 5, 6, xlabel=r'$Q_1$ (pC)', ylabel=r'$Q_2$ (pC)')


def simple_model_to_eval_subthresh_params(stim_list, R_list, C_list, El_list, dt):    
    '''creates a subthreshold voltage trace set for a simple neuron model
    inputs:  Note that although these inputs and output are all in lists in this code
            in practice they are only single floats and arrays floats.
        stim_list: list of arrays
            array of floats contains the waveform of the current injected into model/neuron
        R_list: list of floats
            individual float is resistance of neuron
        C_list: list of floats
            individual float is capacitance of neuron
        El_list: list of floats
            individual float is reversal potential of neuron
        dt: float
            size of time step
    returns:
        voltage_list: list of arrays
            each array is the voltage output corresponding to an individual stim_list array
            
    '''
    voltage_list=[]
    for stim, R, C, El in zip(stim_list, R_list, C_list, El_list): 
        v=El
        v_single_trace=[v]
        for inj in stim:
            v=v + (inj + - 1./R * (v - El)) * dt / C
            v_single_trace.append(v)
        voltage_list.append(np.array(v_single_trace[:-1]))
        
    return voltage_list   

#------------------------------------------------------------
#--open a file from the preprocessor and plot single traces-- 
#------------------------------------------------------------

#set up figure
plt.figure(figsize=(14, 6))
I1_plt=plt.subplot2grid((7,1), (0, 0))
V2_plt=plt.subplot2grid((7,1), (1, 0), rowspan=6)

specimen_id='474637203' #htr3
sub_folder=os.path.join(data_path,os.listdir(data_path)[np.where([specimen_id in fname for fname in os.listdir(data_path)])[0][0]])
file=get_file_path_endswith(sub_folder, '_preprocessor_values.json')
neuron_dict=ju.read(file)
R_NO_asc=neuron_dict['resistance']['R_test_list']['mean']
R_asc=neuron_dict['resistance']['R_fit_ASC_and_R']['mean']
C=neuron_dict['capacitance']['C_test_list']['mean']
El=neuron_dict['El']['El_noise']['measured']['mean']

# get the sweeps
dir_name=os.path.join(relative_path, 'mouse_nwb/specimen_'+ specimen_id)
the_sweeps=ctc.get_ephys_sweeps(int(specimen_id),  os.path.join(dir_name, 'ephys_sweeps.json'))
noise1_sweeps=get_sweep_num_by_name(the_sweeps, 'Noise 1')

# put data in the format required for functions below
n1_s1_data=ctc.get_ephys_data(int(specimen_id), os.path.join(dir_name, 'ephys.nwb')).get_sweep(noise1_sweeps[0])
sr=n1_s1_data['sampling_rate']
dt=1./sr
suthresh_i=n1_s1_data['stimulus'][n1_s1_data['index_range'][0]+int(1.2/dt):int(3./dt)]
suthresh_v=[]
for s in noise1_sweeps:
    data=ctc.get_ephys_data(int(specimen_id), os.path.join(dir_name, 'ephys.nwb')).get_sweep(s)
    suthresh_v.append(data['response'][data['index_range'][0]+int(1.2/dt):int(3./dt)])
V_NO_asc=simple_model_to_eval_subthresh_params([suthresh_i], [R_NO_asc], [C], [El], dt)[0]
#V_asc=simple_model_to_eval_subthresh_params([suthresh_i], [R_asc], [C], [El], dt)[0]

time =np.arange(len(V_NO_asc))*dt*1000
I1_plt.plot(time, suthresh_i*1.e12, 'k')
I1_plt.set_xlim((0,1000))
I1_plt.axis('off')
#note shifting voltage to be on same axis
for v in suthresh_v:
    V2_plt.plot(time, v*1000 +20, 'b', lw=2)
V2_plt.plot(time, V_NO_asc*1000+20, 'r', lw=3)

#---------------------neuron 2----------------

specimen_id='512322162'#ctgf
sub_folder=os.path.join(data_path,os.listdir(data_path)[np.where([specimen_id in fname for fname in os.listdir(data_path)])[0][0]])
file=get_file_path_endswith(sub_folder, '_preprocessor_values.json')
neuron_dict=ju.read(file)
R_NO_asc=neuron_dict['resistance']['R_test_list']['mean']
R_asc=neuron_dict['resistance']['R_fit_ASC_and_R']['mean']
C=neuron_dict['capacitance']['C_test_list']['mean']
El=neuron_dict['El']['El_noise']['measured']['mean']

# get the sweeps
dir_name=os.path.join(relative_path, 'mouse_nwb/specimen_'+ specimen_id)
the_sweeps=ctc.get_ephys_sweeps(int(specimen_id), os.path.join(dir_name, 'ephys_sweeps.nwb'))
noise1_sweeps=get_sweep_num_by_name(the_sweeps, 'Noise 1')

# put data in the format required for functions below
n1_s1_data=ctc.get_ephys_data(int(specimen_id), os.path.join(dir_name, 'ephys.nwb')).get_sweep(noise1_sweeps[0])
sr=n1_s1_data['sampling_rate']
dt=1./sr
suthresh_i=n1_s1_data['stimulus'][n1_s1_data['index_range'][0]+int(1.2/dt):int(3./dt)]
suthresh_v=[]
for s in noise1_sweeps:
    data=ctc.get_ephys_data(int(specimen_id), os.path.join(dir_name, 'ephys.nwb')).get_sweep(s)
    suthresh_v.append(data['response'][data['index_range'][0]+int(1.2/dt):int(3./dt)])
V_NO_asc=simple_model_to_eval_subthresh_params([suthresh_i], [R_NO_asc], [C], [El], dt)[0]
time =np.arange(len(V_NO_asc))*dt*1000

for ii, v in enumerate(suthresh_v):
    if ii==0:
        V2_plt.plot(time, v*1000, 'b', lw=2, label='data')
    else:
        V2_plt.plot(time, v*1000, 'b', lw=2)
V2_plt.plot(time, V_NO_asc*1000, 'r', lw=3, label='model')
V2_plt.set_xlabel('time (ms)')

#add scale bar
V2_plt.plot([800, 800],[-40, -50], 'k', lw=6)
V2_plt.plot([800, 10000],[-50, -50], 'k', lw=6)
V2_plt.annotate('200 ms', xy=(700, -46), fontsize=16)
V2_plt.annotate('10 mV', xy=(850, -56), fontsize=16)

# add legend bar
plt.legend(bbox_to_anchor=(1., .83),frameon=False)
V2_plt.set_xlim((0,1000))
V2_plt.axis('off')
plt.show()


