'''Plots spike rasters of cononical neurons along with their models
'''

import matplotlib.pyplot as plt
import numpy as np
import os
import allensdk.core.json_utilities as ju
import sys
relative_path=os.path.dirname(os.getcwd())
sys.path.append(os.path.join(relative_path, 'libraries'))
from data_library import get_sweep_num_by_name, get_model_spike_times_from_nwb
from allensdk.core.nwb_data_set import NwbDataSet
from allensdk.core.cell_types_cache import CellTypesCache
ctc = CellTypesCache(manifest_file=os.path.join(relative_path,'cell_types_manifest.json'))

#---------------------------------------------------------------
#------------NO SPECIFICATIONS NEEDED---------------------------
#---------------------------------------------------------------


def get_spikes(specimen_id):
    #---get spike times from the data
    dir_name=os.path.join(relative_path, 'mouse_nwb/specimen_'+ str(specimen_id))
    all_sweeps=ctc.get_ephys_sweeps(specimen_id, os.path.join(dir_name, 'ephys_sweeps.json'))
    sweeps=get_sweep_num_by_name(all_sweeps, 'Noise 2')
    nwb=ctc.get_ephys_data(specimen_id, os.path.join(dir_name, 'ephys.json'))
    data=[]
    bio_spike_times=[]
    for s in sweeps:
        bio_spike_times.append(nwb.get_spike_times(s))
        data.append(nwb.get_sweep(s))    
    dt=1./data[0]['sampling_rate']
    stim=data[0]['stimulus']
    
    #---get model spike times
    #find the folder name for specimen
    for dir in folders:
        sp_id=int(os.path.basename(dir)[:9])
        if sp_id == specimen_id:
            folder=dir
    
    out={}
    out['dt']=dt
    out['bio_spike_times']=bio_spike_times
    out['stimulus']=stim
#    out['LIF_spike_times']=get_model_spikes('_GLIF1_neuron_config.json', folder, '(LIF)', sweeps)
#    out['LIFR_spike_times']=get_model_spikes('_GLIF2_neuron_config.json', folder, '(LIF-R)', sweeps)
#    out['LIFASC_spike_times']=get_model_spikes('_GLIF3_neuron_config.json', folder, '(LIF-ASC)', sweeps)
#    out['LIFRASC_spike_times']=get_model_spikes('_GLIF4_neuron_config.json', folder, '(LIF-R-ASC)', sweeps)
#    out['LIFRASCAT_spike_times']=get_model_spikes('_GLIF5_neuron_config.json', folder, '(LIF-R-ASC-A)', sweeps)
    out['LIF_spike_times']=get_model_spike_times_from_nwb('_GLIF1_neuron_config.json', folder, '(LIF)', sweeps)
    out['LIFR_spike_times']=get_model_spike_times_from_nwb('_GLIF2_neuron_config.json', folder, '(LIF-R)', sweeps)
    out['LIFASC_spike_times']=get_model_spike_times_from_nwb('_GLIF3_neuron_config.json', folder, '(LIF-ASC)', sweeps)
    out['LIFRASC_spike_times']=get_model_spike_times_from_nwb('_GLIF4_neuron_config.json', folder, '(LIF-R-ASC)', sweeps)
    out['LIFRASCAT_spike_times']=get_model_spike_times_from_nwb('_GLIF5_neuron_config.json', folder, '(LIF-R-ASC-A)', sweeps)


    return out

#set up figure
plt.figure(figsize=(16, 6))
columns=9
rows=9
I_low=plt.subplot2grid((rows,columns), (0, 1), colspan=4) 
n1_low=plt.subplot2grid((rows,columns), (1, 1), colspan=4, rowspan=4)
n2_low=plt.subplot2grid((rows,columns), (5, 1), colspan=4, rowspan=4)
I_high=plt.subplot2grid((rows,columns), (0, 5), colspan=4)
n1_high=plt.subplot2grid((rows,columns), (1, 5), colspan=4, rowspan=4)
n2_high=plt.subplot2grid((rows,columns), (5, 5), colspan=4, rowspan=4)
ms=12
mew=2

data_path=os.path.join(relative_path,'mouse_struc_data_dir')
folders=np.sort([os.path.join(data_path, f) for f in  os.listdir(data_path)])

#-----------------------------Htr3a------------------------------------------------

out=get_spikes(474637203)

#--------low amplitude-----------------------------------------------------
time=np.arange(len(out['stimulus']))*out['dt']

I_low.plot(time, out['stimulus'], 'k', lw=2)
I_low.set_xlim(9.9, 13.1)
I_low.axis('off')

for ii, bst in enumerate (out['bio_spike_times']):
    n1_low.plot(bst, np.ones(len(bst))*(9-ii), 'k|', ms=ms, mew=mew)

n1_low.plot(out['LIF_spike_times'][0], np.ones(len(out['LIF_spike_times'][0]))*5, '|', ms=ms, mew=mew)
n1_low.plot(out['LIFR_spike_times'][0], np.ones(len(out['LIFR_spike_times'][0]))*4, '|', ms=ms, mew=mew)
n1_low.plot(out['LIFASC_spike_times'][0], np.ones(len(out['LIFASC_spike_times'][0]))*3, '|', ms=ms, mew=mew)
n1_low.plot(out['LIFRASC_spike_times'][0], np.ones(len(out['LIFRASC_spike_times'][0]))*2, '|', ms=ms, mew=mew)
n1_low.plot(out['LIFRASCAT_spike_times'][0], np.ones(len(out['LIFRASCAT_spike_times'][0]))*1, '|', ms=ms, mew=mew)
n1_low.set_xlim(9.9, 13.1)
n1_low.set_ylim(0,10)
n1_low.axis('off')

#---------------high amplitude---------------------------

I_high.plot(time, out['stimulus'], 'k', lw=2)
I_high.set_xlim(17.9, 21.1)
I_high.axis('off')

for ii, bst in enumerate (out['bio_spike_times']):
    n1_high.plot(bst, np.ones(len(bst))*(9-ii), 'k|', ms=ms, mew=mew)

n1_high.plot(out['LIF_spike_times'][0], np.ones(len(out['LIF_spike_times'][0]))*5, '|', ms=ms, mew=mew)
n1_high.plot(out['LIFR_spike_times'][0], np.ones(len(out['LIFR_spike_times'][0]))*4, '|', ms=ms, mew=mew)
n1_high.plot(out['LIFASC_spike_times'][0], np.ones(len(out['LIFASC_spike_times'][0]))*3, '|', ms=ms, mew=mew)
n1_high.plot(out['LIFRASC_spike_times'][0], np.ones(len(out['LIFRASC_spike_times'][0]))*2, '|', ms=ms, mew=mew)
n1_high.plot(out['LIFRASCAT_spike_times'][0], np.ones(len(out['LIFRASCAT_spike_times'][0]))*1, '|', ms=ms, mew=mew)
n1_high.set_xlim(17.9, 21.1)
n1_high.set_ylim(0,10)
n1_high.axis('off')

##------------------------------Ctgf----------------------------------

out=get_spikes(512322162)

for ii, bst in enumerate (out['bio_spike_times']):
    n2_low.plot(bst, np.ones(len(bst))*(9-ii), 'k|', ms=ms, mew=mew)

n2_low.plot(out['LIF_spike_times'][0], np.ones(len(out['LIF_spike_times'][0]))*5, '|', ms=ms, mew=mew)
n2_low.plot(out['LIFR_spike_times'][0], np.ones(len(out['LIFR_spike_times'][0]))*4, '|', ms=ms, mew=mew)
n2_low.plot(out['LIFASC_spike_times'][0], np.ones(len(out['LIFASC_spike_times'][0]))*3, '|', ms=ms, mew=mew)
n2_low.plot(out['LIFRASC_spike_times'][0], np.ones(len(out['LIFRASC_spike_times'][0]))*2, '|', ms=ms, mew=mew)
n2_low.plot(out['LIFRASCAT_spike_times'][0], np.ones(len(out['LIFRASCAT_spike_times'][0]))*1, '|', ms=ms, mew=mew)
n2_low.set_xlim(9.9, 13.1)
n2_low.set_ylim(0,10)
n2_low.axis('off')

#---------------high amplitude---------------------------

for ii, bst in enumerate (out['bio_spike_times']):
    n2_high.plot(bst, np.ones(len(bst))*(9-ii), 'k|', ms=ms, mew=mew)

n2_high.plot(out['LIF_spike_times'][0], np.ones(len(out['LIF_spike_times'][0]))*5, '|', ms=ms, mew=mew)
n2_high.plot(out['LIFR_spike_times'][0], np.ones(len(out['LIFR_spike_times'][0]))*4, '|', ms=ms, mew=mew)
n2_high.plot(out['LIFASC_spike_times'][0], np.ones(len(out['LIFASC_spike_times'][0]))*3, '|', ms=ms, mew=mew)
n2_high.plot(out['LIFRASC_spike_times'][0], np.ones(len(out['LIFRASC_spike_times'][0]))*2, '|', ms=ms, mew=mew)
n2_high.plot(out['LIFRASCAT_spike_times'][0], np.ones(len(out['LIFRASCAT_spike_times'][0]))*1, '|', ms=ms, mew=mew)
n2_high.set_xlim(17.9, 21.1)
n2_high.set_ylim(0,10)
n2_high.axis('off')


plt.annotate('Medium Amplitude', xy=(.3, .95), xycoords='figure fraction', horizontalalignment='center',fontsize=18)
plt.annotate('High Amplitude', xy=(.75, .95), xycoords='figure fraction', horizontalalignment='center',fontsize=18)
plt.annotate('Htr3a', xy=(0, .7), xycoords='figure fraction', horizontalalignment='left',fontsize=18, rotation='vertical')
plt.annotate('data', xy=(.08, .76), xycoords='figure fraction', horizontalalignment='left',fontsize=18, rotation='vertical')
plt.annotate(r'$[$', xy=(.1, .71), xycoords='figure fraction', horizontalalignment='left',fontsize=60)
plt.annotate('model', xy=(.03, .6), xycoords='figure fraction', horizontalalignment='left',fontsize=18, rotation='vertical')
plt.annotate(r'$[$', xy=(.04, .52), xycoords='figure fraction', horizontalalignment='left',fontsize=80)
plt.annotate(r'GLIF$_1$', xy=(.095, .635), xycoords='figure fraction', horizontalalignment='center',fontsize=16)
plt.annotate(r'GLIF$_2$', xy=(.095, .6), xycoords='figure fraction', horizontalalignment='center',fontsize=16)
plt.annotate(r'GLIF$_3$', xy=(.095, .56), xycoords='figure fraction', horizontalalignment='center',fontsize=16)
plt.annotate(r'GLIF$_4$', xy=(.095, .52), xycoords='figure fraction', horizontalalignment='center',fontsize=16)
plt.annotate(r'GLIF$_5$', xy=(.095, .48), xycoords='figure fraction', horizontalalignment='center',fontsize=16)


plt.annotate('Ctgf', xy=(0, .25), xycoords='figure fraction', horizontalalignment='left',fontsize=18, rotation='vertical')
plt.annotate('data', xy=(.08, .38), xycoords='figure fraction', horizontalalignment='left',fontsize=18, rotation='vertical')
plt.annotate(r'$[$', xy=(.1, .33), xycoords='figure fraction', horizontalalignment='left',fontsize=60)
plt.annotate('model', xy=(.03, .22), xycoords='figure fraction', horizontalalignment='left',fontsize=18, rotation='vertical')
plt.annotate(r'$[$', xy=(.04, .14), xycoords='figure fraction', horizontalalignment='left',fontsize=80)
plt.annotate(r'GLIF$_1$', xy=(.095, .25), xycoords='figure fraction', horizontalalignment='center',fontsize=16)
plt.annotate(r'GLIF$_2$', xy=(.095, .21), xycoords='figure fraction', horizontalalignment='center',fontsize=16)
plt.annotate(r'GLIF$_3$', xy=(.095, .17), xycoords='figure fraction', horizontalalignment='center',fontsize=16)
plt.annotate(r'GLIF$_4$', xy=(.095, .13), xycoords='figure fraction', horizontalalignment='center',fontsize=16)
plt.annotate(r'GLIF$_5$', xy=(.095, .09), xycoords='figure fraction', horizontalalignment='center',fontsize=16)


plt.tight_layout(h_pad=-2)
plt.show() 


