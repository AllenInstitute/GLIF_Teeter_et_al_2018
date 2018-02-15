'''This figure plots the zoomed in mechanisms shown in the manuscript by 
loading the data saved in pickle files made via make_and_save_model_data.py.   
'''
import matplotlib.pyplot as plt
import numpy as np
import pickle
import os
import sys
import sys
relative_path=os.path.dirname(os.getcwd())
sys.path.append(os.path.join(relative_path, 'libraries'))
from data_library import get_sweep_num_by_name
from allensdk.core.cell_types_cache import CellTypesCache
ctc = CellTypesCache(manifest_file=os.path.join(relative_path,'cell_types_manifest.json'))

#---------------------------------------------------------------
#------------NO SPECIFICATIONS NEEDED---------------------------
#---------------------------------------------------------------
#set up figure
plt.figure(figsize=(15,14))
current1=plt.subplot2grid((10,2), (0, 0))
data1=plt.subplot2grid((10,2), (1, 0))
LIF_Htr3a=plt.subplot2grid((10,2), (2, 0))
LIFR_Htr3a=plt.subplot2grid((10,2), (3, 0))
LIFASC_Htr3a=plt.subplot2grid((10,2), (4, 0), rowspan=2)
LIFRASC_Htr3a=plt.subplot2grid((10,2), (6, 0), rowspan=2)
LIFRASCAT_Htr3a=plt.subplot2grid((10,2), (8, 0), rowspan=2)

current2=plt.subplot2grid((10,2), (0, 1))
data2=plt.subplot2grid((10,2), (1, 1))
LIF_Ctgf=plt.subplot2grid((10,2), (2, 1))
LIFR_Ctgf=plt.subplot2grid((10,2), (3, 1))
LIFASC_Ctgf=plt.subplot2grid((10,2), (4, 1), rowspan=2)
LIFRASC_Ctgf=plt.subplot2grid((10,2), (6, 1), rowspan=2)
LIFRASCAT_Ctgf=plt.subplot2grid((10,2), (8, 1), rowspan=2)
x_lim=[18, 18.3]

##---474637203Htr3a------

dir_name=os.path.join(relative_path, 'mouse_nwb/specimen_'+ str(474637203))
the_sweeps=ctc.get_ephys_sweeps(474637203, os.path.join(dir_name, 'ephys_sweeps.json'))
nwb=ctc.get_ephys_data(474637203, os.path.join(dir_name, 'ephys.json'))
sweeps=get_sweep_num_by_name(the_sweeps, 'Noise 2')
data=[]
spike_times=[]
for s in sweeps:
    spike_times.append(nwb.get_spike_times(s))
    data.append(nwb.get_sweep(s)) 

print 'loading LIF'
LIF_model=pickle.load(open( "pkl_data/474637203Htr3a-Cre_NO152_LIF_model.pkl", "rb" ))
print 'loading LIFR'
LIFR_model=pickle.load(open( "pkl_data/474637203Htr3a-Cre_NO152_LIFR_model.pkl", "rb" ))
print 'loading LIFASC'
LIFASC_model=pickle.load(open( "pkl_data/474637203Htr3a-Cre_NO152_LIFASC_model.pkl", "rb" ))
print 'loading LIFRASC'
LIFRASC_model=pickle.load(open( "pkl_data/474637203Htr3a-Cre_NO152_LIFRASC_model.pkl", "rb" ))
print 'loading LIFRASCAT'
LIFRASCAT_model=pickle.load(open( "pkl_data/474637203Htr3a-Cre_NO152_LIFRASCAT_model.pkl", "rb" ))

current1.plot(LIF_model['time'], LIF_model['stimulus'], 'k', lw=4)
current1.set_xlim(x_lim)

data_colors=['orange', 'm', 'lightcoral', 'c']

for kk, st in enumerate(spike_times):
    data1.plot(LIF_model['time'], data[kk]['response'], lw=2, color=data_colors[kk])
    data1.plot(st, np.ones(len(st))*-(.065+.01*kk), '.', ms=10, color=data_colors[kk])
data1.set_xlim(x_lim)
  
ylim_v_jump_small=.065
ylim_v_jump_big=2.*ylim_v_jump_small
ylim_ASC_jump=.9

LIF_Htr3a.plot(LIF_model['time'], LIF_model['voltage']+LIF_model['El_reference'], 'b', lw=4, label='model')
LIF_Htr3a.plot(LIF_model['time'], LIF_model['threshold']+LIF_model['El_reference'], '--g', lw=2, label='threshold')
LIF_Htr3a.plot(LIF_model['grid_spike_times'], np.zeros(len(LIF_model['grid_spike_times']))+LIF_model['El_reference']-.01, '.b', ms=16)
LIF_Htr3a.set_xlim(x_lim)
LIF_Htr3a.set_ylim([-.0815, -.0815+ylim_v_jump_small])

LIFR_Htr3a.plot(LIFR_model['time'], LIFR_model['voltage']+LIFR_model['El_reference'], 'b', lw=4, label='model')
LIFR_Htr3a.plot(LIFR_model['time'], LIFR_model['threshold']+LIFR_model['El_reference'], '--g', lw=2, label='threshold')
LIFR_Htr3a.plot(LIFR_model['grid_spike_times'], np.zeros(len(LIFR_model['grid_spike_times']))+LIFR_model['El_reference'], '.b', ms=16)
LIFR_Htr3a.set_xlim(x_lim)
LIFR_Htr3a.set_ylim([-.075, -.075+ylim_v_jump_small])

LIFASC_Htr3a.plot(LIFASC_model['time'], LIFASC_model['voltage']+LIFASC_model['El_reference'], 'b', lw=4, label='model')
LIFASC_Htr3a.plot(LIFASC_model['time'], LIFASC_model['threshold']+LIFASC_model['El_reference'], '--g', lw=2, label='threshold')
LIFASC_Htr3a.plot(LIFASC_model['grid_spike_times'], np.zeros(len(LIFASC_model['grid_spike_times']))+LIFASC_model['El_reference']-.01, '.b', ms=16)
LIFASC_Htr3a.set_xlim(x_lim)
LIFASC_Htr3a.set_ylim([-.09, -.09+ylim_v_jump_big])

LIFASC_Htr3aASC = LIFASC_Htr3a.twinx()
LIFASC_Htr3aASC.plot(LIFASC_model['time'], np.sum(LIFASC_model['AScurrents'], axis=1)*1.e9, 'r', lw=2)
LIFASC_Htr3aASC.set_xlim(x_lim)
LIFASC_Htr3aASC.set_ylim([-.8, -.8+ylim_ASC_jump])

LIFRASC_Htr3a.plot(LIFRASC_model['time'], LIFRASC_model['voltage']+LIFRASC_model['El_reference'], 'b', lw=4, label='model')
LIFRASC_Htr3a.plot(LIFRASC_model['time'], LIFRASC_model['threshold']+LIFRASC_model['El_reference'], '--g', lw=2, label='threshold')
LIFRASC_Htr3a.plot(LIFRASC_model['grid_spike_times'], np.zeros(len(LIFRASC_model['grid_spike_times']))+LIFRASC_model['El_reference']-.01, '.b', ms=16)
LIFRASC_Htr3a.set_xlim(x_lim)
LIFRASC_Htr3a.set_ylim([-.085, -.085+ylim_v_jump_big])

LIFRASC_Htr3aASC = LIFRASC_Htr3a.twinx()
LIFRASC_Htr3aASC.plot(LIFRASC_model['time'], np.sum(LIFRASC_model['AScurrents'],axis=1)*1.e9, 'r', lw=2)
LIFRASC_Htr3aASC.set_xlim(x_lim)
LIFRASC_Htr3aASC.set_ylim([-.8, -.8+ylim_ASC_jump])

LIFRASCAT_Htr3a.plot(LIFRASCAT_model['time'], LIFRASCAT_model['voltage']+LIFRASCAT_model['El_reference'], 'b', lw=4, label='model')
LIFRASCAT_Htr3a.plot(LIFRASCAT_model['time'], LIFRASCAT_model['threshold']+LIFRASCAT_model['El_reference'], '--g', lw=2, label='threshold')
LIFRASCAT_Htr3a.plot(LIFRASCAT_model['grid_spike_times'], np.zeros(len(LIFRASCAT_model['grid_spike_times']))+LIFRASCAT_model['El_reference']-.01, '.b', ms=16)
LIFRASCAT_Htr3a.set_xlim(x_lim)
LIFRASCAT_Htr3a.set_ylim([-.085, -.085+ylim_v_jump_big])

LIFRASCAT_Htr3aASC = LIFRASCAT_Htr3a.twinx()
LIFRASCAT_Htr3aASC.plot(LIFRASCAT_model['time'], np.sum(LIFRASCAT_model['AScurrents'], axis=1)*1.e9, 'r', lw=2)
LIFRASCAT_Htr3aASC.set_xlim(x_lim)
LIFRASCAT_Htr3aASC.set_ylim([-.8, -.8+ylim_ASC_jump])


##---512322162
dir_name=os.path.join(relative_path, 'mouse_nwb/specimen_'+ str(512322162))
the_sweeps=ctc.get_ephys_sweeps(512322162, os.path.join(dir_name, 'ephys_sweeps.json'))
nwb=ctc.get_ephys_data(512322162, os.path.join(dir_name, 'ephys.json'))
sweeps=get_sweep_num_by_name(the_sweeps, 'Noise 2')
data=[]
spike_times=[]
for s in sweeps:
    spike_times.append(nwb.get_spike_times(s))
    data.append(nwb.get_sweep(s)) 

print 'loading LIF'
LIF_model=pickle.load(open( "pkl_data/512322162Ctgf-2A-dgCre_LIF_model.pkl", "rb" ))
print 'loading LIFR'
LIFR_model=pickle.load(open( "pkl_data/512322162Ctgf-2A-dgCre_LIFR_model.pkl", "rb" ))
print 'loading LIFASC'
LIFASC_model=pickle.load(open( "pkl_data/512322162Ctgf-2A-dgCre_LIFASC_model.pkl", "rb" ))
print 'loading LIFRASC'
LIFRASC_model=pickle.load(open( "pkl_data/512322162Ctgf-2A-dgCre_LIFRASC_model.pkl", "rb" ))
print 'loading LIFRASCAT'
LIFRASCAT_model=pickle.load(open( "pkl_data/512322162Ctgf-2A-dgCre_LIFRASCAT_model.pkl", "rb" ))

current2.plot(LIF_model['time'], LIF_model['stimulus'], 'k', lw=4)
current2.set_xlim(x_lim)

for kk, st in enumerate(spike_times):
    data2.plot(LIF_model['time'], data[kk]['response'], lw=2, color=data_colors[kk])
    data2.plot(st, np.ones(len(st))*-(.065+.01*kk), '.', ms=10, color=data_colors[kk])
data2.set_xlim(x_lim)

LIF_Ctgf.plot(LIF_model['time'], LIF_model['voltage']+LIF_model['El_reference'], 'b', lw=4, label='model')
LIF_Ctgf.plot(LIF_model['time'], LIF_model['threshold']+LIF_model['El_reference'], '--g', lw=2, label='threshold')
LIF_Ctgf.plot(LIF_model['grid_spike_times'], np.zeros(len(LIF_model['grid_spike_times']))+LIF_model['El_reference']-.01, '.b', ms=16)
LIF_Ctgf.set_xlim(x_lim)
LIF_Ctgf.set_ylim([-0.10, -0.1+ylim_v_jump_small])

LIFR_Ctgf.plot(LIFR_model['time'], LIFR_model['voltage']+LIFR_model['El_reference'], 'b', lw=4, label='model')
LIFR_Ctgf.plot(LIFR_model['time'], LIFR_model['threshold']+LIFR_model['El_reference'], '--g', lw=2, label='threshold')
LIFR_Ctgf.plot(LIFR_model['grid_spike_times'], np.zeros(len(LIFR_model['grid_spike_times']))+LIFR_model['El_reference']-.01, '.b', ms=16)
LIFR_Ctgf.set_xlim(x_lim)
LIFR_Ctgf.set_ylim([-0.10, -0.1+ylim_v_jump_small])

LIFASC_Ctgf.plot(LIFASC_model['time'], LIFASC_model['voltage']+LIFASC_model['El_reference'], 'b', lw=4, label='model')
LIFASC_Ctgf.plot(LIFASC_model['time'], LIFASC_model['threshold']+LIFASC_model['El_reference'], '--g', lw=2, label='threshold')
LIFASC_Ctgf.plot(LIFASC_model['grid_spike_times'], np.zeros(len(LIFASC_model['grid_spike_times']))+LIFASC_model['El_reference']-.01, '.b', ms=16)
LIFASC_Ctgf.set_xlim(x_lim)
LIFASC_Ctgf.set_ylim([-.1, -.1+ylim_v_jump_big])

LIFASC_CtgfASC = LIFASC_Ctgf.twinx()
LIFASC_CtgfASC.plot(LIFASC_model['time'], np.sum(LIFASC_model['AScurrents'], axis=1)*1e9, 'r', lw=2)
#LIFASC_CtgfASC.plot(LIFASC_model['time'], LIFASC_model['AScurrents'][:,1]*1e9, 'c', lw=2)
LIFASC_CtgfASC.set_xlim(x_lim)
LIFASC_CtgfASC.set_ylim([-.65, -.65+ylim_ASC_jump])

LIFRASC_Ctgf.plot(LIFRASC_model['time'], LIFRASC_model['voltage']+LIFRASC_model['El_reference'], 'b', lw=4, label='model')
LIFRASC_Ctgf.plot(LIFRASC_model['time'], LIFRASC_model['threshold']+LIFRASC_model['El_reference'], '--g', lw=2, label='threshold')
LIFRASC_Ctgf.plot(LIFRASC_model['grid_spike_times'], np.zeros(len(LIFRASC_model['grid_spike_times']))+LIFRASC_model['El_reference']-.01, '.b', ms=16)
LIFRASC_Ctgf.set_xlim(x_lim)
LIFRASC_Ctgf.set_ylim([-.1, -.1+ylim_v_jump_big])

LIFRASC_CtgfASC = LIFRASC_Ctgf.twinx()
LIFRASC_CtgfASC.plot(LIFRASC_model['time'], np.sum(LIFRASC_model['AScurrents'], axis=1)*1e9, 'r', lw=2)
#LIFRASC_CtgfASC.plot(LIFRASC_model['time'], LIFRASC_model['AScurrents'][:,1]*1e9, 'c', lw=2)
LIFRASC_CtgfASC.set_xlim(x_lim)
LIFRASC_CtgfASC.set_ylim([-.65, -.65+ylim_ASC_jump])

LIFRASCAT_Ctgf.plot(LIFRASCAT_model['time'], LIFRASCAT_model['voltage']+LIFRASCAT_model['El_reference'], 'b', lw=4, label='model')
LIFRASCAT_Ctgf.plot(LIFRASCAT_model['time'], LIFRASCAT_model['threshold']+LIFRASCAT_model['El_reference'], '--g', lw=2, label='threshold')
LIFRASCAT_Ctgf.plot(LIFRASCAT_model['grid_spike_times'], np.zeros(len(LIFRASCAT_model['grid_spike_times']))+LIFRASCAT_model['El_reference']-.01, '.b', ms=16)
LIFRASCAT_Ctgf.set_xlim(x_lim)
LIFRASCAT_Ctgf.set_ylim([-.09, -.09+ylim_v_jump_big])

LIFRASCAT_CtgfASC = LIFRASCAT_Ctgf.twinx()
LIFRASCAT_CtgfASC.plot(LIFRASCAT_model['time'], np.sum(LIFRASCAT_model['AScurrents'], axis=1)*1.e9, 'r', lw=2)
#LIFRASCAT_CtgfASC.plot(LIFRASCAT_model['time'], LIFRASCAT_model['AScurrents'][:,1]*1e9, 'c', lw=2)
LIFRASCAT_CtgfASC.set_xlim(x_lim)
LIFRASCAT_CtgfASC.set_ylim([-.65, -.65+ylim_ASC_jump])

##add scale bar
LIFR_CtgfASC = LIFR_Ctgf.twinx() #even though there are not currents for R adding an axis for the current scale bar
LIFR_CtgfASC.set_xlim(x_lim)
LIFR_CtgfASC.set_ylim([-.65, -.65+ylim_ASC_jump/2.])
LIFR_Ctgf.plot([18.2, 18.2],[-.1, -.08], 'k', lw=6) #vertical
LIFR_Ctgf.plot([18.1, 18.2],[-.098, -.098], 'k', lw=6) #horizontal
#LIFR_Ctgf.annotate('100 ms', xy=(18.13, -.07), fontsize=16)
LIFASC_Ctgf.annotate('100 ms', xy=(18.13, .025), fontsize=16)
LIFR_Ctgf.annotate('20 mV', xy=(18.21, -.09), fontsize=16)
LIFR_Ctgf.annotate('100 pA', xy=(18.05, -.095), fontsize=16)
LIFR_CtgfASC.plot([18.1, 18.1],[-.63, -.53], 'k', lw=6) #vertical

plt.tight_layout(w_pad=4., h_pad=0)
#
##--turn off axes
current1.axis('off')
data1.axis('off')
LIF_Htr3a.axis('off')
LIFR_Htr3a.axis('off')
LIFASC_Htr3a.axis('off')
LIFRASC_Htr3a.axis('off')
LIFRASCAT_Htr3a.axis('off')
LIFASC_Htr3aASC.axis('off')
LIFRASC_Htr3aASC.axis('off')
LIFRASCAT_Htr3aASC.axis('off')

current2.axis('off')
data2.axis('off')
LIF_Ctgf.axis('off')
LIFR_Ctgf.axis('off')
LIFR_CtgfASC.axis('off')
LIFASC_Ctgf.axis('off')
LIFRASC_Ctgf.axis('off')
LIFRASCAT_Ctgf.axis('off')
LIFASC_CtgfASC.axis('off')
LIFRASC_CtgfASC.axis('off')
LIFRASCAT_CtgfASC.axis('off')

#--annotate
plt.annotate('Htr3a', xy=(.25, .98), xycoords='figure fraction', horizontalalignment='center',fontsize=18)
plt.annotate('Ctgf', xy=(.75, .98), xycoords='figure fraction', horizontalalignment='center',fontsize=18)
plt.annotate('current', xy=(.5, .93), xycoords='figure fraction', horizontalalignment='center',fontsize=18)
plt.annotate('data', xy=(.5, .83), xycoords='figure fraction', horizontalalignment='center',fontsize=18)
plt.annotate(r'GLIF$_1$', xy=(.5, .74), xycoords='figure fraction', horizontalalignment='center',fontsize=18)
blah=r'GLIF$_2 $'
blah+='\nR'
plt.annotate(blah, xy=(.5, .62), xycoords='figure fraction', horizontalalignment='center', fontsize=18)
blah=r'GLIF$_3 $'
blah+='\nASC'
plt.annotate(blah, xy=(.5, .47), xycoords='figure fraction', horizontalalignment='center', fontsize=18)
blah=r'GLIF$_4 $'
blah+='\nR + ASC'
plt.annotate(blah, xy=(.5, .28), xycoords='figure fraction', horizontalalignment='center', fontsize=18)
blah=r'GLIF$_5 $'
blah+='\nR+ASC\n+AT'
plt.annotate(blah, xy=(.5, .07), xycoords='figure fraction', horizontalalignment='center', fontsize=18)
plt.show()
