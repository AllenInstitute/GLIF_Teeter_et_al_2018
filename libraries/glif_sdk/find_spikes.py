import numpy as np
from allensdk.ephys.feature_extractor import EphysFeatureExtractor
import allensdk.ephys.ephys_extractor as efex
import allensdk.ephys.ephys_features as ft

def find_spikes_list(voltage_list, dt):
    v_set = [ v * 1e3 for v in voltage_list ]  #converts voltage to millivolt for the EphysSweepSetFeatureExtractor
    t_set = [ np.arange(0, len(v)) * dt for v in voltage_list ]
    i_set = [ np.zeros(len(v)) for v in voltage_list ]
    
    ext = efex.EphysSweepSetFeatureExtractor(t_set, v_set, i_set, filter=None)
    ext.process_spikes()
    sweep_spikes = [ s.spikes() for s in ext.sweeps() ]
    
    out_idx = [ np.array([ int(s['threshold_index']) for s in spikes ]) for spikes in sweep_spikes ]
    out_v = [ np.array([ s['threshold_v'] for s in spikes ]) for spikes in sweep_spikes ]
    
    return out_idx, out_v

SHORT_SQUARE_MAX_THRESH_FRAC = 0.1
def find_spikes_ssq_list(voltage_list, dt, dv_cutoff, thresh_frac):
    v_set = [ v * 1e3 for v in voltage_list ]
    t_set = [ np.arange(0, len(v)) * dt for v in voltage_list ]
    i_set = [ np.zeros(len(v)) for v in voltage_list ]

    thresh_frac = max(SHORT_SQUARE_MAX_THRESH_FRAC, thresh_frac)
    
    ext = efex.EphysSweepSetFeatureExtractor(t_set, v_set, i_set, 
                                             dv_cutoff=dv_cutoff, 
                                             thresh_frac=thresh_frac,
                                             filter=None)
    ext.process_spikes()
    sweep_spikes = [ e.spikes() for e in ext.sweeps() ]
    
    out_idx = [ np.array([ int(s['threshold_index']) for s in spikes ]) for spikes in sweep_spikes ]
    out_v = [ np.array([ s['threshold_v'] for s in spikes ]) for spikes in sweep_spikes ]
    
    return out_idx, out_v

ALIGN_CUT_WINDOW = np.array([ 0.002, 0.015 ]) #specifies time before [0] and after [1] spike initiation for window

def align_and_cut_spikes(voltage_list, current_list, dt, spike_window = None):
    ''' This function aligns the spikes to some criteria and returns a current and voltage trace of 
    of the spike over a time window.  Also returns zero crossing,and threshold 
    in reference to the aligned spikes.
    Inputs:
        voltage_list: list of numpy arrays
            each array is a voltage trace over time
        current_list: list of numpy arrays
            each array is the corresponding current injected over time
        dt: float
            size of time step
        spide_window
    
    Returns:
        spike_shapes, current_shapes, aligned_spike_ind, wave_index_of_first_spikes, spike_sweeps
    '''
    if spike_window is None:
        spike_window = ALIGN_CUT_WINDOW

    spike_shapes = [] # will contain voltage in the defined window
    current_shapes = [] # will contain current in the defined window
    index_before_spike = int(spike_window[0] / dt) # index before spike threshold 
    index_after_spike = int(spike_window[1] / dt) # index after spike threshold
    aligned_spike_ind = np.array([]) #indicies of ....
    spike_sweeps = []
    spikes_per_trace = np.array([]) #will contain total number of spikes in each trace
    
    spike_ind_list, _ = find_spikes_list(voltage_list, dt)

    # consider values one sweep at a time
    for jj, voltage_and_current_and_spike in enumerate(zip(voltage_list, current_list, spike_ind_list)):
        voltage, current, whole_trace_spike_ind = voltage_and_current_and_spike 
        
        spikes_per_trace = np.append(spikes_per_trace, len(whole_trace_spike_ind))
        
        alignment_ind = whole_trace_spike_ind                        
        aligned_spike_ind = np.append(aligned_spike_ind, np.ones(len(whole_trace_spike_ind)) * index_before_spike)

        # print 'alignment_ind', alignment_ind
        spike_delimiters = [(ind - index_before_spike, ind + index_after_spike) for ind in alignment_ind]
        for d in spike_delimiters: 
            # this 'if' statement makes sure we don't cause a ValueError
            if min(d) > 0 and max(d) < len(voltage) - 1:
                spike_trace = voltage[d[0]:d[1]]
                current_trace = current[d[0]:d[1]]
                spike_shapes.append(spike_trace)
                current_shapes.append(current_trace)                  
                spike_sweeps.append(jj)

            
    # note: that depending on how things were aligned, all of one of the values will be the same.
    print "spikes_per_trace", spikes_per_trace
    temp = np.append(0, np.cumsum(spikes_per_trace))  
    print 'temp', temp 
    wave_index_of_first_spikes = [int(ii) for ii in list(temp[range(0, len(temp) - 1)])]         
    print "in cut spikes: wave_index_of_first_spikes ", wave_index_of_first_spikes

    return spike_shapes, current_shapes, aligned_spike_ind, wave_index_of_first_spikes, spike_sweeps