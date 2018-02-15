'''Calculates the standard error of the linear regression of the spike cut (Supplimentary Figure 1).
NOTE: THE FIRST COLUMN (THE INDECES) OF THE RESULTING CSV FILE WILL NEED TO BE DELETED IN ORDER
TO USE THE FILE IN OTHER CODE. The values in the resulting csv file may be compared to the saved 
version using compare_two_csv_files.py in the tests folder.
'''

import numpy as np
from scipy import stats
from scipy.optimize import curve_fit, fmin
import logging
import matplotlib.pyplot as plt 
import sys
import os
relative_path=os.path.dirname(os.getcwd())
sys.path.append(os.path.join(relative_path, 'libraries'))
from glif_sdk.find_spikes import align_and_cut_spikes, ALIGN_CUT_WINDOW, find_spikes_list
from data_library import convert_spike_times_to_ind
from allensdk.core.cell_types_cache import CellTypesCache
ctc = CellTypesCache(manifest_file=os.path.join(relative_path,'cell_types_manifest.json'))
from allensdk.core.nwb_data_set import NwbDataSet
import time
from scipy import signal
import pandas as pd

def get_sweep_ind_num_by_name(sweeps, sweep_name):
    '''returns the sweep indicies in the ephys_sweeps.nwb and the nwb sweep numbers for a specific sweep 
    name in the file queried by ctc.get_ephys_sweeps(int(specimen_id))
    '''
    stuff=[ (jj, s['sweep_number']) for jj, s in enumerate(sweeps) if s['stimulus_name'] == sweep_name ]
    ephys_list_indx=[]
    nwb_sweep_indx=[]
    for entry in stuff:
        ephys_list_indx.append(entry[0])
        nwb_sweep_indx.append(entry[1])
    return ephys_list_indx, nwb_sweep_indx
        
def calc_spike_cut_and_v_reset_via_expvar_residuals(all_current_list, 
                                                    all_voltage_list, dt, El_reference, deltaV, 
                                                    max_spike_cut_time=False):
    '''THIS IS A PARTIAL COPY OF THIS FUCTION FOUND IN THE PREPROCESSOR ALTERED TO SAVE THE EXP VARIANCE 
    This function calculates where the spike should be cut based on explained variance.  
    The goal is to find a model where the voltage after a spike maximally explains the
    voltage before a spike.  This will also specify the voltage reset rule. Note that if 
    the slope is greater than one or intercept is greater than zero, they are forced to 
    those values.  Regardless of the required force, the residuals are used.
    inputs: 
        spike_determination_method:  string 
            specifies the method used to find threshold
        all_current_list: list of array like 
            each array is a current traces injected into neuron
        all_voltage_list: list of array-like
            each array is the recorded voltage of the neuron corresponding to the input current trace
    returns: float
        values of the standard error of the best linear regression between the voltage before and
        after spike at different time points. (The best linear regression is used as the model)
    '''
   
    #--find the region of the spike needed for calculation of explained variance
    (temp_v_spike_shape_list, all_i_spike_shape_list, all_thresholdInd, waveIndOfFirstSpikes, spikeFromWhichSweep) \
                            = align_and_cut_spikes(all_voltage_list, all_current_list, dt) 
    
    #--change reference 
    all_v_spike_shape_list=[shape-El_reference-deltaV for shape in temp_v_spike_shape_list]

    # --setting limits to find explained variance
    if max_spike_cut_time and max_spike_cut_time < .010:
        expVarIndRangeAfterSpike = range(int(.001 / dt), int(max_spike_cut_time / dt))  #NOTE: THIS IS USED IN REFERENCE TO SPIKE TIME
    else:
        expVarIndRangeAfterSpike = range(int(.001 / dt), int(.010 / dt))  #NOTE: THIS IS USED IN REFERENCE TO SPIKE TIME

    list_of_endPointArrays = []  # this should end up a list of numpy arrays where each numpy array contains the indices of the v_spike_shape_list that are a certain time after the threshold
    for ii in expVarIndRangeAfterSpike:
        list_of_endPointArrays.append(np.array(all_thresholdInd) + ii)
        
    def line_force_slope_to_1(x,c):
        return x+c
    
    def line_force_int_to_0(x, m):  #TODO: CHANGE THIS TO REST TOD DISCONNECT EVERYTHING.
        return m*x

    # HERE YOU GET THE SLOPE AND INTERCEPT AT EACH POINT
    linRegress_error_4_each_time_end = []
    slope_at_each_time_end=[]
    intercept_at_each_time_end=[]
    varData_4_each_time_end = []
    varModel_4_each_time_end = []
    chi2 = []
    sum_residuals_4_each_time_end=[]
          
    xdata = np.array([v[int(all_thresholdInd[ii])] for ii, v in enumerate(all_v_spike_shape_list)])
    var_of_Vdata_beforeSpike = np.var(xdata)
    for jj, vectorOfIndAcrossWaves in enumerate(list_of_endPointArrays):  # these indices should be in terms of the spike waveforms
        v_at_specificEndPoint = [all_v_spike_shape_list[ii][int(index)] for ii, index in enumerate(vectorOfIndAcrossWaves)]  
        # the model of voltage reset is a linear regression between voltage before the spike and the voltage after the spike but it could be more complicated (for example as a function of current) 
        ydata = np.array(v_at_specificEndPoint)  # this is the voltage at the specified end point
        slope, intercept, r_value, p_value, std_err = stats.linregress(xdata, ydata)              
        linRegress_error_4_each_time_end.append(std_err)

        
    return min(linRegress_error_4_each_time_end)

def extract_data(sweep_name, the_sweeps, nwb):
    '''loads data for desired sweeps from the specified nwb file.
    inputs:
        sweep name: string
            string specifying sweeps to look for in the_sweeps dictionary
        the_sweeps: dictionary
            describes sweeps in experiment returned by the allensdk ctc.get_ephys_sweeps
        nwb: object
            experiment data returned by allensdk ctc.get_ephys_data
    returns:
        sweeps: list of integers
            sweep numbers corresponding to specified sweep_name
        data: list of dictionaries
            data extracted from nwb
        stim: array
            single array of first stimulus array (since all should be the same--note that 
            occasionally they are not the same due to amplitude adjustment during experiment)
        spike_ind: list of arrays
            each array has the indices of the spikes
        spike_times:  list of arrays
            each array has the times of the spikes
        dt: float
            time step of first sweep
    '''
        
    sweeps=get_sweep_num_by_name(the_sweeps, sweep_name)
    data=[]
    spike_times=[]
    for s in sweeps:
        spike_times.append(nwb.get_spike_times(s))
        data.append(nwb.get_sweep(s))    
    dt=1./data[0]['sampling_rate']
    stim=data[0]['stimulus']
    spike_ind=convert_spike_times_to_ind(spike_times, dt)
    
    return sweeps, data, stim, spike_ind, spike_times, dt 

def load_sweep(file_name, sweep_number, desired_dt=None, cut=0, bessel=False):
    '''load a data sweep and do specified data processing.
    Inputs:
        file_name: string
            name of .nwb data file
        sweep_number: 
            number specifying the sweep to be loaded
        desired_dt: 
            the size of the time step the data should be subsampled to
        cut:
            indicie of which to start reporting data (i.e. cut off data before this indicie)
        bessel: dictionary
            contains parameters 'N' and 'Wn' to implement standard python bessel filtering
    Returns:
        dictionary containing
            voltage: array
            current: array
            dt: time step of the returned data
            start_idx: the index at which the first stimulus starts (excluding the test pulse)
    '''
    ds = NwbDataSet(file_name)
    data = ds.get_sweep(sweep_number)

    data["dt"] = 1.0 / data["sampling_rate"]

    if cut > 0:
        data["response"] = data["response"][cut:]
        data["stimulus"] = data["stimulus"][cut:]        

    if bessel:
        sample_freq = 1. / data["dt"]
        filt_coeff = (bessel["freq"]) / (sample_freq / 2.) # filter fraction of Nyquist frequency
        b, a = signal.bessel(bessel["N"], filt_coeff, "low")
        data['response'] = signal.filtfilt(b, a, data['response'], axis=0)

    if desired_dt is not None:
        if data["dt"] != desired_dt:
            data["response"] = subsample_data(data["response"], "mean", data["dt"], desired_dt)
            data["stimulus"] = subsample_data(data["stimulus"], "mean", data["dt"], desired_dt)
            data["start_idx"] = int(data["index_range"][0] / (desired_dt / data["dt"]))
            data["dt"] = desired_dt

    if "start_idx" not in data:
        data["start_idx"] = data["index_range"][0]

    return {
        "voltage": data["response"],
        "current": data["stimulus"],
        "dt": data["dt"],
        "start_idx": data["start_idx"]
        }


def load_sweeps(file_name, sweep_numbers, dt=None, cut=0, bessel=False):
    '''load sweeps and do specified data processing.
    Inputs:
        file_name: string
            name of .nwb data file
        sweep_numbers: 
            sweep numbers to be loaded
        desired_dt: 
            the size of the time step the data should be subsampled to
        cut:
            indicie of which to start reporting data (i.e. cut off data before this indicie)
        bessel: dictionary
            contains parameters 'N' and 'Wn' to implement standard python bessel filtering
    Returns:
        dictionary containing
            voltage: list of voltage trace arrays
            current: list of current trace arrays
            dt: list of time step corresponding to each array of the returned data
            start_idx: list of the indicies at which the first stimulus starts (excluding 
                the test pulse) in each returned sweep
    '''
    data = [ load_sweep(file_name, sweep_number, dt, cut, bessel) for sweep_number in sweep_numbers ]

    return {
        'voltage': [ d['voltage'] for d in data ],
        'current': [ d['current'] for d in data ],
        'dt': [ d['dt'] for d in data ],
        'start_idx': [ d['start_idx'] for d in data ],
        } 

def main(folder):
    '''the error of the spike cutting
    input:
        folder: string
            path to specimen id folder in the data folder
    output:
        Writes explained variance ratio values to a .json file in data folder
    '''

    cell_type_directory=os.path.join(relative_path, 'mouse_nwb') #pointing to directory not the default ctc directory which will download to working directory

    specimen_id=int(os.path.basename(folder)[:9])
    print specimen_id
    cre=os.path.basename(folder)[10:]

    sweeps_file=os.path.join(cell_type_directory,'specimen_'+ str(specimen_id), 'ephys_sweeps.json')
    data_nwb_file_path=os.path.join(cell_type_directory,'specimen_'+ str(specimen_id), 'ephys.nwb')

    #---get spike times from the data
    the_sweeps=ctc.get_ephys_sweeps(specimen_id, sweeps_file)

    bessel = { 'N': 4, 'freq': 10000 }
    cut=0
    
    start_time=time.time()
    
    n1_idx, noise1_sweeps=get_sweep_ind_num_by_name(the_sweeps, 'Noise 1')
    n2_idx, noise2_sweeps=get_sweep_ind_num_by_name(the_sweeps, 'Noise 2')
    noise1_data = load_sweeps(data_nwb_file_path, noise1_sweeps, bessel=bessel)
    noise2_data = load_sweeps(data_nwb_file_path, noise2_sweeps, bessel=bessel)
    
    dt=noise1_data['dt'][0]

    RESTING_POTENTIAL = 'slow_vm_mv'
    subthresh_noise_current_list=[]
    subthresh_noise_voltage_list=[]
    noise_El_list=[]
    for ss in range(0, len(noise1_data['current'])):
        #--subthreshold noise has first epoch of noise with a region of no stimulation before and after (note the selection of end point is hard coded)
        subthresh_noise_current_list.append(noise1_data['current'][ss][noise1_data['start_idx'][ss]:int(6./dt)])
        subthresh_noise_voltage_list.append(noise1_data['voltage'][ss][noise1_data['start_idx'][ss]:int(6./dt)])
        noise_El_list.append(the_sweeps[n1_idx[ss]][RESTING_POTENTIAL]*1e-3)  

    #---------------------------------------------------------------        
    #---------find spiking indicies of spikes in noise--------------
    #---------------------------------------------------------------

#    # note that when using find_spikes_list without removing the testpulse a warning will result from calculating feature_data['base_v'] in the feature extractor (line 375) this not relavent here
    noise1_ind_wo_test_pulse_removed, _ = find_spikes_list(noise1_data['voltage'], dt)
    noise2_ind_wo_test_pulse_removed, _ = find_spikes_list(noise2_data['voltage'], dt)
    #Put all ISI ind in a 
    ISI_length=np.array([])
    for ii in range(len(noise1_ind_wo_test_pulse_removed)):
        ISI_length=np.append(ISI_length,noise1_ind_wo_test_pulse_removed[ii][1:]-noise1_ind_wo_test_pulse_removed[ii][:-1])
    for ii in range(len(noise2_ind_wo_test_pulse_removed)):
        ISI_length=np.append(ISI_length,noise2_ind_wo_test_pulse_removed[ii][1:]-noise2_ind_wo_test_pulse_removed[ii][:-1])
    min_ISI_len=np.min(ISI_length)
    
    # Els calculated from QC
    El_noise=np.mean(noise_El_list)
    std_error= calc_spike_cut_and_v_reset_via_expvar_residuals(noise1_data['current'], noise1_data['voltage'], 
                                                              dt, El_noise, 0, 
                                                              max_spike_cut_time=min_ISI_len*dt)

    return specimen_id, cre, std_error

    

if __name__ == '__main__':
   
    # load data out of configuration files
    folder_path=os.path.join(relative_path,'mouse_struc_data_dir') 
    folders=np.sort([os.path.join(folder_path, f) for f in  os.listdir(folder_path)])

    out=[]
    start_time=time.time()
    for kk, folder in enumerate(folders):
        print kk,'folder', folder
        specimen_id, cre, std_err=main(folder)
        print 'done at', (time.time()-start_time)/60., 'min'
        out.append([specimen_id, cre, std_err])        
        
    df=pd.DataFrame(out, columns=['specimen_id', 'cre', 'standard_err'])
    #df.to_csv('spikecut_standard_err_.csv')
    df.to_csv('out.csv')