'''This is in allensdk.internal.model.glif.threshold_adaptation import calc_spike_component_of_threshold_from_multiblip

'''

import numpy as np
import matplotlib.pyplot as plt
from find_spikes import find_spikes_ssq_list
from scipy.optimize import curve_fit
import logging


def calc_spike_component_of_threshold_from_multiblip(multi_SS, dt, dv_cutoff, thresh_frac,
                                                     MAKE_PLOT=False, SHOW_PLOT=False, BLOCK=False, PUBLICATION_PLOT=False):
    '''Calculate the spike components of the threshold by fitting a decaying exponential function to data to threshold versus time 
    since last spike in the multiblip data. The exponential is forced to decay to the local th_inf (calculated as the mean all of the 
    threshold values of the first spikes in each individual triblip stimulus).  For each multiblip stimulus in a stimulus set if there
    is more than one spike the difference in voltages from the first and second spike are plotted versus the separation in time.  Note that 
    this algorithm should only be implemented on multiblips sweeps where the neuron spike on the first and second blip.  Since there is 
    no easy way to do this, this erroneous data should not be provided to this algorithm (i.e is should be visually checked and eliminated 
    the preprocessor should hold back this data manually for now.)  
    
    #TODO: check to see if this is still true. Notes: The standard SDK spike detection algorithm does not work with the multiblip stimulus 
    due to artifacts when the stimulus turns on and off. Please see the find_multiblip_spikes module for more information.
    
    Input:
    
    multi_SS: dictionary
        contains multiblip information such as current and stimulus
    dt: float
        time step in seconds
    
    Returns:
    
    const_to_add_to_thresh_for_reset: float
        amplitude of the exponential fit otherwise known as a_spike.  Note that this is without any spike cutting
    decay_const: float
        decay constant of exponential. Note the function fit is a negative exponential which will mean this value will 
        either have to be negated when it is used or the functions used will have to have to include the negative.
    thresh_inf: float
    
    '''    
    multi_SS_v=multi_SS['voltage']
    multi_SS_i=multi_SS['current']

    # --get indicies of spikes
    spike_ind, _=find_spikes_ssq_list(multi_SS_v, dt, dv_cutoff, thresh_frac)
#    spike_ind=find_multiblip_spikes(multi_SS_i, multi_SS_v, dt) can depricate find_multiblip_spikes


    # eliminate spurious spikes that may exist
    spike_lt=[np.where(SI<int(2.0/dt))[0] for SI in spike_ind]
    if len(np.concatenate(spike_lt))>0:
        logging.warning('there is a spike before the stimulus in the multiblip')
    spike_ind=[np.delete(SI,ind)for SI, ind in zip(spike_ind, spike_lt)]
    spike_gt=[np.where(SI>int(3.0/dt))[0] for SI in spike_ind]
    if len(np.concatenate(spike_gt))>0:
        logging.warning('there is a spike after the stimulus in the multiblip')    
    spike_ind=[np.delete(SI,ind)for SI, ind in zip(spike_ind, spike_gt)]
    
    # intialize output lists
    time_previous_spike=[]
    threshold=[]
    thresh_first_spike=[]  #will set constant to this

    if MAKE_PLOT:
        plt.figure(figsize=(20,24))
        
    # Loop though each tri blip stimulus in muliblip stimulus
    for k in range(0, len(multi_SS_v)):
        thresh=[multi_SS_v[k][j] for j in spike_ind[k]] # voltage at all spikes in a single tri blip
        if thresh!=[] and len(thresh)>1:# there needs to be more than one spike so that we can find the time difference
            thresh_first_spike.append(thresh[0]) #Note that this finds the first spike (it might not be at the first stimulus blip)
            threshold.append(thresh[1])
            time_previous_spike.append((spike_ind[k][1]-spike_ind[k][0])*dt)
#            Old way when looked at all the spikes instead of just the first two (can be depricated; just here for record keeping)
#            threshold.append(thresh[1:])              
#            time_before_temp=[]
#            for j in range(1,len(thresh)):
#                time_before_temp.append((spike_ind[k][j]-spike_ind[k][j-1])*dt)
#            #for each spike calculate the time from the previous spike   
#            time_previous_spike.append(time_before_temp)    
        if MAKE_PLOT:
            plt.subplot(len(multi_SS_v)+1,1,1)
            plt.plot(np.arange(0, len(multi_SS_i[k]))*dt, multi_SS_i[k]*1e12, lw=2)
            plt.ylabel('current (pA)', fontsize=16)
            plt.xlim([2., 2.12])
            plt.title('Triple Short Square', fontsize=20)
            plt.subplot(len(multi_SS_v)+1,1,k+2)
            plt.plot(np.arange(0, len(multi_SS_v[k]))*dt, multi_SS_v[k], lw=2)
            plt.plot(spike_ind[k]*dt, thresh, '.k', ms=16)
            plt.xlim([2., 2.12])
            
    if MAKE_PLOT:
            plt.ylabel('voltage (V)', fontsize=16)
            plt.xlabel('time (s)', fontsize=16)
    
    if SHOW_PLOT:        
        plt.show(block=False)

    # put numbers into one vector for fitting of exponential function
    thresh_inf=np.mean(thresh_first_spike)  #note this threshold infinity isnt the one coming from single blip
    try: #this try here because sometimes even though have the traces there isnt more than one trace with two spikes
#--these two lines no longer needed because all single values now (depricate with lines up above) instead converst them to arrays
#        threshold=np.concatenate(threshold)
#        time_previous_spike=np.concatenate(time_previous_spike)  #note that this will have nans in it
        threshold=np.array(threshold)
        time_previous_spike=np.array(time_previous_spike)  #note that this will have nans in it
   
        if MAKE_PLOT:
            plt.figure()
            plt.plot(time_previous_spike, threshold, '.k', ms=16)
            plt.ylabel('threshold (mV)')
            plt.xlabel('time since last spike (s)')
    
        # calculate values of exponential function both if force function to local threshold infinity and not forcing to a value
        # (not forcing to a value seems less valid unless a bunch of points are added corresponding to the threshold of the
        # first spike at time equal infinity (because the first spike is a spike that happens where the spike before it was an
        # infinite time away)).  Therefore, the values that are obtained from forcing are the ones that are used.
        p0_force=[.002, -100.]
        p0_fit=[.002, -100., thresh_inf]
        
        #TODO: THIS WOULD BE BETTER IF IT CALLED THE ACTUAL FUNCTION IN THE NEURON METHODS THAT WAY THEY WOULD HAVE TO BE THE SAME
        (popt_force, pcov_force)= curve_fit(exp_force_c, (time_previous_spike, thresh_inf), threshold, p0=p0_force, maxfev=100000)
        (popt_fit, pcov_fit)= curve_fit(exp_fit_c, time_previous_spike, threshold, p0=p0_fit, maxfev=100000)

        # viewing fit functions
        ordered_time_for_plot=np.sort(time_previous_spike) #since time is not in order, making new time vector so that obtained fit curve can be plotted
        fit_force=exp_force_c((ordered_time_for_plot, thresh_inf), popt_force[0], popt_force[1])
        fit_fit=exp_fit_c(ordered_time_for_plot, popt_fit[0], popt_fit[1], popt_fit[2])
        if MAKE_PLOT:
            plt.plot(ordered_time_for_plot, fit_force, 'r', lw=4, label="exp fit (force const to thesh first spike)\n  k=%.3g, amp=%.3g" % (popt_force[1], popt_force[0]))
            plt.plot(ordered_time_for_plot, fit_fit, 'b', lw=4, label="exp fit (fit constant)\n  k=%.3g, amp=%.3g" % (popt_fit[1], popt_fit[0]))
            plt.legend()
            if SHOW_PLOT:        
                plt.show(block=False)
    
            if PUBLICATION_PLOT:       
                plt.figure(figsize=[14, 5])
                ax1=plt.subplot2grid((2, 2), (0,0))
                ax2=plt.subplot2grid((2, 2), (1,0))
                ax3=plt.subplot2grid((2, 2), (0,1), rowspan=2)
                for k in range(0, len(multi_SS_v)):
                    thresh=[multi_SS_v[k][j] for j in spike_ind[k]]
    
                    ax1.plot(np.arange(0, len(multi_SS_i[k]))*dt, multi_SS_i[k]*1.e12, lw=2)
                    ax1.set_ylabel('Current (pA)', fontsize=16)
                    ax1.set_xlim([2., 2.12])
                    ax1.axes.xaxis.set_ticklabels([])
                    #ax1.set_title('Triple Short Square', fontsize=20)
    
                    ax2.plot(np.arange(0, len(multi_SS_v[k]))*dt, multi_SS_v[k]*1.e3, lw=2)
                    ax2.plot(spike_ind[k]*dt, np.array(thresh)*1.e3, '.k', ms=16)
                    ax2.set_ylabel('Voltage (mV)', fontsize=16)
                    ax2.set_xlabel('Time (ms)', fontsize=16)
                    ax2.set_xlim([2., 2.12])
                    
    
                ax3.plot(time_previous_spike, threshold*1.e3, '.k', ms=16)
                ax3.set_ylabel('Threshold (mV)', fontsize=16)
                ax3.set_xlabel('Time since last spike (s)', fontsize=16)
#                ax3.set_title('Spiking component of threshold', fontsize=20)
                ax3.plot(ordered_time_for_plot, fit_force*1.e3, 'r', lw=4)#, label="exp fit: k=%.3g, amp=%.3g" % (popt_force[1], popt_force[0]))
                ax3.legend() 
                plt.tight_layout()
                plt.show()
        
        const_to_add_to_thresh_for_reset=popt_force[0]
        decay_const=popt_force[1]
        
        if decay_const >0:
            logging.critical('This neuron has an increasing decay value for the spike component of the threshold')
        if const_to_add_to_thresh_for_reset<0:
            logging.critical('This neuron has a negative amplitude for the spike component of the threshold') 

        #if the decay constant is positive, or the amplitute is negative set the amplitude to 0 so that there 
        #will be no spike component of the threshold
        if decay_const >0 or const_to_add_to_thresh_for_reset < 0:
            const_to_add_to_thresh_for_reset=0
            decay=-1.0  #note that this number doesnt matter since the amplitude is set to zero
            
        # This decay constant was originally forced to be positive (i.e. decay_const=abs(popt_force[1]))
        # and then it is negated everywhere it is utilized elsewhere in the code.  Now things are forced in 
        # a different way above.  However the decay constant still needs to be negated here for use in the 
        # rest of the code. 
        decay_const=-decay_const
    

    except Exception, e:
        logging.error(e.message)
        const_to_add_to_thresh_for_reset=None 
        decay_const=None
        
    return const_to_add_to_thresh_for_reset, decay_const, thresh_inf

def fit_avoltage_bvoltage(x, v_trace_list, El_list, spike_cut_length, all_spikeInd_list, th_inf, dt, a_spike, 
                     b_spike, fake=False):
    '''This is a version of fit_avoltage_bvoltage_debug that does not require the th_trace, 
    v_component_of_thresh_trace, and spike_component_of_thresh_trace needed for debugging. A
    test should be run to make sure the same output comes out from this and the debug function
    
    This function returns the squared error for the difference between the 'known' voltage
    component of the threshold obtained from the biological neuron and the voltage component 
    of the threshold of the model obtained with the input parameters (so that the minimum can be 
    searched for via fmin). The overall threshold is the sum of threshold infinity the spike component
    of the threshold and the voltage component of the threshold.  Therefore threshold infinity and 
    the spike component of the threshold must be subtracted from the threshold of the neuron in order
    to isolate the voltage component of the threshold.  In the evaluation of the model the actual
    voltage of the neuron is used so that any errors in the other components of the model will not 
    influence the fits here (for example, if a afterspike current was estimated incorrectly)
  
    Notes:
    * The spike component of the threshold is subtracted from the 
        voltage which means that the voltage component of the threshold should only be added to rules.
    * b_spike was fit using a negative value in the function therefore the negative is placed in the 
        equation. 
    * values in this function are in 'real' voltage as opposed to voltage
        relative to resting potential. 
    * current injection during the spike is not taken into account.  This seems reasonable as the 
        ion channels are open during this time and injected current may not greatly influence the neuron.
    
    x: numpy array
        x[0]=a_voltage input, x[1] is b_voltage_input, x[2] is th_inf
    v_trace_list: list of numpy arrays
        voltage traces (v_trace, El, and th_inf must be in the same frame of reference)
    El_list: list of floats
        reversal potential (v_trace, El, and th_inf must be in the same frame of reference)
    spike_cut_length: int
        number of indicies removed after initiation of a spike
    all_spikeInd_list: list of numpy arrays
        indicies of spike trains 
    th_inf: float
        threshold infinity (v_trace, El, and th_inf must be in the same frame of reference)
    dt: float 
        size of time step (SI units)
    a_spike: float
        amplitude of spike component of threshold.
    b_spike: float
        decay constant in spike component of the threshold
    fake: Boolean
        if True makes uses the voltage value of spike step-1 because there is not a voltage value at the spike
        step because it is set to nan in the simulator.
    '''
    a_voltage=x[0]
    b_voltage=x[1]

    total_err=0
    for v_trace, El, all_spikeInd in zip(v_trace_list, El_list, all_spikeInd_list):
        # Calculate values along the whole trace and then take the values at the spike ind
        internal_sp_comp_array=np.zeros(all_spikeInd[0]+spike_cut_length) 
        left_over=0
        #vector of spike component of of the threshold from each spike and previous spikes; note at first spike there is no spike component of threshold so initialized at zero   
        #Note that care has to be taken here to get make sure the right amount of decay is left over
        for spike_number in range(1,len(all_spikeInd)):      #skipping first spike since no residual spike component of threshold
            integration_length=all_spikeInd[spike_number]-all_spikeInd[spike_number-1]+1 #this is the amount of time that needs to be integrated over note that it is one longer than the interval because of the last value is added to the aspike for the next ISI
            local_spike_comp_of_threshold=spike_component_of_threshold_exact(a_spike+left_over, b_spike, np.arange(integration_length)*dt)
            internal_sp_comp_array=np.append(internal_sp_comp_array, local_spike_comp_of_threshold[:-1])
            left_over=local_spike_comp_of_threshold[-1]
        
        # Compute voltage component of threshold at biological spike (subtract th_inf and spike component of threshold
        # from biological voltage values at spike initiation)
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        #     NOTE THAT THERE IS AN ISSUE HERE USING FAKE DATA.  THE -1 IS HERE BECAUSE THE NEURON CROSSES THRESHOLD SOMETIME BETWEEN TWO INDICIES.
        #     FOR THE FAKE DATA THE TIME OF THE SPIKE (THE POINT FOLLOWING WHEN THE VOLTAGE CROSSES THRESHOLD) IS SET TO NAN.
        #     THE INTERPOLATED VOLTAGE CAN BE USED BUT THEN THE INTERPOLATED VOLTAGE MUST BE CALCULATED FOR THE TRUE VOLTAGE 
        #     TRACE AND POSSIBLY IN THE INTEGRATION.
        if fake:
            v_comp_of_th_at_each_spike_via_data=v_trace[all_spikeInd-1]-internal_sp_comp_array[all_spikeInd-1]-th_inf  #USE THIS FOR FAKE DATA
        else:
            v_comp_of_th_at_each_spike_via_data=v_trace[all_spikeInd]-internal_sp_comp_array[all_spikeInd]-th_inf  #USE THIS FOR REAL DATA (although probably not necessary)
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        # For each ISI, calculate the difference between the voltage dependent component of the threshold
        # and the value that would be determined via a model that uses the actual voltage of neuron.  
        sq_err = []    #list to store squared error between the model and biological threshold
        for spike_number in range(1, len(all_spikeInd)):      #loop over all ISI's in data
            v_start_ind = all_spikeInd[spike_number-1]+int(spike_cut_length) #dont want to use voltage during a spike
            end_ind = all_spikeInd[spike_number]
            v_in_ISI=v_trace[v_start_ind:end_ind]
            #voltage component of threshold at the beginning and end of the ISI
            #With fake data if go from one before fake data to fake data this should be exact
            theta0=v_comp_of_th_at_each_spike_via_data[spike_number-1] #this assumes that the voltage component of the threshold does not change over the time period of the spike
            # not sure this makes sense any moretheta0=v_comp_of_th_at_each_spike_via_data[spike_number-1]+sp_comp_of_offset_at_spike_list[spike_number-1] #use this is you want to add the biological component back on
            theta1=v_comp_of_th_at_each_spike_via_data[spike_number]
            tvec=np.arange(len(v_in_ISI))*dt
    
            #analytical solution should be exact with fake data--small differences could be because of the differences in the voltage at 
            #spike indicies are off by one
            model=+theta0*np.exp(-b_voltage*dt*(end_ind-v_start_ind))+a_voltage*np.exp(-b_voltage*tvec[-1])*np.sum(dt*(v_in_ISI-El)*np.exp(b_voltage*tvec))
            err = (theta1-model)**2
            if ~np.isnan(err):
                sq_err.append(err)
        
        total_err+=np.sum(sq_err)
    return total_err

def exp_force_c((t, const), a1, k1):
    return a1*(np.exp(k1*t))+const

def exp_fit_c(t, a1, k1, const):
    return a1*(np.exp(k1*t))+const