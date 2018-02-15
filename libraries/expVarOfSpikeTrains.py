'''Written by Corinne Teeter on 10-6-2014.  This group of functions calculates a 
variety of different explained variance measures of spike trains'''

import numpy as np
from scipy.signal import fftconvolve, gaussian
import warnings

def test_list_array_data_struct(voltage_list):
    '''Test to make sure data is in the correct format
    If voltage_list is not in a list of numpy arrays an exception will be raised.
    '''
    if not isinstance(voltage_list, list):
        raise Exception('data is not in a list')
    else:
        if not all([isinstance(trace, np.ndarray) for trace in voltage_list]):
            raise Exception('the data is not in a list of numpy arrays')
        else:
            pass

def PSTH(convolved_data):
    '''Returns the average trace
    '''
    test_list_array_data_struct(convolved_data)
    return np.mean(convolved_data, axis=0)

def makeBinaryTrainsFromInd(steps, spike_ind_list):
    '''Converts lists of arrays of spike indicies into binary arrays of the length specified 
    where a spike is denoted by a 1 and no spike is zero
    Inputs
        spike_ind_list: list of numpy arrays. Each numpy array contains indicies of spikes
        steps: length of binary array
    Returns list of numpy arrays of binary trains
    Note: right now this requires all the arrays to correspond to sweeps of the same length.'''
    
    test_list_array_data_struct(spike_ind_list)    
    binary_train=[]
    for train in spike_ind_list:
        temp=np.zeros(steps)
        if not train==[]:
            temp[list(train)]=1.
        binary_train.append(temp)
        
    return binary_train  

def makeNormGaussFunc4Convolve(sigma, dt):
    '''Calculates the appropriate gaussian for a given sigma and dt
    input:
        sigma: float
            standard deviation of desired Gaussian
        dt: float
            time step of data
    '''
    convertSigma=sigma/dt #convert the sigma into terms of dt to make window
    windowLength = convertSigma*10 #you want your window to be 10 times larger than your standard deviation
    gaussFunc_unNormalized = gaussian(int(windowLength), convertSigma)
    gaussFunc_norm = gaussFunc_unNormalized/np.sum(gaussFunc_unNormalized)
        
    return gaussFunc_norm 

def makeConvolvedSpikeTrains(spike_ind_list, dt, sigma, steps, convolveType='same', plot=False):
    '''*makes convolved spike trains.
        spike_ind_data    :list of numpy arrays with indicies of spike times
        dt                :sample step size in seconds
        sigma             :standard deviation of gaussian window for convolution
        steps             :length of spike train needed to specify the size of vectors.
                            Note at the moment this assumes the spike trains in the list
                            are the same length.
        Returns convolved: list of numpy arrays of convolved trains.                    
    '''
    test_list_array_data_struct(spike_ind_list)
    
    #--make convolution window--
    gaussFunc_norm=makeNormGaussFunc4Convolve(sigma, dt)
                
    #--convert spike indicies to binary
    binary_train_list=makeBinaryTrainsFromInd(steps, spike_ind_list)    

    #---make sure gaussFunc_norm is smaller than the train lengths otherwise the output of the convolution will be incorrect.
    if all(x==[len(train) for train in binary_train_list][0] for x in [len(train) for train in binary_train_list]): #Check to see if all the trains are the same length
        if len(gaussFunc_norm)>len(binary_train_list[0]):
            halfTrain=len(binary_train_list[0])/2
            halfGauss=len(gaussFunc_norm)/2
            shortenedGauss=gaussFunc_norm[halfGauss-halfTrain:halfGauss+halfTrain]
            #--make sure the new Gaussian is actually shorter
            if len(shortenedGauss)>len(binary_train_list[0]):
                raise Exception('Your Gaussian function is still longer than the train even though you have shortened it')
            #--if the area under the shortened gaussian is almost one, accept it as good enough and renormalize it 
            shortGaussArea=sum(shortenedGauss)
            if shortGaussArea>0.95:
                reNormShortGauss=shortenedGauss+(1.-shortGaussArea)/len(shortenedGauss)
            else:
                raise Exception('Your shortened Gaussian does not have the appropriate area under larger than .95 to enable justification of use')
            gaussFunc_norm=reNormShortGauss
    else:
        raise Exception('write this section to deal with different length trains.')
                   
    #--convolve spike train with a gaussian   
    convolved=[] 
    for i, train in enumerate(binary_train_list):
        #Note using 'full' will make sure the sum of the convolved trains is equal to the number of the 
        #spikes is equal to the number of spikes however the padding will not always be the same and will 
        #have implications when summing across the vectors to calculate the PSTH or taking the dot product
        #(they all have to be the same length)
        convolved.append(fftconvolve(train, gaussFunc_norm, convolveType))
        #print "train original length", len(train), "len new train", len(convolved[i]), 'sigma', sigma, 'len gauss function',len(gaussFunc_norm)     
        if len(train)!=len(convolved[i]) and convolveType=='same':
            print "interation", i, "len(train)", len(train), "len(convolved[i])", len(convolved[i]), "convolveType", convolveType
            raise Exception('The convolved train is not the same length as the original train')

    return convolved

def basic_expVar(trace1, trace2):
    '''This is the fundamental calculation that is used in all different types of explained variation.  
    At a basic level, the explained variance is calculated between two traces.  These traces can be PSTH's
    or single spike trains that have been convolved with a kernel (in this case always a Gaussian)
    Input:
        trace 1 & 2:  1D numpy array containing values of the trace.  (This function requires numpy array
                        to ensure that this is not a multidemensional list.)
    Returns:
        expVar:  float value of explained variance
    '''
    
    var_trace1=np.var(trace1)
    var_trace2=np.var(trace2)
    var_trace1_minus_trace2=np.var(trace1-trace2)

    if var_trace1_minus_trace2 == 0.0:
        return 1.0
    else:
        return (var_trace1+var_trace2-var_trace1_minus_trace2)/(var_trace1+var_trace2)
    

def acrossSet_PWExpVar(convolvedModel_list, convolvedData_list):
    '''returns the mean of pairwise explained variance of a train with a  data set.
    Inputs:
        convolvedModel_list & convolvedData_list: list of numpy arrays. Note
            names are interchangable model and data where just used because of likely
            use cases.
    Returns:
        expVar: float
    '''
    test_list_array_data_struct(convolvedData_list)
    test_list_array_data_struct(convolvedModel_list)

    pw_expVar=[]
    for jj in range(0, len(convolvedModel_list)): 
        for ii in range(0, len(convolvedData_list)):
            pw_expVar.append(basic_expVar(convolvedData_list[ii], convolvedModel_list[jj]))                

#    print '\t pw_expVar in acrossSet_PWExpVar', pw_expVar
    expVar=np.mean(pw_expVar)
    
    return expVar

def fromSpikesToPWExpVar_ofDataSet(spikeInd_list, dt, sigma, steps, plot=False):
    '''makes convolved spike trains.
    spike_ind_data    :list of numpy arrays with indicies of spike times.  Note: for the purposes here
                        these should correspond to the same sweep stimulus since you are finding the  
                        explained variance of the data set.
    dt                :sample step size in seconds
    sigma             :standard deviation of gaussian window for convolution
    steps             :length of spike train needed to specify the size of vectors.
                        Note at the moment this assumes the spike trains in the list
                        are the same length.
    Returns explained variance                        
    '''
    #--make convolved trains
    convolvedTrain_list=makeConvolvedSpikeTrains(spikeInd_list, dt, sigma, steps, convolveType='full', plot=False)   
    avgConvolvedTrain = np.mean(np.array(convolvedTrain_list), axis=0)

    expVar=acrossSet_PWExpVar([avgConvolvedTrain], convolvedTrain_list)
    return expVar
            
def fromSpikesToPWExpVar_ofDataWModel(spikeIndData_list, spikeIndModel_list,  dt, sigma, steps):
    '''makes convolved spike trains.
    spikeIndData_list     :list of numpy arrays with indicies of spike times
    spikeIndModel_list   :list  of numpy arrays of indicies of spike times
    dt                    :sample step size in seconds
    sigma                 :standard deviation of gaussian window for convolution
    steps                 :length of spike train needed to specify the size of vectors.
                            Note at the moment this assumes the spike trains in the list
                            are the same length.
TEST TO MAKE SURE THIS WORKS WHEN YOU HAVE A SINGLE TRAIN
    '''
    test_list_array_data_struct(spikeIndData_list)
    test_list_array_data_struct(spikeIndModel_list)

    convolvedData_list=makeConvolvedSpikeTrains(spikeIndData_list, dt, sigma, steps)
    convolvedModel_list=makeConvolvedSpikeTrains(spikeIndModel_list, dt, sigma, steps)

#    print '*********sigma', sigma
    expVar=acrossSet_PWExpVar(convolvedModel_list, convolvedData_list)   
    return expVar



    
