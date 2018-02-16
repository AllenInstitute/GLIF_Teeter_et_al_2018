import os
import numpy as np
import allensdk.core.json_utilities as ju
import warnings
from allensdk.api.queries.glif_api import GlifApi
import glif_sdk.lims_utilities as lu
from allensdk.core.nwb_data_set import NwbDataSet
from allensdk.api.api import Api

glif_api = GlifApi()

def get_file_path_endswith(path, ends_with):
    '''returns the path of the file ending with the specified string
    inputs:
        path: string
            directory path in which to search for file.  Usually a structured data directory in this analysis
        ends_with: string
            the end of the file name being searched for.  For example "GLIF1_neuron_config.json"
    returns:
        path to file (string)
    '''
    return os.path.join(path, os.listdir(path)[np.where([fname.endswith(ends_with) for fname in os.listdir(path)])[0][0]])

def get_pp_path(specimen_id_directory):
    '''returns the path (string) to the preprocessor file in the specified structured data directory (string)'''
    return get_file_path_endswith(specimen_id_directory, '_preprocessor_values.json')

def model_query(sp_id):
    m_query = lu.query("""SELECT  
            sp.id as specimen_id,
            sp.storage_directory as data_storage_dir,
            nm.id as model_id,
            nm.name as model_name,
            nm.storage_directory as model_storage_dir,
            nm.neuronal_model_template_id as model_template_id,
            nm.created_at as model_created,
            nm.updated_at as model_updated,
            nm.workflow_state as model_WFS,
            err.workflow_state as ephys_WFS,
            nmr.explained_variance_ratio
            FROM specimens sp 
            JOIN neuronal_models nm ON nm.specimen_id = sp.id
            JOIN projects p on p.id = sp.project_id
            JOIN ephys_roi_results err on err.id = sp.ephys_roi_result_id
            JOIN neuronal_model_runs nmr on nmr.neuronal_model_id = nm.id
            WHERE p.code in ('H301')
            AND specimen_id = %d""" % sp_id )   
    return m_query

def query_missing_exp_var(sp_id):
    q=model_query(int(sp_id))
    exp_var=[[m['model_name'], m['explained_variance_ratio']] for m in q]
    # check to make sure there is not more then one of the same model
    names=[ev[0] for ev in exp_var]
    if len(np.unique(names))>len(names):
        print q
        raise Exception('too many models with the same name')
    return exp_var

def get_ev_percent_from_calculated_file(ends_with, specimen_id_directory, opt_state, sweep_name):
    '''get the explained variance ratio from the structured data directory and converts it into 
    a percentage by multiplying it by 100.
    inputs:
        ends_with: string 
            matched with a file in the structured data directory in order to return a value. 
            Must be '_GLIF*_exp_var_ratio_10ms.json' where * can be integer values 1 though 5 
            corresponding to the model level
        specimen_id_directory: string 
            path to the structured data directory containing all the files for a neuron (as labeled 
            with a specimen id). For example: .../mouse_struc_data_dir/474637203_Htr3a-Cre_NO152
        opt_state: string
            corresponds to naming in "*_exp_var_ratio_10ms.json" file.  Can be "before_opt" or "after_opt"
        sweep_name: string
            corresponds to naming in "*_exp_var_ratio_10ms.json" file.  Can be "noise_1" or "noise_2"
    Returns
        ev: float if there is a file with a value, or np.nan if not
    '''
    try:
        file=get_file_path_endswith(specimen_id_directory,ends_with)
        dictionary=ju.read(file)
        ev=dictionary[opt_state][sweep_name]*100.
    except:
        ev=np.nan  
        
    return ev      

def get_ev_from_folder(ends_with, specimen_id_directory, model_string):
    '''get the explained variance ratio from the Allen Institute Cell Type Database
    for the specified model if the requested model exists in the structured data directory.
    i.e. if a LIF-R model is requested but there is not a *_GLIF3_neuron_config.json in 
    the structured data directory then a nan will be returned, regardless of if there is a
    LIF-R in the Database.
    inputs:
        ends_with: string that is to be matched with a file in the structured data directory in order 
            to return a value. i.e. if GLIF2 is being requested but there is not a GLIF2 file in the
            structured data directory, a nan will be returned regardless of whether there is a value 
            of explain variance in the database. For example this would happen if the model was excluded from 
            analysis because of an aberrant parameter.  
        specimen_id_directory: string 
            path to the structured data directory containing all the files for a neuron (as labeled 
            with a specimen id). For example: .../mouse_struc_data_dir/474637203_Htr3a-Cre_NO152
        model_string: string
            string searching for in model name: options '(LIF)', '(LIF-R)', '(LIF-ASC)', '(LIF-R_ASC)', '(LIF-R_ASC_A')
    returns:
        either or nan or the explained variance ratio for the requested model
        
    '''
    sp_id=os.path.basename(specimen_id_directory)[:9]
    print sp_id
    cre=os.path.basename(specimen_id_directory)[10:]
    try:
        nms=glif_api.get_neuronal_models(int(sp_id))[0]['neuronal_models']
        #NOTE THAT THIS COULD BE EMPTY WHICH I HAVE NOT ACCOUNTED FOR AND SINCE IT IS IN A 'TRY' WE WILL NOT KNOW
        exp_var=[[m['name'], m['neuronal_model_runs'][0]['explained_variance_ratio']] for m in nms] #list of madel name, explained variance pairs
    except:
#    if True:
        print sp_id, cre, 'is not in glif api: querying database'
        exp_var=query_missing_exp_var(sp_id)
    if len(exp_var)<0:
        # I have only seen this happen via the glif.api and it was on human data that shouldnt be there
        print sp_id, cre, 'exp var is empty: querying database'
        exp_var=query_missing_exp_var(sp_id)        
    
    if np.any([f.endswith(ends_with) for f in os.listdir(specimen_id_directory)]):            
        LIF_binary=[model_string in ev[0] for ev in exp_var]
        if np.any(LIF_binary):
            return exp_var[np.where(LIF_binary)[0][0]][1]
        else:
            raise Exception('there is %model_string in the directory but not in the api' % model_string)
    else:
        return np.nan 

def get_model_nwb_path_from_folder(ends_with, specimen_id_directory, s):
    '''Get the  path for the specified model if the the
    corresponding model exists in the structured data directory.
    Note: this will only work within the Allen Institute.  Use download_model_nwb
    for use outside the Institute.
    inputs:
        ends_with: string that is to be matched with a file in the structured data directory in 
            order to return a value. i.e. if GLIF2 is being requested but there is not GLIF2 file in the
            structured data directory, a nan will be returned regardless of whether there is a value of explain 
            variance in the database. For example this would happen if the model was excluded from 
            analysis because of an aberrant parameter.  
        specimen_id_directory: path to the structured data directory used in the rest of analysis
        s: search string in the neuronal model name
    returns:
        either a nan or the path of the Allen Institute model .nwb file.
        
    '''
    sp_id=os.path.basename(specimen_id_directory)[:9]
    if sp_id=='580895033':
        print 'skipping 580895033 which is not in the api, returning np.nan'
        return np.nan
    if len(glif_api.get_neuronal_models(int(sp_id)))>1: #basic check that data is the form expected
        raise Exception('why is there more than one list for %d' % sp_id)
    # get models for the neuron
    nms=glif_api.get_neuronal_models(int(sp_id))[0]['neuronal_models']
    paths=[[m['name'], m['neuronal_model_runs'][0]['well_known_files'][0]['path']] for m in nms] #list of model name
    
    if np.any([f.endswith(ends_with) for f in os.listdir(specimen_id_directory)]):            
        LIF_binary=[s in p[0] for p in paths]
        if np.any(LIF_binary):
            return paths[np.where(LIF_binary)[0][0]][1]
        else:
            raise Exception('there is %s in the directory but not in the api' % s)
    else:
        return np.nan 

def download_model_nwb_if_model_exists_in_SDD(ends_with, specimen_id_directory, s):
    '''
    Downloads the .nwb file and returns it's path for the specified model if the 
    corresponding model exists in the structured data directory
    inputs:
        ends_with: string that is to be matched with a file in the structured data directory in 
            order to return a value. i.e. if GLIF2 is being requested but there is not GLIF2 file in the
            structured data directory, a nan will be returned regardless of whether there is a value of explain 
            variance in the database. For example this would happen if the model was excluded from 
            analysis because of an aberrant parameter.  
        specimen_id_directory: path to the structured data directory used in the rest of analysis
        s: search string in the neuronal model name
    returns:
        either a nan or the path of the Allen Institute model .nwb file.
        
    '''
    sp_id=os.path.basename(specimen_id_directory)[:9]
    if sp_id=='580895033':
        print 'skipping 580895033 which is not in the api, returning np.nan'
        return np.nan
    if len(glif_api.get_neuronal_models(int(sp_id)))>1: #basic check that data is the form expected
        raise Exception('why is there more than one list for %d' % sp_id)
    # get models for the neuron
    nms=glif_api.get_neuronal_models(int(sp_id))[0]['neuronal_models']
    wkf_ids=[[m['name'], m['neuronal_model_runs'][0]['well_known_files'][0]['id']] for m in nms] #list
    if np.any([f.endswith(ends_with) for f in os.listdir(specimen_id_directory)]):            
        LIF_binary=[s in p[0] for p in wkf_ids]
        if np.any(LIF_binary):
            base_path=os.path.join(os.path.dirname(os.getcwd()), 'model_nwb_files')
            download_path=os.path.join(base_path, sp_id+'_'+ends_with[1:7]+'model.nwb')
            # if the file does not already exist locally, download it
            try:
                os.stat(base_path)
            except:
                os.mkdir(base_path)
            if not os.path.isfile(download_path):
                wkf_id=wkf_ids[np.where(LIF_binary)[0][0]][1]
                url = Api().construct_well_known_file_download_url(wkf_id)
                Api().retrieve_file_over_http(url, download_path)
            return download_path
        else:
            raise Exception('there is %s in the directory but not in the api' % s)
    else:
        return np.nan 

def get_sweep_num_by_name(sweeps, sweep_name):
    '''returns the sweep numbers for a specific sweep name in the file
    queried by ctc.get_ephys_sweeps(int(specimen_id))d
    '''
    return [ s['sweep_number'] for s in sweeps if s['stimulus_name'] == sweep_name ]

def non_nan_data(neuron_data, index):
    '''if the specified value is a nan it revomes the data from the list.
    inputs:
        neuron_data: list of lists
            each row corresponds to a neuron
        index: integer
            corresponds to the column (index in the sublists) that want to check for nans
    output:
        x: list of lists
            same as neuron_data input where the rows which contained a nan at the specified index are removed.
    '''
    x=[]
    for neuron in neuron_data:
        if not np.isnan(neuron[index]):
            x.append(neuron)
    return x

def float_data(neuron_data, index):
    '''if the specified value is not a float it revomes the data from the list.
    inputs:
        neuron_data: list of lists
            each row corresponds to a neuron
        index: integer
            corresponds to the column (index in the sublists) that want to check for nans
    output:
        x: list of lists
            same as neuron_data input where the rows which contained a nan at the specified index are removed.
    '''
    x=[]
    for neuron in neuron_data:
        try: 
            neuron[index]=float(neuron[index])
            x.append(neuron)
        except:
            pass
    return x


levels=['LIF', 
        'LIF_R', 
        'LIF_ASC',
        'LIF_R_ASC', 
        'LIF_R_ASC_AT']

cononical=[['490205998', 'Scnn1a-Tg2-Cre'],
           ['323834998', 'Scnn1a-Tg3-Cre'],
           ['477490421', 'Pvalb-IRES-Cre'],
           ['467003163', 'Rorb-IRES2-Cre-D'],
           ['512322162', 'Ctgf-2A-dgCre'],
           ['562535995', 'Vip-IRES-Cre'],
           ['569623233', 'Ndnf-IRES2-dgCre'],
           ['490376252', 'Cux2-CreERT2'],
           ['518750800', 'Chat-IRES-Cre-neo'],
           ['490263438', 'Ntsr1-Cre_GN220'],
           ['469704261', 'Nr5a1-Cre'],
           ['474637203', 'Htr3a-Cre_NO152'],
           ['313862134', 'Sst-IRES-Cre'],
           ['488380827', 'Rbp4-Cre_KL100'],
           ['580895033', 'Chrna2-Cre_OE25'],
           ['581058351', 'Nkx2-1-CreERT2']]

excitatory_cre_lines=['Scnn1a-Tg2-Cre', 
#                      'Slc17a6-IRES-Cre',
                      'Nr5a1-Cre',           
                      'Scnn1a-Tg3-Cre',     
                      'Rorb-IRES2-Cre-D', 
                      'Cux2-CreERT2',         
                      'Ntsr1-Cre_GN220',   
                      'Ctgf-2A-dgCre',
#                      'Tlx3-Cre_PL56',
                      'Rbp4-Cre_KL100']   
                      

inhibitory_cre_lines=['Sst-IRES-Cre', 
                      'Pvalb-IRES-Cre',                
                      'Htr3a-Cre_NO152', 
                      'Ndnf-IRES2-dgCre', 
                      'Vip-IRES-Cre', 
                      'Chrna2-Cre_OE25',
                      'Chat-IRES-Cre-neo',
                      'Nkx2-1-CreERT2']

pos_cre_lines=inhibitory_cre_lines+excitatory_cre_lines

# NOTE: Ntsr1-Cre AND Ntsr1-Cre_GN220 ARE THE SAME LINE.  LEAVING BOTH IN HERE SINCE THEY 
# ARE BOTH IN THE DATA SET BUT THE LABEL Ntsr1-Cre IS ELIMINATED IN OTHER DATA LISTS SO SHOULD BE DONE 
# IN THE DATA 


def check_and_organize_data(all_neurons):
    '''does some basic checks of the data.
    Inputs:
        all_neurons: list of lists
            each sublist contains data for a neuron.  Only important aspect for this code is that the 
            specimen ID is in the zeroth position of each sublist and cre-line is in the [1] index
        
    Returns:
        usable_neurons: list of lists
            a subset of all_neurons with the specified specimen IDs removed
        cre_dict: dict
            useable_neuron lists separated into a dictionary of the excitatory, inhibitory, individual cre 
            lines and  cre negative neurons. 
        
    '''
    
    #--check to see if any neurons are not in predefined cre lists 
    cre_in_dir=list(set([neuron[1] for neuron in all_neurons]))
    if np.any([ll not in pos_cre_lines for ll in cre_in_dir]):
        raise Exception('The '+', '.join([cre_in_dir[ii] for ii in np.where([ll not in pos_cre_lines for ll in cre_in_dir])[0]])+' cre lines are not in the list.  Update the list')
    #--check to see if any of the cre lines in the predefined list are not present in the data set.
    if np.any([ll not in cre_in_dir for ll in pos_cre_lines]):
        print 'The '+', '.join([pos_cre_lines[ii] for ii in np.where([ll not in cre_in_dir for ll in pos_cre_lines ])[0]])+' identified cre lines do not exist in this data set.'
        
   
    #-- separating different cre pos lines for distribution plotting
    # creating dictionary
    cre_dict={}
    for cre in pos_cre_lines:
        cre_dict[cre]=[]
    if 'Ntsr1-Cre' in cre_dict.keys():
        raise Exception('Ntsr1-Cre has not yet been changed to Ntsr1-Cre_GN220POS in the data set you are using: output dictionary will incorporate correct name')
        cre_dict.pop('Ntsr1-Cre')     

    # populating dictionary        
    cre_dict['excitatory']=[]
    cre_dict['inhibitory']=[]
    for neuron in all_neurons:
        cre=neuron[1]
        if cre=='Ntsr1-CrePOS':
            cre_dict['Ntsr1-Cre_GN220POS'].append(neuron)
        else:
            cre_dict[cre].append(neuron)
        if cre in excitatory_cre_lines:
            cre_dict['excitatory'].append(neuron)
        elif cre in inhibitory_cre_lines:         
            cre_dict['inhibitory'].append(neuron)
        else:
            raise Exception('cre line is not defined')
            
    return cre_dict

def convert_spike_times_to_ind(spike_times, dt):
    '''converts spike times to the indices of the spikes
    inputs
        spike_times: list of numpy arrays 
            each array contains the indices of the spikes in each sweep
        dt: float
            time step of data
    returns: list of numpy arrays
        spike_ind: list of numpy arrays 
            each array contains the times of the spikes in each sweep
    '''
    spike_ind=[]
    for st in spike_times:
        spike_ind.append((st/dt).astype(int))
    return spike_ind


def get_model_spike_ind_from_nwb(ends_with, specimen_id_directory, model_string, sweeps, dt, where_running):
    ''' Gets the times of spikes from the model nwb file and converts them to indices
    inputs.  If code running outside of the Allen Institute, the model .nwb file will be downloaded.       
        ends_with: string
            end of file searching for:  options "_GLIF1_neuron_config.json","_GLIF2_neuron_config.json' etc."
        specimen_id_directory: string
            path to structured data directory containing neuron_config, preprocessor, etc., files.            
        model_string: string
            string searching for in model name: options '(LIF)', '(LIF-R)', '(LIF-ASC)', '(LIF-R_ASC)', '(LIF-R_ASC_A')
        sweeps: list of integers
            integers refer to the sweep number in the electrophysiology .nwb data file
        dt: float
            time step of data
        where_running: string
            options are 'internal': the code is being run within the Institute and can therefore access the internal file system
                        'external': the code is being run outside the Institute and requires the use of the api to download the model nwb files
        Note that although ends_with and model_string should be appropriately paired, there is no check
        within this module to make sure that they are
    outputs: returns either a 
        nan if the there is not a model in the structured data directory corresponding to what the requested ends_with variable  
        or list of numpy arrays, each array contains the indices of the spikes in each sweep
            '''
    
    model_spike_times=get_model_spike_times_from_nwb(ends_with, specimen_id_directory, model_string, sweeps, where_running)
    if hasattr(model_spike_times, '__len__'): #if there is no model, model_spike_times will be nan and not have a length
        return convert_spike_times_to_ind(model_spike_times, dt)
    else:
        return np.nan


def get_model_spike_times_from_nwb(ends_with, specimen_id_directory, model_string, sweeps, where_running):
    ''' Gets the times of spike from the model nwb file
    inputs       
        ends_with: string
            end of file searching for:  options "_GLIF1_neuron_config.json","_GLIF2_neuron_config.json' etc."
        specimen_id_directory: string
            path to structured data directory containing neuron_config, preprocessor, etc., files.            
        model_string: string
            string searching for in model name: options '(LIF)', '(LIF-R)', '(LIF-ASC)', '(LIF-R_ASC)', '(LIF-R_ASC_A')
        sweeps: list of integers
            integers refer to the sweep number in the electrophysiology .nwb data file
        where_running: string
            options are 'internal': the code is being run within the Institute and can therefore access the internal file system
                        'external': the code is being run outside the Institute and requires the use of the api to download the model nwb files
        Note that although ends_with and model_string should be appropriately paired, there is no check
        within this module to make sure that they are
    outputs: returns either a 
        nan if the there is not a model in the structured data directory corresponding to what the requested ends_with variable  
        or 
        model_spike_times: list of numpy arrays 
            each array contains the times of the spikes in each sweep
        
            '''
    if where_running=='internal':
        path=get_model_nwb_path_from_folder(ends_with, specimen_id_directory, model_string)  #get nwb file path
    elif where_running=='external':
        path=download_model_nwb_if_model_exists_in_SDD(ends_with, specimen_id_directory, model_string)  #get nwb file path
    else:
        raise Exception('specify whether the code is being run internally or externally')
    if isinstance(path, basestring):
        model=NwbDataSet(path)
        model_spike_times=[]
        if sweeps==[]:
            raise Exception ('There are no sweeps to look at')
        for sw in sweeps:
            model_spike_times.append(model.get_spike_times(sw))
        return model_spike_times
    else:
        return np.nan

                
def convert_list_to_numpy(dictionary):
    '''Assumes dictionary input is dictionary with a depth of 1 where values associated with 
    keys are single lists of data.  This converts these lists to numpy arrays.
    '''
    for key in dictionary.keys():
        dictionary[key]=np.array(dictionary[key])
    
    return dictionary
    
def check_spike_times_identical(spike_time_list_of_arrays):
    '''Checks to make sure all of the arrays in the list are the same.
    Here in the analysis code, this function is used to verify that 
    the spike trains in different sweeps from the same stimulus are the
    same.  If they are different it is most likely due to the current injection 
    of all the sweeps not being at the same amplitude (this happens at times
    when a rig operator is adjusting baseline).
    inputs:
        spike_time_list_of_arrays: list of arrays
            each list contains an array of spike times of a sweep
    returns 1 if true, 0 if false and np.nan if the input is a float or nan
    
    '''
    if hasattr(spike_time_list_of_arrays, '__len__'):
        out =np.array([])
        for ii in range(1,len(spike_time_list_of_arrays)):
            if np.array_equal(spike_time_list_of_arrays[ii], spike_time_list_of_arrays[ii-1]):
                out=np.append(out, True)
            else:
                out=np.append(out, False)
        if np.all(out):
            return 1
        else:
            return 0
    else:
        return np.nan
