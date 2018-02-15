'''This code will look at files in the structured folder, identify 
models which should be excluded based on our exclusion criteria and 
eliminate them from structured folder.
'''
import allensdk.core.json_utilities as ju
import os
import numpy as np
import pandas as pd
import shutil
import sys
relative_path=os.path.dirname(os.getcwd())
sys.path.append(os.path.join(relative_path, 'libraries'))
from data_library import get_file_path_endswith, get_pp_path

#-------------------------------------------------
# !!!!!! specify species 'mouse or human' !!!!!!!!
#-------------------------------------------------
species='mouse'
#species='human'
#-------------------------------------------------
#-------------------------------------------------

if species =='mouse':
    structured_data_directory='mouse_struc_data_dir'  #Note that this is accessing data local to this folder not the global 'mouse_data' folder saved one directory up for convenience
elif species=='human':
    structured_data_directory='human_struc_data_dir'
else:
    raise Exception('species not recognized')

def exclude_via_spike_comp_of_thr(folder_path):
    '''return list of specimen ids excluded that should be excluded from the directory
    via spike component of threshold.  These specimens_id should be removed from all 
    levels with reset rules.
    input:
        folder_path: string
            Path to a structured data directory.   
            Inside the path should exist a series of folders with a
            name format: specimen id_cre for each neuron.  Inside those inner folders
            are the neuron configs of the GLIF models and preprocessor files. 
    returns:
        spike_comp_of_thr_exclusion_list: list of strings
        list of specimen ids to be excluded
    '''
    initial_sp_ids=[f[0:9] for f in os.listdir(folder_path)]
    folders=[os.path.join(folder_path, f) for f in  os.listdir(folder_path)]
    
    none_exclusion=[]
    spike_comp_of_thr_exclusion_list=[]
    total_neurons_with_stim_for_reset_rules=0
    exclusion_a_greaterthan_p02=[]
    exclusion_a_lessthan_0=[]
    exclusion_1_over_b_greaterthan_p1=[]
        
    for folder in folders:
        specimen_ID=os.path.basename(folder)[:9]
        pp_file=get_pp_path(folder)
        pp_dict=ju.read(pp_file)
        if pp_dict['threshold_adaptation']['a_spike_component_of_threshold'] is not None and pp_dict['threshold_adaptation']['b_spike_component_of_threshold'] is not None:
            #get overall idea of exclusion reasons
            total_neurons_with_stim_for_reset_rules=total_neurons_with_stim_for_reset_rules+1
            if pp_dict['threshold_adaptation']['a_spike_component_of_threshold']>.02:
                exclusion_a_greaterthan_p02.append(specimen_ID)
            if pp_dict['threshold_adaptation']['a_spike_component_of_threshold']<=0:
                exclusion_a_lessthan_0.append(specimen_ID)
            if 1./pp_dict['threshold_adaptation']['b_spike_component_of_threshold']>.1:
                exclusion_1_over_b_greaterthan_p1.append(specimen_ID)
        else:
            none_exclusion.append(specimen_ID)
    
    spike_comp_of_thr_exclusion_list=list(set(exclusion_a_greaterthan_p02+exclusion_a_lessthan_0+exclusion_1_over_b_greaterthan_p1+none_exclusion))
    print 'of', len(folders), 'neurons,', total_neurons_with_stim_for_reset_rules, 'neurons have the stimuli necessary for reset rules. Of those, excluding ', len(list(set(exclusion_a_greaterthan_p02+exclusion_a_lessthan_0+exclusion_1_over_b_greaterthan_p1))), 'excluded due to a_bad spike or b_spike'
    print '\t this leaves a total of', len(spike_comp_of_thr_exclusion_list), 'excluded due to a_spike or b_spike and a total of ', len(folders)-len(spike_comp_of_thr_exclusion_list), 'for spike component of threshold analysis' 
    return spike_comp_of_thr_exclusion_list

def exclude_via_v_comp_of_th(folder_path):
    '''return list of specimen ids excluded that should be excluded from the directory
    via threshold adaptation.  These specimens_id should be removed from level 5 GLIF. 
    input:
        folder_path: string
            Path to a structured data directory.   
            Inside the path should exist a series of folders with a
            name format: specimen id_cre for each neuron.  Inside those inner folders
            are the neuron configs of the GLIF models and preprocessor files. 
    returns:
        exclusion_list: list of strings
            list of specimen ids to exclude
    '''
    initial_sp_ids=[f[0:9] for f in os.listdir(folder_path)]
    folders=[os.path.join(folder_path, f) for f in  os.listdir(folder_path)]
    none_exclusion=[]
    exclusion_a_lessthan_neg50=[]
    exclusion_b_lessthan_p1=[]
    strange_pp_exclusion=[] #exclusion for preprocessor files that do not have the correct format
    
    for folder in folders:
        specimen_ID=os.path.basename(folder)[:9]
        pp_file=get_pp_path(folder)
        pp_dict=ju.read(pp_file)
        try:
            if pp_dict['threshold_adaptation']['a_voltage_comp_of_thr_from_fitab'] is not None and pp_dict['threshold_adaptation']['b_voltage_comp_of_thr_from_fitab'] is not None:
                if pp_dict['threshold_adaptation']['a_voltage_comp_of_thr_from_fitab']<-50.:
                    exclusion_a_lessthan_neg50.append(specimen_ID)
                if pp_dict['threshold_adaptation']['b_voltage_comp_of_thr_from_fitab']<.1 :
                    exclusion_b_lessthan_p1.append(specimen_ID)
            else: 
                none_exclusion.append(specimen_ID)
        except:
            print folder, 'DOES NOT LOOK LIKE A NORMAL PREPROCESSOR FILE'
            strange_pp_exclusion.append(specimen_ID)   
            

    print len(set(exclusion_a_lessthan_neg50+exclusion_b_lessthan_p1)), 'models of', len(folders), 'have a bad a_volotage or b_voltage. Note that many of these dont have reset rules either.' 
    print '\t', len(none_exclusion), 'of', len(folders), 'total models dont have threshold adaptation'
    print set(strange_pp_exclusion), 'specimen ids have a strange preprocessor file'
    exclusion_list=list(set(exclusion_a_lessthan_neg50+exclusion_b_lessthan_p1+none_exclusion+strange_pp_exclusion))  

    return exclusion_list

def general_exclusions(folder_path, 
                       n_in_cre=5, 
                       resist=True, 
                       th_inf_bad=True, 
                       spike_cut=True, 
                       ev=.2,
                       accidental_exclusion=True): 
    '''Returns a list of specimen ids that will be excluded from all levels and analysis.
    Note that some of these exclusions may be irrelevant for the curated data via the Allen Institute
    Cell Types Database. Nonetheless, I leave these exclusions here for my own use on internal data.
    Inputs:
        folder_path: string
            path to the structured data directory 
        n_in_cre: integer
            specifies number of neurons that have exist in a cre line to include them in the analysis
        resist: boolean
            if True, exclude neurons which have a calculated resistance over 1000 MOhms
        th_inf_bad: boolean
            if True, exclude neurons which have a calculated threshold less than -60 mV
        spike_cut: boolean
            if True, exclude neurons which have an intercept larger than 30 mV after fitting the spike cut length 
        ev: float
            exclude neurons that have an GLIF1 explained variance on noise 1 less than the provided value. 
            Note that noise 1 is used because exclusion criteria are only applied to training data 
        accidental_exclusion: boolean
            one neuron was either accidentally excluded from the analysis or removed for an unknown reason.  
            If True, exclude this neuron. 
    Returns: 
        exclude_me_sp_ids: list of strings
            list of neuron by specimen ids to be eliminated from the structured data directory
    '''
    initial_sp_ids=[f[0:9] for f in os.listdir(folder_path)]
    print 'GENERAL EXCLUSIONS: there will be overlap in numbers below i.e. some models will be excluded for more than one reason'
    print '\tTotal number of preprocessed files:', len(initial_sp_ids)
    folders=[os.path.join(folder_path, f) for f in  os.listdir(folder_path)]
    
    strange_pp_exclusion=[] #exclusion for preprocessor files that do not have the correct format

    # exclude via slope and intercept from spike cutting results
    spike_cutting_exclusions=[]
    if spike_cut:
        for folder in folders:
            specimen_ID=os.path.basename(folder)[:9]
            pp_file=get_pp_path(folder)
            pp_dict=ju.read(pp_file)
            try:
                if pp_dict['spike_cutting']['NOdeltaV']['intercept'] > .03:
                    spike_cutting_exclusions.append(specimen_ID)
            except:
                print folder, 'DOES NOT LOOK LIKE A NORMAL PREPROCESSOR FILE'
                strange_pp_exclusion.append(specimen_ID)   
        print '\t', len(set(spike_cutting_exclusions)), 'neurons were excluded for having an intercept larger than .03'
    
    # exclude based on the measured experimental threshold
    # note that the experimental threshold is the same for all models of the same neuron so just look at GLIF1 file.
    th_inf_exclusion_list=[]
    if th_inf_bad:
        for folder in folders:
            specimen_ID=os.path.basename(folder)[:9]
            pp_file=get_pp_path(folder)
            pp_dict=ju.read(pp_file)
            try:
                if pp_dict['th_inf']['via_Vmeasure']['value']< -.06:
                    th_inf_exclusion_list.append(specimen_ID)
            except:
                print folder, 'DOES NOT LOOK LIKE A NORMAL PREPROCESSOR FILE'
                strange_pp_exclusion.append(specimen_ID)   

        print '\t', len(set(th_inf_exclusion_list)), 'neurons have a th_inf less than -60 mV'                   

    # exclude based on resistance
    resistance_exclusion_list=[]
    if resist:
        for folder in folders:
            specimen_ID=os.path.basename(folder)[:9]
            pp_file=get_pp_path(folder)
            pp_dict=ju.read(pp_file)
            try:
                if pp_dict['resistance']['R_test_list']['mean']>1000.e6:
                    resistance_exclusion_list.append(specimen_ID)
            except:
                print folder, 'DOES NOT LOOK LIKE A NORMAL PREPROCESSOR FILE'
                strange_pp_exclusion.append(specimen_ID)    
        print '\t', len(set(resistance_exclusion_list)), 'neurons have a resistance fit WITHOUT ASC larger than 1000 MOhms.'

        print '\t', len(set(strange_pp_exclusion)), 'neurons have a strange looking preprocessor file.'

    # exclude based on explained variance on training data
    exp_var_exclusion_no_file=[]
    if ev:
        for folder in folders:
            specimen_ID=os.path.basename(folder)[:9]
            try:
                file=get_file_path_endswith(folder, 'GLIF1_exp_var_ratio_10ms.json')
            except:
                exp_var_exclusion_no_file.append(specimen_ID) 
        print '\t', len(set(exp_var_exclusion_no_file)),'neurons have no explained variance file which means they probably had a empty array in a noise 1.  See calc_all_explained_variance.py variable model_GLIF1_n1_after'

    exp_var_exclusion_below=[]
    
    # the following mouse neurons were either accidentally excluded from the analysis or removed for a reason that eludes me now.
    accidental_exclusions=[]
    if accidental_exclusion:
        if os.path.isdir(os.path.join(folder_path,'569739534'+'_Chrna2-Cre_OE25')): #if this directory exists get rid of it.
            accidental_exclusions=['569739534']
        else: pass #neuron already excluded or not in directory
        print '\t', len(set(accidental_exclusions)),'neurons were excluded from the analysis by accident. Set accidental_exclusion flag to False to use it if reprocessing all data.'
 
    
    def check_ev_value(folder,ew):
        '''Checks to see if the explained variance of the training data (noise 1) is below the specified value.
        inputs:
            folder:
                path to folder where files are located
            ew: string
                specifies the unique end of a file name of the file searching for
        returns:
            Nothing.  Appends specimen IDs to be excluded to the 'exp_var_exclusion_below' list
        '''
        specimen_ID=os.path.basename(folder)[:9]
        try:
            file=get_file_path_endswith(folder, ew)  #if file doesnt exist this will fail
            dictionary=ju.read(file)
            if dictionary['after_opt']['noise_1']<ev:
                exp_var_exclusion_below.append(specimen_ID)
        except: 
            print 'cant find a file for', specimen_ID, 'this should not happen if the check_sweeps_and_rm_folders.py was run!'
            pass
        
        
    if ev:
        for folder in folders:
            check_ev_value(folder, 'GLIF1_exp_var_ratio_10ms.json')                
        print '\t', len(set(exp_var_exclusion_below)), 'neurons have a GLIF explained variance on noise 1 training data of less than', ev
                    

    # get the set of all neurons that are still included in analysis after the above exclusions
    init_excluded_id_list=list(set(spike_cutting_exclusions+
                              resistance_exclusion_list+
                              th_inf_exclusion_list+
                              strange_pp_exclusion+
                              exp_var_exclusion_no_file+
                              exp_var_exclusion_below+
                              accidental_exclusions))
    reduced_sp_ids=list(set(initial_sp_ids)-set(init_excluded_id_list)) # specimen ids remaining after above exclusions
    
    # remove data that does not have at least a specified number (n_in_cre) of neurons in a cre line    
    if n_in_cre is not False:
        small_cre_line_exclusion=np.array([])
        cre_list=[]
        for folder in folders:
            specimen_ID=os.path.basename(folder)[:9]
            if specimen_ID in reduced_sp_ids:
                cre_list.append({'sp':specimen_ID, 'cre': os.path.basename(folder)[10:]})
                
        df=pd.DataFrame(cre_list)
        for cre in df['cre'].unique():
            if len(df[df['cre']==cre])<n_in_cre:
                small_cre_line_exclusion=np.append(small_cre_line_exclusion, (df[df['cre']==cre]['sp'].values))

    # create list of specimen IDs whose folder should be completely eliminated
    exclude_me_sp_ids=list(set(small_cre_line_exclusion.tolist()+
                               init_excluded_id_list))

    print 'A total of',len(exclude_me_sp_ids), 'out of', len(folders), 'neurons are excluded via general exclusion criteria leaving',len(folders)-len(exclude_me_sp_ids), 'for this analysis'
    
    return exclude_me_sp_ids

def count(path, ew):
    '''count the number of files in the structured data directory that 
    ends with the specified input
    Inputs:
        path: string
            path to the structured data directory
        ew: string
            string the file name ends with
    Returns:
        n: integer
            number of files in the structured data directory that end 
            with the specified input string
        
    '''
    folders=[os.path.join(path, f) for f in  os.listdir(structured_data_directory)]
    n=0
    for folder in folders:
        try:
            get_file_path_endswith(folder, ew)
            n=n+1
        except:
            pass 
    return n

if __name__ == '__main__':
    '''Remove files from the structured data directory based on exclusion criteria.
    '''
    
    # find and remove neurons that should be generally excluded from the directory
    exclude_specimen_ids=general_exclusions(structured_data_directory)
    folders=[os.path.join(structured_data_directory, f) for f in  os.listdir(structured_data_directory)]
    for folder in folders:
        sp_id=os.path.basename(folder)[0:9]
        if sp_id in exclude_specimen_ids:
            shutil.rmtree(folder)
            
    # find and remove bad spike component of threshold from the directory 
    bad_spike_comp_of_th=exclude_via_spike_comp_of_thr(structured_data_directory)
    folders=[os.path.join(structured_data_directory, f) for f in  os.listdir(structured_data_directory)]
    for folder in folders:
        sp_id=os.path.basename(folder)[0:9]
        if (sp_id in bad_spike_comp_of_th):
            for f in os.listdir(folder):
                if ('GLIF2' in f )  or ('GLIF4' in f) or ('GLIF5' in f):
                    os.remove(os.path.join(folder,f))
    
    # find and remove bad spike component of threshold from the directory                 
    bad_voltage_comp_of_th=exclude_via_v_comp_of_th(structured_data_directory)   
    for folder in folders:
        sp_id=os.path.basename(folder)[0:9]                 
        if sp_id in bad_voltage_comp_of_th:
            for f in os.listdir(folder):
                if ('GLIF5' in f):
                    os.remove(os.path.join(folder,f))

    # count up the files in the directory
    for folder in folders:
        if not get_file_path_endswith(folder,'GLIF1_neuron_config.json'):
            print 'nope'

    
    cre_list=[]
    folders=[os.path.join(structured_data_directory, f) for f in  os.listdir(structured_data_directory)]
    for folder in folders:
        specimen_ID=os.path.basename(folder)[:9]
        cre_list.append({'sp':specimen_ID, 'cre': os.path.basename(folder)[10:]})
    df=pd.DataFrame(cre_list)
    print df.groupby('cre').size()    
    
    # count up the files in the directory
    print 'TOTALS'
    print 'GLIF1 has', count(structured_data_directory,'_GLIF1_neuron_config.json')
    print 'GLIF2 has', count(structured_data_directory,'_GLIF2_neuron_config.json')
    print 'GLIF3 has', count(structured_data_directory,'_GLIF3_neuron_config.json')
    print 'GLIF4 has', count(structured_data_directory,'_GLIF4_neuron_config.json')
    print 'GLIF5 has', count(structured_data_directory,'_GLIF5_neuron_config.json')
    print 'DID YOU REMEMBER TO DOWNLOAD THE 580895033_Chrna2-Cre_OE25 NEURON? PLEASE SEE README.txt FOR MORE INFORMATION.'
