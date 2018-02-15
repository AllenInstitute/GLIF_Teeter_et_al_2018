'''
Grab GLIF configuration (one for each GLIF model, *_neuron_config.json) and preprocessor 
files (one for each neuron *_preprocessor_values.json) from the Allen Institute Cell Types
database and arranges them for use in the analysis code.  A data folder (mouse_struc_data_dir 
or human_struc_data_dir) should appear within this directory.  Note that that cleaned up 
versions of these directories are provided in this repository and are available in the root 
directory.'''

import os
import sys
relative_path=os.path.dirname(os.getcwd())
from allensdk.api.queries.glif_api import GlifApi
from allensdk.core.cell_types_cache import CellTypesCache
ctc = CellTypesCache(manifest_file=os.path.join(relative_path,'cell_types_manifest.json'))
import pandas as pd
from shutil import copyfile
import numpy as np
from allensdk.config import enable_console_log 
enable_console_log()

#---------------------------------------------------------------------------
#----------------------SPECIFY SPECIES--------------------------------------
#---------------------------------------------------------------------------
species="mouse"
#species="human"
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------

model_template_ids=[395310469, 395310479, 395310475, 471355161, 395310498] #model template ids associated with data in database
model_names=['GLIF1', 'GLIF2', 'GLIF3', 'GLIF4', 'GLIF5']

def get_files_from_LIMS_public(output_path, glif_sp_ids=None, type='mouse'):
    '''This will grab cre positive data config files from LIMS and sort them and put them in 
    the specified output folder.  
    input:
        output_path: string
            specifies path for files to be placed in
        glif_sp_ids: list of strings or integers
            specimen ids of cells specifically want to grab.  If none it will get all available on the 
            Allen Institue Cell Types Database.
        type: string
            can be 'mouse' or 'human'. Note that if mouse is specified is will only grab cre positive mouse cells 
            (code can be altered to get cre negative cells).
    output:
        Does not return values but creates the specified 'output_path' folder.  
        Inside the folder a series of folders are created with the name format:
        specimenid_cre.  Inside those inner folders are the neuron configs of 
        the available GLIF models along with the preprocessor files.
    '''

    glif_api = GlifApi()     
    ctc = CellTypesCache(manifest_file=os.path.join(relative_path,'cell_types_manifest.json'))

    # select the specimen ids to grab from the data base (cre positive or human which have at least 1 GLIF model)
    if glif_sp_ids==None: #if no specimen id's are specified grab all data in the cell types manifest
        specimen_id_list = []
        if type=='mouse':
            for c in ctc.get_cells():
                if c['reporter_status']=='cre reporter positive':
                    specimen_id_list.append(c['id'])
        elif type=='human':
            print 'getting human'
            for c in ctc.get_cells(species=['Homo Sapiens']):
                #print c
                specimen_id_list.append(c['id'])
            print specimen_id_list
        # reduce list to cells that have a GLIF model
        glif_sp_ids=[]
        for sp in specimen_id_list:
            models=glif_api.get_neuronal_models(sp)[0]
            for m in models['neuronal_models']:
                if 'LIF' in m['name']:
                    glif_sp_ids.append(m['specimen_id'])
                    
        glif_sp_ids=list(set(glif_sp_ids))
        print len(glif_sp_ids), 'cre positive specimens with at least 1 LIF model'

    # create the overall output directory if it doesn't exist
    try:
        os.makedirs(output_path)
    except:
        pass
    
    # go get the files corresponding to the specimen ids from the Allen Cell Types Database 
    # and put them into a specified output directory 
    for id in glif_sp_ids:
        model_query=glif_api.get_neuronal_models(id)[0]['neuronal_models']
        df=pd.DataFrame(model_query)
        for mt_id, short_name in zip(model_template_ids, model_names):
            dff=df[df['neuronal_model_template_id']==mt_id]
            if len(dff)>=2:
                print dff
                raise Exception("This is public data, there should not be more than 1 model")
            elif len(dff)==1:
                use_me=dff
                #go get the file 
                path=use_me['well_known_files'].iloc[0][0]['path'] 
                if type=='mouse':
                    cre=(str(use_me['name'].values).split(')_'))[1].split(';')[0]
                elif type=='human':
                    cre='human'
                else:
                    raise Exception('specified species not known')
                # convert old non complete cre names
                if 'Ntsr1-Cre' in cre:
                    cre='Ntsr1-Cre_GN220'
                if 'Chat-IRES-Cre' in cre:
                    cre='Chat-IRES-Cre-neo'
                dir_name=os.path.join(output_path, str(id)+'_'+cre)
                try:    
                    os.makedirs(dir_name)
                except:
                    pass
                if path.endswith('_neuron_config.json'):
                    pass
                else:
                    print path
                    raise Exception('the file doesnt end with _neuron_config.json')       
                try:   
                    copyfile(path, os.path.join(dir_name, str(id)+'_'+cre+'_'+short_name+'_neuron_config.json'))
                except:
                    print 'couldnt make ', os.path.join(dir_name, str(id)+'_'+cre+'_'+short_name+'_neuron_config.json')
                if mt_id==model_template_ids[0]:
                    model_path=os.path.dirname(path)
                    pp_path=os.path.join(model_path,
                        os.listdir(model_path)[np.where([fname.endswith('_preprocessor_values.json') for fname in os.listdir(model_path)])[0][0]]) 
                    try:   
                        copyfile(pp_path, os.path.join(dir_name, str(id)+'_'+cre+'_preprocessor_values.json'))
                    except:
                        print 'couldnt make ', os.path.join(dir_name, str(id)+'_'+cre+'_preprocessor_values.json')
                        raise Exception('there should be a preprocessed file')
            elif len(dff)<1:
                use_me=pd.DataFrame()
                path=None
    


if __name__ == '__main__':
    
    if species=='mouse':
        get_files_from_LIMS_public('mouse_struc_data_dir', type='mouse') #this will create a directory called 'mouse_struc_data_dir' and put the mouse files inside
    elif species=='human':
        get_files_from_LIMS_public('human_struc_data_dir', type='human') #this will create a directory called 'human_struc_data_dir' and put the human files inside
    else:
        raise Exception('the desired species is not recognized')