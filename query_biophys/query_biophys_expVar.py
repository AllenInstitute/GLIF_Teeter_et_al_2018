'''Written by Corinne Teeter. Grab the explained variance for the published biophys models'''

import numpy as np
from allensdk.core.cell_types_cache import CellTypesCache
import allensdk.internal.core.lims_utilities as lu
import pandas as pd
from allensdk.api.queries.glif_api import GlifApi
import os
import sys
relative_path=os.path.dirname(os.getcwd())
sys.path.append(os.path.join(relative_path, 'libraries'))

# Find all mouse cells with models
glif_api = GlifApi()
ctc = CellTypesCache(manifest_file=os.path.join(relative_path,'cell_types_manifest.json'))
specimen_id_list=[]
temp = ctc.get_cells()
for c in temp:
    if c['species']=='Mus musculus':
        specimen_id_list.append(c['id'])

print len(specimen_id_list), 'mouse specimens in public database'

def get_expVar(specimen_id_list, keyword):
    '''Grab explained variance value of specimen id list in public database
    Inputs:
        specimen_id_list: list of integers
            desired specimen ids of data in AIBS public database
        keyword: string
            string to search for in the 'name' entry of the  AIBS public database 
    Outputs:
        ev: list of floats
            explained variance for input specimen ids
        
    '''
    ev=[]
    for sp in specimen_id_list:
        models=glif_api.get_neuronal_models(sp)[0]
        for m in models['neuronal_models']:
        #    print m['name']
            if keyword in m['name']:
                if len(m['neuronal_model_runs'])>1:
                    print ('there is more than 1 model run')
                else:
                    ev.append(m['neuronal_model_runs'][0]['explained_variance_ratio'])
    return ev

def find_model_spid(specimen_id_list, keyword):
    '''Given a list of specimen ids, returns list that have a specified model
    inputs:
        specimen_id_list: list of integers
            desired specimen ids of data in AIBS public database 
        keyword: string
            string to search for in the 'name' entry of the  AIBS public database 
    outputs:
        ev: list of integers
            specimen ids out of list that contain keywork in the 'name'
        
    '''
    reduced_ids=[]
    for sp in specimen_id_list:
        models=glif_api.get_neuronal_models(sp)[0]
        for m in models['neuronal_models']:
        #    print m['name']
            if keyword in m['name']:
                reduced_ids.append(m['specimen_id'])
            
    return list(set(reduced_ids))


# get biophys perisomatic models 
biophys_ps_sp_ids=find_model_spid(specimen_id_list, 'Biophysical - perisomatic')               
print len(biophys_ps_sp_ids), 'specimens with Biophysical - perisomatic'
biophys_ps_expVar=get_expVar(biophys_ps_sp_ids, 'Biophysical - perisomatic')
biophys_ps_eV_values=np.array(biophys_ps_expVar)[np.array([b!=None for b in biophys_ps_expVar])] #get rid of 'None'
print 'total Biophysical - perisomatic with explained variance value:', len(biophys_ps_eV_values), ', mean:',np.mean(biophys_ps_eV_values), ', median:', np.median(biophys_ps_eV_values)

# get biophys all active
biophys_aa_sp_ids=find_model_spid(specimen_id_list, 'Biophysical - all active')               
print len(biophys_aa_sp_ids), 'specimens with Biophysical - all active'
biophys_aa_expVar=get_expVar(biophys_aa_sp_ids, 'Biophysical - all active')
biophys_aa_eV_values=np.array(biophys_aa_expVar)[np.array([b!=None for b in biophys_aa_expVar])] #get rid of 'None'
print 'total Biophysical - all active with explained variance value:', len(biophys_aa_eV_values), ', mean:',np.mean(biophys_aa_eV_values), ', median:', np.median(biophys_aa_eV_values)

