'''Runs calc_all_explained_variance.py either locally or via a cluster to 
calculate explained variance of noise 1 and noise 2 before and after 
optimization. These values are written into the data folders as json files.
Note that if you are running on your own cluster, all code that involves
cluster computing will need to be altered (including paths).  Specify species 
desired in the __main__ code at bottom of this file. IT IS HIGHLY RECOMMENDED 
THAT "compare_calc_and_db_EV.py "IN THE "sanity_checks" DIRECTORY IS RUN AFTER
 THE COMPLETION OF THIS CODE.
'''

import os
import sys
relative_path=os.path.dirname(os.getcwd())
sys.path = [os.path.join(relative_path, 'libraries')] + sys.path
import numpy as np
import subprocess


#-------------------------------------------------
#---------SPECIFY THESE --------------------------
#-------------------------------------------------
'''THE INPUT TO THIS SCRIPT SHOULD BE TWO INTEGER VALUES THAT SPECIFY 
THE START AND END INDICIES OF THE FILE ONE WOULD LIKE TO PROCESS.  THE 
NUMBER OF FILES DOES NOT NEED TO BE KNOWN AS THE CODE WILL STOP WHEN
IT RUNS OUT OF NEURONS TO PROCESS.

Example command line run:
>> python calc_all_exp_var_RUN.py 0 1000
'''

# specify species 'mouse or human'
species='mouse'
#species='human'
#-------------------------------------------------
#-------------------------------------------------

if species =='mouse':
    structured_data_directory='mouse_struc_data_dir'  #Note that this is accessing data local to this folder not the global directory provided in this repository 
    nwb_directory=os.path.join(relative_path,'mouse_nwb')
elif species=='human':
    structured_data_directory='human_struc_data_dir' #Note that this is accessing data local to this folder not the global directory provided in this repository 
    nwb_directory=os.path.join(relative_path,'human_nwb')
else:
    raise Exception('species not recognized')


# files for cluster submission from within the Allen Institute
QSUB_HEADERS = ['#!/bin/bash',
                '#PBS -l walltime=3:00:00',
                '#PBS -m a',
                '#PBS -r y']

QSUB_BODY="""
ANALYSIS_PATH=/allen/aibs/mat/Corinne/GLIF_paper_analysis/create_data_dir
LIBRARY_PATH=/allen/aibs/mat/Corinne/GLIF_paper_analysis/libraries
PYTHON_BINARY=/home/corinnet/anaconda/bin/python2.7
SDK_PATH=/allen/aibs/mat/Corinne/allensdk

export PYTHONPATH=$GLIF_PATH:$SDK_PATH:$LIBRARY_PATH
"""

def inititiate_list_of_qsub_on_cluster(qsub_file_list, priority='med'):
    '''specifies just priority and submits job to the cluster
    qsub_file_list: list of strings 
        strings specify path of qsub files to be run
     priority: string
         sets priority on cluster
         defalt is 'med', options are 'high', 'med', and 'low'
    '''

    for qsub_file in qsub_file_list:
        if priority =='high':
            qsub_args = ['qsub', '-W x=QOS:high', qsub_file]
        elif priority == 'med':
            qsub_args = ['qsub', '-W x=QOS:med', qsub_file]
        elif priority =='low':        
            qsub_args = ['qsub', '-W x=QOS:low', qsub_file]
        else:
            raise Exception('cluster priority is not defined')
        
        subprocess.call(qsub_args)


def generate_qsub_single_files(specimen_id_directory, nwb_directory, cluster_name='uno'):
    '''generate and saves the qsub files for running on the cluster
    input:
        specimen_id_directory: string
            path to specimen id folder in the mouse_struc_data_dir (or human_struc_data_dir) folder
        nwb_directory: string
            path where experiment data in nwb format is saved
        cluster_name: string
            name of cluster to submit job to
    returns:
        qsub_file_name: string
            full path to saved qsub file
        
    '''
    specimen_id=os.path.basename(specimen_id_directory)[:9]
    cre=os.path.basename(specimen_id_directory)[10:]
    name=specimen_id+'_'+cre
    
    cluster_file_path='/allen/aibs/mat/Corinne/expVar_cluster_output' #specifies a directory to write the files used and generated via the cluster
    try:
        os.makedirs(cluster_file_path)
    except: pass
    
    cmd_lines = [ ( '$PYTHON_BINARY $ANALYSIS_PATH/calc_all_explained_variance.py '+specimen_id_directory+' '+nwb_directory)]  #this will need to be more complicated
    
    #--specify names of files of cluster run
    qsub_file_name = os.path.join(cluster_file_path, name+'_expVar.qsub')
    err_file_name = os.path.join(cluster_file_path, name+'_expVar.err')
    out_file_name = os.path.join(cluster_file_path, name+'_expVar.out')
    job_name = specimen_id + 'ev'

    with open(qsub_file_name, 'wb') as f:

        qsub_headers = QSUB_HEADERS + [ '#PBS -q %s' % cluster_name,
                                        '#PBS -l ncpus=%d' % len(cmd_lines), 
                                        '#PBS -e ' + err_file_name, 
                                        '#PBS -o ' + out_file_name, 
                                        '#PBS -N ' + job_name ]

        job_cmds = cmd_lines 
        lines = qsub_headers + [QSUB_BODY] + job_cmds
        f.write('\n'.join(lines) + '\n')

    return qsub_file_name

    
def run(specimen_id_directory, nwb_directory, run_location='local'):
    '''
    Runs explained variance calculation either locally or on a cluster.
    input:
        specimen_id_directory: string
            path to specimen id folder in the mouse_struc_data_dir (or human_struc_data_dir) folder
        nwb_directory: string
            path where experiment data in nwb format is saved
        run_location: string
            options: 'local' or 'cluster'
            runs code either locally or on a cluster
    ''' 
    if run_location=='local':
        import calc_all_explained_variance
        calc_all_explained_variance.main(specimen_id_directory, nwb_directory)
        
    elif run_location =='cluster':
        inititiate_list_of_qsub_on_cluster([generate_qsub_single_files(specimen_id_directory, nwb_directory)], priority='med')
    else:
        raise Exception('Where should this be run?')

        
if __name__ == '__main__':

    # specify subsections of data to be processed.  This is done because some times it is easier to just run
    # the code in blocks on a local machine.
    start=int(sys.argv[1]) #folders are processed in numerical order.  Specify data number integer at which to start (order 0-through total number of data files; not the specimen id number) 
    finish=int(sys.argv[2]) #folders are processed in numerical order.  Specify data number integer at which to end (order 0-through total number of data files; not the specimen id number) 
    
    # sort the data so that specifying start and end integers works
    folders=np.sort([os.path.join(structured_data_directory, f) for f in  os.listdir(structured_data_directory)])
    
    # initiate code run either locally or on a cluster for all specified files
    for specimen_id_directory, ii in zip(folders[start:finish], range(start,finish)):
        print 'running', ii, 'of', finish-1  #prints the sequential file number being run
        run(specimen_id_directory, nwb_directory, run_location='local') #run the calculation
    
    
        
