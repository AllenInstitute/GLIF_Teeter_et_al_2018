'''Not perfectly implemented.  spikecut_cal_err.py seemed to run fast enough 
without cluster implimentation.  However if it is needed in the future this 
and that code can be easily implimented
'''

import os
import sys
import numpy as np

QSUB_HEADERS = ['#!/bin/bash',
                '#PBS -l walltime=1:00:00',
                '#PBS -m a',
                '#PBS -r y']

QSUB_BODY="""
ANALYSIS_PATH=/home/corinnet/workspace/GLIF_analysis/internal
PYTHON_BINARY=/home/corinnet/anaconda/bin/python2.7
SDK_PATH=/home/corinnet/workspace/allensdk/

export PYTHONPATH=$GLIF_PATH:$SDK_PATH:
"""

import subprocess

def inititiate_list_of_qsub_on_cluster(qsub_file_list, priority='med'):

    for qsub_file in qsub_file_list:
        if priority =='high':
            qsub_args = ['qsub', '-W x=QOS:high', qsub_file]
        elif priority == 'med':
            qsub_args = ['qsub', '-W x=QOS:med', qsub_file]
        elif priority =='low':        
            qsub_args = ['qsub', '-W x=QOS:low', qsub_file]
        else:
            raise Exception('cluster priority is not defined')
        
        call_simple(qsub_args)

def call(cmd, pipe_input=None):

    # Open a pipe to the qsub command.
    p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, close_fds=True)
    (input, output) = (p.stdin, p.stdout)

    # Send job_string to qsub
    if pipe_input != None:
        input.write(pipe_input)
    input.close()
    return output.read()

def call_simple(cmd):
    subprocess.call(cmd)

def generate_qsub_single_files(folder, cluster_name='uno'):
    '''generate the qsub files for running on the cluster
    run_list: list of configuration files to be run via the qsub file
    '''
    specimen_id=os.path.basename(folder)[:9]
    cre=os.path.basename(folder)[10:]
    name=specimen_id+'_'+cre
    
    output_path='/allen/aibs/mat/Corinne/spikecut_cluster_output'
    try:
        os.makedirs(output_path)
    except: pass
    
    cmd_lines = [ ( '$PYTHON_BINARY $ANALYSIS_PATH/spike_width_calc_err.py '+folder)]  #this will need to be more complicated
    
    #--specify names of files of cluster run
    qsub_file_name = os.path.join(output_path, name+'_expVar.qsub')
    err_file_name = os.path.join(output_path, name+'_expVar.err')
    out_file_name = os.path.join(output_path, name+'_expVar.out')
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

    
def run(folder, run_location='local'):
    '''
    '''
        
        
    if run_location=='local':
        import spikecut_calc_err
        spikecut_calc_err.main(folder)
        
    elif run_location =='cluster':
        inititiate_list_of_qsub_on_cluster([generate_qsub_single_files(folder)], priority='high')
    else:
        raise Exception('Where should this be run?')

        
if __name__ == '__main__':

    
    # load data out of configuration files
    folder_path='/home/corinnet/workspace/GLIF_analysis/data'
    folders=np.sort([os.path.join(folder_path, f) for f in  os.listdir(folder_path)])

    for folder in ['/home/corinnet/workspace/GLIF_analysis/data/313861828_Scnn1a-Tg3-Cre']: #in folders:
        print 'folder', folder
        run(folder, run_location='local')
    
    
