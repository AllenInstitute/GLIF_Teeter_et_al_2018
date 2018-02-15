'''Written by Corinne Teeter 8-25-2016.  This library will have global definitions and functions for all the 
publication plots.
'''
import os
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 16})
from data_library import pos_cre_lines, cononical, non_nan_data, float_data
                    
cre_neg_color='lime'
color_dict={'all_neurons':'k',
            
            'excitatory': 'b',                          #RGB [0,0,255]         
            'Cux2-CreERT2':'Teal',                   #RGB [0,128,128]
            'Nr5a1-Cre':'CadetBlue',                 #RGB [95,158,160]
            'Scnn1a-Tg2-Cre':'DarkSlateBlue',        #RGB [72,61,139]
            'Rorb-IRES2-Cre-D':'MidnightBlue',         #RGB [25,25,112]
            'Scnn1a-Tg3-Cre':'Cyan',                 #RGB [0,255,255]
            'Rbp4-Cre_KL100':'DarkBlue',             #RGB [0,0,139]
            'Ctgf-2A-dgCre':'MediumSlateBlue',       #RGB [123,104,238]
            'Ntsr1-Cre_GN220':'DodgerBlue',          #RGB [30,144,255]

            'inhibitory': 'r',                          #RGB [255, 0,0]
            'Sst-IRES-Cre':'FireBrick',              #RGB [178,34,34]
            'Pvalb-IRES-Cre':'DeepPink',             #RGB [255,20,147]       
            'Htr3a-Cre_NO152':'IndianRed',           #RGB [205,92,92]
            'Ndnf-IRES2-dgCre':'DarkRed',            #RGB [139,0,0]
            'Chat-IRES-Cre-neo':'Crimson',               #RGB [220,20,60]
            'Vip-IRES-Cre':'LightCoral',             #RGB [240,128,128]
            'Chrna2-Cre_OE25': 'OrangeRed',          #RGB [255,69,0]
            'Nkx2-1-CreERT2':'HotPink',               #RGB [255,105,180]            
            } 

def draw_unity():
    axes=plt.gca()
    ends=axes.get_xlim()
    plt.plot(ends, ends, 'k')

def distribution_plot(cre_dict, x_index, y_index, xlabel='', ylabel=''):
    '''creates a scatter plot and data distributions of CRE lines containing 10 or more points.
    input:
        cre_dict: dictionary containing list of lists
            keys corresponds to of the excitatory, inhibitory, individual cre 
            lines and cre negative neurons. Each key contains a list of lists of neuron data
        x_index: integer
            index (column index) of data in lists of lists to be plotted on the x-axis
        y_index: integer
            index (column index) of data in lists of lists to be plotted on the y-axis
        xlabel: string
            label for xaxis    
        ylabel: string
            label for yaxis
    returns:
        percentile_dict: dictionary
            contains percentiles for the distributions of the keys in cre_dict
    '''
    plt.figure()
    ax_main = plt.subplot2grid((4,4), (1, 0), colspan=3, rowspan=3)
    ax_x = plt.subplot2grid((4,4), (0, 0), colspan=3)
    ax_y = plt.subplot2grid((4,4), (1,3), rowspan=3)
    ax_y.get_yaxis().set_ticks([])
    ax_y.get_xaxis().set_ticks([])
    
    usable_neurons=cre_dict['excitatory']+cre_dict['inhibitory']
    
    # make scatter plot

    for key in pos_cre_lines:
        for neuron in cre_dict[key]:
            ax_main.plot(neuron[x_index], neuron[y_index], '.',color=color_dict[key], ms=12)
    
    #note this is a hack and assumes that specimen id is always the first number in the list        
    for con in cononical:
        #find the values
        neuron =[neuron for neuron in cre_dict[con[1]] if neuron[0]==con[0]][0]
        ax_main.plot(neuron[x_index], neuron[y_index], '*',color=color_dict[con[1]], markeredgewidth='3', ms=24)
        
    
    x_lims=ax_main.get_xlim()
    y_lims=ax_main.get_ylim()
    ax_main.set_ylabel(ylabel)
    ax_main.set_xlabel(xlabel)
    ax_x.get_yaxis().set_ticks([])
    ax_x.get_xaxis().set_ticks([])
    
    #--------------------------------------------------------------------------
    #----------------------------distributions---------------------------------
    #--------------------------------------------------------------------------

    #----------------------xaxis
    #--total
    x=[neuron[x_index] for neuron in usable_neurons]
    v, edges=np.histogram(x)
    bincenters = 0.5*(edges[1:]+edges[:-1])
    ax_x.plot(bincenters, v/np.float(len(x)), 'k', lw=5)
    ax_x.set_xlim(x_lims)
    
    percentile_dict={}
    #-- excitatory
    percentile_dict['excitatory']={'percentile':{xlabel:{'25':np.NAN, '50':np.NAN, '75':np.NAN}}}
    x=[neuron[x_index] for neuron in cre_dict['excitatory']]
    v, edges=np.histogram(x)
    bincenters = 0.5*(edges[1:]+edges[:-1])
    ax_x.plot(bincenters, v/np.float(len(x)), 'b', lw=5)
    percentile_dict['excitatory']['percentile'][xlabel]['25']=np.percentile(x, 25)
    percentile_dict['excitatory']['percentile'][xlabel]['50']=np.percentile(x, 50)
    percentile_dict['excitatory']['percentile'][xlabel]['75']=np.percentile(x, 75)
    
    #-- inhibitory
    percentile_dict['inhibitory']={'percentile':{xlabel:{'25':np.NAN, '50':np.NAN, '75':np.NAN}}}
    x=[neuron[x_index] for neuron in cre_dict['inhibitory']]
    v, edges=np.histogram(x)
    bincenters = 0.5*(edges[1:]+edges[:-1])
    ax_x.plot(bincenters, v/np.float(len(x)), 'r', lw=5)
    percentile_dict['inhibitory']['percentile'][xlabel]['25']=np.percentile(x, 25)
    percentile_dict['inhibitory']['percentile'][xlabel]['50']=np.percentile(x, 50)
    percentile_dict['inhibitory']['percentile'][xlabel]['75']=np.percentile(x, 75)
    
    #--plot individual distributions
    for cre in pos_cre_lines:
        percentile_dict[cre]={'percentile':{xlabel:{'25':np.NAN, '50':np.NAN, '75':np.NAN}}}
        x=[neuron[x_index] for neuron in cre_dict[cre]]
        if len(x)>=10:
            v, edges=np.histogram(x)
            bincenters = 0.5*(edges[1:]+edges[:-1])
            ax_x.plot(bincenters, v/np.float(len(x)), color=color_dict[cre], lw=2)
            percentile_dict[cre]['percentile'][xlabel]['25']=np.percentile(x, 25)
            percentile_dict[cre]['percentile'][xlabel]['50']=np.percentile(x, 50)
            percentile_dict[cre]['percentile'][xlabel]['75']=np.percentile(x, 75)
        else:
            print cre, 'has less than 10 points'
    
    #--------------------yaxis------------------------------------------------        
    #--total
    y=[neuron[y_index] for neuron in usable_neurons]
    v, edges=np.histogram(y)
    bincenters = 0.5*(edges[1:]+edges[:-1])
    ax_y.plot(v/np.float(len(y)), bincenters, 'k', lw=5)
    ax_y.set_ylim(y_lims)
    
    #-- excitatory
    percentile_dict['excitatory']['percentile'][ylabel]={'25':np.NAN, '50':np.NAN, '75':np.NAN}
    x=[neuron[y_index] for neuron in cre_dict['excitatory']]
    v, edges=np.histogram(x)
    bincenters = 0.5*(edges[1:]+edges[:-1])
    ax_y.plot(v/np.float(len(x)), bincenters, 'b', lw=5)
    percentile_dict['excitatory']['percentile'][ylabel]['25']=np.percentile(x, 25)
    percentile_dict['excitatory']['percentile'][ylabel]['50']=np.percentile(x, 50)
    percentile_dict['excitatory']['percentile'][ylabel]['75']=np.percentile(x, 75)
    
    #-- inhibitory
    percentile_dict['inhibitory']['percentile'][ylabel]={'25':np.NAN, '50':np.NAN, '75':np.NAN}
    x=[neuron[y_index] for neuron in cre_dict['inhibitory']]
    v, edges=np.histogram(x)
    bincenters = 0.5*(edges[1:]+edges[:-1])
    ax_y.plot(v/np.float(len(x)), bincenters, 'r', lw=5)
    percentile_dict['inhibitory']['percentile'][ylabel]['25']=np.percentile(x, 25)
    percentile_dict['inhibitory']['percentile'][ylabel]['50']=np.percentile(x, 50)
    percentile_dict['inhibitory']['percentile'][ylabel]['75']=np.percentile(x, 75)
    
    #--individual cre lines
    for cre in pos_cre_lines:
        percentile_dict[cre]['percentile'][ylabel]={'25':np.NAN, '50':np.NAN, '75':np.NAN}
        x=[neuron[y_index] for neuron in cre_dict[cre]]
        if len(x)>=10:
            v, edges=np.histogram(x)
            bincenters = 0.5*(edges[1:]+edges[:-1])
            ax_y.plot(v/np.float(len(x)),bincenters, color=color_dict[cre], lw=2)
            percentile_dict[cre]['percentile'][ylabel]['25']=np.percentile(x, 25)
            percentile_dict[cre]['percentile'][ylabel]['50']=np.percentile(x, 50)
            percentile_dict[cre]['percentile'][ylabel]['75']=np.percentile(x, 75)
        else:
            print cre, 'has less than 10 points'
    
    return percentile_dict

def distribution_analysis(data, index, percentile_boundries=[25, 75]): 
    '''Return the specimen ids that are within the percentile boundries of the distribution. 
    Also returns mean and quartiles and data.  Note that this is a super 
    sloppy function as it has integer division that doesnt really get
    a specific number of indicies, n is just a ball park number you want. Also
    specimen id index is hard coded so this function greatly depends on the list order 
    being correct 
    input:
        data: list of lists
            each list is data of a neuron
        index: integer
            index of the column of the inner lists that want
        percentile_boundries: list with two numbers 
            first number is lower percentile used to find specimen IDs, 
            second number is upper percentile used to find specimen IDs.
    '''
    #TODO: FIGURE OUT WHY: non_nan_data doesnt always work sometimes needs to be changed to float_data like for running feature data
    reduced_data=non_nan_data(data, index)  #make sure no nans are in the data column interested in
    data_column=[neuron[index] for neuron in reduced_data]
    if len(reduced_data)>=5:
        median=np.median(data_column)
        quartiles=[np.percentile(data_column, 25), np.percentile(data_column, 75)]
        tightness=np.abs(quartiles[1]-quartiles[0])
        #---get specimen IDs of data within values
        sorted_ind=np.argsort(data_column)
        #convert precentiles in indexes
        min_index=np.ceil((percentile_boundries[0]/100.)*len(sorted_ind))  #using cieling here because we want the index to be within the boundries
        max_index=np.floor((percentile_boundries[1]/100.)*len(sorted_ind))  #using floor here because we want the index to be within the boundries
        chosen_ind=sorted_ind[int(min_index):int(max_index)]
        nearest_specimen_ids=np.array([neuron[5] for neuron in reduced_data])[chosen_ind]
    else:
        nearest_specimen_ids=[]
        median=np.nan 
        quartiles=np.nan
        tightness=np.nan
    
    return nearest_specimen_ids, median, quartiles, data_column, tightness
