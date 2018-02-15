'''plots explained variance of different levels of example neurons
'''
import matplotlib.pyplot as plt
import numpy as np
import os
import allensdk.core.json_utilities as ju
import sys
relative_path=os.path.dirname(os.getcwd())
sys.path.append(os.path.join(relative_path, 'libraries'))


sigma=np.array([.0001, .001, .004, .01, .020, .1, 1.0])

#-----------plots of explained variance----------------
#set up figure
plt.figure(figsize=(16, 6))

#--474637203 #Ht3ra

ev_data=ju.read('json_data/474637203htr3_ev_data.json') #
for ii, ev in enumerate(ev_data):
    ax=plt.subplot(2,5, ii+1)
    ax.plot(sigma*1000., ev[0], 'k', lw=2)
    ax.plot(sigma*1000., ev[1], 'b', lw=2)
    ax.plot(sigma*1000., ev[2], 'r', lw=5)
    ax.set_xticklabels([])
        
    ax.set_ylim(.5, 1)
    ax.set_xlim(0, 20)
    if ii==0:
        ax.set_yticks(np.linspace(.5, 1, 6))
        ax.set_yticklabels([50, 60, 70, 80, 90, 100])
    else:
        ax.set_yticklabels([])

##--512322162 ctgf

ev_data=ju.read('json_data/512322162ctgf_ev_data.json')
for ii, ev in enumerate(ev_data):
    ax=plt.subplot(2,5, ii+6)
    ax.plot(sigma*1000., ev[0], 'k', lw=2)
    ax.plot(sigma*1000., ev[1], 'b', lw=2)
    ax.plot(sigma*1000., ev[2], 'r', lw=5)
    
    ax.set_ylim(.5, 1)
    ax.set_xlim(0, 20)
    
    if ii==0:
        ax.set_yticks(np.linspace(.5, 1, 6))
        ax.set_yticklabels([50, 60, 70, 80, 90, 100])
    else:
        ax.set_yticklabels([])

plt.annotate(r'$GLIF_1$', xy=(.2, .95), xycoords='figure fraction', horizontalalignment='center',fontsize=18)
plt.annotate(r'$GLIF_2$', xy=(.36, .95), xycoords='figure fraction', horizontalalignment='center',fontsize=18)
plt.annotate(r'$GLIF_3$', xy=(.52, .95), xycoords='figure fraction', horizontalalignment='center',fontsize=18)
plt.annotate(r'$GLIF_4$', xy=(.68, .95), xycoords='figure fraction', horizontalalignment='center',fontsize=18)
plt.annotate(r'$GLIF_5$', xy=(.84, .95), xycoords='figure fraction', horizontalalignment='center',fontsize=18)
plt.annotate('% Explained variance', xy=(.09, .5), xycoords='figure fraction', verticalalignment='center', horizontalalignment='center',fontsize=18, rotation='vertical')
plt.annotate('time (ms)', xy=(.5, .02), xycoords='figure fraction', horizontalalignment='center',fontsize=18)
plt.annotate('Htr3a', xy=(.05, .75), xycoords='figure fraction', horizontalalignment='left',fontsize=18, rotation='vertical')
plt.annotate('Ctgf', xy=(.05, .31), xycoords='figure fraction', horizontalalignment='left',fontsize=18, rotation='vertical')
    

plt.show()