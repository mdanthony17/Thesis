#!/usr/bin/python
import sys, array, os
sys.path.insert(0, '/Users/Matt/Desktop/Xenon/xenon1t/fit_nr_band')

import matplotlib
matplotlib.use('QT4Agg')

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 16

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.patches as mpatches

import config_xe1t
import numpy as np
import time, tqdm
import cPickle as pickle

import pandas as pd
import scipy.interpolate

alpha = 0.3
l_file_types = ['ambe', 'rn']
d_acceptances = {}

for file_type in l_file_types:
    if file_type == 'ambe':
        d_acceptances[file_type] = {}
        d_acceptances[file_type]['file_path'] = '/Users/Matt/Desktop/Xenon/xenon1t/fit_nr_band/resources/bias_smearing_ambe.p'
        d_acceptances[file_type]['color_bias'] = '#0082c8'
        d_acceptances[file_type]['color_smearing'] = '#911eb4'
        d_acceptances[file_type]['df'] = pd.DataFrame(pickle.load(open(d_acceptances[file_type]['file_path'], 'rb')))
        #d_acceptances[file_type]['df'] = d_acceptances[file_type]['df'][(d_acceptances[file_type]['df']['s1s'] <= config_xe1t.l_s1_settings[2])]
        d_acceptances[file_type]['label'] = 'NR Calibration'
    elif file_type == 'rn':
        d_acceptances[file_type] = {}
        d_acceptances[file_type]['file_path'] = '/Users/Matt/Desktop/Xenon/xenon1t/fit_nr_band/resources/bias_smearing_wimps.p'
        d_acceptances[file_type]['color_bias'] = '#f58231'
        d_acceptances[file_type]['color_smearing'] = '#e6194b'
        d_acceptances[file_type]['df'] = pd.DataFrame(pickle.load(open(d_acceptances[file_type]['file_path'], 'rb')))
        #d_acceptances[file_type]['df'] = d_acceptances[file_type]['df'][(d_acceptances[file_type]['df']['s1s'] <= config_xe1t.l_s1_settings[2])]
        d_acceptances[file_type]['label'] = 'ER Calibration and First Science Run'



s_path_to_plots = '../images/'

l_handlers_s1_bias = []
l_handlers_s1_smearing = []
l_handlers_s2_bias = []
l_handlers_s2_smearing = []

fig_s1_bias, ax_s1_bias = plt.subplots(1)
#fig_s1_smearing, ax_s1_smearing = plt.subplots(1)
fig_s2_bias, ax_s2_bias = plt.subplots(1)
#fig_s2_smearing, ax_s2_smearing = plt.subplots(1)


ax_s1_smearing = ax_s1_bias.twinx()
ax_s2_smearing = ax_s2_bias.twinx()

for file_type in l_file_types:
    d_acceptances[file_type]['s1'] = {}
    d_acceptances[file_type]['s1']['x_values'] = np.asarray(d_acceptances[file_type]['df']['s1s'], dtype=np.float32)
    mask_s1 = d_acceptances[file_type]['s1']['x_values'] <= config_xe1t.l_s1_settings[2]
    d_acceptances[file_type]['s1']['x_values'] = d_acceptances[file_type]['s1']['x_values'][mask_s1]

    d_acceptances[file_type]['s1']['bias_values_lower'] = np.asarray(d_acceptances[file_type]['df']['s1bias_mins'], dtype=np.float32)[mask_s1]
    d_acceptances[file_type]['s1']['bias_values_upper'] = np.asarray(d_acceptances[file_type]['df']['s1bias_maxs'], dtype=np.float32)[mask_s1]
    d_acceptances[file_type]['s1']['smear_values_lower'] = np.asarray(d_acceptances[file_type]['df']['s1smearing_mins'], dtype=np.float32)[mask_s1]
    d_acceptances[file_type]['s1']['smear_values_upper'] = np.asarray(d_acceptances[file_type]['df']['s1smearing_maxs'], dtype=np.float32)[mask_s1]
    
    #print d_acceptances[file_type]['s1']['bias_values_lower']
    #print d_acceptances[file_type]['s1']['bias_values_upper']

    d_acceptances[file_type]['handler_bias'] = ax_s1_bias.fill_between(d_acceptances[file_type]['s1']['x_values'], d_acceptances[file_type]['s1']['bias_values_lower'], d_acceptances[file_type]['s1']['bias_values_upper'], color=d_acceptances[file_type]['color_bias'], alpha=alpha, label='Bias: %s' % d_acceptances[file_type]['label'])
    d_acceptances[file_type]['handler_smearing'] = ax_s1_smearing.fill_between(d_acceptances[file_type]['s1']['x_values'], d_acceptances[file_type]['s1']['smear_values_lower'], d_acceptances[file_type]['s1']['smear_values_upper'], color=d_acceptances[file_type]['color_smearing'], alpha=alpha, label='Smearing: %s' % d_acceptances[file_type]['label'])
    
    l_handlers_s1_bias.append(d_acceptances[file_type]['handler_bias'])
    l_handlers_s1_bias.append(d_acceptances[file_type]['handler_smearing'])
    #l_handlers_s1_smearing.append(d_acceptances[file_type]['handler_smearing'])






    # now S2
    d_acceptances[file_type]['s2'] = {}
    d_acceptances[file_type]['s2']['x_values'] = np.asarray(d_acceptances[file_type]['df']['s2s'], dtype=np.float32)
    mask_s2 = d_acceptances[file_type]['s2']['x_values'] <= 5000
    d_acceptances[file_type]['s2']['x_values'] = d_acceptances[file_type]['s2']['x_values'][mask_s2]
    
    d_acceptances[file_type]['s2']['bias_values_lower'] = np.asarray(d_acceptances[file_type]['df']['s2bias_mins'], dtype=np.float32)[mask_s2]
    d_acceptances[file_type]['s2']['bias_values_upper'] = np.asarray(d_acceptances[file_type]['df']['s2bias_maxs'], dtype=np.float32)[mask_s2]
    d_acceptances[file_type]['s2']['smear_values_lower'] = np.asarray(d_acceptances[file_type]['df']['s2smearing_mins'], dtype=np.float32)[mask_s2]
    d_acceptances[file_type]['s2']['smear_values_upper'] = np.asarray(d_acceptances[file_type]['df']['s2smearing_maxs'], dtype=np.float32)[mask_s2]
    
    print d_acceptances[file_type]['s2']['x_values']
    #print d_acceptances[file_type]['s2']['bias_values_lower']
    #print d_acceptances[file_type]['s2']['bias_values_upper']

    d_acceptances[file_type]['handler_bias'] = ax_s2_bias.fill_between(d_acceptances[file_type]['s2']['x_values'], d_acceptances[file_type]['s2']['bias_values_lower'], d_acceptances[file_type]['s2']['bias_values_upper'], color=d_acceptances[file_type]['color_bias'], alpha=alpha, label='Bias: %s' % d_acceptances[file_type]['label'])
    d_acceptances[file_type]['handler_smearing'] = ax_s2_smearing.fill_between(d_acceptances[file_type]['s2']['x_values'], d_acceptances[file_type]['s2']['smear_values_lower'], d_acceptances[file_type]['s2']['smear_values_upper'], color=d_acceptances[file_type]['color_smearing'], alpha=alpha, label='Smearing: %s' % d_acceptances[file_type]['label'])
    
    l_handlers_s2_bias.append(d_acceptances[file_type]['handler_bias'])
    l_handlers_s2_bias.append(d_acceptances[file_type]['handler_smearing'])
    #l_handlers_s2_smearing.append(d_acceptances[file_type]['handler_smearing'])





#ax_s1_bias.set_xlim(0, max_x)
ax_s1_bias.set_xlabel('Photoelectrons')
ax_s1_bias.set_ylabel('Bias')

ax_s1_smearing.set_xlabel('Photons Detected by PMTs')
ax_s1_smearing.set_ylabel('Smearing')

ax_s1_bias.legend(handles=l_handlers_s1_bias, loc='best', prop={'size': 16}, frameon=False)
ax_s1_smearing.legend(handles=l_handlers_s1_smearing, loc='best', prop={'size': 16}, frameon=False)



ax_s2_bias.set_xlabel('Photoelectrons')
ax_s2_bias.set_ylabel('Bias')

ax_s2_smearing.set_xlabel('Photons Detected by PMTs')
ax_s2_smearing.set_ylabel('Smearing')

ax_s2_bias.legend(handles=l_handlers_s2_bias, loc='best', prop={'size': 16}, frameon=False)
ax_s2_smearing.legend(handles=l_handlers_s2_smearing, loc='best', prop={'size': 16}, frameon=False)


fig_s1_bias.tight_layout()
fig_s2_bias.tight_layout()

fig_s1_bias.savefig('%sxe1t_pax_s1_bias.png' % (s_path_to_plots))
#fig_s1_smearing.savefig('%sxe1t_pax_s1_smearing.png' % (s_path_to_plots))
fig_s2_bias.savefig('%sxe1t_pax_s2_bias.png' % (s_path_to_plots))
#fig_s2_smearing.savefig('%sxe1t_pax_s2_smearing.png' % (s_path_to_plots))


plt.show()
