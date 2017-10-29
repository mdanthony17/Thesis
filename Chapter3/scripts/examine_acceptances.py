#!/usr/bin/python
import sys, array, os

import matplotlib
matplotlib.use('QT4Agg')

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 16

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.patches as mpatches

import numpy as np
import time, tqdm
import cPickle as pickle

import pandas as pd
import scipy.interpolate

alpha = 0.3
l_file_types = ['ambe', 'rn', 'wimp']
d_acceptances = {}

for file_type in l_file_types:
    if file_type == 'ambe':
        d_acceptances[file_type] = {}
        d_acceptances[file_type]['file_path'] = '/Users/Matt/Desktop/Xenon/xenon1t/fit_nr_band/resources/sr0_efficiency_ambe_conditions_goodnoise.csv'
        d_acceptances[file_type]['color'] = '#0082c8'
        d_acceptances[file_type]['df'] = pd.read_csv(d_acceptances[file_type]['file_path'])
        d_acceptances[file_type]['label'] = 'NR Calibration'
    elif file_type == 'rn':
        d_acceptances[file_type] = {}
        d_acceptances[file_type]['file_path'] = '/Users/Matt/Desktop/Xenon/xenon1t/fit_nr_band/resources/sr0_efficiency_good_conditions.csv'
        d_acceptances[file_type]['color'] = '#e6194b'
        d_acceptances[file_type]['df'] = pd.read_csv(d_acceptances[file_type]['file_path'])
        d_acceptances[file_type]['label'] = 'ER Calibration'
    elif file_type == 'wimp':
        d_acceptances[file_type] = {}
        d_acceptances[file_type]['file_path'] = '/Users/Matt/Desktop/Xenon/xenon1t/fit_nr_band/resources/sr0_efficiency_weightedaverage_conditions.csv'
        d_acceptances[file_type]['color'] = '#3cb44b'
        d_acceptances[file_type]['df'] = pd.read_csv(d_acceptances[file_type]['file_path'])
        d_acceptances[file_type]['label'] = 'Science Run'



s_path_to_plots = '../images/'
l_handlers = []
fig_pf_s1, ax_pf_s1 = plt.subplots(1)

max_x = -1


for file_type in l_file_types:
    d_acceptances[file_type]['pf_s1'] = {}
    d_acceptances[file_type]['pf_s1']['x_values'] = np.asarray(d_acceptances[file_type]['df']['photons_detected'], dtype=np.float32)
    d_acceptances[file_type]['pf_s1']['y_values_lower'] = np.asarray(d_acceptances[file_type]['df']['minimum'], dtype=np.float32)
    d_acceptances[file_type]['pf_s1']['y_values_mean'] = np.asarray(d_acceptances[file_type]['df']['median'], dtype=np.float32)
    d_acceptances[file_type]['pf_s1']['y_values_upper'] = np.asarray(d_acceptances[file_type]['df']['maximum'], dtype=np.float32)
    
    if max(d_acceptances[file_type]['pf_s1']['x_values']) > max_x:
        max_x = max(d_acceptances[file_type]['pf_s1']['x_values'])


    d_acceptances[file_type]['handler'] = ax_pf_s1.fill_between(d_acceptances[file_type]['pf_s1']['x_values'], d_acceptances[file_type]['pf_s1']['y_values_lower'], d_acceptances[file_type]['pf_s1']['y_values_upper'], color=d_acceptances[file_type]['color'], alpha=alpha, label=d_acceptances[file_type]['label'])
    #ax_pf_s1.plot(d_acceptances[file_type]['pf_s1']['x_values'], d_acceptances[file_type]['pf_s1']['y_values_mean'], color=d_acceptances[file_type]['color'], linestyle='--')
    l_handlers.append(d_acceptances[file_type]['handler'])


ax_pf_s1.set_ylim(0, 1.03)
ax_pf_s1.axhline(1, color='black', linestyle='--')

ax_pf_s1.set_xlim(0, max_x)
ax_pf_s1.set_xlabel('Photons Detected by PMTs')
ax_pf_s1.set_ylabel('Efficiency')

ax_pf_s1.legend(handles=l_handlers, loc='lower right', prop={'size': 16}, frameon=False)

fig_pf_s1.tight_layout()

fig_pf_s1.savefig('%sxe1t_pax_acceptance.png' % (s_path_to_plots))


plt.show()
