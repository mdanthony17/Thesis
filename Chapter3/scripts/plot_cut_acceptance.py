import numpy as np
import matplotlib
matplotlib.use('QT4Agg')

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 16

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.patches as mpatches
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
from matplotlib.ticker import LinearLocator

import pandas as pd

from scipy.interpolate import UnivariateSpline
from scipy.special import gamma

figure_size = [7,4]

fig, ax_cut = plt.subplots(1, figsize=figure_size)

cut_acceptance_s1_intercept = 0.96517
cut_acceptance_s1_intercept_uncertainty = np.abs(0.0421*cut_acceptance_s1_intercept)

cut_acceptance_s1_slope = -0.000466585
cut_acceptance_s1_slope_uncertainty = np.abs(0.0421*cut_acceptance_s1_slope)

cut_acceptance_s2_intercept = 0.851629
cut_acceptance_s2_intercept_uncertainty = np.abs(0.0176*cut_acceptance_s2_intercept)

cut_acceptance_s2_slope = 7.20775e-06
cut_acceptance_s2_slope_uncertainty = np.abs(0.0176*cut_acceptance_s2_slope)


num_points = 100
num_mc_samples = 10000

a_s1_lb = np.zeros(num_points)
a_s1_med = np.zeros(num_points)
a_s1_ub = np.zeros(num_points)
a_s2_lb = np.zeros(num_points)
a_s2_med = np.zeros(num_points)
a_s2_ub = np.zeros(num_points)

a_holder_s1 = np.zeros(num_mc_samples)
a_holder_s1 = np.zeros(num_mc_samples)

low_s1 = 3
high_s1 = 200
low_s2 = 200
high_s2 = 15000

a_s1 = np.logspace(np.log10(low_s1), np.log10(high_s1), num_points)
a_s2 = np.logspace(np.log10(low_s2), np.log10(high_s2), num_points)

for current_point in xrange(num_points):
    current_s1 = a_s1[current_point]
    current_s2 = a_s2[current_point]

    #for i in xrange(num_mc_samples):
    a_holder_s1 = np.random.normal(cut_acceptance_s1_intercept, cut_acceptance_s1_intercept_uncertainty, num_mc_samples) + current_s1*np.random.normal(cut_acceptance_s1_slope, cut_acceptance_s1_slope_uncertainty, num_mc_samples)

    a_s1_lb[current_point], a_s1_med[current_point], a_s1_ub[current_point] = np.percentile(a_holder_s1, [16, 50, 84])

    if (a_s1_lb[current_point] > 1):
        a_s1_lb[current_point] = 1

    if (a_s1_med[current_point] > 1):
        a_s1_med[current_point] = 1
        
    if (a_s1_ub[current_point] > 1):
        a_s1_ub[current_point] = 1


    a_holder_s2 = np.random.normal(cut_acceptance_s2_intercept, cut_acceptance_s2_intercept_uncertainty, num_mc_samples) + current_s2*np.random.normal(cut_acceptance_s2_slope, cut_acceptance_s2_slope_uncertainty, num_mc_samples)

    a_s2_lb[current_point], a_s2_med[current_point], a_s2_ub[current_point] = np.percentile(a_holder_s2, [16, 50, 84])

    if (a_s2_lb[current_point] > 1):
        a_s2_lb[current_point] = 1

    if (a_s2_med[current_point] > 1):
        a_s2_med[current_point] = 1
        
    if (a_s2_ub[current_point] > 1):
        a_s2_ub[current_point] = 1




handle_s1, = ax_cut.plot(a_s1, a_s1_med, color='#0082c8', linestyle='--', linewidth=1, label='S1 Cut Acceptance')
ax_cut.fill_between(a_s1, a_s1_lb, a_s1_ub, facecolor='#0082c8', alpha=0.3)


ax_cut.set_xlim(low_s1, high_s1)
ax_cut.set_xscale('log')
ax_cut.set_ylim(0, 1.)

ax_cut.set_xlabel('S1 [PE]')
ax_cut.set_ylabel('Cut Acceptance')

ax_s2_cut = ax_cut.twiny()

handle_s2, = ax_s2_cut.plot(a_s2, a_s2_med, color='#f58231', linestyle='--', linewidth=1, label='S2 Cut Acceptance')
ax_s2_cut.fill_between(a_s2, a_s2_lb, a_s2_ub, facecolor='#f58231', alpha=0.3)


ax_s2_cut.set_ylim(0, 1.)
ax_s2_cut.set_xscale('log')
ax_s2_cut.set_xlim(low_s2, high_s2)
ax_s2_cut.set_xlabel('S2 [PE]')

#ax_cut.legend([handle_s1, handle_s2], ['S1 Cut Acceptance', 'S2 Cut Acceptance'], loc='best', prop={'size': 16}, frameon=False)
ax_s2_cut.legend([handle_s1, handle_s2], ['S1 Cut Acceptance', 'S2 Cut Acceptance'], loc='best', prop={'size': 16}, frameon=False)

fig.tight_layout()

fig.savefig('../images/xe1t_cut_acceptances.png')

plt.show()



