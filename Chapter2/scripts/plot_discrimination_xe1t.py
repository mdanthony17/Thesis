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

figure_size = (10, 7)

df_ambe = pd.read_table('../data/ambe_nr.csv', sep=',')
df_rn = pd.read_table('../data/rn_er.csv', sep=',')
# cs1_pe, cs2_bottom_pe

fig, ax_disc = plt.subplots(1, figsize=figure_size)

marker_size = 12
alpha_level = 1

ax_disc.scatter(df_ambe['cs1_pe'], df_ambe['cs2_bottom_pe'], c='#0082c8', marker='o', linewidth=0, alpha=alpha_level, s=marker_size, label='Nuclear Recoils (AmBe)')
ax_disc.scatter(df_rn['cs1_pe'], df_rn['cs2_bottom_pe'], c='#e6194b', marker='o', linewidth=0, alpha=alpha_level, s=marker_size, label=r'Electronic Recoils ($^{222}\mathrm{Rn}$)')

ax_disc.set_xlim(0, 70)
#ax_disc.set_xscale('log')
ax_disc.set_ylim(100, 8000)
ax_disc.set_yscale('log')


ax_disc.set_xlabel('S1 [PE]')
ax_disc.set_ylabel('S2 [PE]')

ax_disc.legend(loc='upper left', prop={'size': 16}, frameon=False,scatterpoints=1)

fig.savefig('../images/xe1t_disc.png')

plt.show()