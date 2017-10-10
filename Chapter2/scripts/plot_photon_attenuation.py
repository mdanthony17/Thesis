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

figure_size = (9, 6)

df_attenuation = pd.read_table('./xcom_photon_cross_sections.txt', sep=' ')
df_attenuation['pp'] = df_attenuation['pp_nf'] + df_attenuation['pp_ef']

fig, ax_attenuation = plt.subplots(1, figsize=figure_size)



handle_compton,  = ax_attenuation.plot(df_attenuation['energy'], df_attenuation['compton'], color='#f58231', linestyle='--', linewidth=2, label='Compton Scattering')
handle_pp,  = ax_attenuation.plot(df_attenuation['energy'], df_attenuation['pp'], color='#911eb4', linestyle='--', linewidth=2, label='Pair Production')
handle_pe,  = ax_attenuation.plot(df_attenuation['energy'], df_attenuation['pe_abs'], color='#008080', linestyle='--', linewidth=2, label='Photoelectric Absorption')
handle_tot,  = ax_attenuation.plot(df_attenuation['energy'], df_attenuation['tot_nc'], color='#000000', linestyle='-', linewidth=2, label='Total')


ax2 = ax_attenuation.twinx()
density_xenon = 2.84 # g /cm^3
ax2.plot(df_attenuation['energy'], 1. / (df_attenuation['tot_nc'] * density_xenon), color='#000000', linestyle='-', linewidth=2, label='Total')

ax_attenuation.set_xlim(1e-3, 10)
ax_attenuation.set_xscale('log')
ax_attenuation.set_ylim(1e-4, 1e4)
ax_attenuation.set_yscale('log')

ax2.set_ylim(1./(1e-4*density_xenon), 1./(1e4*density_xenon))
ax2.set_yscale('log')
ax2.set_ylabel('Attenuation Length [cm]')

ax_attenuation.set_xlabel('Photon Energy [MeV]')
ax_attenuation.set_ylabel('Mass Attenuation Coefficiency [$\mathrm{cm}^2 \, \, \mathrm{g}^{-1}$]')

ax_attenuation.legend(handles=[handle_tot, handle_pe, handle_compton, handle_pp], loc='best', prop={'size': 16}, frameon=False)

fig.savefig('../images/photon_attenuation.png')

plt.show()
