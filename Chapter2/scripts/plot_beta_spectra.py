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

#figure_size = (9, 5)

df_beta = pd.read_table('./kr85_beta_spec.txt')
df_beta['x'] *= 1000

a_energies = np.linspace(0, 687, 300)
func_beta = UnivariateSpline(df_beta['x'], df_beta['y'], k=5)

"""
# plot out theoretical
alpha = 1./137.
speed_of_light = 3e8
e_mass = 511 / (speed_of_light)**2
h_bar = 1.05e-34 # J*s


def fermi_func(atomic_number, atomic_mass, e_kin):
    s = (1. - alpha**2. * atomic_number**2.)**0.5
    energy = e_kin + e_mass*speed_of_light**2
    momentum = ((energy/speed_of_light)**2 - (e_mass*speed_of_light)**2)**0.5
    radius_of_nucleus = 1.25e-15*atomic_mass**(1./3.)
    eta = alpha * atomic_number * energy / (momentum*speed_of_light)
    rho = radius_of_nucleus / h_bar

    f = 2 * (1 + s) / gamma(1 + 2*s)**2. * (2*momentum*rho)**(2*s - 2) * np.exp(np.pi*eta) * np.absolute(gamma(s + np.complex(0,eta)))**2

    return f * momentum * energy * (687 - e_kin)**2


a_rate = np.zeros(len(a_energies))
for i, energy in enumerate(a_energies):
    a_rate[i] = fermi_func(36, 85, energy)

"""

fig, ax_beta = plt.subplots(1)

integral = func_beta.integral(0, 687)
a_y = func_beta(a_energies)/integral

roi_max = 15

roi_integral = func_beta.integral(0, roi_max)
print roi_integral/integral

mask_lt30 = (a_energies < roi_max)

s_plot = r'$^{85} \mathrm{Kr} \, \, \beta^-$ Spectrum'
s_plot += '\n'
s_plot += '$E_{\mathrm{max}} = 687 \, \, \mathrm{keV}$'

ax_beta.plot(a_energies, a_y, color='#0082c8', linestyle='-', linewidth=4)
ax_beta.fill_between(a_energies[mask_lt30], [0 for i in xrange(len(a_energies[mask_lt30]))], a_y[mask_lt30], facecolor='red', alpha=0.3)
ax_beta.text(470, 0.9*max(a_y), s_plot, size='medium', color='k')

ax_beta.set_xlim(1, 687)
#ax_beta.set_xscale('log')
ax_beta.set_ylim(0, 1.1*max(a_y))
#ax_beta.get_yaxis().set_visible(False)
ax_beta.set_ylabel(r'$p(E)$')
ax_beta.set_yticks([])

ax_beta.set_xlabel('$\mathrm{e}^-$ Kinetic Energy [keV]')
#ax_beta.set_ylabel('PDF')

fig.tight_layout()

fig.savefig('../images/kr85_beta_rates.png')

plt.show()
