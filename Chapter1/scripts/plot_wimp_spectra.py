import wimp_functions

import numpy as np
import matplotlib
matplotlib.use('QT4Agg')
import matplotlib.style
matplotlib.style.use('classic')


matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 16

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.patches as mpatches
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
from matplotlib.ticker import LinearLocator

figure_size = (9, 5)

a_energies = np.linspace(0.01, 80, 200)
fixed_wimp_mass = 100 # GeV/c^2
fixed_cross_section = 1e-47 #cm^2

earth_velocity = 245
wimp_velocity_dispersion = 220.
escaping_velocity = 544.
wimp_local_density = 0.3


#print wimp_functions._wimp_recoil_spectrum(70, mass=100, sigma=1e-47, A=131.293, form_factor='Helm')


d_materials = {}
d_materials['Xe'] = {'A':131.2, 'color':'#0082c8'}
d_materials['Ge'] = {'A':72.64, 'color':'#e6194b'}
d_materials['Ar'] = {'A':39.948, 'color':'#3cb44b'}

l_handles_materials = []
l_handles_masses = []


fig, (ax_materials, ax_masses) = plt.subplots(1, 2, figsize=figure_size)
fig_mat_only, ax_mat_only = plt.subplots()
fig_mass_only, ax_mass_only = plt.subplots()

for mat in d_materials:
    a_rates = np.zeros(len(a_energies))
    for i, energy in enumerate(a_energies):
        a_rates[i] = wimp_functions.getWIMPRateIsotope(d_materials[mat]['A'], fixed_wimp_mass, fixed_cross_section, energy, earth_velocity, wimp_velocity_dispersion, escaping_velocity, wimp_local_density)
    d_materials[mat]['a_rates'] = a_rates

    d_materials[mat]['handle'], = ax_materials.plot(a_energies, d_materials[mat]['a_rates'], color=d_materials[mat]['color'], linestyle='-', label=mat)
    d_materials[mat]['handle_mat'], = ax_mat_only.plot(a_energies, d_materials[mat]['a_rates'], color=d_materials[mat]['color'], linestyle='-', label=mat)


#handle_xe = ax_materials.plot(a_energies, d_materials[mat]['a_rates'], color='blue', linestyle='-')

ax_materials.set_ylim(1e-9, 1e-6)
ax_materials.set_yscale('log')
ax_materials.set_ylabel(r'$\frac{dR}{dE} \, \, [\mathrm{events} \, \, \mathrm{kg}^{-1} \, \, \mathrm{keV}^{-1} \, \, \mathrm{day}^{-1}]$')
ax_materials.set_xlabel(r'$E \, \, [\mathrm{keV}]$')
ax_materials.legend(handles=[d_materials['Xe']['handle'], d_materials['Ge']['handle'], d_materials['Ar']['handle']], loc='best', frameon=False, prop={'size': 11})

ax_mat_only.set_xlim(0, 80)
ax_mat_only.set_ylim(1e-9, 1e-6)
ax_mat_only.set_yscale('log')
ax_mat_only.set_ylabel(r'$\frac{dR}{dE} \, \, [\mathrm{events} \, \, \mathrm{kg}^{-1} \, \, \mathrm{keV}^{-1} \, \, \mathrm{day}^{-1}]$')
ax_mat_only.set_xlabel(r'$E \, \, [\mathrm{keV}]$')
ax_mat_only.legend(handles=[d_materials['Xe']['handle_mat'], d_materials['Ge']['handle_mat'], d_materials['Ar']['handle_mat']], loc='best', frameon=False, prop={'size': 11})
ax_mat_only.annotate(r'100 GeV/$\mathrm{c}^2$', xy=(0.5, 0.9), xycoords='axes fraction', ha='center', va='center')
ax_mat_only.annotate(r'$\sigma = 10^{-47} \mathrm{cm}^2$', xy=(0.5, 0.82), xycoords='axes fraction', ha='center', va='center')


l_wimp_masses = [10, 30, 50, 100, 500, 1000]
l_colors = plt.get_cmap('jet')(np.linspace(0, 1.0, len(l_wimp_masses)))
l_handles_masses_only = []

for current_color, wimp_mass in zip(l_colors, l_wimp_masses):
    a_rates = np.zeros(len(a_energies))
    for i, energy in enumerate(a_energies):
        a_rates[i] = wimp_functions.getWIMPRateIsotope(d_materials['Xe']['A'], wimp_mass, fixed_cross_section, energy, earth_velocity, wimp_velocity_dispersion, escaping_velocity, wimp_local_density)

    current_handle, = ax_masses.plot(a_energies, a_rates, color=current_color, linestyle='-', label=r'$%d \, \, \frac{\mathrm{GeV}}{\mathrm{c}^2}$' % (wimp_mass))
    l_handles_masses.append(current_handle)

    current_handle, = ax_mass_only.plot(a_energies, a_rates, color=current_color, linestyle='-', label=r'$%d \, \, \frac{\mathrm{GeV}}{\mathrm{c}^2}$' % (wimp_mass))
    l_handles_masses_only.append(current_handle)


ax_masses.set_ylim(1e-11, 1e-5)
ax_masses.set_yscale('log')
ax_masses.legend(ncol=2, handles=l_handles_masses, loc='best', frameon=False, prop={'size': 11})
ax_masses.set_ylabel(r'$\frac{dR}{dE} \, \, [\mathrm{events} \, \, \mathrm{kg}^{-1} \, \, \mathrm{keV}^{-1} \, \, \mathrm{day}^{-1}]$')
ax_masses.set_xlabel(r'$E \, \, [\mathrm{keV}]$')


ax_mass_only.annotate(r'$\sigma = 10^{-47} \mathrm{cm}^2$', xy=(0.5, 0.85), xycoords='axes fraction', ha='center', va='center')

ax_mass_only.set_ylim(1e-11, 1e-5)
ax_mass_only.set_yscale('log')
ax_mass_only.legend(ncol=2, handles=l_handles_masses_only, loc='best', frameon=False, prop={'size': 11})
ax_mass_only.set_ylabel(r'$\frac{dR}{dE} \, \, [\mathrm{events} \, \, \mathrm{kg}^{-1} \, \, \mathrm{keV}^{-1} \, \, \mathrm{day}^{-1}]$')
ax_mass_only.set_xlabel(r'$E \, \, [\mathrm{keV}]$')


fig.tight_layout()
fig_mat_only.tight_layout()
fig.savefig('../images/wimp_recoil_rates.png')
fig_mat_only.savefig('../images/wimp_recoil_rates_mat_only.png')
fig_mass_only.savefig('../images/wimp_recoil_rates_mass_only.png')

plt.show()


