import wimp_functions

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

a_energies = np.linspace(0.01, 80, 200)
fixed_wimp_mass = 100 # GeV/c^2
fixed_cross_section = 1e-44 #cm^2

#print wimp_functions._wimp_recoil_spectrum(70, mass=100, sigma=1e-47, A=131.293, form_factor='Helm')


d_materials = {}
d_materials['xe'] = {'A':131.2}
d_materials['ge'] = {'A':72.64}
d_materials['ar'] = {'A':39.948}

fig, (ax_materials, ax_masses) = plt.subplots(1, 2)

for mat in d_materials:
    a_rates = np.zeros(len(a_energies))
    a_rates = wimp_functions.wimp_recoil_spectrum(a_energies, mass=fixed_wimp_mass, sigma=fixed_cross_section, A=d_materials[mat]['A'], form_factor='Helm')
    d_materials[mat]['a_rates'] = a_rates

    ax_materials.plot(a_energies, d_materials[mat]['a_rates'], color='blue', linestyle='-')



#handle_xe = ax_materials.plot(a_energies, d_materials[mat]['a_rates'], color='blue', linestyle='-')



ax_materials.set_yscale('log')


plt.show()


