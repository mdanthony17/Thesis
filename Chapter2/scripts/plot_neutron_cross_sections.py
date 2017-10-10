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

from scipy.interpolate import UnivariateSpline, interp1d
from scipy.special import gamma

figure_size = (9, 6)

d_neutrons = {}
d_neutrons['el'] = {'filename':'nr_cross_section_nel_all.txt', 'color':'#0082c8', 'linestyle':'-', 'data':{}, 'label':'Elastic Scatters'}
d_neutrons['g'] = {'filename':'nr_cross_section_ng_all.txt', 'color':'#3cb44b', 'linestyle':'-', 'data':{}, 'label':'Radiative Absorption'}
d_neutrons['inl'] = {'filename':'nr_cross_section_ninl_all.txt', 'color':'#ffe119', 'linestyle':'-', 'data':{}, 'label':'Inelastic Scatters'}

# keep track of when you need to increase atomic number
# special: #, //
b_last_special_character = True

l_isotopes = [128, 129, 130, 131, 132, 134, 136]
d_abundances = {128:.00190, 129:.164, 130:.04071, 131:.21232, 132:.26909, 134:.10436, 136:.26909}


for scatter_type in d_neutrons:
    print scatter_type
    f_data = open('../data/%s' % d_neutrons[scatter_type]['filename'], 'r')
    
    counter_mass_number = 0
    l_current_energies = []
    l_current_cross_sections = []
    d_neutrons[scatter_type]['data'][l_isotopes[counter_mass_number]] = {}
    for line in f_data:
        if line[0] == '\n':
            continue
        #print line
        if line[0] == '#' or line[0] == '/':
            if b_last_special_character == False and counter_mass_number < len(l_isotopes):
            
                # save arrays
                print l_isotopes[counter_mass_number]
                d_neutrons[scatter_type]['data'][l_isotopes[counter_mass_number]]['energy'] = np.asarray([0, l_current_energies[0]-1e-10] + l_current_energies + [l_current_energies[-1]+1e-10, 1e10])
                d_neutrons[scatter_type]['data'][l_isotopes[counter_mass_number]]['cross_section'] = np.asarray([0, 0] + l_current_cross_sections + [0, 0])
                d_neutrons[scatter_type]['data'][l_isotopes[counter_mass_number]]['cross_section_spline'] = interp1d(d_neutrons[scatter_type]['data'][l_isotopes[counter_mass_number]]['energy'], d_neutrons[scatter_type]['data'][l_isotopes[counter_mass_number]]['cross_section'])
                #print d_neutrons[scatter_type]['data'][l_isotopes[counter_mass_number]]
            
                counter_mass_number += 1
                if counter_mass_number < len(l_isotopes):
                    d_neutrons[scatter_type]['data'][l_isotopes[counter_mass_number]] = {}
                
                # reset lists
                l_current_energies = []
                l_current_cross_sections = []
            
            b_last_special_character = True
        else:
           l_current_line = line.split()
           l_current_energies.append(float(l_current_line[0]))
           l_current_cross_sections.append(float(l_current_line[1]))
           b_last_special_character = False


    f_data.close()


a_energies = np.logspace(-3, 3, 100000)

# combine datasets relative to abundance

for scatter_type in d_neutrons:
    d_neutrons[scatter_type]['nat_energy'] = a_energies
    d_neutrons[scatter_type]['nat_cross_section'] = np.zeros(len(d_neutrons[scatter_type]['nat_energy']))
    
    for mass_number in l_isotopes:
        #print d_abundances[mass_number]*d_neutrons[scatter_type]['data'][mass_number]['cross_section']
        d_neutrons[scatter_type]['nat_cross_section'] += d_abundances[mass_number]*d_neutrons[scatter_type]['data'][mass_number]['cross_section_spline'](a_energies)




fig, ax_neutron = plt.subplots(1, figsize=figure_size)

#s_plot = r'$^{85} \mathrm{Kr} \, \, \beta^-$ Spectrum'
#s_plot += '\n'
#s_plot += '$E_{\mathrm{max}} = 687 \, \, \mathrm{keV}$'

for scatter_type in d_neutrons:
    ax_neutron.plot(d_neutrons[scatter_type]['nat_energy'], d_neutrons[scatter_type]['nat_cross_section'], color=d_neutrons[scatter_type]['color'], linestyle='-', linewidth=1, label=d_neutrons[scatter_type]['label'])


#ax_beta.text(470, 0.9*max(a_y), s_plot, size='medium', color='k')

ax_neutron.set_xlim(1e-3, 10)
ax_neutron.set_xscale('log')
ax_neutron.set_ylim(1e-2, 1e3)
ax_neutron.set_yscale('log')

ax_neutron.legend(loc='best', prop={'size': 16}, frameon=False)

ax_neutron.set_xlabel('Kinetic Energy [keV]')
ax_neutron.set_ylabel('Cross Section [bn]')

fig.savefig('../images/neutron_cross_sections.png')

plt.show()

