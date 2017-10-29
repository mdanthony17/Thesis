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
from rootpy.plotting import Hist, Hist2D, Canvas, Legend
from rootpy.io import File

import pandas as pd

from scipy.interpolate import UnivariateSpline
from scipy.special import gamma

s_path_to_file = '/Users/Matt/Desktop/Xenon/xenon1t/fit_nr_band/mcmc_analysis/radiogenic_neutrons_spectra/RadiogenicNeutrons_SR0.root'

f_rad = File(s_path_to_file, 'read')
h_energy = f_rad.hEdTotal_1tFVcylSR0

#print h_energy[:].value
#print h_energy.xedges()

c1 = Canvas()
h_energy.Draw()

h_energy.SetStats(0)
h_energy.SetTitle('')
h_energy.GetXaxis().SetTitle('Nuclear Recoil Energy [keV]')
h_energy.GetYaxis().SetTitle(r'NR Rate [s^{-1}]')


c1.Update()


raw_input('Press enter to continue...')






