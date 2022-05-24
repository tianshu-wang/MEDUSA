#!/usr/bin/env python
import numpy as np
import fornax
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from runinfo import *
from styles import *
from matplotlib import rcdefaults
import argparse

# Parse command line
parser = argparse.ArgumentParser(description="Plot the angular means of several variables vs. radius and mass.")
parser.add_argument("-s","--skip",help="number of output dumps to skip by [5]",type=int,default=5)
parser.add_argument("-t","--tstop",help="max number of output dumps to plot [200]",type=int,default=200)
parser.add_argument("--rmax",help="max radius in km of plot window [50.0]",type=float,default=50.0)
parser.add_argument("--Mmax",help="max mass in Msun of plot window [1.0]",type=float,default=1.0)
args = parser.parse_args()

dump_stop = dump_start + args.tstop
grid = fornax.GeomInterface(dump_dir+'grid.h5')

# Set/normalize colormap
cNorm = colors.Normalize(vmin=dump_start, vmax=dump_stop+1)
scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cm.get_cmap('jet'))

# Plot parameters for different variables
plots = [
{'var':'rho', 'name':'rho'             , 'logy':False, 'ymin':0.0   , 'ymax':3.5e14, 'ylabel':r'$\rho$ [g cm$^{-3}$]'},
{'var':'vr' , 'name':'u1'              , 'logy':False, 'ymin':-1.0e8, 'ymax':1.0e8 , 'ylabel':r'$v_r$ [km s$^{-1}$]'},
{'var':'e'  , 'name':'u'               , 'logy':False, 'ymin':0.0   , 'ymax':2.0e34, 'ylabel':r'$e$ [erg cm$^{-3}$]'},
{'var':'Ye' , 'name':'comp0'           , 'logy':False, 'ymin':0.0   , 'ymax':0.5   , 'ylabel':r'$Y_e$ [g cm$^{-3}$]'},
{'var':'Er0', 'name':'Er0'             , 'logy':False, 'ymin':0.0   , 'ymax':1.8e33, 'ylabel':r'$\mathcal{E}_{\nu_e}$ [erg cm$^{-3}$]'},
{'var':'Fr0', 'name':'Fr0'             , 'logy':False, 'ymin':0.0   , 'ymax':2.0e39, 'ylabel':r'$F_{r,\nu_e}$ [erg s$^{-1}$ cm$^{-2}$]'},
{'var':'L0' , 'name':'L0'              , 'logy':False, 'ymin':0.0   , 'ymax':2.0e50, 'ylabel':r'$L_{r,\nu_e}$ [erg s$^{-1}$]'},
{'var':'S'  , 'name':'eos3'            , 'logy':False, 'ymin':0.0   , 'ymax':20.0  , 'ylabel':r'$S$ [$k_\mathrm{B}$ baryon$^{-1}$]'},
{'var':'T'  , 'name':'eos2'            , 'logy':False, 'ymin':0.0   , 'ymax':30.0  , 'ylabel':r'$T$ [MeV]'},
{'var':'P'  , 'name':'eos0'            , 'logy':False, 'ymin':0.0   , 'ymax':1.4e34, 'ylabel':r'$P$ [erg cm$^{-3}$]'}
]

# Loop over all variables and plot figure
for p in plots:
    p['fig'], (p['ax1'], p['ax2']) = plt.subplots(1,2,sharey=True)

for d in range(dump_start,dump_stop+1,args.skip):
    dump = fornax.DataInterface(dump_dir+"dump_%05d.h5" % d)
    r = grid.rc
    M = dump.get_Mass(grid)
  
    for p in plots:
        var = dump.get_mean_4pi(p['name'],grid)
        line = p['ax1'].plot(r,var)
        colorVal = scalarMap.to_rgba(d)
        plt.setp(line,color=colorVal,antialiased=True)
        line = p['ax2'].plot(M,var)
        colorVal = scalarMap.to_rgba(d)
        plt.setp(line,color=colorVal,antialiased=True)

    dump.close()

# Set plot parameters & save figure
for p in plots:
    p['ax1'].set_xlim(0.0,args.rmax)
    p['ax1'].set_xlabel(r'$r$ [km]')
    if (p['ymax']>p['ymin']):
        p['ax1'].set_ylim(p['ymin'],p['ymax'])
    p['ax1'].set_ylabel(p['ylabel'])

    p['ax2'].set_xlim(0.0,args.Mmax)
    p['ax2'].set_xlabel(r'$M$ [$M_\odot$]')

    if p['logy']:
        p['ax1'].set_yscale("log")
        p['ax2'].set_yscale("log")

    plt.figure(p['fig'].number)
    p['fig'].suptitle(title_text)
    savefile = 'rad_ccsn_{0:s}_{1:s}_4pi'.format(save_text,p['var'])
    plt.savefig(savefile+'.eps')
    plt.savefig(savefile+'.png')


print 'done'
