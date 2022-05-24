#!/usr/bin/env python
import numpy as np
import fornax
import matplotlib.pyplot as plt
from runinfo import *
from styles import *
from matplotlib import rcdefaults
import argparse

U_KM  = 1.0e5
U_DEG = np.pi/180.0

# Parse command line
parser = argparse.ArgumentParser(description="Plot the mean/polar luminosity at a given radius.")
parser.add_argument("-r","--radius",help="radius in km at which to compute mean luminosity [100.0]",type=float,default=100.0)
parser.add_argument("--theta",help="angle in degrees over which to compute polar luminosity [20.0]",type=float,default=20.0)
parser.add_argument("--recompute",help="recompute values; otherwise, read from file [False]",action="store_true")
parser.add_argument("--xmin",help="xmin of plot window [0.0]",type=float,default=0.0)
parser.add_argument("--xmax",help="xmax of plot window [1.0]",type=float,default=1.0)
parser.add_argument("--ymin",help="ymin of plot window [0.0]",type=float,default=0.0)
parser.add_argument("--ymax",help="ymax of plot window [15.0]",type=float,default=15.0)
parser.add_argument("--dt",help="window width for time averaging [0.010]",metavar="DT_MEAN",type=float,default=0.010)
args = parser.parse_args()

# Open grid object
grid = fornax.GeomInterface(dump_dir+'grid.h5')
datafile = 'Luminosity_r{:d}.npz'.format(int(args.radius))

if args.recompute:
    # Recompute and save data to file
    N = dump_stop - dump_start + 1
    time = np.zeros(N)
    L0   = np.zeros(N)
    L1   = np.zeros(N)
    L2   = np.zeros(N)
    L0_pole = np.zeros(N)
    L1_pole = np.zeros(N)
    L2_pole = np.zeros(N)
    jmax = np.where(grid.th <= args.theta*U_DEG)[0][-1]
    dOmega = grid.dOmega()
    dOmega_pole = dOmega[0:jmax,:]
    A = 4.0*np.pi*(args.radius*U_KM)**2

    dd   = 0
    for d in range(dump_start,dump_stop+1):
        dump = fornax.DataInterface(dump_dir+"dump_%05d.h5" % d)
        # Get the radial flux at radius r, average over solid angle,
        # then multiply by surface area
        F = dump.F_at_r(0,args.radius,grid)
        L0[dd] = A*np.sum(F*dOmega)/np.sum(dOmega)
        L0_pole[dd] = A*np.sum(F[0:jmax,:]*dOmega_pole)/np.sum(dOmega_pole)
        F = dump.F_at_r(1,args.radius,grid)
        L1[dd] = A*np.sum(F*dOmega)/np.sum(dOmega)
        L1_pole[dd] = A*np.sum(F[0:jmax,:]*dOmega_pole)/np.sum(dOmega_pole)
        F = dump.F_at_r(2,args.radius,grid)
        L2[dd] = A*np.sum(F*dOmega)/np.sum(dOmega)
        L2_pole[dd] = A*np.sum(F[0:jmax,:]*dOmega_pole)/np.sum(dOmega_pole)
        time[dd] = dump.get_dataset('Time')[0] - tbounce
        dump.close()
        print 'd={0:3d}, t={1:5.3f}, L0={2:1.3e}, L0_pole={3:1.3e}, L1={4:1.3e}, L1_pole={5:1.3e}, L2={6:1.3e}, L2_pole={7:1.3e}'.format(d,time[dd],L0[dd],L0_pole[dd],L1[dd],L1_pole[dd],L2[dd],L2_pole[dd])
        dd += 1
    np.savez(datafile,time=time,L0=L0,L1=L1,L2=L2,L0_pole=L0_pole,L1_pole=L1_pole,L2_pole=L2_pole)
else:
    # Load data from file
    data = np.load(datafile)
    time = data['time']
    L0   = data['L0']
    L1   = data['L1']
    L2   = data['L2']
    L0_pole = data['L0_pole']
    L1_pole = data['L1_pole']
    L2_pole = data['L2_pole']

# Take the mean over a window of width dt_mean
L0,t = fornax.tmean(L0,time,args.dt)
L1,t = fornax.tmean(L1,time,args.dt)
L2,t = fornax.tmean(L2,time,args.dt)
L0_pole,t = fornax.tmean(L0_pole,time,args.dt)
L1_pole,t = fornax.tmean(L1_pole,time,args.dt)
L2_pole,t = fornax.tmean(L2_pole,time,args.dt)
time = t

# Reset and load plotting style
rcdefaults()
set_style()

# Plot & save figure
plt.figure()
plt.plot(time,L0/1.0e52,'b-',label=r'$\nu_e$')
plt.plot(time,L1/1.0e52,'g-',label=r'$\nu_\bar{e}$')
plt.plot(time,L2/1.0e52,'r-',label=r'$\nu_\mu$')
plt.xlim(args.xmin,args.xmax)
plt.ylim(args.ymin,args.ymax)
plt.xlabel('Time after bounce [s]')
plt.ylabel(r'Luminosity [10$^{52}$ erg s$^{-1}$]')
plt.title(title_text)
plt.legend(loc='upper right')
savefile = 'rad_ccsn_{0:s}_Luminosity_mean_r{1:d}'.format(save_text,int(args.radius))
plt.savefig(savefile+'.eps')
plt.savefig(savefile+'.png')

plt.figure()
plt.plot(time,L0_pole/1.0e52,'b-',label=r'$\nu_e$')
plt.plot(time,L1_pole/1.0e52,'g-',label=r'$\nu_\bar{e}$')
plt.plot(time,L2_pole/1.0e52,'r-',label=r'$\nu_\mu$')
plt.xlim(args.xmin,args.xmax)
plt.ylim(args.ymin,args.ymax)
plt.xlabel('Time after bounce [s]')
plt.ylabel(r'Luminosity at $\theta=0$ [10$^{52}$ erg s$^{-1}$]')
plt.title(title_text)
plt.legend(loc='upper right')
savefile = 'rad_ccsn_{0:s}_Luminosity_pole_r{1:d}'.format(save_text,int(args.radius))
plt.savefig(savefile+'.eps')
plt.savefig(savefile+'.png')


print 'done'
