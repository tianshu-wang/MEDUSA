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
parser = argparse.ArgumentParser(description="Plot the net heating in the gain region.")
parser.add_argument("--recompute",help="recompute values; otherwise, read from file [False]",action="store_true")
parser.add_argument("--theta",help="angle in degrees over which to compute polar luminosity [20.0]",type=float,default=20.0)
parser.add_argument("--xmin",help="xmin of plot window [0.0]",type=float,default=0.0)
parser.add_argument("--xmax",help="xmax of plot window [1.0]",type=float,default=1.0)
parser.add_argument("--ymin",help="ymin of plot window [0.0]",type=float,default=0.0)
parser.add_argument("--ymax",help="ymax of plot window [8.0]",type=float,default=8.0)
parser.add_argument("--dt",help="window width for time averaging [0.010]",metavar="DT_MEAN",type=float,default=0.010)
args = parser.parse_args()

# Open grid object
grid = fornax.GeomInterface(dump_dir+'grid.h5')
datafile = 'netheat_gain.npz'

if args.recompute:
    # Recompute and save data to file
    N        = dump_stop - dump_start + 1
    time     = np.zeros(N)
    netheat0 = np.zeros(N)
    netheat1 = np.zeros(N)
    netheat2 = np.zeros(N)
    netheat_pole0 = np.zeros(N)
    netheat_pole1 = np.zeros(N)
    netheat_pole2 = np.zeros(N)
    ir50 = np.where(grid.r >= 50.0)[0][0]
    jmax = np.where(grid.th <= args.theta*U_DEG)[0][-1]

    dd       = 0
    for d in range(dump_start,dump_stop+1):
        dump = fornax.DataInterface(dump_dir+"dump_%05d.h5" % d)
        # Get the entropy, use it to detect the max shock radius
        S = dump.get_var('S',grid)
        rsmax = 0.0
        for k in range(grid.nphi):
            for j in range(grid.nth):
                rsmax = max(rsmax,fornax.findshock(grid.rc,S[:,j,k]))

        # Get the mass and calculate the net heating rate from r=50km to rsmax
        # Ignore regions where the cooling rate dominates (max with 0)
        irsmax = np.where(grid.r >= rsmax)[0][0]
        mass = dump.get_var('rho',grid)*grid.dV()
        netheat = mass*(dump.get_var('heat0',grid) - dump.get_var('cool0',grid))
        netheat = np.maximum(netheat,0.0)
        netheat0[dd] = np.sum(np.sum(netheat[ir50:irsmax,:,:],axis=0))
        netheat_pole0[dd] = np.sum(np.sum(netheat[ir50:irsmax,:jmax,:],axis=0))

        netheat = mass*(dump.get_var('heat1',grid) - dump.get_var('cool1',grid))
        netheat = np.maximum(netheat,0.0)
        netheat1[dd] = np.sum(np.sum(netheat[ir50:irsmax,:,:],axis=0))
        netheat_pole1[dd] = np.sum(np.sum(netheat[ir50:irsmax,:jmax,:],axis=0))
 
        netheat = mass*(dump.get_var('heat2',grid) - dump.get_var('cool2',grid))
        netheat = np.maximum(netheat,0.0)
        netheat2[dd] = np.sum(np.sum(netheat[ir50:irsmax,:,:],axis=0))
        netheat_pole2[dd] = np.sum(np.sum(netheat[ir50:irsmax,:jmax,:],axis=0))

        time[dd] = dump.get_dataset('Time')[0] - tbounce
        dump.close()
        print 'd={0:3d}, t={1:5.3f}, netheat0={2:1.3e}, netheat_pole0={3:1.3e}, netheat1={4:1.3e}, netheat_pole1={5:1.3e}, netheat2={6:1.3e}, netheat_pole2={7:1.3e}'.format(d,time[dd],netheat0[dd],netheat_pole0[dd],netheat1[dd],netheat_pole1[dd],netheat2[dd],netheat_pole2[dd])
        dd += 1
    np.savez(datafile,time=time,netheat0=netheat0,netheat1=netheat1,netheat2=netheat2,
        netheat_pole0=netheat_pole0,netheat_pole1=netheat_pole1,netheat_pole2=netheat_pole2)
else:
    data     = np.load(datafile)
    time     = data['time']
    netheat0 = data['netheat0']
    netheat1 = data['netheat1']
    netheat2 = data['netheat2']
    netheat_pole0 = data['netheat_pole0']
    netheat_pole1 = data['netheat_pole1']
    netheat_pole2 = data['netheat_pole2']

# Take the mean over a window of width dt_mean
nh0,t  = fornax.tmean(netheat0,time,args.dt)
nhp0,t = fornax.tmean(netheat_pole0,time,args.dt)

# Plot & save figure
plt.figure()
plt.plot(t,nh0/1.0e51,label=r'$\nu_e$')
plt.xlabel('Time after bounce [s]')
plt.ylabel(r'$\int_{4\pi} (\mathcal{H}-\mathcal{C})_{\nu_e} \,d\Omega$ [$10^{51}$ erg g$^{-1}$ s$^{-1}$]')
plt.xlim(args.xmin,args.xmax)
plt.ylim(args.ymin,args.ymax)
savefile = 'rad_ccsn_{0:s}_netheat_gain'.format(save_text)
plt.savefig(savefile+'.eps')
plt.savefig(savefile+'.png')

plt.figure()
plt.plot(t,nhp0/1.0e51,label=r'$\nu_e$')
plt.xlabel('Time after bounce [s]')
plt.ylabel(r'$\int_{\pi/8} (\mathcal{H}-\mathcal{C})_{\nu_e} \,d\Omega$ [$10^{51}$ erg g$^{-1}$ s$^{-1}$]')
plt.xlim(args.xmin,args.xmax)
scl = 0.5*(1.0-np.cos(args.theta*U_DEG))
plt.ylim(args.ymin*scl,args.ymax*scl)
savefile = 'rad_ccsn_{0:s}_netheat_gain_pole'.format(save_text)
plt.savefig(savefile+'.eps')
plt.savefig(savefile+'.png')


print 'done'
