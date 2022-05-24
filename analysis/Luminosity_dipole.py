#!/usr/bin/env python
import numpy as np
import fornax
import matplotlib.pyplot as plt
from scipy.special import sph_harm
from runinfo import *
from styles import *
from matplotlib import rcdefaults
import argparse

U_KM  = 1.0e5
U_DEG = np.pi/180.0

# Parse command line
parser = argparse.ArgumentParser(description="Plot the normalized dipoles of the shock radius and luminosity (at a given radius).")
parser.add_argument("-r","--radius",help="radius in km at which to compute mean luminosity [100.0]",type=float,default=100.0)
parser.add_argument("--recompute",help="recompute values; otherwise, read from file [False]",action="store_true")
parser.add_argument("--xmin",help="xmin of plot window [0.0]",type=float,default=0.0)
parser.add_argument("--xmax",help="xmax of plot window [1.0]",type=float,default=1.0)
parser.add_argument("--ymin",help="ymin of plot window [-0.25]",type=float,default=-0.25)
parser.add_argument("--ymax",help="ymax of plot window [0.25]",type=float,default=0.25)
parser.add_argument("--dt",help="window width for time averaging [0.010]",metavar="DT_MEAN",type=float,default=0.010)
args = parser.parse_args()

# Open grid object
grid = fornax.GeomInterface(dump_dir+'grid.h5')
datafile = 'Luminosity_dipole_r{:d}.npz'.format(int(args.radius))

if args.recompute:
    # Recompute and save data to file
    N      = dump_stop - dump_start + 1
    time   = np.zeros(N)
    a1rs   = np.zeros(N)
    a1L0   = np.zeros(N)
    a1L1   = np.zeros(N)
    a1L2   = np.zeros(N)
    a1Ltot = np.zeros(N)
    dOmega = grid.dOmega()
    # Compute spherical harmonics (note scipy switches roles of theta & phi)
    PH,TH  = np.meshgrid(grid.phic,grid.thc)
    y00    = np.real(sph_harm(0,0,PH,TH))
    y01    = np.real(sph_harm(0,1,PH,TH))
    eps    = 1.0e-10

    dd     = 0
    for d in range(dump_start,dump_stop+1):
        dump = fornax.DataInterface(dump_dir+"dump_%05d.h5" % d)
        # Get the entropy, use it to detect the shock radius
        S  = dump.get_var('S',grid)
        rs = np.zeros((grid.nth,grid.nphi))
        for k in range(grid.nphi):
            for j in range(grid.nth):
                rs[j,k] = fornax.findshock(grid.rc,S[:,j,k])
        a0 = np.sum(rs*y00*dOmega)
        a1 = np.sum(rs*y01*dOmega)
        a1rs[dd] = a1/max(a0,eps)
        L0 = dump.L_at_r(0,args.radius,grid)
        a0 = np.sum(L0*y00*dOmega)
        a1 = np.sum(L0*y01*dOmega)
        a1L0[dd] = a1/max(a0,eps)
        L1 = dump.L_at_r(1,args.radius,grid)
        a0 = np.sum(L1*y00*dOmega)
        a1 = np.sum(L1*y01*dOmega)
        a1L1[dd] = a1/max(a0,eps)
        L2 = dump.L_at_r(2,args.radius,grid)
        a0 = np.sum(L2*y00*dOmega)
        a1 = np.sum(L2*y01*dOmega)
        a1L2[dd] = a1/max(a0,eps)
        Ltot = L0 + L1 + L2
        a0 = np.sum(Ltot*y00*dOmega)
        a1 = np.sum(Ltot*y01*dOmega)
        a1Ltot[dd] = a1/max(a0,eps)
        time[dd] = dump.get_dataset('Time')[0] - tbounce
        dump.close()
        print 'd={0:3d}, t={1:5.3f}, a1rs={2:6.3f}, a1L0={3:6.3f}, a1L1={4:6.3f}, a1L2={5:6.3f}, a1Ltot={6:6.3f}'.format(d,time[dd],a1rs[dd],a1L0[dd],a1L1[dd],a1L2[dd],a1Ltot[dd])
        dd += 1
    np.savez(datafile,time=time,a1rs=a1rs,a1L0=a1L0,a1L1=a1L1,a1L2=a1L2,a1Ltot=a1Ltot)
else:
    # Load data from file
    data   = np.load(datafile)
    time   = data['time']
    a1rs   = data['a1rs']
    a1L0   = data['a1L0']
    a1L1   = data['a1L1']
    a1L2   = data['a1L2']
    a1Ltot = data['a1Ltot']

# Reset and load plotting style
rcdefaults()
set_style()

# Plot & save figure
plt.figure()
plt.plot(time,a1rs,'k-',label=r'$r_\mathrm{shock}$')
plt.plot(time,a1L0,'r-',label=r'$L_{\nu_e,100}$')
plt.xlim(args.xmin,args.xmax)
plt.ylim(args.ymin,args.ymax)
plt.xlabel('Time after bounce [s]')
plt.ylabel(r'$a_1/a_0$')
plt.figtext(0.15,0.82,r'${0:d} M_\odot$'.format(model))
plt.figtext(0.15,0.75,'LS220')
plt.figtext(0.15,0.68,method)
plt.legend(loc='lower right')
savefile = 'rad_ccsn_{0:s}_Ltot_dipole_r{1:d}'.format(save_text,int(args.radius))
plt.savefig(savefile+'.eps')
plt.savefig(savefile+'.png')


print 'done'
