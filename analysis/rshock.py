#!/usr/bin/env python
import numpy as np
import fornax
import matplotlib.pyplot as plt
from runinfo import *
from styles import *
import argparse
from matplotlib import rcdefaults

U_KM  = 1.0e5

# Parse command line
parser = argparse.ArgumentParser(description="Plot the min/mean/max shock radii.")
parser.add_argument("--recompute",help="recompute values; otherwise, read from file [False]",action="store_true")
parser.add_argument("--xmin",help="xmin of plot window [0.0]",type=float,default=0.0)
parser.add_argument("--xmax",help="xmax of plot window [1.0]",type=float,default=1.0)
parser.add_argument("--ymin",help="ymin of plot window [0.0]",type=float,default=0.0)
parser.add_argument("--ymax",help="ymax of plot window [200.0]",type=float,default=200.0)
parser.add_argument("--dt",help="window width for time averaging [0.010]",metavar="DT_MEAN",type=float,default=0.010)
args = parser.parse_args()

# Open grid object
grid = fornax.GeomInterface(dump_dir+'grid.h5')
datafile = 'rshock.npz'

if args.recompute:
    # Recompute and save data to file
    N      = dump_stop - dump_start + 1
    time   = np.zeros(N)
    rs     = np.zeros((grid.nth,grid.nphi))
    rsmin  = np.zeros(N)
    rsmax  = np.zeros(N)
    rsmean = np.zeros(N)
    dOmega = grid.dOmega()
    dd     = 0
    for d in range(dump_start,dump_stop+1):
        dump = fornax.DataInterface(dump_dir+"dump_%05d.h5" % d)
        # Get the entropy, use it to detect the shock radius
        S  = dump.get_var('S',grid)
        for k in range(grid.nphi):
            for j in range(grid.nth):
                rs[j,k] = fornax.findshock(grid.rc,S[:,j,k])
        rsmin[dd]  = min(rs)
        rsmax[dd]  = max(rs)
        rsmean[dd] = np.sum(rs*dOmega)/np.sum(dOmega)
        time[dd]   = dump.get_dataset('Time')[0] - tbounce
        dump.close()
        print 'd={0:3d}, t={1:5.3f}, rsmin={2:5.1f}, rsmean={3:5.1f}, rsmax={4:5.1f}'.format(d,time[dd],rsmin[dd],rsmean[dd],rsmax[dd])
        dd += 1
    np.savez(datafile,time=time,rsmin=rsmin,rsmean=rsmean,rsmax=rsmax)
else:
    # Load data from file
    data   = np.load(datafile)
    time   = data['time']
    rsmin  = data['rsmin']
    rsmean = data['rsmean']
    rsmax  = data['rsmax']

# Take the mean over a window of width dt_mean
rsmin,t  = fornax.tmean(rsmin,time,args.dt)
rsmean,t = fornax.tmean(rsmean,time,args.dt)
rsmax,t  = fornax.tmean(rsmax,time,args.dt)
time     = t

# Reset and load plotting style
rcdefaults()
set_style()

# Plot & save figure
plt.figure()
plt.plot(time,rsmin ,'b',label='Min')
plt.plot(time,rsmax ,'r',label='Max')
plt.plot(time,rsmean,'g',label='Avg')
plt.xlim(args.xmin,args.xmax)
plt.ylim(args.ymin,args.ymax)
plt.xlabel('Time after bounce [s]')
plt.ylabel('Shock radius [km]')
plt.title(title_text)
plt.figtext(0.75,0.82,r'${0:d} M_\odot$'.format(model))
plt.figtext(0.75,0.75,eos)
plt.figtext(0.75,0.68,method)
plt.legend(loc='upper left')
savefile = 'rad_ccsn_{0:s}_rshock'.format(save_text)
plt.savefig(savefile+'.eps')
plt.savefig(savefile+'.png')


print 'done'
