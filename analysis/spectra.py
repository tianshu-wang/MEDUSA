#!/usr/bin/env python
import numpy as np
import fornax
import matplotlib.pyplot as plt
from runinfo import *
#from styles import *
from matplotlib import rcdefaults
import argparse
from matplotlib import animation
import matplotlib

U_KM  = 1.0e5
U_DEG = np.pi/180.0

# Parse command line
parser = argparse.ArgumentParser(description="Plot the mean spectra at a given radius.")
parser.add_argument("-r","--radius",help="radius in km at which to compute mean spectra [100.0]",type=float,default=100.0)
parser.add_argument("--recompute",help="recompute values; otherwise, read from file [False]",action="store_true")
parser.add_argument("--xmin",help="xmin of plot window [0.0]",type=float,default=0.0)
parser.add_argument("--xmax",help="xmax of plot window [60.0]",type=float,default=60.0)
parser.add_argument("--ymin",help="ymin of plot window [0.0]",type=float,default=0.0)
parser.add_argument("--ymax",help="ymax of plot window [2.0]",type=float,default=2.0)
parser.add_argument("--dt",help="window width for time averaging [0.010]",metavar="DT_MEAN",type=float,default=0.010)
parser.add_argument("--fps",help="fps [15]",type=int,default=15)
parser.add_argument("--dumps_pad",help="start this many milliseconds before bounce [40]",type=int,default=40)
parser.add_argument("--dumps_fin",help="end this many milliseconds after bounce [100]",type=int,default=100)
args = parser.parse_args()

# Open grid object
grid = fornax.GeomInterface(dump_dir+'grid.h5')
datafile = 'Spectra_r{:d}.npz'.format(int(args.radius))

if args.recompute:
    # Recompute and save data to file
    N = args.dumps_fin +args.dumps_pad + 1
    ngroup0=grid.nr1
    ngroup1=grid.nr2
    ngroup2=grid.nr3
    time = np.zeros(N)
    eg0=np.zeros(ngroup0)
    eg1=np.zeros(ngroup1)
    eg2=np.zeros(ngroup2)
    deg0=np.zeros(ngroup0)
    deg1=np.zeros(ngroup1)
    deg2=np.zeros(ngroup2)
    Lspec0   = np.zeros((N,ngroup0))
    Lspec1   = np.zeros((N,ngroup1))
    Lspec2   = np.zeros((N,ngroup2))
    dOmega = grid.dOmega()

    dd   = 0
    for d in range(dump_start-args.dumps_pad,args.dumps_fin+dump_start+1):
        dump = fornax.DataInterface(dump_dir+"dump_%05d.h5" % d)
        # Get the luminosity at radius r
        eg0,deg0,L = dump.Lspec_at_r(0,args.radius,grid)
        Lspec0[dd,:] = np.sum(L*dOmega,0)/np.sum(dOmega)
        eg1,deg1,L = dump.Lspec_at_r(1,args.radius,grid)
        Lspec1[dd,:] = np.sum(L*dOmega,0)/np.sum(dOmega)
        eg2,deg2,L = dump.Lspec_at_r(2,args.radius,grid)
        Lspec2[dd,:] = np.sum(L*dOmega,0)/np.sum(dOmega)
        time[dd] = dump.get_dataset('Time')[0] - tbounce
        dump.close()
        print 'd={0:3d}, t={1:5.3f}'.format(d,time[dd])
        dd += 1
    np.savez(datafile,time=time,eg0=eg0,eg1=eg1,eg2=eg2,Lspec0=Lspec0,Lspec1=Lspec1,Lspec2=Lspec2)
else:
    # Load data from file
    data = np.load(datafile)
    time = data['time']
    eg0 = data['eg0']
    eg1 = data['eg1']
    eg2 = data['eg2']
    Lspec0   = data['Lspec0']
    Lspec1   = data['Lspec1']
    Lspec2   = data['Lspec2']

# Take the mean over a window of width dt_mean
#Lspec0,t = fornax.tmean(Lspec0,time,args.dt)
#Lspec1,t = fornax.tmean(Lspec1,time,args.dt)
#Lspec2,t = fornax.tmean(Lspec2,time,args.dt)
#time = t

# Reset and load plotting style
rcdefaults() #bug with animation garbling
#set_style()

# Plot & save figure
#plt.figure()
#plt.plot(eg0,Lspec0[-1],'b-',label=r'$\nu_e$')
#plt.plot(eg1,Lspec1[-1],'g-',label=r'$\nu_\bar{e}$')
#plt.plot(eg2,Lspec2[-1],'r-',label=r'$\nu_\mu$')
#plt.xlim(args.xmin,args.xmax)
#plt.ylim(args.ymin,args.ymax)
#plt.xlabel('Energy [MeV]')
#plt.ylabel(r'Spectra [10$^{52}$ erg s$^{-1}$ MeV$^{-1}$]')
#plt.title(title_text)
#plt.legend(loc='upper right')
#savefile = 'rad_ccsn_{0:s}_Spectra{1:d}'.format(save_text,int(args.radius))
#plt.savefig(savefile+'.eps')
#plt.savefig(savefile+'.png')

fig = plt.figure()
ax = plt.axes(xlim=(0,60),ylim=(0,4))
line0, = ax.plot([], [], lw=2)
line1, = ax.plot([], [], lw=2)
line2, = ax.plot([], [], lw=2)
#time_text = ax.text(5, 3,'')
annotation=ax.annotate('',xy=(0.2,0.8),xycoords="axes fraction", size=15,color='k', fontweight='bold')
ax.annotate(title_text,xy=(0.7,0.3),xycoords="axes fraction", size=15,color='k', fontweight='bold')
line0.set_label(r'$\nu_e$')
line1.set_label(r'$\bar{\nu}_e$')
line2.set_label(r'"$\nu_\mu$"')
plt.legend(frameon=False,bbox_to_anchor=[0.9,0.9])

def init():
    line0.set_data([], [])
    line1.set_data([], [])
    line2.set_data([], [])
    #time_text.set_text('initial')
    annotation.set_text('initial')
    #annotation.set_animated(True)
    #ax.annotate('',xy=(0.2,0.8),xycoords="axes fraction", size=15,color='k', fontweight='bold')
    return line0, line1, line2,annotation#, time_text

def animate(i):
    line0.set_data(eg0,Lspec0[i])
    line1.set_data(eg1,Lspec1[i])
    line2.set_data(eg2,Lspec2[i])
    time_ann = '%0.3f %s' %(time[i],'s')
    #time_text.set_text(time_ann )
    annotation.set_text(time_ann)
    print 'creating frame %03d' %i
    #ax.annotate(time_ann,xy=(0.2,0.8),xycoords="axes fraction", size=15,color='k', fontweight='bold')
    return line0, line1, line2,annotation

plt.xlabel(r'Energy (MeV)')
plt.ylabel(r'$L_\epsilon\,(10^{52}\,\mathrm{erg/s/MeV})$')

anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=args.dumps_fin+args.dumps_pad+1, interval=150, blit=False)
writer = animation.writers['ffmpeg'](fps=args.fps)

anim.save('rad_ccsn_{0:s}_spectra_r{1:d}.mp4'.format(save_text,int(args.radius)),writer=writer, extra_args=['-vcodec', 'libx264'])

