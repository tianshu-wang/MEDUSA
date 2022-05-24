import numpy as np
import fornax
import matplotlib.pyplot as plt
from runinfo import *
from styles import *
from matplotlib import rcdefaults
import argparse

# Parse command line
parser = argparse.ArgumentParser(description="Plot the average energy")
parser.add_argument("--recompute",help="recompute values; otherwise, read from file [False]",action="store_true")
parser.add_argument("-r","--radius",help="radius in km at which to compute average energy [100.0]",type=float,default=100.0)
parser.add_argument("--xmin",help="xmin of plot window [0.0]",type=float,default=0.0)
parser.add_argument("--xmax",help="xmax of plot window [1.0]",type=float,default=1.0)
parser.add_argument("--ymin",help="ymin of plot window [0.0]",type=float,default=0.0)
parser.add_argument("--ymax",help="ymax of plot window [30.0]",type=float,default=30.0)
args = parser.parse_args()

# Open grid object
grid = fornax.GeomInterface(dump_dir+'grid.h5')
datafile = 'avenergy_r{:d}.npz'.format(int(args.radius))

if args.recompute:
    # Recompute and save data to file
    N        = dump_stop - dump_start + 1
    time     = np.zeros(N)
    Frad0_angavg = np.zeros(len(range(grid.nr1)))
    Frad1_angavg = np.zeros(len(range(grid.nr2)))
    Frad2_angavg = np.zeros(len(range(grid.nr3)))
    nres = len(range(grid.nr3))
    avenergyrms0 = np.zeros(N)
    avenergyrms1 = np.zeros(N)
    avenergyrms2 = np.zeros(N)
    avenergy0 = np.zeros(N)
    avenergy1 = np.zeros(N)
    avenergy2 = np.zeros(N)
    dOmega = grid.dOmega()
    ir = np.where(grid.r < args.radius)[0][-1]
    idx = ["%.2d" % i for i in range(grid.nr1)]

    dd = 0
    for d in range(dump_start,dump_stop+1):
        dump = fornax.DataInterface(dump_dir+"dump_%05d.h5" % d)

        for i in idx:

          Frad0=dump.get_dataset3d('Frad0/dir0/''%s%02s' %('g',i) ,grid)[ir]
          Frad0_angavg[int(i)]=np.sum(Frad0*dOmega)/np.sum(dOmega)

          Frad1=dump.get_dataset3d('Frad1/dir0/''%s%02s' %('g',i) ,grid)[ir]
          Frad1_angavg[int(i)]=np.sum(Frad1*dOmega)/np.sum(dOmega)

          Frad2=dump.get_dataset3d('Frad2/dir0/''%s%02s' %('g',i) ,grid)[ir]
          Frad2_angavg[int(i)]=np.sum(Frad2*dOmega)/np.sum(dOmega)

        avenergy0[dd]=np.sum(Frad0_angavg)/np.sum(Frad0_angavg/grid.egroup[:nres])
        avenergy1[dd]=np.sum(Frad1_angavg)/np.sum(Frad1_angavg/grid.egroup[nres:2*nres])
        avenergy2[dd]=np.sum(Frad2_angavg)/np.sum(Frad2_angavg/grid.egroup[2*nres:])

        avenergyrms0[dd]=(np.sum(Frad0_angavg*grid.egroup[:nres]**2)/np.sum(Frad0_angavg))**0.5
        avenergyrms1[dd]=(np.sum(Frad1_angavg*grid.egroup[nres:2*nres]**2)/np.sum(Frad1_angavg))**0.5
        avenergyrms2[dd]=(np.sum(Frad2_angavg*grid.egroup[2*nres:]**2)/np.sum(Frad2_angavg))**0.5

        time[dd] = dump.get_dataset('Time')[0] - tbounce
        dump.close()
        print 'd={0:3d}, t={1:5.3f}'.format(d,time[dd])
        dd += 1
    np.savez(datafile,time=time,avenergy0=avenergy0,avenergy1=avenergy1,avenergy2=avenergy2,avenergyrms0=avenergyrms0,avenergyrms1=avenergyrms1,avenergyrms2=avenergyrms2)

else:
    # Load data from file
    data = np.load(datafile)
    time = data['time']
    avenergy0 = data['avenergy0']
    avenergy1 = data['avenergy1']
    avenergy2 = data['avenergy2']
    avenergyrms0 = data['avenergyrms0']
    avenergyrms1 = data['avenergyrms1']
    avenergyrms2 = data['avenergyrms2']

rcdefaults()
set_style()

#Plot & save figure
plt.figure()
plt.plot(time,avenergy0,'b-',label=r'$\nu_e$')
plt.plot(time,avenergy1,'g-',label=r'$\nu_\bar{e}$')
plt.plot(time,avenergy2,'r-',label=r'$\nu_\mu$')
plt.xlim(args.xmin,args.xmax)
plt.ylim(args.ymin,args.ymax)
plt.xlabel('Time after bounce [s]')
plt.ylabel(r'Average Energy [MeV]')
plt.title(title_text)
plt.legend(loc='upper right')
savefile = 'rad_ccsn_{0:s}_avenergy{1:d}'.format(save_text,int(args.radius))
plt.savefig(savefile+'.eps')
plt.savefig(savefile+'.png')

plt.figure()
plt.plot(time,avenergyrms0,'b-',label=r'$\nu_e$')
plt.plot(time,avenergyrms1,'g-',label=r'$\nu_\bar{e}$')
plt.plot(time,avenergyrms2,'r-',label=r'$\nu_\mu$')
plt.xlim(args.xmin,args.xmax)
plt.ylim(args.ymin,args.ymax)
plt.xlabel('Time after bounce [s]')
plt.ylabel(r'rms Energy [MeV]')
plt.title(title_text)
plt.legend(loc='upper right')
savefile = 'rad_ccsn_{0:s}_avenergyrms{1:d}'.format(save_text,int(args.radius))
plt.savefig(savefile+'.eps')
plt.savefig(savefile+'.png')

print 'done'
