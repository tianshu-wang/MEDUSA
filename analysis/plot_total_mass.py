import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('total_mass.txt')
t = data[:,0]/1.0e-3
Mtot = data[:,1]/1.99e33

plt.plot(t,Mtot)
plt.xlabel(r'$t$ [ms]')
plt.ylabel(r'$M_\mathrm{tot}$ [$M_\odot$]')
plt.ylim((0,2))
plt.savefig('total_mass.png')

print 'done'
