import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('total_energy.txt')
t = data[:,0]/1.0e-3
Etot = data[:,1]
Ugrav = data[:,2]
total_energy = data[:,3]

plt.plot(t,Etot/1.0e53,label=r'$E_\mathrm{tot}$')
plt.plot(t,Ugrav/1.0e53,label=r'$U_\mathrm{grav}$')
plt.plot(t,total_energy/1.0e53,label=r'$E_\mathrm{tot}+U_\mathrm{grav}$')
plt.xlabel(r'$t$ [ms]')
plt.ylabel(r'Energy [$10^{53}$ erg]')
plt.legend(loc='best')
plt.savefig('energy_components.png')

plt.figure()
plt.plot(t,total_energy/1.0e49)
plt.xlabel(r'$t$ [ms]')
plt.ylabel(r'$E_\mathrm{tot}$ [$10^{49}$ erg]')
plt.savefig('total_energy.png')

print 'done'
plt.show()
