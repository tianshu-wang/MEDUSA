#!/usr/bin/env python
import numpy as np
import fornax
import matplotlib.pyplot as plt
from matplotlib import rcdefaults
from styles import *

dump_dir = '../dumps/'
dump_start = 0
dump_stop = 4

U_KM  = 1.0e5
U_DEG = np.pi/180.0

# Open grid object
grid = fornax.GeomInterface(dump_dir+'grid.h5')
x = 0.5*(grid.x[1:] + grid.x[:-1])

dump = fornax.DataInterface(dump_dir+"dump_%05d.h5" % dump_stop)
rho = dump.get_dataset('rho')
P = dump.get_dataset('eos0')
v = dump.get_dataset('u1')
shock_flag = dump.get_dataset('shock_flag')
dump.close()

# Reset and load plotting style
rcdefaults()
set_style()

# Plot & save figure
plt.figure()
plt.plot(x,rho,'kx')
plt.xlabel(r'$x$')
plt.ylabel(r'$\rho$')
plt.savefig('sod_rho.png')

plt.figure()
plt.plot(x,P,'kx')
plt.xlabel(r'$x$')
plt.ylabel(r'$P$')
plt.savefig('sod_P.png')

plt.figure()
plt.plot(x,P/rho,'kx')
plt.xlabel(r'$x$')
plt.ylabel(r'$P/\rho$')
plt.savefig('sod_u.png')

plt.figure()
plt.plot(x,v,'kx')
plt.xlabel(r'$x$')
plt.ylabel(r'$v$')
plt.savefig('sod_v.png')

plt.figure()
plt.plot(x,shock_flag,'kx')
plt.xlabel(r'$x$')
plt.ylabel(r'Shock Flag')

print 'done'
plt.show()
