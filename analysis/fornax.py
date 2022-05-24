import h5py
import numpy as np
from scipy.integrate import romberg
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

U_KM    = 1.0e5	
U_C     = 2.99792458e10
U_MEV   = 1.60217657e-6
U_NA    = 6.0221409e+23
U_MSUN  = 1.9891e33
U_G     = 6.67384e-8
U_KB    = 1.3806488e-16
U_MP    = 1.67262158e-24
U_K2MEV = 8.617342337686e-11

one_six = 1./6.
cdict =    {    'red':    ((0.0, 0.0, 0.0),
                     (one_six, 0.0, 0.0),
                     (2*one_six, 50./255., 50./255.),
                     #(0.5, 240./255., 240./255.),
                     (0.5, 200./255., 200./255.),
                     (4*one_six, 255./255., 255./255.),
                     (5*one_six, 204./255., 204./255.),
                     (1., 127./255., 127./255.)),
            'green':((0.0, 0.0, 0.0),
                     (one_six, 68./255., 68./255.),
                     (2*one_six, 186./255., 186./255.),
                     (0.5, 240./255., 240./255.),
                     (4*one_six, 170./255., 170./255.),
                     (5*one_six, 67./255, 67./255.),
                     (1., 0.0, 0.0)),
            'blue':    ((0., 127./255., 127/255.),
                     (one_six, 204./255., 204./255.),
                     #(2*one_six, 150./255., 150./255.),
                     (2*one_six, 1., 1.),
                     #(0.5, 240./255., 240./255.),
                     (0.5, 200./255., 200./255.),
                     (4*one_six, 0.0, 0.0),
                     (5*one_six, 0.0, 0.0),
                     (1., 0.0, 0.0))
        }

#cdict['alpha'] = ((0.0, 0.0, 0.0),
#                     (0.25, 1./16., 1./16.),
#                     (0.5, 0.25, 0.25),
#                     (0.75, 9./16., 9./16.),
#                     (1., 1., 1.))

cdict['alpha'] = ((0.,1.,1.),(1.,1.,1.))

#fornax_cmap = colors.LinearSegmentedColormap('fornax_cmap', cdict)
plt.register_cmap(name='fornax_cmap', data=cdict)
cmap = plt.get_cmap('fornax_cmap')
#plt.set_cmap(cmap)

def tmean(x,t,dt_mean):
    """Compute the mean over a window of width dt_mean"""
    dt = np.diff(t)
    dt = np.append(dt,dt[-1])
    imin = np.where(t>=t.min()+0.5*dt_mean)[0][0]
    imax = np.where(t<=t.max()-0.5*dt_mean)[0][-1]
    tt = t[imin:imax+1]
    xx = np.zeros(imax-imin+1)
    i = 0
    for time in tt:
        imin = np.where(t>=time-0.5*dt_mean)[0][0]
        imax = np.where(t<=time+0.5*dt_mean)[0][-1]
        mask = slice(imin,imax+1)
        xx[i] = np.sum(x[mask]*dt[mask])/np.sum(dt[mask])
        i += 1
    return xx,tt

class DataInterface:

    def __init__(self, name):
        self.f = h5py.File(name,"r")
        self.d = {}        

    def close(self):
        self.f.close()

    def get_dataset(self,name):
        self.d[name] = self.f[name].value
        return self.d[name]

    def get_dataset3d(self,name,geom):
        self.d[name] = self.f[name].value
        return self.d[name].reshape(geom.nr,geom.nth,geom.nphi)

    def get_var(self,name,geom):
        if name in ('rho','u1','u2','u3','u','comp0','eos0','eos1','eos2','eos3','velocity0','velocity1'):
            return self.get_dataset3d(name,geom)
        elif name=='Er0':
            return self.get_dataset3d('Erad0/total',geom)
        elif name=='Fr0':
            return self.get_dataset3d('Frad0/dir0/total',geom)
        elif name=='Er1':
            return self.get_dataset3d('Erad1/total',geom)
        elif name=='Fr1':
            return self.get_dataset3d('Frad1/dir0/total',geom)
        elif name=='Er2':
            return self.get_dataset3d('Erad2/total',geom)
        elif name=='Fr2':
            return self.get_dataset3d('Frad2/dir0/total',geom)
        elif name=='Frtot':
            return self.get_var('Fr0',geom) + self.get_var('Fr1',geom)
        elif name=='L0':
            return self.get_var('Fr0',geom)*geom.dA()
        elif name=='L1':
            return self.get_var('Fr1',geom)*geom.dA()
        elif name=='Ltot':
            return self.get_var('L0',geom) + self.get_var('L1',geom)
        elif name=='heat0':
            return self.get_dataset3d('Erad0/heat',geom)
        elif name=='heat1':
            return self.get_dataset3d('Erad1/heat',geom)
        elif name=='heat2':
            return self.get_dataset3d('Erad2/heat',geom)
        elif name=='cool0':
            return self.get_dataset3d('Erad0/cool',geom)
        elif name=='cool1':
            return self.get_dataset3d('Erad1/cool',geom)
        elif name=='cool2':
            return self.get_dataset3d('Erad2/cool',geom)
        elif name=='P':
            return self.get_var('eos0',geom)
        elif name=='cs':
            return self.get_var('eos1',geom)
        elif name=='T':
            return self.get_var('eos2',geom)
        elif name=='S':
            return self.get_var('eos3',geom)
        elif name=='Gamma1':
            return self.get_var('eos4',geom)
        elif name=='Ye':
            return self.get_var('comp0',geom)
        else:
            raise ValueError('Unknown variable name!')

    def get_Mass(self,geom):
        rho = self.get_dataset3d('rho',geom)
        return np.cumsum(np.sum(np.sum(rho*geom.dphi,axis=2)*geom.dmu,axis=1)*np.diff((geom.r*U_KM)**3)/3.0,axis=0)/U_MSUN

    def get_mean_4pi(self,name,geom):
        dOmega = geom.dOmega()
        var = self.get_var(name,geom)
        return np.sum(np.sum(var*dOmega[np.newaxis,:,:],axis=2),axis=1)/np.sum(dOmega)

    def Espec_at_r(self,species,r,geom):
        ir = np.where(geom.r < r)[0][-1]
        dr = (r-geom.r[ir])/(geom.r[ir+1]-geom.r[ir])
        if geom.NDIM > 1:
            print "NOT YET IMPLEMENTED"
            return None, None

        if species==0:
            e0 = 0
            e1 = geom.nr1
        elif species==1:
            e0 = geom.nr1
            e1 = e0+geom.nr2
        else :
            e0 = geom.nr1+geom.nr2
            e1 = e0+geom.nr3

        ng = e1-e0

        E = np.empty(ng,dtype=np.double)

        for i in range(ng):
            name = "Erad%d/g%02d" % (species, i)
            E[i] = (1-dr)*self.f[name][ir] + dr*self.f[name][ir+1]

        return geom.egroup[e0:e1], E

    def Fspec_at_r(self,species,r,geom):
        ir = np.where(geom.r < r)[0][-1]
        dr = (r-geom.r[ir])/(geom.r[ir+1]-geom.r[ir])
        if geom.NDIM > 2:
            print "NOT YET IMPLEMENTED"
            return None, None

        if species==0:
            e0 = 0
            e1 = geom.nr1
        elif species==1:
            e0 = geom.nr1
            e1 = e0+geom.nr2
        else :
            e0 = geom.nr1+geom.nr2
            e1 = e0+geom.nr3

        ng = e1-e0
        if   geom.NDIM == 1:
            F = np.empty(ng,dtype=np.double)
        elif geom.NDIM == 2:
            F = np.empty((geom.nth,ng),dtype=np.double)

        for i in range(ng):
            name = "Frad%d/dir0/g%02d" % (species, i)
            if   geom.NDIM == 1:
                F[i]   = (1-dr)*self.f[name][ir] + dr*self.f[name][ir+1]
            elif geom.NDIM == 2:
                F[:,i] = (1-dr)*self.f[name][ir,:] + dr*self.f[name][ir+1,:]

        return geom.egroup[e0:e1], geom.degroup[e0:e1]/U_MEV, F

    def F_at_r(self,species,r,geom):
        ir = np.where(geom.r < r)[0][-1]
        dr = (r-geom.r[ir])/(geom.r[ir+1]-geom.r[ir])
        if geom.NDIM > 2:
            print "NOT YET IMPLEMENTED"
            return None

        name = "Fr%d" % (species)
        F = self.get_var(name,geom)

        return (1-dr)*F[ir,:,:] + dr*F[ir+1,:,:]

    def Lspec_at_r(self,species,r,geom):
        eg, deg, L = self.Fspec_at_r(species,r,geom)
        L *= 4.0*np.pi*(r*U_KM)**2/1.0e52

        return eg, deg, L

    def L_at_r(self,species,r,geom):
        F = self.F_at_r(species,r,geom)  # flux [erg cm^-2 s^-1]
        L = F*geom.dA_at_r(r)  # luminosity [erg s^-1]

        return L

    def get_Mdot_at_r(self,r,geom):
        ir = np.where(geom.r < r)[0][-1]
        dr = (r-geom.r[ir])/(geom.r[ir+1]-geom.r[ir])

        rho = self.get_dataset3d('rho',geom)
        vr  = -self.get_dataset3d('u1',geom)
#        Mdot = np.sum(np.sum(rho*vr*geom.dphi,axis=2)*geom.dmu,axis=1)*(geom.rc*1.0e5)**2
        Mdot = (1-dr)*(rho*vr)[ir,:,:] + dr*(rho*vr)[ir+1,:,:]
        return np.sum(np.sum(Mdot*geom.dA_at_r(r),axis=1),axis=0)

    def netheat_at_r(self,species,r,geom):
        ir = np.where(geom.r < r)[0][-1]
        dr = (r-geom.r[ir])/(geom.r[ir+1]-geom.r[ir])

        hname = 'Erad%d/heat' % species
        cname = 'Erad%d/cool' % species
        if   geom.NDIM == 1:
            heat = (1-dr)*self.f[hname][ir    ]*self.f['rho'][ir    ] + dr*self.f[hname][ir+1    ]*self.f['rho'][ir+1    ]
            cool = (1-dr)*self.f[cname][ir    ]*self.f['rho'][ir    ] + dr*self.f[cname][ir+1    ]*self.f['rho'][ir+1    ]
        elif geom.NDIM == 2:
            heat = (1-dr)*self.f[hname][ir,:  ]*self.f['rho'][ir,:  ] + dr*self.f[hname][ir+1,:  ]*self.f['rho'][ir+1,:  ]
            cool = (1-dr)*self.f[cname][ir,:  ]*self.f['rho'][ir,:  ] + dr*self.f[cname][ir+1,:  ]*self.f['rho'][ir+1,:  ]
        elif geom.NDIM == 3:
            heat = (1-dr)*self.f[hname][ir,:,:]*self.f['rho'][ir,:,:] + dr*self.f[hname][ir+1,:,:]*self.f['rho'][ir+1,:,:]
            cool = (1-dr)*self.f[cname][ir,:,:]*self.f['rho'][ir,:,:] + dr*self.f[cname][ir+1,:,:]*self.f['rho'][ir+1,:,:]

        return heat - cool

    def get_Qdot(self,geom):
        rho = self.get_dataset3d('rho',geom)
        S   = self.get_dataset3d('eos3',geom)
        T   = self.get_dataset3d('eos2',geom)
        q0  = self.get_dataset3d('Erad0/heat',geom) - self.get_dataset3d('Erad0/cool',geom)
        q1  = self.get_dataset3d('Erad1/heat',geom) - self.get_dataset3d('Erad1/cool',geom)
#        v1  = self.get_dataset3d('u1',geom)
#        v2  = self.get_dataset3d('u2',geom)
#        v3  = self.get_dataset3d('u3',geom)
#        Phi = -U_G*self.get_Mass(geom)/geom.rc**2  # This could be improved...
#        Phi = np.tile(Phi[:,np.newaxis,np.newaxis],(1,geom.nth,geom.nphi))
        e   = self.get_dataset3d('u',geom)
#        e   = 1.5*(rho/mp)*kB*(T/K_to_MeV)
        q   = rho*(q0+q1)
#        E   = 0.5*rho*(v1**2+v2**2+v3**2) + e + rho*Phi
        dV  = geom.dV()
        Qdot  = 0.0
        Mgain = 0.0
        egain = 0.0
        Sgain = 0.0
        for k in range(geom.nphi):
            for j in range(geom.nth):
                rshock = findshock(geom.rc,S[:,j,k])
                mask   = findgain(geom.rc,rshock,q[:,j,k])
                if (mask.size>0):
                    Qdot  += np.sum((q*dV)[mask,j,k])
                    Mgain += np.sum((rho*dV)[mask,j,k])
#                    Etotg += np.sum(   E[mask,j,k]*dV[mask,j,k])
                    egain += np.sum((e*dV)[mask,j,k])
                    Sgain += np.sum((rho*S*dV)[mask,j,k])
        return Qdot,Mgain,egain,Sgain
        
    def slice_dataset_at_r_orig(self,name,r,geom):
        ir = np.where(geom.rc < r)[0][-1]
        dr = (r-geom.rc[ir])/(geom.rc[ir+1]-geom.rc[ir])
        if geom.NDIM==1:
            return (1-dr)*self.f[name][ir] + dr*self.f[name][ir+1]
        elif geom.NDIM==2:
            return (1-dr)*self.f[name][ir,:]+dr*self.f[name][ir+1,:]
        elif geom.NDIM==3:
            return np.meshgrid(geom.th, geom.phi, indexing='ij'), self.f[name][ir,:,:]

    def slice_dataset_at_r_damn(self,name,r,geom):
        for i in range(geom.rc.shape[0]):
            if geom.rc[i] > r:
                ir = i
                break
        rg = geom.rc[ir-2:ir+2]
        fg = self.f[name][ir-2:ir+2]
        rs,fs = mspline(rg,fg,50)
        for i in range(50):
            if rs[i] > r:
                ir = i
                break
        dr = (r-rs[ir-1])/(rs[ir]-rs[ir-1])
        return (1-dr)*fs[ir-1] + dr*fs[ir]

    def slice_dataset_at_r(self,name,r,geom,plaw=False):
        ir = np.where(geom.rc <= r)[0][-1]
        dr = (r-geom.rc[ir])/(geom.rc[ir+1]-geom.rc[ir])
        if plaw==False:
            return (1-dr)*self.f[name][ir] + dr*self.f[name][ir+1]
        alpha = np.log(self.f[name][ir+1]/self.f[name][ir])/np.log(geom.rc[ir+1]/geom.rc[ir])
        return self.f[name][ir]*np.power(r/geom.rc[ir],alpha)

    def slice_dataset_at_th(self,name,th,geom):
        ith = np.where(geom.th[:-1] <= th)[0][-1]
        if geom.NDIM==2:
            return self.f[name][:,ith]
        if geom.NDIM==3:
            return self.f[name][:,ith,:]

    def slice_dataset_at_ph(self,name,phi,geom):
        iph_temp = np.where(geom.phi[:-1] <= phi)
        if iph_temp[0].shape[0] == 0:
            iph = 0
        else:
            iph = iph_temp[0][-1]
        return np.outer(geom.r,np.sin(geom.th)), np.outer(geom.r, np.cos(geom.th)), self.f[name][:,:,iph]

    def gradient(self,name,geom,normalize=False):
        if geom.NDIM==3:
            grad = np.empty(self.d[name].shape + (3,), dtype=float)
        else:
            grad = np.empty(self.d[name].shape + (128,3), dtype=float)

        #rhat = np.empty(self.d[name].shape, dtype=float)
        #thhat = np.empty(self.d[name].shape, dtype=float)
        #phhat = np.empty(self.d[name].shape, dtype=float)
        for i in range(grad.shape[0]):
            im = max(0,i-1)
            ip = min(grad.shape[0]-1,i+1)
            for j in range(grad.shape[1]):
                jm = max(0,j-1)
                jp = min(grad.shape[1]-1,j+1)
                sth = np.sin(geom.thc[j])
                cth = np.cos(geom.thc[j])
                if geom.NDIM == 3:
                    for k in range(grad.shape[2]):
                        km = max(0,k-1)
                        kp = min(grad.shape[2]-1,k+1)
                        sph = np.sin(geom.phic[k])
                        cph = np.cos(gepm.phic[k])
                        rhat = np.array([sth*cph, sth*sph, cth])
                        thhat = np.array([cth*cph, cth*sph, -sth])
                        phhat = np.array([-sph, cph, 0])
                        grad[i,j,k] = rhat*(self.d[name][ip,j,k] - self.d[name][im,j,k])/(geom.rc[ip] - geom.rc[im])
                        grad[i,j,k] = thhat*(self.d[name][i,jp,k] - self.d[name][i,jm,k])/(geom.rc[i]*(geom.thc[jp] - geom.thc[jm]))
                        grad[i,j,k] = phhat*(self.d[name][i,j,kp] - self.d[name][i,j,km])/(geom.rc[i]*sth*(geom.phic[kp] - geom.phic[km]))    
                else:
                    for k in range(128):
                        sph = np.sin(2*np.pi*(k+0.5)/128.)
                        cph = np.cos(2*np.pi*(k+0.5)/128.)
                        rhat = np.array([sth*cph, sth*sph, cth])
                        thhat = np.array([cth*cph, cth*sph, -sth])
                        grad[i,j,k] = rhat*(self.d[name][ip,j] - self.d[name][im,j])/(geom.rc[ip] - geom.rc[im])
                        grad[i,j,k] = thhat*(self.d[name][i,jp] - self.d[name][i,jm])/(geom.rc[i]*(geom.thc[jp] - geom.thc[jm]))

        if normalize:
            mag = np.sqrt(np.sum(np.power(grad,2),axis=3))
            grad[:,:,:,0] /= mag
            grad[:,:,:,1] /= mag
            grad[:,:,:,2] /= mag
            #max_grad = np.amax(np.sqrt(np.sum(np.power(grad,2),axis=3)))
            #grad /= max_grad

        return grad
        

class GeomInterface:

    def __init__(self, name):
        self.f = h5py.File(name,"r")
        self.x = self.f["X"].value
        self.NDIM = 1
        if "Y" in self.f:
            self.y = self.f["Y"].value
            self.NDIM = 2
        if "Z" in self.f:
            self.z = self.f["Z"].value
            self.NDIM = 3

        # Set up the 1-d coordinate arrays
        if self.NDIM==1:
            self.r    = self.x
            self.th   = np.array((0,np.pi))
            self.phi  = np.array((0,2*np.pi))
        if self.NDIM==2:
            self.r    = np.sqrt(self.x[:,0]**2 + self.y[:,0]**2)
            self.th   = np.arctan2(self.y[:,0],self.x[:,0])
            self.phi  = np.array((0,2*np.pi))
        if self.NDIM==3:
            self.r    = np.sqrt(self.x[:,0,0]**2 + self.y[:,0,0]**2 + self.z[:,0,0]**2)
            Rcyl      = np.sqrt(self.x[:,0,0]**2 + self.y[:,0,0]**2)
            self.th   = np.arctan2(Rcyl,self.z[:,0,0])
            self.phi  = np.arctan2(self.y[:,0,0],self.x[:,0,0])

        self.rc   = 0.5*(self.r[1:] + self.r[:-1])
        self.dr   = np.diff(self.r)
        self.nr   = len(self.rc)

        self.thc  = 0.5*(self.th[1:] + self.th[:-1])
        self.dth  = np.diff(self.th)
        self.nth  = len(self.thc)

        self.phic = 0.5*(self.phi[1:] + self.phi[:-1])
        self.dphi = np.diff(self.phi)
        self.nphi = len(self.phic)

        self.dmu = -np.diff(np.cos(self.th))

        if "nr1" in self.f:
            self.nr1 = self.f["nr1"].value
            self.egroup = self.f["egroup"].value
            self.degroup = self.f["degroup"].value
        if "nr2" in self.f:
            self.nr2 = self.f["nr2"].value
            self.nr3 = self.f["nr3"].value

    # The following are always 3-d arrays, even in 1-d and 2-d
    def dV(self):
        r2dr = np.diff((self.r*U_KM)**3)/3.0
        return r2dr[:,np.newaxis,np.newaxis]*self.dOmega()

    def dOmega(self):
        return self.dmu[np.newaxis,:,np.newaxis]*self.dphi[np.newaxis,np.newaxis,:]

    def dA(self):
        return (self.rc[:,np.newaxis,np.newaxis]*U_KM)**2*self.dOmega()

    def dA_at_r(self,r):
        return (r*U_KM)**2*self.dOmega()

    def get_mask(self):
        return self.f["mask"].value

    def get_scale(self,d):
        return np.sqrt(self.f["gcov%d" % d].value)

    def ijk_given_x(self, x):
        if self.NDIM==1:
            r = x
        if self.NDIM==2:
            r  = np.sqrt(x[0]**2 + x[1]**2)
            th = np.arctan2(x[1],x[0])
        elif self.NDIM==3:
            r   = np.sqrt(x[0]**2 + x[1]**2 + x[2]**2)
            R   = np.sqrt(x[0]**2 + x[1]**2)
            th  = np.arctan2(R/x[2])
            phi = np.arctan2(x[1],x[0])
            iphi = phi/self.dphi
        ir = np.where(self.r < r)[0][-1]
        if self.NDIM >= 2:
            ith = np.where(self.th < th)[0][-1]

        if self.NDIM==1:
            ind = (ir,)
        elif self.NDIM==2:
            ind = (ir,ith)
        elif self.NDIM==3:
            ind = (ir,ith,iphi)

        return ind


class VolRend:

    def __init__(self,data,grad,geom,rclip,nsamp,colormap,crange):
        self.data = data
        self.grad = grad
        self.grid = geom
        self.rclip = rclip
        self.nsamp = nsamp
        self.cmin = crange[0]
        self.cmax = crange[1]
        th0 = 30.*np.pi/180.
        phi0 = 60.*np.pi/180.
        self.light = np.array([np.cos(th0)*np.cos(phi0), np.cos(th0)*np.sin(phi0), -np.sin(th0)])
        theCM = colormap
        theCM._init()
        #alphas = np.zeros(theCM.N)
        #ia = np.arange(0,theCM.N)
        #c = [0, 64, 128, 192, 255]
        #m = [0.5, 0.7, 0.9, 1., 0.9]
        #wsq = [64,64,64,64,64]
        #for i in range(len(c)):
        #    alphas += m[i]*np.exp(-np.power(ia-c[i],2)/wsq[i])
        alphas = np.linspace(0,1.,theCM.N)
        betas = np.linspace(1.,0.5,theCM.N)
        alphas = np.minimum(alphas,betas)
        
        #print theCM.N
        #alphas[0:theCM.N/2] = 0.
        theCM._lut[:-3,-1] = alphas
        cNorm = colors.Normalize(vmin=crange[0], vmax=crange[1], clip=True)
        self.scalarMap = cmx.ScalarMappable(cNorm, cmap=theCM)
        #cmap.set_under()
        print self.scalarMap.to_rgba(crange[0]/2.), self.scalarMap.to_rgba(crange[0]), self.scalarMap.to_rgba(crange[1])

    def build_ray(self,origin,n):
        """
        n is assumed to be a unit vector pointing away from the data
        """
        r = np.sqrt(np.sum(np.power(origin,2)))
        dr = self.rclip - r
        disc = np.power(np.sum(n*origin),2)
        disc = disc + dr*(dr + 2*r)
        if disc <= 0:
            return np.array([]), np.array([]), np.array([]), np.array([])
        l0 = -np.sum(n*origin) + np.sqrt(disc)
        l1 = -np.sum(n*origin) - np.sqrt(disc)
        x0 = origin + l0*n
        x1 = origin + l1*n

        dx = (x1-x0)/(self.nsamp-1.)

        xc = x0 + np.outer(np.arange(self.nsamp),dx)
        rc = np.sqrt(np.sum(xc*xc,1))
        thc = np.arccos(xc[:,2]/rc)
        phc = np.fmod(np.arctan2(xc[:,1],xc[:,0])+2*np.pi,2*np.pi)

        ic = np.abs(np.subtract.outer(self.grid.rc,rc)).argmin(0)
        jc = np.abs(np.subtract.outer(self.grid.thc,thc)).argmin(0)
        #if self.grid.NDIM==3:
        kc = np.abs(np.subtract.outer(self.grid.phic,phc)).argmin(0)
        #else:
        #    kc = np.array([])

        return ic, jc, kc, xc

    def integrate_pixel(self, vals, gvals, x):

        # black background
        rgb = np.array([0., 0., 0.])
        if len(vals) < 2:
            return rgb

        dl = np.sqrt(np.power(x[1:]-x[:-1],2))
        ltot = np.sum(dl)
        #dl = self.nsamp/5.*dl/ltot
        #dl = np.sqrt(np.sum(np.power(x[-1] - x[0],2)))
        wgt = 1.    #ltot/(2*self.rclip)


        alpha_prod = 1.;
        c0 = np.array(self.scalarMap.to_rgba(vals[0]))
        l = (x[1] - x[0])/dl[0]
        ldotg0 = np.dot(self.light,gvals[0])

        for i in range(len(vals)-1):
            #c0 = np.array(self.scalarMap.to_rgba(vals[i]))
            c1 = np.array(self.scalarMap.to_rgba(vals[i+1]))
            l = (x[i+1] - x[i])/dl[i]
            ldotg1 = np.dot(self.light,gvals[i+1])
            ldotg = 0.5*(ldotg0 + ldotg1)
            #print l, gvals
            alpha = 0.5*(c0[3] + c1[3])*(10./self.nsamp)*wgt*0.5*(1+ldotg)
            col = 0.5*(c0[0:3] + c1[0:3])*alpha
            c0 = np.copy(c1)
            ldotg0 = ldotg1
            #print "alpha = ", alpha, vals[i]
            #alpha_prod = alpha_prod*(1 - alpha)
            rgb = rgb*(1 - alpha) + col
            #rgb = rgb + col*alpha_prod

        return rgb

    def tau_func(self,dtau):
        return np.exp(dtau)
        

    def tau_inv_func(self,dtau):
        return np.exp(-dtau)
        #etau = (np.power(tau-2,2)-1)/3.
        #etau = np.maximum(etau,0.)
        #return etau

    def opac(self,val):
        b0 = 60.
        b1 = 1.e-6
        g0 = 1.e-6
        g1 = 60.
        g2 = 1.e-6
        r0 = 1.e-6
        r1 = 6.
        if val < self.cmin:
            return np.array([1.e-6,1.e-6,1.e-6])
        if val > self.cmax:
            return np.array([1.e-6,1.e-6,1.e-6])

        cmid = 0.5*(self.cmin+self.cmax)
        if val < cmid:
            bd = (val - self.cmin)/(cmid-self.cmin)
            b = (1-bd)*b0 + bd*b1
            gd = bd
            g = (1-gd)*g0 + gd*g1
            r = r0
        else:
            b = b1
            gd = (val - cmid)/(self.cmax - cmid)
            g = (1-gd)*g1 + gd*g2
            rd = gd
            r = (1-rd)*r0 + rd*r1
        
        return np.array([g,g,g])

    def integrate_pixel_huh(self, vals, gvals, x):

        # black background
        rgb = np.array([0., 0., 0.])
        if len(vals) < 2:
            return rgb

        tau = np.array([0., 0., 0.])

        dl = np.sqrt(np.power(x[1:]-x[:-1],2))
        dl /= (2*self.rclip)
        #ltot = np.sum(dl)
        #dl /= ltot

        c0 = np.array(self.scalarMap.to_rgba(vals[0]))[0:3]
        a0 = self.opac(vals[0])
        for i in range(len(vals)-1):
            c1 = np.array(self.scalarMap.to_rgba(vals[i+1]))[0:3]
            a1 = self.opac(vals[i+1])
            dtau = dl[i]*0.5*(a0 + a1)
            f = (1 - 1/dtau)*c0 + c1/dtau
            g = c0/dtau - c1*(1 + 1/dtau)
            tau += dtau

            for j in range(3):
                if dtau[j] < 1e-6:
                    rgb[j] += 0.5*(2*rgb[j] + c0[j] + c1[j])*dtau[j] + (3*rgb[j] + 2*c0[j] + c1[j])*dtau[j]*dtau[j]/6.
                else :
                    etau = self.tau_func(dtau[j])
                    rgb[j] = (rgb[j] + f[j])*etau + g[j]

            if np.all(tau > 5):
                break

            c0 = c1
            a0 = a1

        rgb *= self.tau_inv_func(tau)

        return rgb

        

    def build_image(self, origin, npix):
        image = np.empty([npix,npix,3],dtype=float)
        r0 = np.sqrt(np.sum(origin*origin))
        th0 = np.arccos(origin[2]/r0)
        phi0 = np.fmod(np.arctan2(origin[1],origin[0])+2*np.pi,2*np.pi)
        n2 = origin/r0    # radial unit vector, points along z in local coordinates
        n0 = np.array([np.cos(th0)*np.cos(phi0), np.cos(th0)*np.sin(phi0), -np.sin(th0)])    # theta unit vector, points along x in local coordinates
        n1 = np.array([-np.sin(phi0), np.cos(phi0), 0]) # phi unit vector, points along y in local coordinates
        R = np.empty([3,3],dtype=float)
        ix = np.array([1,0,0])
        iy = np.array([0,1,0])
        iz = np.array([0,0,1])
        R[0] = np.array([np.sum(ix*n0), np.sum(ix*n1), np.sum(ix*n2)])
        R[1] = np.array([np.sum(iy*n0), np.sum(iy*n1), np.sum(iy*n2)])
        R[2] = np.array([np.sum(iz*n0), np.sum(iz*n1), np.sum(iz*n2)])
        angle = 2.1*np.arctan(self.rclip/r0)
        for i in range(npix):
            thl = 0.5*angle - (i+0.5)*angle/npix
            #sthl = np.sin(thl)
            #cthl = np.cos(thl)
            for j in range(npix):
                phl = 0.5*angle - (j+0.5)*angle/npix
                #sphl = np.sin(phl)
                #cphl = np.cos(phl)
                #nl = np.array([cphl*sthl, sphl*sthl, cthl])
                nl = np.array([thl, phl, 0.])
                nl[2] = np.sqrt(1 - thl*thl - phl*phl)
                n = np.dot(R,nl)

                ic, jc, kc, xc = self.build_ray(origin, n)
                #print ic
                #print jc
                #print i, j, xc, n, np.sum(n*n)

                if self.grid.NDIM==2:
                    if len(ic) != 0:
                        vals = self.data[[ic,jc]]
                        gvals =    np.array([])    #self.grad[[ic,jc,kc]]
                    else:
                        vals = np.array([])
                        gvals = np.array([])
                elif self.grid.NDIM==3:
                    if len(ic) != 0:
                        vals = self.data[[ic,jc,kc]]
                        gvals = np.array([])    #self.data[[ic,jc,kc]]
                    else:
                        vals = np.array([])
                        gvals = np.array([])

                image[i,j] = self.integrate_pixel(vals[::-1], -gvals[::-1], xc[::-1])

                del ic, jc, kc, xc
                #print ret.shape, image[i,j,:].shape

        return image




#################################
#    Some useful utilities        #
#################################

def mspline(xin,yin,n,xlog=False,ylog=False):

    if xlog:
        x = np.log(xin)
    else:
        x = xin
    if ylog:
        y = np.log(yin+1.e-20)
    else:
        y = yin

    d = np.diff(y)/np.diff(x)
    m = np.empty(x.shape, dtype=np.double)
    for i in range(1,m.shape[0]-1):
        m[i] = 0.5*(d[i]+d[i-1])
    m[0] = d[0]
    m[x.shape[0]-1]

    a = m[:-1]/d
    b = m[1:]/d

    p = a*a + b*b
    for i in range(a.shape[0]):
        if p[i] > 9:
            t = 3./np.sqrt(p[i])
            m[i] = t*a[i]*d[i]
            m[i+1] = t*b[i]*d[i]

    xs = np.linspace(x[0],x[x.shape[0]-1],num=n,endpoint=True)
    ys = np.empty(xs.shape[0], dtype=np.double)
    xl = x[0]
    dx = x[1]-x[0]
    ix = 0
    for si,xt in enumerate(xs):
        if xt > xl+dx:
            ix += 1
            xl = x[ix]
            dx = x[ix+1]-x[ix]
        t = (xt-xl)/dx
        h00 = (1+2*t)*np.power(1-t,2)
        h10 = t*np.power(1-t,2)
        h01 = t*t*(3-2*t)
        h11 = t*t*(t-1)
        ys[si] = y[ix]*h00 + dx*m[ix]*h10 + y[ix+1]*h01 + dx*m[ix+1]*h11

    if xlog:
        xs = np.exp(xs)
    if ylog:
        ys = np.exp(ys)

    return xs,ys

def integrate(f,x0,x1):
    return romberg(f,x0,x1)

def findgain(r, rshock, q):
    """
    given 1-d radius, shock radius, and net heating, this routine returns the gain radius
    """
    return np.where(np.logical_and(np.logical_and(r<rshock,r>50.0),q>0.0))[0]
    
def findshock(r, S):
    """
    given 1-d radius and entropy, this routine returns the shock radius
    """
    dSdr = np.diff(S)/np.diff(r)
    # If it looks like rshock is jumping around a lot, try increasing the
    # the following threshhold slightly...
    dSdrth = max(-0.1,dSdr.min())
    ishock = np.where(dSdr<=dSdrth)[0][-1]
    if (ishock == 0 or ishock == len(dSdr)-1):
        return 0.0

    # return r at the min cell-centered radius
    rshock = r[ishock]

#    # return r at the min of the quadratic interpolant
#    rdp = 0.5*(r[imin+2]+r[imin+1])
#    rdm = 0.5*(r[imin+1]+r[imin  ])
#    drs = rdp-rdm
#    dp = (dS[imin+1] - dS[imin  ])/(r[imin+2]-r[imin+1])
#    dm = (dS[imin  ] - dS[imin-1])/(r[imin+1]-r[imin  ])
#    A = 0.5*(dp - dm)/drs
#    B = dp - 2*A*rdp
#    rshock = -B/(2*A)

    return rshock

def interp(x,y,x0):
    i = np.where(x < x0)[0][-1]
    dx = (x0-x[i])/(x[i+1]-x[i])
    return (1-dx)*y[i] + dx*y[i+1]

def a0(X,geom):
    if   geom.NDIM==1:
        return None
    elif geom.NDIM==2:
        return 2.0*np.pi*np.sum(X*geom.y00*geom.dmu)
    elif geom.NDIM==3:
        return np.sum(np.sum(X*geom.dphi,axis=2)*geom.y00*geom.dmu)

def a1(X,geom):
    if   geom.NDIM==1:
        return None
    elif geom.NDIM==2:
        return 2.0*np.pi*np.sum(X*geom.y01*geom.dmu)
    elif geom.NDIM==3:
        return np.sum(np.sum(X*geom.dphi,axis=2)*geom.y01*geom.dmu)

