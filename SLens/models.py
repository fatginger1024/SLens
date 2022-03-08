from math import *
import numpy as np
from scipy.integrate import quad,quad_vec
from scipy.special import hyp2f1,gamma,gammainc
from scipy.optimize import bisect

from .load_data import tables


class Distances():
    
    def __init__(self,z0=0,z1=.3,z2=1.5,Om=.25,Ol=.75,Or=8.4e-5,H0=70):
        
        self.Om = Om
        self.Or = Or
        self.Ol = Ol
        self._z0 = z0
        self._z1 = z1
        self._z2 = z2
        self.H0 = H0
        h = (self.H0/3.086e19)/100
        G  = 6.67e-11 # m^3  kg^-1 s^-2
        Ez = np.sqrt(self.Or*(1+self._z1)**4+self.Om*(1+self._z1)**3+self.Ol)
        self.Hz = (self.H0/3.086e19)*Ez
        self._Da1 = self.angular_diameter_distance(0,self._z1)
        self._rhoz = 3*self.Hz**2/(8*pi*G)/1.989e30/(3.24e-23)**3  # M_sun * Mpc^(-3)
        self._Sigmacr = self.get_Sigmacr()  #Msun Mpc^-2
        
        
    def comoving_distance(self,z1,z2):
        invE = lambda z: 1./np.sqrt(self.Or*(1+z)**4+self.Om*(1+z)**3+self.Ol)
        res = quad(invE,z1,z2)[0]*3e5/self.H0
        return res

    def angular_diameter_distance(self,z1,z2):
        invE = lambda z: 1/np.sqrt(self.Or*(1+z)**4+self.Om*(1+z)**3+self.Ol)
        res = quad(invE,z1,z2)[0]*3e5/self.H0
        return res/(1+z2)
    
    def DA(self,OR=0.):
        DH = 3e5/self.H0
        Dm1 = self.comoving_distance(self._z0,self._z1)
        Dm2 = self.comoving_distance(self._z0,self._z2)
        return (Dm2*np.sqrt(1+OR*Dm1**2/DH**2)-Dm1*np.sqrt(1+OR*Dm2**2/DH**2))/(1+self._z2)
    
    def get_Sigmacr(self):
        G = 6.67e-11*3.24e-23*1.989e30*1e-6  #Mpc Msun^-1 (km/s)^2
        c = 3e5 #km/s
        Da1 = self.angular_diameter_distance(self._z0,self._z1)
        Da2 = self.angular_diameter_distance(self._z0,self._z2)
        Da12 = self.angular_diameter_distance(self._z1,self._z2)
        return c**2/(4*pi*G) * Da2/(Da1*Da12)
        
class concentration():
    
    def __init__(self,):
        pass
    
    def OL(self,z,OL0=.75,OM0=.25):
        return OL0/(OL0+OM0*(1+z)**3)

    def OM(self,z):
        return 1-self.OL(z)

    def Phi(self,z):
        return self.OM(z)**(4/7.)-self.OL(z)+(1+self.OM(z)/2)*(1+self.OL(z)/70)

    def D(self,z,OM0=.25):
        return self.OM(z)/OM0 * self.Phi(0)/self.Phi(z) /(1+z)
    
    def sigma(self,M,z):
            xi = 1/(M/1e10)
            return self.D(z)*22.26*xi**.292/(1+1.53*xi**.275+3.36*xi**.198)
    
    def nu(self,M,z,delta_sc=1.686):
        return delta_sc/self.sigma(M,z)

    def concentration(self,M,z):
        def c0(z):
            return 3.395*(1+z)**(-.215)
        def beta(z):
            return .307*(1+z)**(.540)
        def gamma1(z):
            return .628*(1+z)**(-.047)
        def gamma2(z):
            return .317*(1+z)**(-.893)
        def nu0(z):
            a = 1/(1+z)
            return (4.135-.564/a-.210/a**2+.0557/a**3-.00348/a**4)/self.D(z)
        return c0(z)*(self.nu(M,z)/nu0(z))**(-gamma1(z))*(1+(self.nu(M,z)/nu0(z))**(1/beta(z)))**(-beta(z)*(gamma1(z)+gamma2(z)))
    
        
    
class gnfw(tables):
    
    def __init__(self,z1=.3,z2=1.5,M200=1e13,c=5,gamma=1):
        
        tables.__init__(self,)
        self.gamma = gamma # power law index of the gnfw profile
        self.z1 = z1
        self.z2 = z2
        self.c = c #concentration parameter
        self.rhoz = self.rhoz_dict[str(np.round(self.z1,2))+'+'+str(np.round(self.z2,2))]
        self.rhos = 200/3*self.rhoz*(3-gamma)*c**(gamma)/hyp2f1(3-gamma,3-gamma,4-gamma,-c)
        self.M200 = M200
        self.r200 = (self.M200/(4/3*pi*200*self.rhoz))**(1/3) #Mpc
        self.rs = self.r200/self.c #Mpc
        self.Da1 = self.Da1_dict[str(np.round(self.z1,2))+'+'+str(np.round(self.z2,2))]
        self.thetas = self.rs/self.Da1 #in radians
        self.Sigmacr = self.Sigmacr_dict[str(np.round(self.z1,2))+'+'+str(np.round(self.z2,2))]
        self.kappas =  self.rhos*self.rs/self.Sigmacr
        self.b = 4*self.rhos*self.rs/self.Sigmacr
        
     
    def gnfw_Sigma(self,x,f):
        return 2*self.rhos*self.rs*f(x,self.gamma)
    
    def gnfw_alpha(self,x,g):
        
        return self.b*g(x,self.gamma)
    
    def gnfw_kappa(self,x,f):
        
        return self.b*f(x,self.gamma)/2
    
    def gnfw_gamma(self,x,f,g):
        
        return self.gnfw_alpha(x,g)/x - self.gnfw_kappa(x,f)
    
    
    def gnfw_detA(self,x,f,g):
        x = np.atleast_1d(x)
        detA = x*0.
        detA = (1-self.gnfw_kappa(x,f))**2-self.gnfw_gamma(x,f,g)**2
        
        return detA
        
    def gnfw_beta(self,x,g):
        x = np.atleast_1d(x)
        beta = x*0.
        beta = x-self.gnfw_alpha(x,g)
        
        return beta

    @staticmethod
    def f(x,gamma):
        def quad_func(y):
            return (y+x)**(gamma-4)*(1-np.sqrt(1-y**2))
            
        x = np.atleast_1d(x)
        
        return x**(1-gamma)*((1+x)**(gamma-3)+(3-gamma)*quad_vec(quad_func,0,1)[0])
    
    @staticmethod
    def g(x,gamma):
        
        def hypF(a,b,c,z):
            return hyp2f1(a,b,c,z)
        
        def quad_func(y):
            return (y+x)**(gamma-3)*(1-np.sqrt(1-y**2))/y

    
        x = np.atleast_1d(x)
        out = x**(2-gamma)*(hypF(3-gamma,3-gamma,4-gamma,-x)/(3-gamma)+quad_vec(quad_func,0,1)[0])
        
        return out
       
class mass_size():
    
    #@staticmethod
    def get_Re_from_Mstar(self,Mstar,sigmaR=.147):
        
        def mu(Mstar,mu0R=.855,betaR=1.366):
            return mu0R+betaR*(np.log10(Mstar)-11.4)
        
        loc = mu(Mstar)
        logRe =  np.random.normal(loc=loc,scale=sigmaR,size=1)
        
        return 10**logRe
    
class Sersic():
    """
    =======================================================
    Galaxy parameters:
    =======================================================
    - _Mstar: the mass of the galaxy
    - _Re: effective radius of the galaxy, i.e. the projected 
    radius encircling half of the total luminosity
    - _m: Sersic law index, default set to be 4. Normally 
    1<= m <= 4
    =======================================================
    """
    
    def __init__(self,Mstar=1e12,Re=3,m=4):
        self._Mstar =  Mstar#Msun
        self._Re = Re/1e3 #Mpc
        self._m = m
        
        
    def b_func(self):
        m = self._m
        return 2*m-1/3.+4/(405*m)+46/(25515*m**2)+131/(1148175*m**3)-2194697/(30690717750*m**4)

    def L(self):
        Re = self._Re
        m =  self._m
        return Re**2*2*pi*m/self.b_func()**(2*m)*gamma(2*m)

    def Sigma(self,R):
        m = self._m
        Re = self._Re
        return self._Mstar*np.exp(-self.b_func()*(R/Re)**(1./m))/self.L()
   
    
    def M2d(self,R):
        Re = self._Re
        m = self._m
        b = self.b_func()
        return 2*pi*self._Mstar/self.L()*m*Re**2/b**(2*m)*gammainc(2*m,b*R**(1/m)/Re**(1/m))*gamma(2*m)
    
    def galaxy_alpha(self,R,Sigmacr,x):
        
        alpha = self.M2d(R)/(pi*R**2)/Sigmacr*x
        
        return alpha
    
    def galaxy_beta(self,R,Sigmacr,x):
        
        beta = x - self.galaxy_alpha(R,Sigmacr,x)
        
        return beta
    
    def galaxy_kappa(self,R,Sigmacr):
        
        kappa = self.Sigma(R)/Sigmacr
        
        return kappa
    
    def galaxy_gamma(self,R,Sigmacr,x):
        
        gamma = self.galaxy_alpha(R,Sigmacr,x)/x - self.galaxy_kappa(R,Sigmacr)
        
        return gamma
    
    def galaxy_detA(self,R,Sigmacr,x):
        
        detA = (1-self.galaxy_kappa(R,Sigmacr))**2-self.galaxy_gamma(R,Sigmacr,x)**2
        
        return detA
        
    
