from SLens.models import Sersic, gnfw
from SLens.load_data import fxgx

import numpy as np
from math import *
from scipy.special import hyp2f1
from scipy.interpolate import interp2d


class gnfwSersic(gnfw,Sersic,fxgx):
    
    def __init__(self,z1=.3,z2=1.5,M200=1e13,Mstar=10**11.5,c=5,Re=3,alpha=1,gamma=1): 
        gnfw.__init__(self,z1=z1,z2=z2,M200=1e13,c=5,gamma=1)
        Sersic.__init__(self,Mstar=alpha*Mstar,Re=Re)
        fxgx.__init__(self,)
   
        self.fx_approx = interp2d(self._x,self._gamma,self.fvals)
        self.gx_approx = interp2d(self._x,self._gamma,self.gvals)
        
        
      
        
    def lens_alpha(self,x):
        r = x*self.rs
        gnfw = self.gnfw_alpha(x,self.gx_approx)
        galaxy = self.galaxy_alpha(r,self.Sigmacr,x)
        
        return gnfw + galaxy
    
    def lens_kappa(self,x):
        r = x*self.rs
        gnfw = self.gnfw_kappa(x,self.fx_approx)
        galaxy = self.galaxy_kappa(r,self.Sigmacr)
        
        return gnfw + galaxy
    
    def lens_gamma(self,x):
        r = x*self.rs
        gnfw = self.gnfw_gamma(x,self.fx_approx,self.gx_approx)
        galaxy = self.galaxy_gamma(r,self.Sigmacr,x)
        
        return gnfw + galaxy
    
    def lens_detA(self,x):
        x = np.atleast_1d(x)
        detA = x*0.
        detA = (1-self.lens_kappa(x))**2-self.lens_gamma(x)**2
        
        return detA
        
    def lens_beta(self,x):
        x = np.atleast_1d(x)
        beta = x*0.
        beta = x-self.lens_alpha(x)
        
        return beta
        
        
   
        
        

        
