from SLens.models import Sersic
from SLens.load_data import tables, fxgx
from scipy.interpolate import interp2d

import numpy as np
from math import *
from scipy.special import hyp2f1
from scipy.optimize import fmin,bisect



class gnfwSersic(Sersic,fxgx,tables):
    
    def __init__(self,z1=.3,z2=1.5,M200=1e13,Mstar=10**11.5,c=5,Re=3,alpha=1,gamma=1,source_mag=25.,dist=1): 
        Sersic.__init__(self,)
        fxgx.__init__(self,)
        tables.__init__(self,)
        self.fx_approx = interp2d(self._x,self._gamma,self.fvals)
        self.gx_approx = interp2d(self._x,self._gamma,self.gvals)
        self.dist = dist
        self.gamma = gamma # power law index of the gnfw profile
        self.z1 = z1
        self.z2 = z2
        self.mag_unlensed = source_mag
        self.c = c #concentration parameter
        self.rhoz = self.rhoz_dict[str(np.round(self.z1,2))+'+'+str(np.round(self.z2,2))]
        self.rhos = 200/3*self.rhoz*(3-gamma)*self.c**(gamma)/hyp2f1(3-gamma,3-gamma,4-gamma,-self.c)
        self.M200 = M200
        self.Mstar = Mstar
        self.r200 = (self.M200/(4/3*pi*200*self.rhoz))**(1/3) #Mpc
        self.rs = self.r200/self.c #Mpc
        
        self.Da1 = self.Da1_dict[str(np.round(self.z1,2))+'+'+str(np.round(self.z2,2))]
        self.thetas = self.rs/self.Da1 #in radians
        self.Sigmacr = self.Sigmacr_dict[str(np.round(self.z1,2))+'+'+str(np.round(self.z2,2))]
        self.kappas =  self.rhos*self.rs/self.Sigmacr
        self.b = 4*self.rhos*self.rs/self.Sigmacr
      
        self.galaxy = Sersic(Mstar=alpha*Mstar,Re=Re)
        
        
      
        
    def lens_alpha(self,x):
        r = x*self.rs
        alpha_galaxy = self.galaxy.M2d(r)/(pi*r**2)/self.Sigmacr*x
        
        return self.b*self.gx_approx(x,self.gamma)+alpha_galaxy
    
    def lens_kappa(self,x):
        r = x*self.rs
        kappa_galaxy = self.galaxy.Sigma(r)/self.Sigmacr
        
        return self.b*self.fx_approx(x,self.gamma)/2+kappa_galaxy
    
    def lens_gamma(self,x):
        r = x*self.rs
        gamma_galaxy = (self.galaxy.M2d(r)/(pi*r**2)-self.galaxy.Sigma(r))/self.Sigmacr
        
        return self.b*self.gx_approx(x,self.gamma)/x-self.b*self.fx_approx(x,self.gamma)/2+gamma_galaxy
    
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
        
        
    def get_search_range(self,tol=26.3):
        xmin,xmax = -6,0
        xval = np.logspace(xmin,xmax,10000)
        thetas = xval*self.thetas*206265
        betas = self.lens_beta(xval)*self.thetas*206265
        #check if singularity is in range(0,2*rt2)
        ind = np.argmin(betas)
        xval_min = xval[ind]
        caus1 = betas[ind]
        rt_eins_ind = np.argmin(np.abs(self.lens_alpha(xval)-xval))
        eins = xval[rt_eins_ind]*self.thetas*206265
        msk = (xval>xval_min*1.1)*(betas<0)
        mag = np.abs(1/self.lens_detA(xval[msk]))
        app_m = self.get_lensed_mag(mag)
        mg_msk = app_m<tol
        
        try:
            beta_mag = betas[msk][mg_msk][0]
            
        except ValueError or IndexError:
            print(betas[msk][mg_msk])
            beta_mag = 0
        
        return np.abs(caus1),beta_mag
    
    def get_cross_section(self,tol=26.3):
        xmin,xmax = -4,0
        xval = np.logspace(xmin,xmax,10000)
        thetas = xval*self.thetas*206265
        alphas = lambda x:-self.lens_alpha(x)
        xmax_alpha = fmin(alphas,1e-2,disp=0)
        betas = self.lens_beta(xval)*self.thetas*206265
        ind = np.argmin(betas)
        xval_min = xval[ind]
        caus1 = betas[ind]
        rt_eins_ind = np.argmin(np.abs(self.lens_alpha(xval)-xval))
        eins = xval[rt_eins_ind]*self.thetas*206265
        msk = (xval>xval_min*1.1)*(betas<0)
        mag = np.abs(1/self.lens_detA(xval[msk]))
        app_m = self.get_lensed_mag(mag)
        mg_msk = app_m<tol
        
        try:
            beta_mag = betas[msk][mg_msk][0]
        except ValueError or IndexError:
            beta_mag = 0
            
       
            
        dist_beta = (self.dist/(self.thetas*206265))
        
        if dist_beta <= np.abs(self.lens_beta(xval_min)):
            dist_image1 = xval[np.argmin(np.abs(self.lens_beta(xval) - dist_beta))]
            try:
                func_image2 = lambda x: self.lens_beta(x) + dist_beta
                dist_image2 = bisect(func_image2,xmax_alpha,xval[-1])
                
            except ValueError:
                dist_image2 = xval[np.argmin(np.abs(self.lens_beta(xval) + dist_beta))]
            
            pos_image1 = dist_image1 * (self.thetas*206265)
            pos_image2 = -dist_image2 * (self.thetas*206265)
            mu_image1 = np.abs(1/self.lens_detA(dist_image1))
            mu_image2 = np.abs(1/self.lens_detA(dist_image2))
        else:
            
            pos_image1 = 0
            pos_image2 = 0
            mu_image1 = 0
            mu_image2 = 0
            
        return np.abs(beta_mag),eins,[pos_image1,pos_image2],[mu_image1,mu_image2]

    def get_lensed_mag(self,m):
        return self.mag_unlensed - 2.5*np.log10(m)
    
        
        
        

        
