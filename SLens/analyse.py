from SLens.cp_model import gnfwSersic

from math import *
import numpy as np
from scipy.optimize import fmin,bisect



class analyser(gnfwSersic):
    
    def __init__(self,z1=.3,z2=1.5,M200=1e13,Mstar=10**11.5,c=5,Re=3,alpha=1,gamma=1,source_mag=25.,dist=1):
        
        gnfwSersic.__init__(self,z1=z1,z2=z2,M200=M200,Mstar=Mstar,c=c,Re=Re,alpha=alpha,gamma=gamma)
        self.dist = dist
        self.mag_unlensed = source_mag
        
    def get_search_range(self,tol=26.3,buffer=.1):
        xmin,xmax = -6,0
        xval = np.logspace(xmin,xmax,10000)
        thetas = xval*self.thetas*206265
        betas = self.lens_beta(xval)*self.thetas*206265
        ind = np.argmin(betas)
        xval_min = xval[ind]
        caus1 = betas[ind]
        rt_eins_ind = np.argmin(np.abs(self.lens_alpha(xval)-xval))
        eins = xval[rt_eins_ind]*self.thetas*206265
        msk = (xval>xval_min*(1+buffer))*(betas<0)
        mag = np.abs(1/self.lens_detA(xval[msk]))
        app_m = self.get_lensed_mag(mag)
        mg_msk = app_m<tol
        
        try:
            beta_mag = betas[msk][mg_msk][0]
            
        except ValueError or IndexError:
            
            beta_mag = 0
        
        return np.abs(caus1),np.abs(beta_mag)
        
    
    def get_cross_section(self,tol=26.3,buffer=.1):
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
        msk = (xval>xval_min*(1+buffer))*(betas<0)
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
            
        return np.abs(beta_mag),eins,pos_image1,pos_image2,mu_image1,mu_image2

    def get_lensed_mag(self,m):
        return self.mag_unlensed - 2.5*np.log10(m)
    
        
        
        

        

        