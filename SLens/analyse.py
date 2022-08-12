from .cp_model import gnfwSersic

from math import *
import numpy as np
from scipy.optimize import minimize,fmin,fmin_slsqp,bisect,brentq




class analyser(gnfwSersic):
    
    def __init__(self,z1=.3,z2=1.5,M200=1e13,Mstar=10**11.5,c=5,Re=3,alpha=1,gamma=1,source_mag=25.,dist=1,m_tol=26.3):
        
        gnfwSersic.__init__(self,z1=z1,z2=z2,M200=M200,Mstar=Mstar,c=c,Re=Re,alpha=alpha,gamma=gamma)
        self.dist = dist
        self.mag_unlensed = source_mag
        self.m_tol = m_tol
        self.attr = "attr"
    
 
    
    @property    
    def attr(self):
        return self._attr

    @attr.setter
    def attr(self,val):
        """
        Function that sets the lens statistics.
        Returns: 
        - Einstein radius (tangential critical curve), 
        - radial critical curve,
        - tangential caustic,
        - radial caustic,
        - strong lensing limit (in the source plane) 
          with applied apparent magnitude threshold.
        - strong lensing limit (in the lens plane) 
          with applied apparent magnitude threshold.  
        """
        minsearch = fmin(self.lens_beta,1e-2,disp=0)
        xval_min = minsearch[0]
        radial_critical = xval_min*self.thetas*206265
        sol = brentq(self.lens_beta,self._x[0],self._x[-1])
        eins = sol*self.thetas*206265
        caus_r = self.lens_beta(xval_min)[0]*self.thetas*206265
        caus_t = self.lens_beta(sol)[0]*self.thetas*206265
        x_try = np.linspace(1e-4,2e-2,100)
        ind_try = np.argmin(self.magnitude_difference(x_try))
        x0 = .5*sol
        sol = fmin_slsqp(self.magnitude_difference,x0,bounds=((xval_min,sol),),disp=0)[0]
        beta_mag = self.lens_beta(sol)[0]*self.thetas*206265
        self._attr = (
                eins,
                radial_critical,
                np.abs(caus_t),
                np.abs(caus_r),
                np.abs(beta_mag),
                sol
                )
    
    def get_search_range(self):
        
        return self.attr[3],self.attr[4]
        
      
    def get_cross_section(self,buffer=.1):
        
        dist_beta = (self.dist/(self.thetas*206265))
        func1 = lambda x: np.abs(self.lens_beta(x)-dist_beta)
        func2 = lambda x: np.abs(self.lens_beta(x)+dist_beta)
        dist_image1 = fmin(func1,.1,xtol=1e-6,disp=0)[0]
        pos_image1 = dist_image1 * (self.thetas*206265)
        mu_image1 = np.abs(1/self.lens_detA(dist_image1))[0]
        dist_image2 = fmin(func2,.5*dist_image1,xtol=1e-6,disp=0)[0]
        pos_image2 = -dist_image2 * (self.thetas*206265)
        mu_image2 = np.abs(1/self.lens_detA(dist_image2))[0]
        
        if self.dist <= self.attr[4]:
            pass
             
        else:
            pos_image1 = 0
            mu_image1 = 0
            pos_image2 = 0
            mu_image2 = 0
        
        
        return self.attr[4],self.attr[0],pos_image1,pos_image2,mu_image1,mu_image2
    
    def plotter(self,):
        const = self.thetas*206265 
        x = np.linspace(self.attr[1]/const,self.attr[0]/const,100)
        beta = self.lens_beta(x)
        y = self.magnitude_difference(x)
        
        import matplotlib.pyplot as plt
        
        plt.figure(figsize=(4,4))
        plt.plot(x,y)
        plt.xlabel("x")
        plt.ylabel(r"$m_{lensed}-m_{tol}$")
        plt.show()
        
          

    def get_lensed_mag(self,m):
        return self.mag_unlensed - 2.5*np.log10(m)
    
    def magnitude_difference(self,x):
        m_tol = self.m_tol
        return np.abs(self.get_lensed_mag(1./np.abs(self.lens_detA(x)))-m_tol)
    
    
    
    
    
    
if __name__ == "__main__":
    stat = analyser(M200=10**13.5,dist=.5)
    stat.plotter()
    stat.get_search_range()
    
    
        

        

        
