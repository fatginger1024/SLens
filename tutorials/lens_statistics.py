from SLens import gnfwSersic

import numpy as np
import matplotlib.pyplot as plt



class lensplot(gnfwSersic):
    
    def __init__(self,z1=.3,z2=1.5,M200=1e13,Mstar=10**11.5,c=5,Re=3,alpha=1,gamma=1):
        
        gnfwSersic.__init__(self,z1=z1,z2=z2,M200=M200,Mstar=Mstar,c=c,Re=Re,alpha=alpha,gamma=gamma)
        
        
    
 
    def plot_fxgx(self,write=False):
        fig,ax = plt.subplots(1,2,figsize=(10,5))
        [ax[0].loglog(self._x,self.gvals[i],
                      ls='-',label=r'$\gamma={}$'.format(round(self._gamma[i],2))) for i in range(self._gamma_num)]
        [ax[1].loglog(self._x,self.fvals[i],
                      ls='-',label=r'$\gamma={}$'.format(round(self._gamma[i],2))) for i in range(self._gamma_num)]
        ax[0].set_title("Interpolated g(x)")
        ax[0].set_xlabel('x')
        ax[0].set_ylabel('g(x)')
        ax[0].legend()
        ax[1].set_title("Interpolated f(x)")
        ax[1].set_xlabel('x')
        ax[1].set_ylabel('f(x)')
        ax[1].legend()
        plt.show()
        if write:
            fig.savefig("../plots/fxgx.eps",format="eps")

    
    def plot_Sigma(self,write=False):
        Sigmas = self.lens_Sigma(self._x)
        Sigma_gnfw = self.gnfw_Sigma(self._x,self.fx_approx)
        Sigma_galaxy = self.Sigma(self._x*self.rs)
        fig,ax = plt.subplots(1,1,figsize=(5,5))
        ax.loglog(self._x,Sigmas,ls='-',color="darkred",label="total",zorder=2) 
        ax.loglog(self._x,Sigma_gnfw,ls='-',color='rebeccapurple',label="dark matter",zorder=1) 
        ax.loglog(self._x,Sigma_galaxy,ls='-',color='orangered',label="baryons",zorder=1) 
        ax.set_title("Projected mass density")
        ax.set_xlabel('x')
        ax.set_ylabel(r"$\Sigma(x)$")
        ax.legend()

        plt.show()
        if write:
            fig.savefig("../plots/Sigma.eps",format="eps")
            
    def plot_alpha(self,write=False):
        alphas = self.lens_alpha(self._x)
        alpha_gnfw = self.gnfw_alpha(self._x,self.gx_approx)
        alpha_galaxy = self.galaxy_alpha(self._x*self.rs,self.Sigmacr,self._x)
        fig,ax = plt.subplots(1,1,figsize=(5,5))
        ax.loglog(self._x,alphas,ls='-',color="darkred",label="total",zorder=2) 
        ax.loglog(self._x,alpha_gnfw,ls='-',color='rebeccapurple',label="dark matter",zorder=1) 
        ax.loglog(self._x,alpha_galaxy,ls='-',color='orangered',label="baryons",zorder=1) 
        ax.set_title("Deflection angle")
        ax.set_xlabel('x')
        ax.set_ylabel(r"$\alpha(x)$")
        ax.legend()

        plt.show()
        if write:
            fig.savefig("../plots/alpha.eps",format="eps")
            
    def plot_beta(self,write=False):
        betas = self.lens_beta(self._x)
        beta_gnfw = self.gnfw_beta(self._x,self.gx_approx)
        beta_galaxy = self.galaxy_beta(self._x*self.rs,self.Sigmacr,self._x)
        fig,ax = plt.subplots(1,1,figsize=(5,5))
        ax.plot(np.concatenate((-self._x[::-1],0,self._x),axis=None),
                np.concatenate((-betas[::-1],0,betas),axis=None),
                ls='-',color="darkred",label="total",zorder=2) 
        ax.plot(np.concatenate((-self._x[::-1],0,self._x),axis=None),
                np.concatenate((-beta_gnfw[::-1],0,beta_gnfw),axis=None),
                ls='-',color='rebeccapurple',label="dark matter",zorder=1) 
        ax.plot(np.concatenate((-self._x[::-1],0,self._x),axis=None),
                np.concatenate((-beta_galaxy[::-1],0,beta_galaxy),axis=None),
                ls='-',color='orangered',label="baryons",zorder=1) 
        ax.set_title("Lens equation")
        ax.set_xlabel('x')
        ax.set_ylabel(r"$\beta(x)$")
        ax.set_xlim([-self._x[-1]*.1,self._x[-1]*.1])
        ax.set_ylim([-.15,.15])
        ax.legend()

        plt.show()
        if write:
            fig.savefig("../plots/beta.eps",format="eps")
            
            
    def plot_kappa(self,write=False):
        kappas = self.lens_kappa(self._x)
        kappa_gnfw = self.gnfw_kappa(self._x,self.fx_approx)
        kappa_galaxy = self.galaxy_kappa(self._x*self.rs,self.Sigmacr)
        fig,ax = plt.subplots(1,1,figsize=(5,5))
        ax.loglog(self._x,kappas,ls='-',color="darkred",label="total",zorder=2) 
        ax.loglog(self._x,kappa_gnfw,ls='-',color='rebeccapurple',label="dark matter",zorder=1) 
        ax.loglog(self._x,kappa_galaxy,ls='-',color='orangered',label="baryons",zorder=1) 
        ax.set_title("Convergence")
        ax.set_xlabel('x')
        ax.set_ylabel(r"$\kappa(x)$")
        ax.legend()

        plt.show()
        if write:
            fig.savefig("../plots/kappa.eps",format="eps")
            
            
    def plot_gamma(self,write=False):
        gammas = self.lens_gamma(self._x)
        gamma_gnfw = self.gnfw_gamma(self._x,self.fx_approx,self.gx_approx)
        gamma_galaxy = self.galaxy_gamma(self._x*self.rs,self.Sigmacr,self._x)
        fig,ax = plt.subplots(1,1,figsize=(5,5))
        ax.loglog(self._x,gammas,ls='-',color="darkred",label="total",zorder=2) 
        ax.loglog(self._x,gamma_gnfw,ls='-',color='rebeccapurple',label="dark matter",zorder=1) 
        ax.loglog(self._x,gamma_galaxy,ls='-',color='orangered',label="baryons",zorder=1) 
        ax.set_title("Shear")
        ax.set_xlabel('x')
        ax.set_ylabel(r"$\gamma(x)$")
        ax.legend()

        plt.show()
        if write:
            fig.savefig("../plots/gamma.eps",format="eps")
            
   
            




if __name__=="__main__":
    plot = lensplot()
    plot.plot_fxgx()
    plot.plot_alpha()
    plot.plot_beta()
    plot.plot_kappa()
    plot.plot_gamma()
