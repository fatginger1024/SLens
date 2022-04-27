from SLens import analyser,load_MICE, load_COSMOS, ReConc_loader

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import interp2d, griddata
from scipy.stats import poisson



class truth_properties():
    
    def __init__(self,gamma_truth:float=1.3,alpha_truth:float=1.2):
        self.gamma_truth = gamma_truth
        self.alpha_truth = alpha_truth
        
class pairplots(truth_properties,analyser,load_MICE,load_COSMOS):
    
    def __init__(self,gamma_truth:float=1.3,alpha_truth:float=1.2,sample_num:int=10000):
        truth_properties.__init__(self,gamma_truth=gamma_truth,alpha_truth=alpha_truth)
        analyser.__init__(self,)
        load_MICE.__init__(self,)
        load_COSMOS.__init__(self,)
        ReConc_loader.__init__(self,)
        self.sample_num = sample_num
        
    def get_MICE(self,h=.7):
        sample_num = self.sample_num
        ind = np.arange(len(self.Mstar_arr))
        np.random.shuffle(ind)
        Mh = self.Mh_arr-np.log10(h) 
        Mstar = self.Mstar_arr-np.log10(h)
        zlens = self.zlens_arr
        Re = self.Re_arr
        Conc = self.Conc_arr
        Mh = Mh[ind][:sample_num]
        Mstar = Mstar[ind][:sample_num]
        zlens = zlens[ind][:sample_num]
        Re = Re[ind][:sample_num]
        c = Conc[ind][:sample_num]

        
        Re = np.log10(Re)
        dat = np.fromfile("/Users/tardis/Downloads/GLens1/GLens_para/lenses_all/r{}_alpha{}.bin".format(self.alpha_truth,self.gamma_truth),dtype='f8').reshape(-1,8)
        ind = np.arange(len(dat))
        np.random.shuffle(ind)
        dat = dat[ind][:sample_num]
        
        
        
        return dat, Mh, Mstar, zlens, c, Re
    
    def get_COSMOS(self,):
        sample_num = self.sample_num
        r_mag = self.COSMOS_rmag
        redshift = self.COSMOS_redshift
        mask = (redshift>=0)  *  (r_mag>0)
        ind = np.random.choice(np.arange(len(r_mag[mask])),sample_num)
        r_mag = np.array(r_mag[mask])[ind]
        redshift = np.array(redshift[mask])[ind]
        
        dat = np.fromfile("/Users/tardis/Downloads/GLens1/GLens_para/lenses_all/r{}_alpha{}_Msource.bin".format(self.alpha_truth,self.gamma_truth),dtype='f8').reshape(-1,3)
        ind = np.arange(len(dat))
        np.random.shuffle(ind)
        dat = dat[ind][:sample_num]
        
        return dat, redshift, r_mag

      
        
    def get_dataframe1(self,h=.7):
        
        dat, Mh, Mstar, zlens, c, Re = self.get_MICE()
        ind = np.arange(len(dat))
        kind = np.zeros(len(dat)+len(Mh))
        kind[-len(c):] = np.ones(len(c))
        d = {r'$\log\,M_{200}$':np.concatenate((dat[:,0]-np.log10(h),Mh),axis=None),r'$\log\,M_{*}$':np.concatenate((dat[:,1]-np.log10(h),Mstar),axis=None),r'$z_{L}$':np.concatenate((dat[:,2],zlens),axis=None),r'$c$':np.concatenate((dat[:,3],c),axis=None),r'$\log\,R_e$':np.concatenate((dat[:,4],Re),axis=None),'kind':kind}

        df = pd.DataFrame(data=d)
        df = df.replace(1,'MICE')
        df = df.replace(0,'Lenses')
        
        return df
    
    def get_dataframe2(self):
        
        dat, redshift,r_mag = self.get_COSMOS()
        kind = np.zeros(len(dat)+len(r_mag))
        kind[-len(r_mag):] = np.ones(len(r_mag))
        d = {r'$z_{S}$':np.concatenate((dat[:,2],redshift),axis=None),r'$m_0$':np.concatenate((dat[:,1],r_mag),axis=None),'kind':kind}
        df = pd.DataFrame(data=d)
        df = df.replace(1,'COSMOS')
        df = df.replace(0,'Sources')
        
        return df
        
        
    def get_pairplot1(self,):
        palette = {"Lenses":'lightsteelblue','MICE':'bisque'}    
        df = self.get_dataframe1()
        g = sns.PairGrid(df, hue="kind",corner=True)
        g.map_lower(sns.kdeplot,levels=[.003,.05,.34])
        g.map_diag(sns.histplot,kde=True,bins=20,stat='count')
        g.map_diag(sns.histplot,palette=palette,bins=20,stat='count')
        g.add_legend()
        
        g.savefig("./plots/diag1.eps",format='eps',transparent=True)
        
    def get_pairplot2(self,):
        
        palette = {"Sources":'lightsteelblue','COSMOS':'bisque'}    
        df = self.get_dataframe2()
        g = sns.PairGrid(df, hue="kind",corner=True)
        g.map_lower(sns.kdeplot,levels=[.003,.05,.34])
        g.map_diag(sns.histplot,kde=True,bins=20)
        g.map_diag(sns.histplot,palette=palette,bins=20)
        g.add_legend()
       
        g.savefig("./plots/diag2.eps",format='eps',transparent=True)
        
        
class interpolations(truth_properties):
    
    def __init__(self,gamma_truth:float=1.3,alpha_truth:float=1.2):
        truth_properties.__init__(self,gamma_truth=gamma_truth,alpha_truth=alpha_truth)
        self.num = 5
        self.alpha = np.linspace(1.0,1.8,self.num)
        self.gamma = np.linspace(.8,1.8,self.num)
        self.Alpha = np.round(np.repeat(self.alpha,self.num),2)
        self.Gamma = np.round(np.tile(self.gamma,self.num),2)
        self.totnum = np.zeros(len(self.Alpha)) 
        self.Args1 = np.zeros(len(self.Alpha))       
        self.Args2 = np.zeros(len(self.Alpha)) 
        self.Args3 = np.zeros(len(self.Alpha)) 
        self.Args4 = np.zeros(len(self.Alpha)) 
        self.load_data()
        self.interp_func()
        
        
    
    def load_data(self,):
        
        
        for i in range(len(self.Alpha)):  
            dat = np.fromfile("../Downloads/GLens1/GLens_para/lenses_all/r{}_alpha{}.bin".format(self.Alpha[i],self.Gamma[i]),dtype='f8').reshape(-1,8)    
            zlens = dat[:,2]
            einsrad = dat[:,5]
            ind = einsrad > 0
            exp_e = np.log10(einsrad[ind])
            mo1 = np.mean(exp_e)
            mo2 = np.var(exp_e)
            N = len(einsrad[ind])
            self.totnum[i] = len(zlens[zlens>0])
            self.Args1[i] = mo1
            self.Args2[i] = mo2
            self.Args3[i] = mo2/N
            self.Args4[i] = np.mean((exp_e-mo1)**4)/N- mo2**2 * (N-3)/(N*(N-1))
        
        self.totnum = np.asarray(self.totnum).reshape(len(self.alpha),len(self.gamma))
        self.Args1 = np.asarray(self.Args1).reshape(len(self.alpha),len(self.gamma))
        self.Args2 = np.asarray(self.Args2).reshape(len(self.alpha),len(self.gamma))
        self.Args3 = np.asarray(np.sqrt(self.Args3)).reshape(len(self.alpha),len(self.gamma))
        self.Args4 = np.asarray(np.sqrt(self.Args4)).reshape(len(self.alpha),len(self.gamma))

    def interp_func(self,): 
        # create a meshgrid
        self.gg,self.aa = np.meshgrid(self.gamma,self.alpha)
        self.Num_interp = 1000
        self.new_Gamma = np.linspace(.8,1.8,self.Num_interp)
        self.new_Alpha = np.linspace(1.,1.8,self.Num_interp)
        
        self.func2d = interp2d(self.gg.flatten(), self.aa.flatten(), self.totnum.flatten(),'cubic')
        # create function of mu
        self.func2d_mu = interp2d(self.gg.flatten(), self.aa.flatten(), self.Args1.flatten(),'cubic')
        # create function of sigma^2
        self.func2d_sig = interp2d(self.gg.flatten(), self.aa.flatten(), self.Args2.flatten(),'cubic')
        self.func2d_meanvar = interp2d(self.gg.flatten(), self.aa.flatten(), self.Args3.flatten(),'cubic')
        self.func2d_varvar = interp2d(self.gg.flatten(), self.aa.flatten(), self.Args4.flatten(),'cubic')
        self.new_totnum = self.func2d(self.new_Gamma,self.new_Alpha)
        self.new_mu = self.func2d_mu(self.new_Gamma,self.new_Alpha)
        self.new_sig = self.func2d_sig(self.new_Gamma,self.new_Alpha)
        self.new_meanvar = self.func2d_meanvar(self.new_Gamma,self.new_Alpha)
        self.new_varvar = self.func2d_varvar(self.new_Gamma,self.new_Alpha)

        self.NSL_obs = int(self.func2d(self.gamma_truth,self.alpha_truth))
        self.mu_obs = self.func2d_mu(self.gamma_truth,self.alpha_truth)
        self.sig_obs = self.func2d_sig(self.gamma_truth,self.alpha_truth)
        
        
            
    def get_interpolations(self,):
        
        Vals = [self.totnum,self.Args1,self.Args2,self.Args3**2,self.Args4**2]
        Vals_interp = [self.new_totnum,self.new_mu,self.new_sig,self.new_meanvar**2,self.new_varvar**2]
        Vals_title = [r"$\mathcal{N}_{\rm{SL}}$",r"$\mu_{\log\,\theta_E}$",r"$\sigma^2_{\log\,\theta_E}$",r"${Var[\mu_{\log\,\theta_E}]}$",r"${Var[\sigma^2_{\log\,\theta_E}]}$"]

        Vals_interp_title =[r"$\mathcal{N}_{\rm{SL},spl}$",r"$\mu_{\log\,\theta_E,spl}$",r"$\sigma^2_{\log\,\theta_E,spl}$",r"$Var[\mu_{\log\,\theta_E}]_{spl}$",r"$Var[\sigma^2_{\log\,\theta_E}]_{spl}$"]

        fig,ax = plt.subplots(2,5,figsize=(15,6))
        for i in range(len(Vals)):
            im1 = ax[0,i].imshow(Vals[i],extent=[-1,1,-1,1],origin='lower')
            ax[0,i].set_xticks(np.linspace(-1,1,6))
            ax[0,i].set_xticklabels(np.round(np.linspace(.8,1.8,6),2))
            ax[0,i].set_yticks(np.linspace(-1,1,6))
            ax[0,i].set_yticklabels(np.round(np.linspace(1.,1.8,6),2))
            divider = make_axes_locatable(ax[0,i])
            cax1 = divider.append_axes('right', size='5%', pad=0.05)
            cbar = fig.colorbar(im1, cax=cax1, orientation='vertical')
            if i == 0:
                cbar.set_ticks([100000,200000,300000])
                cbar.set_ticklabels(['100k','200k','300k'])
            if i in [3,4]:
                cbar.formatter.set_powerlimits((0, 0))
                cbar.ax.yaxis.set_offset_position('left')
                cbar.update_ticks()

            ax[0,i].set_ylabel(r'$\alpha_{\rm{sps}}$')
            ax[0,i].set_xlabel(r'$\gamma_{DM}$')
            ax[0,i].set_title(Vals_title[i],fontsize=15)

            im2 = ax[1,i].imshow(Vals_interp[i],extent=[-1,1,-1,1],origin='lower')
            ax[1,i].set_xticks(np.linspace(-1,1,6))
            ax[1,i].set_xticklabels(np.round(np.linspace(.8,1.8,6),2))
            ax[1,i].set_yticks(np.linspace(-1,1,6))
            ax[1,i].set_yticklabels(np.round(np.linspace(1.,1.8,6),2))
            divider = make_axes_locatable(ax[1,i])
            cax2 = divider.append_axes('right', size='5%', pad=0.05)
            cbar = fig.colorbar(im2, cax=cax2, orientation='vertical')
            if i == 0:
                cbar.set_ticks([100000,200000,300000])
                cbar.set_ticklabels(['100k','200k','300k'])
            if i in [3,4]:
                cbar.formatter.set_powerlimits((0, 0))
                cbar.ax.yaxis.set_offset_position('left')
                cbar.update_ticks()

            ax[1,i].set_ylabel(r'$\alpha_{\rm{sps}}$')
            ax[1,i].set_xlabel(r'$\gamma_{DM}$')
            ax[1,i].set_title(Vals_interp_title[i],fontsize=15)


        #fig.tight_layout(pad=.5,h_pad=-1.3,w_pad=.1)
        plt.subplots_adjust(wspace=.8, hspace=-.2)
       
        plt.show()
        fig.savefig("./plots/interp_3.eps",format='eps')
        
        
    @staticmethod  
    def probN(x,mu=1):
        return poisson.pmf(mu,x)
    
    @staticmethod  
    def probmu(x,mu=1,sig=.05):
        return np.exp(-(x-mu)**2/(2*sig**2))/(sig*np.sqrt(2*np.pi))

    @staticmethod  
    def probsig(x,mu=1,sig=.05):
        return np.exp(-(x-mu)**2/(2*sig**2))/(sig*np.sqrt(2*np.pi))
    
    @staticmethod  
    def get_levels(Z):
        Z = Z.flatten()
        r1 = np.linspace(Z.min(),Z.max(),10000)
        r_arr = np.asarray([np.sum(Z[Z>r1[i]])/Z.sum() for i in range(len(r1))])
        lv1 = r1[np.argmin(np.abs(r_arr-.68))]
        lv2 = r1[np.argmin(np.abs(r_arr-.95))]
        lv3 = r1[np.argmin(np.abs(r_arr-.997))]

        return [lv3,lv2,lv1]
    
    def get_Z(self,):
        
        Z1 = interpolations.probN(x=self.new_totnum,mu=self.NSL_obs)
        Z2 = interpolations.probmu(x=self.new_mu,mu=self.mu_obs,sig=self.new_meanvar)
        Z3 = interpolations.probsig(x=self.new_sig,mu=self.sig_obs,sig=self.new_varvar)
        Z4 = Z1*Z2*Z3
        lvs1 = interpolations.get_levels(Z1)
        lvs2 = interpolations.get_levels(Z2)
        lvs3 = interpolations.get_levels(Z3)
        lvs4 = interpolations.get_levels(Z4)
        
        return lvs1, lvs2, lvs3, lvs4, Z1, Z2, Z3, Z4
    
    def get_new_Z4(self,):

        Num_Z4 = 500
        Gamma_Z4 = np.linspace(1.294,1.306,Num_Z4)
        Alpha_Z4 = np.linspace(1.182,1.218,Num_Z4)
        
        totnum_Z4 = self.func2d(Gamma_Z4,Alpha_Z4)
        mu_Z4 = self.func2d_mu(Gamma_Z4,Alpha_Z4)
        sig_Z4 = self.func2d_sig(Gamma_Z4,Alpha_Z4)
        meanvar_Z4 = self.func2d_meanvar(Gamma_Z4,Alpha_Z4)
        varvar_Z4 = self.func2d_varvar(Gamma_Z4,Alpha_Z4)

        Z1_new = interpolations.probN(x=totnum_Z4,mu=self.NSL_obs)
        Z2_new = interpolations.probmu(x=mu_Z4,mu=self.mu_obs,sig=meanvar_Z4)
        Z3_new = interpolations.probsig(x=sig_Z4,mu=self.sig_obs,sig=varvar_Z4)
        Z4_new = Z1_new*Z2_new*Z3_new
        lvs1_new = interpolations.get_levels(Z1_new)
        lvs2_new = interpolations.get_levels(Z2_new)
        lvs3_new = interpolations.get_levels(Z3_new)
        lvs4_new = interpolations.get_levels(Z4_new)

        return Gamma_Z4, Alpha_Z4, lvs4_new, Z4_new.reshape(Num_Z4,Num_Z4)

    
    
    def get_contour(self,):
        
        lvs1, lvs2, lvs3, lvs4, Z1, Z2, Z3, Z4 = self.get_Z()
        Gamma_Z4, Alpha_Z4, lvs4, Z4 = self.get_new_Z4()
        x = self.new_Gamma
        y = self.new_Alpha
        
        fig,ax = plt.subplots(1,1,figsize=(4,4))
        ax.contourf(x,y,Z1,levels=np.concatenate((lvs1,Z1.max()),axis=None),cmap='twilight',zorder=1)
        ax.contourf(x,y,Z2,levels=np.concatenate((lvs2,Z2.max()),axis=None),cmap='Blues',zorder=2)
        ax.contourf(x,y,Z3,levels=np.concatenate((lvs3,Z3.max()),axis=None),cmap='Greens',zorder=2)
        #ax.contourf(x,y,Z4,levels=np.concatenate((lvs4,Z4.max()),axis=None),cmap='Reds',zorder=4)
        ax.contourf(Gamma_Z4,Alpha_Z4,Z4,levels=np.concatenate((lvs4,Z4.max()),axis=None),cmap='Reds',zorder=4)
        ax.text(1.09,1.76,r'$\mathcal{N}_{\rm{SL}}$',fontsize=14, style='oblique', ha='center',va='top', wrap=True)
        ax.text(1.3,1.75,r'$\mu_{\log\,\theta_E}$',fontsize=14, style='oblique', ha='center',va='top', wrap=True)
        ax.text(1.49,1.77,r'$\sigma^2_{\log\,\theta_E}$',fontsize=14, style='oblique', ha='center',va='top', wrap=True)
        ax.set_ylabel(r'$\alpha_{\rm{sps}}$')
        ax.set_xlabel(r'$\gamma_{\rm{DM}}$')

        axins = ax.inset_axes([0.58, 0.3, 0.35, 0.35])
        #axins.contourf(x,y,Z4,levels=np.concatenate((lvs4,Z4.max()),axis=None),cmap='Reds')
        axins.contourf(Gamma_Z4,Alpha_Z4,Z4,levels=np.concatenate((lvs4,Z4.max()),axis=None),cmap='Reds')
        axins.text(1.303,1.215,'joint',fontsize=10, family='serif',style='italic', ha='center',va='top', wrap=True)
        axins.plot(self.gamma_truth,self.alpha_truth,lw=0,marker='h',color="khaki",mec="k",markersize=5,zorder=2,label="truth")
        axins.tick_params(axis='x', labelsize=8) 
        axins.tick_params(axis='y', labelsize=8) 
        axins.set_xlim([1.294,1.306])
        axins.set_ylim([1.182,1.218])
        axins.legend(loc=3,fontsize="x-small")

        ax.indicate_inset_zoom(axins, edgecolor="black")
        fig.tight_layout()
        plt.show()
        fig.savefig("./plots/constrain_by_3.eps",format='eps')
        
        
    def get_limits(self,):
        
        Gamma_Z4, Alpha_Z4, lvs4, Z4 = self.get_new_Z4()
        
        cs=plt.contour(Gamma_Z4,Alpha_Z4,Z4,levels=np.concatenate((lvs4,Z4.max()),axis=None),cmap='Reds')
        path1 = cs.collections[0].get_paths()[0]._vertices
        path2 = cs.collections[1].get_paths()[0]._vertices
        path3 = cs.collections[2].get_paths()[0]._vertices

        xmin1 = path3[:,0].min();xmax1 = path3[:,0].max()
        ymin1 = path3[:,1].min();ymax1 = path3[:,1].max()
        xmin3 = path1[:,0].min();xmax3 = path1[:,0].max()
        ymin3 = path1[:,1].min();ymax3 = path1[:,1].max()
        xmean1 = .5*(xmin1+xmax1)
        ymean1 = .5*(ymin1+ymax1)
        xmean3 = .5*(xmin3+xmax3)
        ymean3 = .5*(ymin3+ymax3)
        
        print("1 sigma: ",xmean1,xmean1-xmin1,xmax1-xmean1,ymean1,ymean1-ymin1,ymax1-ymean1)
        print("3 sigma: ",xmean3,xmean3-xmin3,xmax3-xmean3,ymean3,ymean3-ymin3,ymax3-ymean3)




        
        
        
if __name__=="__main__":
    
    #plots = pairplots()
    #plots.get_pairplot1()
    #plots.get_pairplot2()
    plots = interpolations()
    plots.get_interpolations()
    plots.get_contour()
    plots.get_limits()



    
    