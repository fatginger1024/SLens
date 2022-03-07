from .analyse import analyser
from .load_data import load_MICE,load_COSMOS,ReConc_loader

import os
import sys
import numpy as np
from tqdm import tqdm

class SimRun(analyser,load_MICE,load_COSMOS,ReConc_loader):
    
    def __init__(self,z1=.3,z2=1.5,M200=1e13,Mstar=10**11.5,c=5,Re=3,alpha=1,gamma=1,source_mag=25.,dist=1):
        
        
        load_MICE.__init__(self,)
        load_COSMOS.__init__(self,)
        ReConc_loader.__init__(self,)
        analyser.__init__(self,z1=z1,z2=z2,M200=M200,Mstar=Mstar,c=c,Re=Re,
                          alpha=alpha,gamma=gamma,source_mag=source_mag,dist=dist)
        
        
    def run_one(self,i,h=.7,buffer=.1,alpha=1,gamma=1):
        
        zlens = self.zlens_arr[i]
        Mh_lens = self.Mh_arr[i]
        Mstar_lens = self.Mstar_arr[i]
        ra_lens = self.ra_arr[i]
        dec_lens = self.dec_arr[i]
        scale_rad =  self.Re_arr[i]
        conc = self.Conc_arr[i]
        z_bins = np.linspace(0,6,21)
        z_digit = np.digitize(self.source_generator[:,0],z_bins)
        zlens_digit = np.digitize(zlens,z_bins)
        rmag_max = np.min(self.source_generator[:,1][z_digit==zlens_digit])
        analyse = analyser(z1=zlens,z2=5.9,M200=10**Mh_lens*h,Mstar=10**Mstar_lens*h,c=conc,Re=scale_rad,
                          alpha=alpha,gamma=gamma,source_mag=rmag_max)
       
        search_lim = analyse.get_search_range()[1]
        lambda_rate = self.TOTnum * (np.pi*search_lim**2) / (5000*3600**2)
        INT_num = np.random.poisson(lambda_rate,size=1)
        Lens_arr = [0.]
        Source_arr = [0.]

        if INT_num !=0:
            ind_rand = np.random.choice(np.arange(len(self.source_generator)),size=INT_num)
            source_choice = self.source_generator[ind_rand]
            source_mask = source_choice[:,0] > (zlens+buffer)
            ind_source = ind_rand[source_mask]
            if np.any(source_mask) == 1:

                for j in range(len(ind_source)):

                    z2_source,rmag_source = self.source_generator[ind_source][j]
                    dist = search_lim*np.sqrt(np.random.uniform(size=1))

                    try:
                        analyse_einsrad = analyser(z1=zlens,z2=z2_source,M200=10**Mh_lens*h,Mstar=10**Mstar_lens*h,
                                                        c=conc,Re=scale_rad,alpha=alpha,gamma=gamma,
                                                   source_mag=rmag_source,dist=dist)
                        stat = analyse_einsrad.get_cross_section()
                        dist_real = stat[1]


                    except RuntimeError or KeyError:
                        print("Runtime Error.")
                        
                    if dist < dist_real:

                        eins_rad = stat[0]
                        pos_img1 = stat[2]
                        pos_img2 = stat[3]
                        mu_img1 = stat[4]
                        mu_img2 = stat[5]
                        Lens_arr = np.concatenate(np.array([Mh_lens,Mstar_lens,zlens,conc,np.log10(scale_rad),
                                             eins_rad,pos_img1,pos_img2,mu_img1,mu_img2]),axis=None)
                        Source_arr = np.array([eins_rad,rmag_source,z2_source])
                        print("found one!")
                        print(Lens_arr)
                        print(Source_arr)

                        break


        return Lens_arr,Source_arr
    
    
def run_mocks(sim_num):
    gammas = np.linspace(.8,1.8,5)
    alphas = np.linspace(1.,1.8,5)
    Gamma = np.repeat(gammas,5)
    Alpha = np.tile(alphas,5)
    cat = load_MICE()
    num_data = len(cat.Mstar_arr)
    Lens_arr = []
    Source_arr = []
    """
    num_bus = 64
    bussize = int(num_data/num_bus)+1
    count = count#int(sys.argv[1])

    if count == num_bus-1:
        index = np.arange(count*bussize,num_data)

    else:
        index = np.arange(count*bussize,(count+1)*bussize,1)
    ind_cut = ind_cut#int(sys.argv[2])
    Lens_arr = []
    Source_arr = []
    
    for item in np.arange(num_data)[index]:
        
        gamma = Gamma[ind_cut]
        alpha = Alpha[ind_cut]
        sim = SimRun(gamma=gamma,alpha=alpha)
        P,M = sim.run_one(i=item,alpha=alpha,gamma=gamma)    
        if P[0] != 0.:
            Lens_arr.append(P)
            Source_arr.append(M)
    if len(Lens_arr)>0:
        dat1 = np.concatenate(Lens_arr,axis=None).reshape(-1,10)
        dat2 = np.concatenate(Source_arr,axis=None).reshape(-1,3)
        dirbase = "SLens/mocks/"+str(count)
        if not os.path.exists(dirname):
            os.mkdir(dirname)


        #dat1.tofile(dirbase+'/gamma{}_alpha{}.bin'.format(gamma,alpha),format='f8')
        #dat2.tofile(dirbase+'/gamma{}_alpha{}_Msource.bin'.format(gamma,alpha),format='f8')
    """
    sim = SimRun()
    for i in tqdm(range(num_data)):
       
        
        P,M = sim.run_one(i,alpha=Alpha[sim_num],gamma=Gamma[sim_num])
        if P[0] != 0.:
            Lens_arr.append(P)
            Source_arr.append(M)
    if len(Lens_arr)>0:
        dat1 = np.concatenate(Lens_arr,axis=None).reshape(-1,10)
        dat2 = np.concatenate(Source_arr,axis=None).reshape(-1,3)
        dirbase = "./mocks/"
        #if not os.path.exists(dirbase):
        #    os.mkdir(dirbase)
        dat1.tofile(os.path.join(dirbase,'gamma{}_alpha{}.bin').format(Gamma[sim_num],Alpha[sim_num]),format='f8')
        dat2.tofile(os.path.join(dirbase,'gamma{}_alpha{}_Msource.bin').format(Gamma[sim_num],Alpha[sim_num]),format='f8')

if __name__ == "__main__":
    start = int(sys.argv[1])
    end = int(sys.argv[2])
    for i in range(start,end+1):
        run_mocks(i)




