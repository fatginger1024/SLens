from .analyse import analyser
from .load_data import load_MICE,load_COSMOS,ReConc_loader

import os
import sys
import time
import numpy as np
from tqdm import tqdm
from multiprocessing import Process, Queue, Value, Array

class SimRun(analyser,load_MICE,load_COSMOS,ReConc_loader):
    
    def __init__(self,z1=.3,z2=1.5,M200=1e13,Mstar=10**11.5,c=5,Re=3,alpha=1,gamma=1,source_mag=25.,dist=1):
        
        
        load_MICE.__init__(self,data="128")
        load_COSMOS.__init__(self,)
        ReConc_loader.__init__(self,)
        analyser.__init__(self,z1=z1,z2=z2,M200=M200,Mstar=Mstar,c=c,Re=Re,
                          alpha=alpha,gamma=gamma,source_mag=source_mag,dist=dist)
        self.z_bins = np.linspace(0,6,21)
        self.z_digit = np.digitize(self.source_generator[:,0],self.z_bins)
        
    @property
    def pops(self):
        return self._pops
        
    @pops.setter
    def pops(self,i):  
        zlens = self.zlens_arr[i]
        Mh_lens = self.Mh_arr[i]
        Mstar_lens = self.Mstar_arr[i]
        ra_lens = self.ra_arr[i]
        dec_lens = self.dec_arr[i]
        scale_rad =  self.Re_arr[i]
        conc = self.Conc_arr[i]
        self._pops = (zlens,Mh_lens,Mstar_lens,ra_lens,dec_lens,scale_rad,conc)
        
    
    def count_source(self,i,h=.7,alpha=1,gamma=1):
        
        self.pops = i
        zlens,Mh_lens,Mstar_lens,ra_lens,dec_lens,scale_rad,conc = self.pops
        zlens_digit = np.digitize(zlens,self.z_bins)
        rmag_max = np.min(self.source_generator[:,1][self.z_digit==zlens_digit])
        analyse = analyser(z1=zlens,z2=5.9,M200=10**Mh_lens*h,Mstar=10**Mstar_lens*h,c=conc,Re=scale_rad,
                          alpha=alpha,gamma=gamma,source_mag=rmag_max)
       
        search_lim = analyse.get_search_range()[1]
        lambda_rate = self.TOTnum * (np.pi*search_lim**2) / (5000*3600**2)
        INT_num = np.random.poisson(lambda_rate,size=1)
        
        return self.pops, search_lim, INT_num
        
          
        
    def run_one(self,i,h=.7,buffer=.1,alpha=1,gamma=1):
        
        pops, search_lim, INT_num = self.count_source(i,alpha=alpha,gamma=gamma)
        zlens,Mh_lens,Mstar_lens,ra_lens,dec_lens,scale_rad,conc = pops
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
                        dist_real = stat[0]


                    except RuntimeError or KeyError:
                        print("Runtime Error.")
                        
                    if dist < dist_real and stat[2] !=0:

                        eins_rad = stat[1]
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

def qinit(q,index):
    for i in index:
        q.put(i)

def map_func(q,count,sim,arr1,arr2,alpha,gamma):
    while not q.empty():
        ind = q.get()
        P, M = sim.run_one(ind,alpha=alpha,gamma=gamma)
        if P[0] != 0.:
            arr1[ind*10:(ind+1)*10] = P
            arr2[ind*3:(ind+1)*3] = M
            count.value += 1
            num = count.value

def main(sim_num,num_proc=8):
    gammas = np.linspace(.8,1.8,5)
    alphas = np.linspace(1.,1.8,5)
    Gamma = np.repeat(gammas,5)
    Alpha = np.tile(alphas,5)
    cat = load_MICE(data="128")
    num_data = len(cat.Mstar_arr)
    index = np.arange(num_data)
    Lens_arr = Array('d', np.zeros(num_data*10))
    Source_arr = Array('d', np.zeros(num_data*3))
    process_list = []
    q = Queue(len(index))
    num = Value('i', 0)
    qinit(q,index)
    sim = SimRun()
    dirbase = "./mocks/"
    if not os.path.exists(dirbase):
        os.mkdir(dirbase)
    for i in range(num_proc):
        p = Process(target=map_func, 
                    args=(q,
                          num,
                          sim,
                          Lens_arr,
                          Source_arr,
                          Alpha[sim_num],
                          Gamma[sim_num]
                         )
                   )
        p.start()
        process_list.append(p)
    for i in process_list:
        p.join() 
    fp = dirbase+'gamma{}_alpha{}.txt'.format(Gamma[sim_num],Alpha[sim_num])
    time.sleep(60)
    np.savetxt(fp,Lens_arr[:])
    np.savetxt(fp[:-4]+"_source.txt",Source_arr[:])
    
    
if __name__ == "__main__":
    sim_num = int(sys.argv[1])
    main(sim_num)
   




