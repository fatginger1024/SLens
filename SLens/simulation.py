from SLens.analyse import analyser
from SLens.load_data import load_MICE,load_COSMOS, ReConc_loader

import numpy as np

class SimRun(analyser,load_MICE,load_COSMOS,ReConc_loader):
    
    def __init__(self,z1=.3,z2=1.5,M200=1e13,Mstar=10**11.5,c=5,Re=3,alpha=1,gamma=1,source_mag=25.,dist=1):
        
        analyser.__init__(self,z1=z1,z2=z2,M200=M200,Mstar=Mstar,c=c,Re=Re,
                          alpha=alpha,gamma=gamma,source_mag=source_mag,dist=dist)
        load_MICE.__init__(self,)
        load_COSMOS.__init__(self,)
        ReConc_loader.__init__(self,)
        
        
    def run_one(self,i,h=.7,buffer=.1):
        
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

       
        search_lim = self.get_search_range()[1]
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
                        stat = self.get_cross_section()
                        dist_real = stat[1]


                    except RuntimeError or KeyError:
                        print("Runtime Error.")
                        
                    if dist < dist_real:

                        eins_rad = stat[0]
                        pos_img1 = stat[2]
                        pos_img2 = stat[3]
                        mu_img1 = stat[4]
                        mu_img2 = stat[5]
                        Lens_arr = np.array([Mh_lens,Mstar_lens,zlens,conc,np.log10(scale_rad),eins_rad])
                        Source_arr = np.array([eins_rad,rmag_source,z2_source,pos_img1,pos_img2,mu_img1,mu_img2])
                        print("found one!")
                        print(Lens_arr)
                        print(Source_arr)

                        break


        return Lens_arr,Source_arr
    
    
    def run_mocks(self,):
        parr = np.linspace(1.,1.8,5)
        parc = [0.4,0.8,1,1.2,1.6,2]
        alph = np.linspace(.8,1.8,5)
        Alpha = np.repeat(alph,5)
        Beta = np.tile(parr,5)
        num_data = len(dat['z_cgal'])
        num_bus = 64
        bussize = int(num_data/num_bus)+1
        count = int(sys.argv[1])
        
        if count == num_bus-1:
            index = np.arange(count*bussize,num_data)
        else:
            index = np.arange(count*bussize,(count+1)*bussize,1)
        ind_cut = int(sys.argv[2])
        Lens_arr = []
        Source_arr = []
        
        for item in np.arange(num_data)[index]:
            alpha = Alpha[ind_cut]
            ratio = Beta[ind_cut]
            P,M = run_one(i=item,ratio=ratio,alpha=alpha)    
            if P[0] != 0.:
                Lens_arr.append(P)
                Source_arr.append(M)
        if len(Lens_arr)>0:
            dat1 = np.concatenate(Lens_arr,axis=None).reshape(-1,8)
            dat2 = np.concatenate(Source_arr,axis=None).reshape(-1,3)
            dirname = "/net/vdesk/data2/qzhou/GLens/lens_output4/bus"+str(count)
            if not os.path.exists(dirname):
                os.mkdir(dirname)


            dat1.tofile(dirname+'/r{}_alpha{}.bin'.format(alpha,ratio),format='f8')
            dat2.tofile(dirname+'/r{}_alpha{}_Msource.bin'.format(alpha,ratio),format='f8')

    
    

