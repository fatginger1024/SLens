import numpy as np
from ..models import gnfw, mass_size, concentration, Distances
from ..load_data import load_MICE

class _fxgx(gnfw):

    def __init__(self,N=100,gamma_num=5,write=False):
        super(gnfw,self).__init__()
        self._N = N
        self._x = np.logspace(-6,0,self._N)
        self._gamma_num = gamma_num
        self._gamma = np.linspace(.8,1.8,self._gamma_num)
        self._xx,self._gg = np.meshgrid(self._x,self._gamma)
        if write == True:
            dirbase = "../test_data/"
            self.fvals = self.fval()
            self.gvals = self.gval()
            self.fvals.tofile(dirbase+"fx.bin",format='f8')
            self.gvals.tofile(dirbase+"gx.bin",format='f8')
            
        else:
            dirbase = "../test_data/"
            self.fvals = np.fromfile(dirbase+"fx.bin",dtype="f8").reshape(self._gamma_num,self._N)
            self.gvals = np.fromfile(dirbase+"gx.bin",dtype="f8").reshape(self._gamma_num,self._N)
        

    def fval(self):
        fout = np.array([self.f(xxi,ggi) for xxi,ggi in zip(self._xx.flatten(),self._gg.flatten())]).reshape(self._gamma_num,self._N)
        return fout
        
    
    def gval(self):
        gout = np.array([self.g(xxi,ggi) for xxi,ggi in zip(self._xx.flatten(),self._gg.flatten())]).reshape(self._gamma_num,self._N)
        return gout
    
    


class _ReConc_generator(mass_size,concentration,load_MICE):
    
    def __init__(self,h=.7,write=False):
        load_MICE.__init__(self,)
        mass_size.__init__(self,)
        concentration.__init__(self,)
        self.h = h
        if write == True:
            self.data_gen()
        
      
    
    def data_gen(self,):
        
        scale_rad = np.zeros(len(self.Mstar_arr))
        conc_arr = np.zeros(len(self.Mstar_arr))

        for i in range(len(self.Mstar_arr)):
            scale_rad[i] =  self.get_Re_from_Mstar(10**self.Mstar_arr[i]/self.h)
            conc_arr[i] = self.concentration(10**self.Mh_arr[i],self.zlens_arr[i])
            
        dirbase = "SLens/test_data/"
        np.c_[scale_rad,conc_arr].tofile(dirbase+"RadConc.bin",format="f8")

        return scale_rad,conc
    
    
class _tables(Distances):
    
    def __init__(self,write=False):
        Distances.__init__(self,)
        if write == True:
            self.Da1_dict,self.Sigmacr_dict,self.rhoz_dict = self.get_tables()
            self.write()
            
        else:
            self.Da1_dict,self.Sigmacr_dict,self.rhoz_dict = self.load()

                
    def write(self,):
        dirbase = "../test_data/"
        f1 = open(dirbase+"Da1.pkl","wb")
        pickle.dump(self.Da1_dict,f1)
        f1.close()

        f2 = open(dirbase+"Sigmacr.pkl","wb")
        pickle.dump(self.Sigmacr_dict,f2)
        f2.close()

        f3 = open(dirbase+"rhoz.pkl","wb")
        pickle.dump(self.rhoz_dict,f3)
        f3.close()
        
    def load(self,):
        dirbase = "../test_data/"
        with open(dirbase+"Da1.pkl", "rb") as f1:
            Da1_dict = pickle.load(f1)
            f1.close()
        with open(dirbase+"Sigmacr.pkl", "rb") as f2:
            Sigmacr_dict = pickle.load(f2)
            f2.close()
        with open(dirbase+"rhoz.pkl", "rb") as f3:
            rhoz_dict = pickle.load(f3)
            f3.close()  
            
        return Da1_dict,Sigmacr_dict,rhoz_dict

        
                    
    def get_tables(self): 
        Da1_dict = {}
        Sigmacr_dict = {}
        rhoz_dict = {}
        z1 = np.round(np.arange(0.01,1+0.01,.01),2)
        z2 = np.round(np.arange(0.01,6+0.01,.01),2)
        for i,z1i in enumerate(z1):
            for j,z2j in enumerate(z2):
                if z1i < z2j:
                    d = Distances(z1=z1i,z2=z2j)
                    Da1 = d._Da1
                    Sigma_cr = d._Sigmacr  #Msun Mpc^-2
                    rhoz = d._rhoz
                    Da1_dict[str(z1i)+'+'+str(z2j)] = Da1
                    Sigmacr_dict[str(z1i)+'+'+str(z2j)] = Sigma_cr
                    rhoz_dict[str(z1i)+'+'+str(z2j)] = rhoz

        return Da1_dict, Sigmacr_dict, rhoz_dict
    
