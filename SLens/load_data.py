import numpy as np
import pickle
from astropy.table import Table 
    
class tables():

    def __init__(self,):

        self.Da1_dict,self.Sigmacr_dict,self.rhoz_dict = self.load()

    def load(self,):
        dirbase = "SLens/test_data/"
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


class fxgx():

    def __init__(self,N=100,gamma_num=5):
        
        self._N = N
        self._x = np.logspace(-6,0,self._N)
        self._gamma_num = gamma_num
        self._gamma = np.linspace(.8,1.8,self._gamma_num)
        self._xx,self._gg = np.meshgrid(self._x,self._gamma)
        
        dirbase = "SLens/test_data/"
        self.fvals = np.fromfile(dirbase+"fx.bin",dtype="f8").reshape(self._gamma_num,self._N)
        self.gvals = np.fromfile(dirbase+"gx.bin",dtype="f8").reshape(self._gamma_num,self._N)
        
        
        
class ReConc_loader():
    
    def __init__(self,):
        dirbase = "SLens/test_data/"
        RadConc = np.fromfile(dirbase+"RadConc.bin",dtype="f8").reshape(-1,2)
        self.Re_arr = RadConc[:,0]
        self.Conc_arr = RadConc[:,1]
        
    

class load_MICE():
    def __init__(self,data="128"):    
        dirbase = "SLens/test_data/"
        if data == "128":
            dat = Table.read(dirbase+"9981.fits",format="fits")
        else:
            dat = Table.read(dirbase+"11601.fits",format="fits")
        self.Mstar_arr = np.array(dat['lmstellar'])
        self.Mh_arr = np.array(dat['lmhalo'])
        self.zlens_arr = np.array(dat['z_cgal'])
        self.ra_arr = np.array(dat['ra_gal'])
        self.dec_arr = np.array(dat['dec_gal'])


class load_COSMOS():
    def __init__(self,):
        dirbase = "SLens/test_data/"
        COSMOS2015 = Table.read(dirbase+'COSMOS2015_Laigle+_v1.1.fits',format='fits')
        self.COSMOS_rmag = np.array(COSMOS2015['r_MAG_APER2'])
        self.COSMOS_redshift = np.array(COSMOS2015['ZPDF'])
        mask = (self.COSMOS_redshift>=0)  *  (self.COSMOS_rmag>0)
        self.TOTnum = len(self.COSMOS_redshift[mask])*2500
        self.source_generator = np.c_[self.COSMOS_redshift[mask],self.COSMOS_rmag[mask]]


       
       

