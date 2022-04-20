from SLens import Sersic
from scipy.integrate import quad
from math import *
import unittest


class tests(unittest.TestCase):
        
    def test_mass1(self):
        galaxy = Sersic()
        self.assertAlmostEqual(galaxy.M2d(galaxy._Re)/(galaxy._Mstar/2),1,places=2)

    def test_mass2(self):
        galaxy = Sersic()
        inte_func = lambda R: galaxy.Sigma(R)*R
        inte_out = quad(inte_func,0,galaxy._Re)[0]
        self.assertAlmostEqual(2*pi*inte_out/(galaxy._Mstar/2),1,places=2)
 
    
if __name__ == '__main__':
    unittest.main()
