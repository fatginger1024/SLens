from SLens import analyser
import numpy as np
from math import *
import unittest
from scipy.optimize import fmin


class tests(unittest.TestCase):
        
    def test_attr1(self):
        method = analyser()
        lens_eqn = method.lens_beta
        x_eins = method.attr[0]/(method.thetas*206265)
        self.assertAlmostEqual(method.lens_alpha(x_eins)[0]/x_eins,1,places=5)
        self.assertAlmostEqual(method.attr[2]/(lens_eqn(x_eins)*method.thetas*206265),1,places=5)
        
    def test_attr2(self):
        method = analyser()
        x_radial = method.attr[1]/(method.thetas*206265)
        lens_eqn = method.lens_beta
        xtol=1e-5
        x = np.array([x_radial-xtol,x_radial,x_radial+xtol])
        y = lens_eqn(x)
        self.assertTrue(np.argmin(y),1)
        self.assertAlmostEqual(method.attr[3]/abs(lens_eqn(x_radial)[0]*method.thetas*206265),1,places=5)
        
 
        
if __name__ == '__main__':
    unittest.main()
