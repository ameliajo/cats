import unittest
import sys
sys.path.append('../')
from ftrans import *
 
class Test_FTRANS(unittest.TestCase):
 
    def test1(self):
        T = 0.0255
        t = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10] 
        X = [ 0.1, 0.2, 0.3, 0.5, 0.8 ]
        Q = [ 0.1, 0.3, 0.7, 0.4, 0.0 ]
        correctPC = [ 0.0, 0.752123, 0.752096, 0.752059, 0.752011, 0.751953, \
                      0.751883, 0.751804, 0.751713, 0.751612]
        correctPS = [ 0.0, 0.00519993, 0.00779978, 0.0103995, 0.0129990, \
                      0.0155982, 0.0181972, 0.0207958, 0.0233940, 0.0259917]
        PC, PS = ftrans(T, X, Q, t)
        np.testing.assert_almost_equal(PC,correctPC,6)
        np.testing.assert_almost_equal(PS,correctPS,6)

    def test2(self):
        T = 0.0255
        t = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10] 
        X = [ 0.1, 0.7, 0.10, 0.11, 0.20 ]
        Q = [ 0.0, 0.5, 0.0,  0.4,  0.0  ]
        correctPC = [ 0.0, 0.186750, 0.186749, 0.186748, 0.186746, 0.186745, \
                      0.186742, 0.186740, 0.186736, 0.186733 ]
        correctPS = [ 0.0, 0.496969E-3, 0.745452E-3, 0.993934E-3, 0.00124241,\
                      0.00149089, 0.00173937, 0.00198784, 0.0022363, 0.00248477]
        PC, PS = ftrans(T, X, Q, t);
        np.testing.assert_almost_equal(PC,correctPC,6)
        np.testing.assert_almost_equal(PS,correctPS,6)



        #np.testing.assert_almost_equal( 5.61770521E-2,sterp(-0.01, betas, sLog),6)



    
if __name__ == '__main__':
    unittest.main()
