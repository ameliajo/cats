import unittest
import sys
sys.path.append('../')
from bessl import *
 
class Test_BESSL(unittest.TestCase):
 
    def test1(self):
        B = [0.0]*10

        # given a small x value (xVal < 0.05)
        # when results of one is fed into the next (round 1)
        xVal = 7.4491602E-5;
        nMax = bessl(xVal,B,10);
        correctB = [ 1, 3.72458E-5, 6.93624E-10, 8.61153E-15, 8.01859E-20, 0, 0, 0, 0, 0 ]
        np.testing.assert_almost_equal(B,correctB,6)
        np.testing.assert_equal(nMax,5)

        xVal = 2.895559E-7
        nMax = bessl(xVal,B,10);
        correctB = [1.0, 0.144778E-06, 0.104803E-13, 0, 0, 0, 0, 0, 0, 0]
        np.testing.assert_almost_equal(B,correctB,6)
        np.testing.assert_equal(nMax,3)

    def test2(self):

        B = [1.0, 0.143330E-4, 0.102718E-9, 0.490751E-15, 0, 0, 0, 0, 0, 0]
        xVal = 7.4491602E-3;
        nMax = bessl(xVal,B,10);
        correctB = [ 1.0, 0.372458E-2, 0.693625E-5, 0.861154E-8, 0.801859E-11,\
                   0.597318E-14, 0.370793E-17, 0, 0, 0 ]
        np.testing.assert_almost_equal(B,correctB,6)
        np.testing.assert_equal(nMax,7)

        xVal = 2.895559E-5
        nMax = bessl(xVal,B,10);
        correctB = [1.0, 0.144778E-4, 0.104803E-9, 0.505773E-15, 0, 0, 0, 0, 0, 0]
        np.testing.assert_almost_equal(B,correctB,6)
        np.testing.assert_equal(nMax,4)


    def test3(self):
        B = [0.0]*10
        xVal = 0.05
        nMax = bessl(xVal,B,10)
        correctB = [1.000625, 2.500781E-2, 3.125651E-4, 2.604574E-6, 1.627808E-8,\
           8.138869E-11, 3.391145E-13, 1.211110E-15, 3.784685E-18, 1.051294E-20]
        np.testing.assert_almost_equal(B,correctB,6)
        np.testing.assert_equal(nMax,0)

        nMax = bessl(xVal,B,10)
        correctB = [1.000625, 2.500781E-2, 3.125651E-4, 2.604574E-6, 1.627808E-8,\
            8.138869E-11, 3.391145E-13, 1.211110E-15, 3.784685E-18, 1.051294E-20]
        np.testing.assert_almost_equal(B,correctB,6)
        np.testing.assert_equal(nMax,0)


    
if __name__ == '__main__':
    unittest.main()
