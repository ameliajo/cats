from helper import *
from convolution import *


class TestConvol(unittest.TestCase):

    def test6(self):
        t1 = [0.2, 0.6, 0.8, 2.0, 6.0, 8.0]
        t2 = [0.2, 0.6, 0.8, 2.0, 6.0, 8.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        delta = 0.03
        nn = 18
        t3 = convol(t1, t2, delta, nn)
        approxEqualVec(t3, \
          [3.8459762, 2.6993367, 1.0195307, 0.53364442, 0.37281623, \
        0.384, 0.624, 1.008, 1.8, 2.16, 0.96, 0, 0, 0, 0, 0, 0, 0], 1e-6)




if __name__ == '__main__':
    
    unittest.main()
