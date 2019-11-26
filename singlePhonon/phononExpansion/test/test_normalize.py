from helper import *
from normalize import *


class TestNormalize(unittest.TestCase):

    def test4(self):
        P = [0.02*i for i in range(1,7)]
        betas = [0.5*i for i in range(len(P))]
        wgt = 0.8
        delta = 0.5
        P = normalize(betas,P,wgt)
        approxEqualVec(\
            [2.62155807e-2, 5.24311614e-2, 7.86467421e-2, 0.10486322, \
             0.1310779, 0.15729348], P, 1e-5)

    

if __name__ == '__main__':
    
    unittest.main()
