from helper import *
from start  import *


class TestStart(unittest.TestCase):

    def test5(self):
        rho   = [0.1, 0.2, 0.3, 0.5, 0.8, 1.3]
        betas = [0.0001/0.001 * i for i in range(len(rho))]
        lambda_s, T1 = start(betas,rho,1.0)
        approxEqual(41.517752, lambda_s,1e-6)
        approxEqualVec(T1,\
          [ 1.9662112, 2.06616, 0.8135182, 0.632185, 0.5964, 0.649624 ],1e-5)




if __name__ == '__main__':
    
    unittest.main()
