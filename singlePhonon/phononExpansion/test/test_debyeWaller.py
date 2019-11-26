from helper import *
from getDebyeWaller import *

class TestDebyeWaller(unittest.TestCase):

    def test2(self):
        P = [1,2,3,4,5,6,7,8,9,10,11,12]
        betas = [2.0*x for x in range(len(P))]
        approxEqual(1444532.840,getDebyeWaller(betas,P),1e-6)
        
    def test3(self):
        P = [0.01,0.02,0.03,0.04,0.05,0.06]
        betas = [0.1*x for x in range(len(P))]
        approxEqual(0.035514341, getDebyeWaller(betas,P),1e-6)



if __name__ == '__main__':
    
    unittest.main()
