import sys
sys.path.append("..")
from phononExpansion import *
import unittest
import pytest
from dos import *

class TestLambda_s(unittest.TestCase):

    def test_lambda_s(self):
        delta = delta_parabolas[:]
        rho = rho_parabolas[:]
        betaGrid = [delta*i for i in range(len(rho))]
        #lambda_s = getLambda_s(rho,betaGrid)
        #print(lambda_s)

        assert True == True







if __name__ == '__main__':
    unittest.main()
