import sys
sys.path.append("..")
from phononExpansion import *
import unittest
import pytest

class TestLambda_s(unittest.TestCase):

    def test_lambda_s(self):
        delta = 0.03
        rho = [-(i-20)**2+400 for i in range(40)] + [10 for i in range(40,47)] + [-(i-70)**2+500 for i in range(48,93)] + [0.0]
        betaGrid = [delta*i for i in range(len(rho))]
        lambda_s = getLambda_s(rho,betaGrid)
        print(lambda_s)

        assert True == True







if __name__ == '__main__':
    unittest.main()
