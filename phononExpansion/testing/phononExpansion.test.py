import sys
sys.path.append("../")
from phononExpansion import *
import numpy as np
import unittest
import pytest

class TestConvolution(unittest.TestCase):

    def test_convolution(self):

        nphon = 20
        delta = 0.0255
        rho = [ 0, .0005, .001, .002, .0035, .005, .0075, .01, .013, .0165,    \
               .02, .0245, .029, .034, .0395, .045, .0506, .0562, .0622, .0686,\
               .075, .083, .091, .099, .107, .115, .1197, .1214, .1218, .1195, \
               .1125, .1065, .1005, .09542, .09126, .0871, .0839, .0807,       \
               .07798, .07574, .0735, .07162, .06974, .06804, .06652, .065,    \
               .0634, .0618, .06022, .05866, .0571, .05586, .05462, .0535,     \
               .0525, .0515, .05042, .04934, .04822, .04706, .0459, .04478,    \
               .04366, .04288, .04244, .042, 0. ]
        betas = list(np.linspace(0,15,16))
        alphas = [1e-3, 1e-1] + list(np.linspace(1,10,10))

        sab = contin(nphon, delta, rho, alphas, betas)
        print(sab[:16])

        #for i in range(len(t3)):
        #    assert t3[i] == pytest.approx(correct_t3[i])






if __name__ == '__main__':
    unittest.main()
