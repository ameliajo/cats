import unittest
import sys
sys.path.append('../')
from gtg import *
 
class Test_GTG(unittest.TestCase):
 
    def test1(self):
        #np.testing.assert_almost_equal( 5.61770521E-2,sterp(-0.01, betas, sLog),6)
        wgt = 0.444444
        T = 0.0255
        AM = 0.991409;
        X = [0.0025, 0.0050, 0.0075, 0.0100, 0.0125, 0.0150, 0.0175, 0.0200, 0.0225, \
         0.0250, 0.0275, 0.0300, 0.0325, 0.0350, 0.0375, 0.0400, 0.0425, 0.0450, \
         0.0475, 0.0500, 0.0525, 0.0550, 0.0575, 0.0600, 0.0625, 0.0650, 0.0675, \
         0.0700, 0.0725, 0.0750, 0.0775, 0.0800, 0.0825, 0.0850, 0.0875, 0.0900, \
         0.0925, 0.0950, 0.0975, 0.1000, 0.1025, 0.1050, 0.1075, 0.1100, 0.1125, \
         0.1150, 0.1175, 0.1200, 0.1225, 0.1250, 0.1275, 0.1300, 0.1325, 0.1350, \
         0.1375, 0.1400, 0.1425, 0.1450, 0.1475, 0.1500, 0.1525, 0.1550, 0.1575, \
         0.1600, 0.1625, 0.1650, 0.16750]
        Q = [4.9E-4, 9.8E-4, 1.91E-3, 3.38E-3, 4.85E-3, 7.21E-3, 9.66E-3, 1.245E-2, \
         1.588E-2, 1.931E-2, 2.353E-2, 2.794E-2, 3.26E-2, 3.799E-2, 4.338E-2, \
         4.884E-2, 5.433E-2, 5.994E-2, 6.622E-2, 7.249E-2, 7.971E-2, 8.755E-2, \
         9.539E-2, 1.0324E-1, 1.1108E-1, 1.173E-1, 1.205E-1, 1.2158E-1, \
         1.2081E-1, 1.1662E-1, 1.1015E-1, 1.0426E-1, 9.838E-2, 9.403E-2, \
         8.995E-2, 8.616E-2, 8.302E-2, 7.988E-2, 7.745E-2, 7.526E-2, 7.313E-2,\
         7.129E-2, 6.945E-2, 6.783E-2, 6.634E-2, 6.484E-2, 6.327E-2, 6.171E-2, \
         6.016E-2, 5.863E-2, 0.E+00, 5.588E-2, 5.467E-2, 5.356E-2, 5.258E-2, \
         5.16E-2, 5.055E-2, 4.949E-2, 4.84E-2, 4.726E-2, 4.613E-2, 4.502E-2, \
         4.392E-2, 4.299E-2, 4.256E-2, 4.213E-2, 2.471E-2]
        t = [i*0.1 for i in range(80)]
        correctPC = [9.39557, 9.39535, 9.39469, 9.39359, 9.39205, 9.39006, \
        9.38764, 9.38478, 9.38148, 9.37774, 9.37357, 9.36895, 9.36390, 9.35841, \
        9.35249, 9.34612, 9.33933, 9.33210, 9.32444, 9.31634, 9.30782, 9.29886, \
        9.28947, 9.27966, 9.26942, 9.25875, 9.24766, 9.23615, 9.22421, 9.21185, \
        9.19907, 9.18588, 9.17227, 9.15824, 9.14380, 9.12895, 9.11369, 9.09803, \
        9.08195, 9.06547, 9.04859, 9.03131, 9.01363, 8.99555, 8.97708, 8.95822, \
        8.93896, 8.91932, 8.89929, 8.87888, 8.85809, 8.83692, 8.81537, 8.79345, \
        8.77115, 8.74849, 8.72546, 8.70207, 8.67832, 8.65420, 8.62973, 8.60491, \
        8.57974, 8.55422, 8.52835, 8.50215, 8.47560, 8.44872, 8.42151, 8.39396, \
        8.36609, 8.33790, 8.30938, 8.28055, 8.25141, 8.22195, 8.19218, 8.16211, \
        8.13174, 8.10107]
        correctPS = [ 0, 4.48291E-2, 8.96539E-2, 0.134470, 0.179274, 0.2240610, \
        0.2688268, 0.3135672, 0.3582781, 0.4029552, 0.4475944, 0.4921915, \
        0.5367422, 0.5812424, 0.6256880, 0.6700746, 0.7143983, 0.7586548, \
        0.8028399, 0.8469497, 0.8909798, 0.9349264, 0.9787851, 1.022552, \
        1.066223, 1.109794, 1.153261, 1.196620, 1.239867, 1.282998, 1.326009, \
        1.368896, 1.411655, 1.454283, 1.496774, 1.539126, 1.581335, 1.623396, \
        1.665306, 1.707061, 1.748658, 1.790092, 1.831360, 1.872458, 1.913383, \
        1.954130, 1.994696, 2.035079, 2.075273, 2.115275, 2.155083, 2.194692, \
        2.234099, 2.273300, 2.312293, 2.351074, 2.389639, 2.427985, 2.466109, \
        2.504007, 2.541677, 2.579116, 2.616319, 2.653284, 2.690008, 2.726488, \
        2.762721, 2.798704, 2.834433, 2.869907, 2.905122, 2.940075, 2.974763, \
        3.009185, 3.043336, 3.077215, 3.110819, 3.144144, 3.177190, 3.209952]

        tbar, PC, PS = GTG( wgt, T, AM, X, Q, t )

        np.testing.assert_almost_equal(tbar,4.9117341E-2,6)
        np.testing.assert_almost_equal(PC,correctPC,5)
        np.testing.assert_almost_equal(PS,correctPS,5)





    
if __name__ == '__main__':
    unittest.main()