import unittest
from scint import *
 
class Test_SCINT(unittest.TestCase):
 
    def test_numbers_3_4(self):
        T = 0.0255
        t = [i*0.1 for i in range(80)]
        GC = [ 9.3955738, 9.3953535, 9.3946929, 9.3935919, 9.3920506, 9.3900693, \
         9.3876482, 9.3847876, 9.3814877, 9.3777490, 9.3735719, 9.3689570, \
         9.3639046, 9.3584154, 9.3524901, 9.3461293, 9.3393338, 9.3321043, \
         9.3244417, 9.3163469, 9.3078208, 9.2988643, 9.2894786, 9.2796647, \
         9.2694237, 9.2587568, 9.2476652, 9.2361503, 9.2242133, 9.2118557, \
         9.1990787, 9.1858840, 9.1722730, 9.1582473, 9.1438085, 9.1289582, \
         9.1136981, 9.0980301, 9.0819559, 9.0654773, 9.0485962, 9.0313145, \
         9.0136343, 8.9955575, 8.9770861, 8.9582224, 8.9389685, 8.9193264, \
         8.8992985, 8.8788871, 8.8580944, 8.8369229, 8.8153749, 8.7934528, \
         8.7711591, 8.7484964, 8.7254672, 8.7020741, 8.6783197, 8.6542067, \
         8.6297379, 8.6049159, 8.5797436, 8.5542238, 8.5283592, 8.5021529, \
         8.4756078, 8.4487267, 8.4215128, 8.3939689, 8.3660983, 8.3379039, \
         8.3093890, 8.2805566, 8.2514099, 8.2219522, 8.1921867, 8.1621168, \
         8.1317456, 8.1010765 ]
        GS = [0.0, 0.0448291, 0.0896540, 0.1344704, 0.1792741, 0.2240610, 0.2688268, \
         0.3135672, 0.3582781, 0.4029552, 0.4475944, 0.4921915, 0.5367422, \
         0.5812424, 0.6256880, 0.6700746, 0.7143983, 0.7586548, 0.8028399, \
         0.8469497, 0.8909798, 0.9349264, 0.9787851, 1.0225520, 1.0662230, \
         1.1097941, 1.1532611, 1.1966202, 1.2398672, 1.2829981, 1.3260091, \
         1.3688961, 1.4116552, 1.4542826, 1.4967742, 1.5391262, 1.5813347, \
         1.6233959, 1.6653061, 1.7070613, 1.7486578, 1.7900919, 1.8313599, \
         1.8724580, 1.9133825, 1.9541299, 1.9946964, 2.0350785, 2.0752727, \
         2.1152753, 2.1550828, 2.1946918, 2.2340988, 2.2733004, 2.3122931, \
         2.3510736, 2.3896386, 2.4279847, 2.4661087, 2.5040073, 2.5416773, \
         2.5791155, 2.6163188, 2.6532841, 2.6900082, 2.7264882, 2.7627210, \
         2.7987037, 2.8344333, 2.8699069, 2.9051217, 2.9400748, 2.9747634, \
         3.0091848, 3.0433362, 3.0772151, 3.1108187, 3.1441444, 3.1771898, \
         3.2099522]
        EPS = [0.0, 0.00204, 0.00408, 0.00612, 0.00816, 0.01020, 0.01224, 0.01428, \
         0.01632, 0.01836, 0.02040, 0.02244, 0.02448, 0.02652, 0.02856, 0.03060, \
         0.03264, 0.03468, 0.03672, 0.03876, 0.04080, 0.04284, 0.04488, 0.04692, \
         0.04896, 0.05100, 0.05304, 0.05508, 0.05712, 0.05916, 0.06120, 0.06324, \
         0.06528, 0.06732, 0.06936, 0.07140, 0.07344, 0.07548, 0.07752, 0.07956,\
         0.08160, 0.08364, 0.08568, 0.08772, 0.08976, 0.09180, 0.09384, 0.09588,\
         0.09792, 0.09996, 0.10200, 0.10404, 0.10608, 0.10812, 0.11016, 0.11220, \
         0.11424, 0.11628, 0.11832, 0.12036, 0.12240, 0.12444, 0.12648, 0.12852, \
         0.13056, 0.13260, 0.13464, 0.13668, 0.13872, 0.14076, 0.14280, 0.14484, \
         0.14688, 0.14892, 0.15096, 0.15300, 0.15504, 0.15708, 0.15912, 0.16116, \
         0.16320, 0.16524, 0.16728, 0.16932, 0.17136, 0.17340, 0.17544, 0.17748, \
         0.17952, 0.18156, 0.18360, 0.18564, 0.18768, 0.18972, 0.19176, 0.19380, \
         0.19584, 0.19788, 0.19992, 0.20196]
        A = 1.4166678E-4
        B = 2.5280929E-3
        F = 0.3108374;
        S = SCINT(t,GC,GS,EPS,T,A,B,F)
        correctS = [ 0.0561771, 0.0562363, 0.0563434, 0.0564458, 0.0565435, \
        0.0566365, 0.0567247, 0.0568082, 0.0568869, 0.0569608, 0.0570299, \
        0.0570942, 0.0571537, 0.0572083, 0.0572581, 0.0573031, 0.0573432, \
        0.0573785, 0.0574089, 0.0574344, 0.0574551, 0.0574709, 0.0574818, \
        0.0574878, 0.0574889, 0.0574852, 0.0574766, 0.0574631, 0.0574447, \
        0.0574215, 0.0573934, 0.0573604, 0.0573226, 0.0572799, 0.0572324, \
        0.0571800, 0.0571228, 0.0570608, 0.0569939, 0.0569223, 0.0568458, \
        0.0567646, 0.0566786, 0.0565879, 0.0564924, 0.0563922, 0.0562873, \
        0.0561777, 0.0560634, 0.0559444, 0.0558208, 0.0556926, 0.0555598, \
        0.0554224, 0.0552805, 0.0551340, 0.0549829, 0.0548274, 0.0546674, \
        0.0545030, 0.0543342, 0.0541609, 0.0539833, 0.0538014, 0.0536151, \
        0.0534246, 0.0532298, 0.0530308, 0.0528275, 0.0526202, 0.0524087, \
        0.0521930, 0.0519734, 0.0517497, 0.0515219, 0.0512903, 0.0510547, \
        0.0508152, 0.0505718, 0.0503247, 0.0500737, 0.0498190, 0.0495606, \
        0.0492985, 0.0490328, 0.0487636, 0.0484907, 0.0482144, 0.0479346, \
        0.0476514, 0.0473648, 0.0470749, 0.0467817, 0.0464853, 0.0461856, \
        0.0458828, 0.0455769, 0.0452679, 0.0449559, 0.0446409 ]


        np.testing.assert_almost_equal(S,correctS, 6)
 
if __name__ == '__main__':
    unittest.main()
