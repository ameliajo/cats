import numpy as np
rho = [0.0000,3.8183e-01,1.0796,1.6625,1.9735,2.0196,1.8515,1.5652,1.3384,1.2414,1.2130,1.2339,1.2963,1.3604,1.4245,1.5091,1.5999,1.6795,1.7703,1.8850,2.0024,2.1078,2.2058,2.2969,2.3736,2.4315,2.4631,2.4630,2.4420,2.4283,2.4497,2.5256,2.6468,2.8022,2.9994,3.2493,3.5338,3.8369,4.1585,4.5100,4.8879,5.2840,5.6880,6.0939,6.5062,6.9240,7.3146,7.6278,7.8765,8.1019,8.3148,8.4660,8.5106,8.4430,8.3176,8.2258,8.1672,8.0634,7.8788,7.7084,7.5812,7.4298,7.2166,7.0165,6.8529,6.6832,6.4839,6.2897,6.1278,5.9917,5.8309,5.6189,5.3965,5.2295,5.0961,4.9383,4.7503,4.5634,4.3927,4.2340,4.0757,3.9060,3.7421,3.5991,3.4636,3.3199,3.1826,3.0533,2.9071,2.7299,2.5375,2.3434,2.1513,1.9532,1.7492,1.5582,1.3939,1.2399,1.1004,9.9699e-01,9.2203e-01,8.4542e-01,7.7743e-01,7.3029e-01,6.8805e-01,6.4727e-01,6.1932e-01,5.9475e-01,5.6481e-01,5.4473e-01,5.3324e-01,5.1494e-01,4.9853e-01,4.9162e-01,4.8152e-01,4.6689e-01,4.5969e-01,4.5248e-01,4.4006e-01]
spacing = 0.001265
spacing = 0.00255
energies = [i*spacing for i in range(len(rho))]
invArea = 1.0/np.trapz(rho,x=energies)
rho = [x*invArea for x in rho]

fullE = [-x for x in energies[::-1]]+energies[1:]
fullRho = rho[::-1]+rho[1:]


def getVal(energies,val,vec):
    for i in range(len(energies)-1):
        if energies[i] <= val and val <= energies[i+1]:
            m = (vec[i+1]-vec[i])/(energies[i+1]-energies[i])
            b = vec[i]
            return m*(val-energies[i])+b
    return 0.0
 