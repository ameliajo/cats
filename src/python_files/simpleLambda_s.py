import numpy as np
from testRho import *

def getVal(energies,val,vec):
    for i in range(len(energies)-1):
        if energies[i] <= val and val <= energies[i+1]:
            m = (vec[i+1]-vec[i])/(energies[i+1]-energies[i])
            b = vec[i]
            return m*(val-energies[i])+b
    return 0.0
 
waterRho = rho[:]
rho = [0, 0.09, 0.16, 0.21, 0.24, 0.25, 0.24, 0.21, 0.16, 0.09, 0]
rho = [0,1.0,1.0,2.0,3.0,5.0,3.0,2.0,1.0,1.0,0.0]
rho = [0.0,2.0,4.0,6.0,8.0,9.0,7.0,5.0,3.0,1.0,0.0]
rho = waterRho[:]
betas = [0.1*i for i in range(len(rho))]
# simple water
rho = [0, .0005, .001, .002, .0035, .005, .0075, .01, .013, .0165, .02, .0245, .029, .034, .0395, .045, .0506, .0562, .0622, .0686, .075, .083, .091, .099, .107, .115, .1197, .1214, .1218, .1195, .1125, .1065, .1005, .09542, .09126, .0871, .0839, .0807, .07798, .07574, .0735, .07162, .06974, .06804, .06652, .065, .0634, .0618, .06022, .05866, .0571, .05586, .05462, .0535, .0525, .0515, .05042, .04934, .04822, .04706, .0459, .04478, .04366, .04288, .04244, .042, 0.]
# Be in BeO
rho = [ 3.100000E-6, 6.899985E-6, 1.579987E-5, 3.249975E-5, 5.624949E-5, 8.334943E-5, 1.141990E-4, 1.542986E-4, 1.998983E-4, 2.561974E-4, 3.232470E-4, 3.888969E-4, 4.642953E-4, 5.549946E-4, 6.520937E-4, 7.426447E-4, 8.467409E-4, 9.801895E-4, 1.118239E-3, 1.260638E-3, 1.427084E-3, 1.598535E-3, 1.781880E-3, 1.997277E-3, 2.222176E-3, 2.447475E-3, 2.693269E-3, 2.959668E-3, 3.249061E-3, 3.548913E-3, 3.871502E-3, 4.231850E-3, 4.647983E-3, 5.072545E-3, 5.508525E-3, 5.996425E-3, 6.571393E-3, 7.270280E-3, 8.078150E-3, 9.242100E-3, 1.058603E-2, 1.170882E-2, 1.260280E-2, 1.355699E-2, 1.489420E-2, 1.686302E-2, 1.985052E-2, 2.108993E-2, 1.985079E-2, 1.877167E-2, 1.770130E-2, 1.647055E-2, 1.539779E-2, 1.454566E-2, 1.401839E-2, 1.371271E-2, 1.347211E-2, 1.325255E-2, 1.314101E-2, 1.313275E-2, 1.303920E-2, 1.263506E-2, 1.216104E-2, 1.174264E-2, 1.108204E-2, 1.048161E-2, 1.010476E-2, 9.748349E-3, 9.503098E-3, 9.594495E-3, 9.700141E-3, 9.128369E-3, 8.656028E-3, 9.131014E-3, 9.605228E-3, 9.105205E-3, 8.259261E-3, 8.786475E-3, 1.086055E-2, 1.277622E-2, 1.621424E-2, 1.946884E-2, 2.976416E-2, 9.697986E-2, 1.552680E-1, 1.328978E-1, 1.235346E-1, 1.214221E-1, 8.868644E-2, 5.501578E-2, 3.047213E-2, 1.783501E-2, 8.844040E-3, 4.163439E-3, 2.958511E-3, 2.121391E-3, 2.761060E-3, 5.546849E-3, 9.687513E-3, 1.308667E-2, 1.468507E-2, 1.572126E-2, 1.673388E-2, 1.808334E-2, 2.032280E-2, 2.189195E-2, 2.199840E-2, 2.159049E-2, 2.082654E-2, 2.048958E-2, 2.079742E-2, 2.099254E-2, 2.091494E-2, 1.985058E-2, 1.765808E-2, 1.589056E-2, 1.491489E-2, 1.428881E-2, 1.395833E-2, 1.452532E-2, 1.367085E-2, 1.145167E-2, 1.047933E-2, 1.045589E-2, 1.133254E-2, 1.124611E-2, 9.129051E-3, 7.375213E-3, 6.719129E-3, 6.179936E-3, 5.740469E-3, 5.386990E-3, 5.107686E-3, 4.892018E-3, 4.627346E-3, 4.353081E-3, 4.173784E-3, 4.024547E-3, 3.690213E-3, 3.332776E-3, 3.125731E-3, 2.892161E-3, 2.661229E-3, 2.460275E-3, 2.270017E-3, 1.943553E-3, 1.456127E-3, 9.878817E-4, 3.999894E-4, 1.361800E-5 ]


kb = 8.6173332e-5

betas = [0.0999714*i for i in range(len(rho))]
betas = [0.0493192*i for i in range(len(rho))]
betas = [0.001*i/(kb*296.3) for i in range(len(rho))]
betas = [0.001*i/(kb*900) for i in range(len(rho))]

#rho = [0,0.526532,0.526532,1.05306,1.5796,2.63266,1.5796,1.05306,0.526532,0.526532,0]

#rho = [0,0.545441,0.969672,1.27269,1.45451,1.51511,1.45451,1.27269,0.969672,0.545441, 0]

doQuad = True
doQuad = False

if doQuad:
    integrationVal = 0.0
    N = 100
    points,weights = np.polynomial.legendre.leggauss(N)
    xs = [2*b/betas[-1]-1 for b in betas]
    for i in range(N):
        x, w = points[i], weights[i]
        b = betas[-1]*(x+1)*0.5
        integrationThing = betas[-1]*0.5
        rhoVal = getVal(xs,x,rho)
        integrationVal += w*rhoVal*integrationThing
    rho = [x/integrationVal for x in rho]
    print("\nintegration val",integrationVal,"\n")
else:
    invArea = 1.0/np.trapz(rho,x=betas)
    rho = [invArea*x for x in rho]

P = [rho[1]/(betas[1]*betas[1])] + [rho[i]/(2.0*betas[i]*np.sinh(betas[i]*0.5)) for i in range(1,len(rho))]
fullP = P[::-1] + P[1:]


fullBetas = [-b for b in betas][::-1] + betas[1:]

if doQuad:
    xs = [b/betas[-1] for b in fullBetas]
    integrationVal = 0.0
    N = 10
    points,weights = np.polynomial.legendre.leggauss(N)
    for i in range(N):
        x, w = points[i], weights[i]
        b = betas[-1]*x
        integrationThing = betas[-1]
        pVal = getVal(xs,x,fullP)
        print("------- ",pVal*np.exp(-b*0.5)*integrationThing)
        integrationVal += w*pVal*np.exp(-b*0.5)*integrationThing
    print(integrationVal)
else:
    integrand = [fullP[i]*np.exp(-fullBetas[i]*0.5) for i in range(len(fullP))]
    print(np.trapz(integrand,x=fullBetas))







