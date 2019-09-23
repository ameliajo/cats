import sys
sys.path.append('../../python/')
from gasket import *

alphas = [0.001*i  for i in range(0,15)]
alphas[0] = 0.0005
betas  = [0.02*i for i in range(0,51)]
print(alphas)
print()
print(betas)
X = [0.0025*i for i in range(69)]
Q = [0.0, 0.00049, 0.00098, 0.00191, 0.00338, 0.00485, 0.00721, 0.00966, 0.01245, 0.01588, 0.01931, 0.02353, 0.02794, 0.0326, 0.03799, 0.04338, 0.04884, 0.05433, 0.05994, 0.06622, 0.07249, 0.07971, 0.08755, 0.09539, 0.10324, 0.11108, 0.1173, 0.1205, 0.12158, 0.12081, 0.11662, 0.11015, 0.10426, 0.09838, 0.09403, 0.08995, 0.08616, 0.08302, 0.07988, 0.07745, 0.07526, 0.07313, 0.07129, 0.06945, 0.06783, 0.06634, 0.06484, 0.06327, 0.06171, 0.06016, 0.05863, 0.0, 0.05588, 0.05467, 0.05356, 0.05258, 0.0516, 0.05055, 0.04949, 0.0484, 0.04726, 0.04613, 0.04502, 0.04392, 0.04299, 0.04256, 0.04213, 0.02471, 0.0]
t = [0.1 * i for i in range(300)]
T = 296*8.617333e-5
continWgt = 0.444444
freeGasWgt = 0.055556
oscWgt = 1.0-continWgt-freeGasWgt
oscEnergies = [0.205   , 0.48    ]
oscWeights  = [0.166667, 0.333333]

X = X[1:]
Q = Q[1:]
sab = gasket(X,Q,t,alphas,betas,T,continWgt,freeGasWgt,oscWgt,oscEnergies,oscWeights)

f = open("gasket_output.py",'w')
f.write("gasket_alphas = ")
f.write(str(alphas))
f.write("\n")
f.write("gasket_betas = ")
f.write(str(betas))
f.write("\n")
f.write("gasket_sab = ")
f.write(str(sab))
f.close()

