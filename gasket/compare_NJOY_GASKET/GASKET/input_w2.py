import sys
import numpy as np
sys.path.append('../../python/');              from gasket import *
sys.path.append('../../phononDistributions/'); from waterDataContinuous import *
sys.path.append('../');                        from writeOutput import *

alphas = [0.01,0.1,1.0,5.0,10.0,20.0]
betas = list(np.linspace(0.0,10,101))
t      = list(np.linspace(0,3000,3001))

T = 296.6*8.617333e-5

sab,H,F = gasket(X,Q,t,alphas,betas,T,continWgt,freeGasWgt,oscWgt,oscEnergies,oscWeights)
fileName = 'gasket_'+sys.argv[0].split('_')[-1]
writeOutput(fileName,alphas,betas,sab,H,F,T,t,title)
