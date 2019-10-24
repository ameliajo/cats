import sys
import numpy as np
sys.path.append('../../newGasket/');           from main import *
sys.path.append('../../phononDistributions/'); from waterDataContinuous import *
sys.path.append('../');                        from writeOutput import *

alphas = [0.001,0.01,0.5,1.0,5.0,10.0,20.0]
betas = list(np.linspace(0.0,20,401))
t      = list(np.linspace(0,96.3574,8001))

T = 296.6*8.617333e-5; invT = 1.0/T
rhoBetas = [x*invT for x in X]

sab,H,F = simpleGASKET(rhoBetas,Q,t,alphas,betas) 
fileName = 'simple_'+sys.argv[0].split('_')[-1]
writeOutput(fileName,alphas,betas,sab,H,F,T,t,title)
