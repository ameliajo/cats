import sys
import numpy as np
sys.path.append('../../quadGasket/');          from main import *
#sys.path.append('../../quadGasket/');          from quadGasket import *
sys.path.append('../../phononDistributions/'); from waterDataContinuous import *
sys.path.append('../');                        from writeOutput import *

alphas = [1e-5,0.001,0.01]#,0.5,1.0,5.0,10.0,20.0]
#alphas = [0.01]
betas = list(np.linspace(0.0,5,101))
t     = list(np.linspace(0,96.3574,1001))

T = 296.6*8.617333e-5; invT = 1.0/T
rhoBetas = [x*invT for x in X]

sab,H,F = quadGASKET(rhoBetas,Q,t,alphas,betas) 
#sab,H,F = quadGASKET(alphas,betas,T,rhoBetas,Q,t) 
fileName = 'quad_'+sys.argv[0].split('_')[-1]
writeOutput(fileName,alphas,betas,sab,list(H),list(F),T,t,title)

