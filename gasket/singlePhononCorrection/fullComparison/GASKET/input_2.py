import sys
import numpy as np
sys.path.append('../../../newGasket/');              from main import *
sys.path.append('../../../../phononDistributions/'); from waterData import *
sys.path.append('../');                        from writeOutput import *

t      = list(np.linspace(0,96.2,3001))
alphas = [0.1]
betas = list(np.linspace(0.0,20,201))

T = 296.0*8.617333e-5; invT = 1.0/T
rhoBetas = [x*invT for x in X]

sab,H,F = simpleGASKET(rhoBetas,Q,t,alphas,betas) 
fileName = 'gasket_'+sys.argv[0].split('_')[-1]

writeOutput(fileName,alphas,betas,sab,H,F,T,t,title,'None','None')
