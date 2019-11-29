import sys
import numpy as np
sys.path.append('../../../../singlePhonon')
sys.path.append('../../../../phononExpansion/help')
sys.path.append('../../../newGasket')
sys.path.append('../../');              from plot_correction_effects import *
sys.path.append('../../../../phononDistributions/'); from waterData import *
sys.path.append('../');                        from writeOutput import *

t      = list(np.linspace(0,96.2,3001))
alphas = [0.0001]
betas = list(np.linspace(0.0,20,201))

T = 296.0*8.617333e-5; invT = 1.0/T
rhoBetas = [x*invT for x in X]

sab,H,F,oscBegin,correctionBegin = simpleGASKET(rhoBetas,Q,t,alphas,betas,True) 
fileName = 'gasket_corrected_'+sys.argv[0].split('_')[-1]

oscBegin        = 'None' if oscBegin        == None else str(float(oscBegin))
correctionBegin = 'None' if correctionBegin == None else str(float(correctionBegin))

writeOutput(fileName,alphas,betas,sab,H,F,T,t,title,oscBegin,correctionBegin)
