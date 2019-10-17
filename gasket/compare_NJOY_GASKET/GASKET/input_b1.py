import sys
import numpy as np
sys.path.append('../../python/');              from gasket import *
sys.path.append('../../phononDistributions/'); from beoData import *
sys.path.append('../');                        from writeOutput import *

alphas = [0.1,0.5,1.0,5.0,10.0,15.0]
betas = list(np.linspace(0,10,801))
t      = list(np.linspace(0,6280,6280))

T = 296.6*8.617333e-5

sab,H,F = gasket(X,Q,t,alphas,betas,T,continWgt,freeGasWgt,oscWgt,oscEnergies,oscWeights)

fileName = 'gasket_'+sys.argv[0].split('_')[-1]
writeOutput(fileName,alphas,betas,sab,H,F,T,t,title)
"""
nFile = sys.argv[0].split('_')
f = open("gasket_output_"+nFile[-1],'w')
f.write("alphas = "+str(alphas)   +"\n")
f.write("betas  = "+str(betas)    +"\n")
f.write("sab    = "+str(sab)      +"\n")
f.write("H      = "+str(H)        +"\n")
f.write("F      = "+str(F)        +"\n")
f.write("temp   = "+str(T)        +"\n")
f.write("num_t  = "+str(len(t))   +"\n")
f.write("delta_t= "+str(t[1]-t[0])+"\n")
f.write("name   = '"+title+"'\n")
f.close()

"""
