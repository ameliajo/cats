import sys
import numpy as np
sys.path.append('../newGasket/')
from testData import *
from main import *

# This shows us that the envelope of F(t) and H(t) closes at 2pi/db where db
# is the spacing in beta.

alphas = [1.0];betas = [1.0]
t = np.linspace(0,150,10001)
T = 0.025; invT = 1.0/T

rhoBetas1 = [x*invT for x in X1]
rhoBetas2 = [x*invT for x in X2]
rhoBetas3 = [x*invT for x in X3]
rhoBetas4 = [x*invT for x in X4]

plt.subplot(2,1,1); plt.plot(X1,Q1)
plt.subplot(2,1,2); plt.plot(X2,Q2)
#plt.subplot(4,1,3); plt.plot(X3,Q3)
#plt.subplot(4,1,4); plt.plot(X4,Q4)
plt.show()


plt.subplot(2,1,1)
sab,H,F = simpleGASKET(rhoBetas1,Q1,t,alphas,betas) 
plt.plot(t,F,label='F(t)'); plt.plot(t,H,label='H(t)')
plt.ylabel('time');         plt.legend(loc='best')
plt.plot([3.14159*2/(rhoBetas1[1]-rhoBetas1[0]),\
          3.14159*2/(rhoBetas1[1]-rhoBetas1[0])],[-max(F),max(F)])
plt.subplot(2,1,2)
sab,H,F = simpleGASKET(rhoBetas2,Q2,t,alphas,betas) 
plt.plot(t,F,label='F(t)'); plt.plot(t,H,label='H(t)')
plt.ylabel('time');         plt.legend(loc='best')
plt.plot([3.14159*2/(rhoBetas2[1]-rhoBetas2[0]),\
          3.14159*2/(rhoBetas2[1]-rhoBetas2[0])],[-max(F),max(F)])
plt.legend(loc='best')
#plt.subplot(3,1,3)
#sab,H,F = simpleGASKET(rhoBetas3,Q3,t,alphas,betas) 
#plt.plot(t,F,label='F(t)')
#plt.plot(t,H,label='H(t)')
#plt.legend(loc='best')
#plt.subplot(4,1,4)
#sab,H,F = simpleGASKET(rhoBetas4,Q4,t,alphas,betas) 
#plt.plot(t,F,label='F(t)')
#plt.plot(t,H,label='H(t)')
#plt.legend(loc='best')


plt.show()


