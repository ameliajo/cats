from phononDists import *
import numpy as np
import matplotlib.pyplot as plt
from math import pi
import matplotlib.colors as colors
import matplotlib.cm as cmx



def prepPlot(vec):
    cnorm = colors.Normalize(vmin=0,vmax=len(vec)+1)
    scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('tab10')) #hot autumn tab10
    #scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('tab10')) #hot autumn tab10
    mymap = colors.LinearSegmentedColormap.from_list('funTestColors',\
            [scalarMap.to_rgba(a) for a in range(len(vec))])
    colorBar = plt.contourf([[0,0],[0,0]], vec, cmap=mymap)
    plt.clf()
    return scalarMap, colorBar



def finishPlotting(colorBar,title):
    ax = plt.gca()
    plt.colorbar(colorBar).ax.set_ylabel('time t')
    #plt.title(title)
    #ax.set_facecolor('xkcd:light grey blue') # off white
    #plt.yscale('log')
    plt.show()





P = [rho[1]/(betas[1]*betas[1])] + \
    [rho[i]/(2.0*betas[i]*np.sinh(betas[i]*0.5)) for i in range(1,len(rho))]
fullP = P[::-1]+P[1:]

fullBetas = [-x for x in betas][::-1] + betas[1:]

def getF(t,fullBetas,fullP):
    F_integrand = [fullP[i]*np.exp(-fullBetas[i]*0.5)*np.cos(fullBetas[i]*t) \
                   for i in range(len(fullP))]
    return np.trapz(F_integrand,x=fullBetas)

def getG(t,fullBetas,fullP):
    G_integrand = [fullP[i]*np.exp(-fullBetas[i]*0.5)*np.sin(fullBetas[i]*t) \
                   for i in range(len(fullP))]
    return np.trapz(G_integrand,x=fullBetas)

#plt.plot(fullBetas,rho[::-1]+rho[1:],label='rho(b)')
#plt.plot(fullBetas,fullP,label='P(b)')
#plt.plot(fullBetas,[fullP[i]*np.exp(-fullBetas[i]*0.5) for i in range(len(fullBetas))],label='P(b)*exp(-b/2)')
#times = [0.0,1.0,2.0,3.0,10.0,30.0,60.0,61.0,62.85]
times = [0.0,1.0,2.0,3.0,10.0]
scalarMap,colorBar = prepPlot(times)
for i,t in enumerate(times):
    plt.plot(fullBetas,[np.cos(fullBetas[i]*t)*fullP[i]*np.exp(-fullBetas[i]*0.5) for i in range(len(fullBetas))],color=scalarMap.to_rgba(i))
    plt.fill(fullBetas,[np.cos(fullBetas[i]*t)*fullP[i]*np.exp(-fullBetas[i]*0.5) for i in range(len(fullBetas))],color=scalarMap.to_rgba(i),alpha=0.2)

plt.xlabel('beta')
finishPlotting(colorBar,'')

exit()

alpha = 1.01
beta = 0.006
beta = 0.00255*3/(kb*296)
print(beta,2*pi/beta)

F_contrib, G_contrib, sabContrib = [], [], []
partial1, partial2, partial3, integral = [], [], [], []
#t_vec = np.linspace(0,50000,10000)
#t_vec = np.linspace(0,8*2*pi/betas[1],1e4)
#t_vec = np.linspace(0,40*2*pi/betas[1],1e5)
t_vec = np.linspace(0,20*2*pi/betas[1],1e5)
for i,t in enumerate(t_vec):
       
    F = getF(t,fullBetas,fullP)
    G = getG(t,fullBetas,fullP)
    F_contrib.append(F)
    G_contrib.append(G)
    sabContrib.append( np.exp(alpha*F) * np.cos(beta*t - alpha*G) )
    partial1.append( np.cos(beta*t - alpha*G) )
    partial2.append( np.exp(alpha*F) )
    integral.append(np.trapz(sabContrib,x=t_vec[:i+1]))
    #partial3.append( np.exp(alpha*F) * np.cos(beta*t) )
    if (int(i)%1000==0):
        print((i*100/len(t_vec)))
 

plt.plot(t_vec,F_contrib,label='F')
#plt.plot(t_vec,G_contrib,label='G')
#plt.plot(t_vec,partial2,label='e^(a*F)')
#plt.plot(t_vec,partial1,label='cos(b*t-a*G)')
#plt.plot(t_vec,sabContrib,label='e^(a*F) * cos(b*t-a*G)')
#plt.plot(t_vec,integral,label='int e^(a*F) * cos(b*t-a*G)')
plt.legend(loc='best')
plt.show()
"""
numSpaces = int(1e5/2)
t_vec = np.linspace(0,4*2*pi/betas[1],2*numSpaces)
for i,t in enumerate(t_vec):
    if (int(i)%1000==0):
        print(i)
    F_contrib.append(getF(t,fullBetas,fullP))
    G_contrib.append(getG(t,fullBetas,fullP))
plt.plot(t_vec[:numSpaces],F_contrib[numSpaces:2*numSpaces])
plt.plot(t_vec[1:numSpaces+1],F_contrib[0:numSpaces])
plt.plot(t_vec[:numSpaces],G_contrib[0:numSpaces])
plt.plot(t_vec[1:numSpaces+1],G_contrib[numSpaces:2*numSpaces])
plt.show()
"""


