import matplotlib.pyplot as plt
import numpy as np
from scipy import special 
from math import pi
from testRho import *

       
alphas = list(np.linspace(1e-3,5,5))
betas  = list(np.linspace(0.001,5,40))

N = 10
weights = special.roots_laguerre(N, mu=False)[1]
zeros   = special.roots_laguerre(N, mu=False)[0]

alpha = 0.2
beta = 0.5



def get_g(alpha,t):
    N = 10
    weights = special.roots_laguerre(N, mu=False)[1]
    zeros   = special.roots_laguerre(N, mu=False)[0]
    g = 0.0
    for i in range(N):
        xi = zeros[i]
        rhoVal = getRho(2.0*xi)
        G1 = rhoVal*np.cos(2*xi*t)/(2*xi*np.sinh(xi))*np.exp(2*xi)
        G2 = rhoVal*np.cos(2*xi*t)/(2*xi*np.sinh(xi))
        g += weights[i]*(G1+G2)
    return g*alpha

def get_g_2(alpha,t):
    fullBetas = [-b for b in betas[::-1]] + betas[1:]
    P = [getRho(abs(beta))/(2.0*beta*np.sinh(beta*0.5)) for beta in fullBetas]
    argument = [P[b]*np.exp(-fullBetas[b]*0.5)*np.cos(fullBetas[b]*t) \
                for b in range(len(fullBetas))]
    return alpha * np.trapz(argument,x=fullBetas)

print(get_g(2.0,5.0))
print(get_g_2(2.0,5.0))
    

"""
def get_f(alpha,beta,t):
    N = 10
    weights = special.roots_laguerre(N, mu=False)[1]
    zeros   = special.roots_laguerre(N, mu=False)[0]
    f = 0.0
    for i in range(N):
        xi = zeros[i]
        rhoVal = getRho(energies,rho,2.0*xi)
        F1 = rhoVal*np.sin(2*xi*t)/(2*xi*np.sinh(xi))*np.exp(2*xi)
        F2 = rhoVal*np.sin(2*xi*t)/(2*xi*np.sinh(xi))
        f += weights[i]*(-F1+F2)
    return beta*t - f*alpha


def get_lambda():
    N = 10
    weights = special.roots_laguerre(N, mu=False)[1]
    zeros   = special.roots_laguerre(N, mu=False)[0]
    lambda_s = 0.0
    for i in range(N):
        bi = zeros[i]
        rhoVal = getRho(energies,rho,2.0*bi)
        F1 = rhoVal/(2*bi*np.sinh(bi))*np.exp(-bi*0.5)
        lambda_s += weights[i]*F1
    return lambda_s 

def getSAB(alpha,beta):
    N = 10
    weights = special.roots_laguerre(N, mu=False)[1]
    zeros   = special.roots_laguerre(N, mu=False)[0]
    sabVal = 0.0

    for i in range(N):
        ti = zeros[i]
        sabVal += np.exp(get_g(alpha,ti)*ti)   * \
                  np.cos(get_f(alpha,beta,ti)) * weights[i]
    sabVal *= np.exp(-get_lambda()*alpha)/pi
    return sabVal


for alpha in alphas:
    print(alpha)
    sabVals = []
    for beta in betas:
        sabVals.append(getSAB(alpha,beta))

    plt.plot(betas,sabVals)
plt.show()



"""

    


