import matplotlib.pyplot as plt
import numpy as np
from numpy import sin,cos,exp

def func(x):
    return 5.0*x*exp(-x)

def closeEnough(xL,xR,tol=1e-5):
    xM = (xL+xR)*0.5
    yL = func(xL)
    yR = func(xR)
    guess = (yL+yR)*0.5
    real  = func((xL+xR)*0.5)
    return abs(real-guess) < 1e-4

spacing = 4.0
xL = 0.0
xR = spacing
final_X = [xL]
final_Y = [func(xL)]

yL = func(xL)
yR = func(xR)

while xR < 25:
    xM = (xL+xR)*0.5
    yM = func(xM)
    guess = (yL+yR)*0.5
    if abs(guess-yM) < 1e-3:
        final_X.append(xR)
        final_Y.append(yR)
        xL = xR
        yL = yR
        xR = xL + spacing
        yR = func(xR)
    else:
        xR = xM
        yR = yM


plt.plot(final_X,final_Y,'ro',label='Function on adaptive grid')
plt.plot(final_X,[0.0]*len(final_X),'go',label='Timestep locations')
plt.legend(loc='best')
plt.show()



