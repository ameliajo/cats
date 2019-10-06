import matplotlib.pyplot as plt
import numpy as np
from numpy import sin,cos,exp


def func(x):
    return 5.0*x*exp(-x)

spacing = 4.0

def closeEnough(xL,xR,tol=1e-5):
    xM = (xL+xR)*0.5
    yL = func(xL)
    yR = func(xR)
    guess = (yL+yR)*0.5
    real  = func((xL+xR)*0.5)
    return abs(real-guess) < 1e-4
    #return abs((real-guess)/real) < tol

xL = 0.0;        yL = func(xL)
xR = xL+spacing; yR = func(xR)

final_X = [xL]
final_Y = [yL]

while xR < 25:
    yL = func(xL)
    yR = func(xR)
    if not closeEnough(xL,xR):
        xR = (xL+xR)*0.5
        yR = func(xR)
    else:
        final_X.append(xR)
        final_Y.append(yR)
        xL = xR
        xR = xL + spacing

plt.plot(final_X,final_Y,'ro')
plt.plot(final_X,[0.0]*len(final_X),'go')
plt.show()





#plt.plot(x,y)
#plt.plot(x,y,'ro')
#plt.plot(halfway_X,halfway_guess,'go')
#plt.plot(halfway_X,halfway_real,'bo')
#plt.show()





