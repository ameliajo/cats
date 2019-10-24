import matplotlib.pyplot as plt
import numpy as np

xVec = [6.375E-3, 1.275E-2, 1.9125E-2, 2.55E-2, 3.1875E-2, 3.825E-2, 4.4625E-2, \
     5.1E-2, 5.7375E-2, 6.375E-2, 6.63E-2, 6.885E-2, 7.14E-2, 7.395E-2, 7.65E-2,\
     8.2875E-2, 8.925E-2, 9.5625E-2, 1.02E-1, 1.08375E-1, 1.1475E-1, 1.21125E-1,\
     1.275E-1, 1.33875E-1, 1.4025E-1, 1.46625E-1, 1.53E-1, 1.59375E-1, 1.6575E-1,\
     1.7e-1]

fVec = [1.25E-3, 5.0E-3, 1.125E-2, 2.0E-2, 3.125E-2, 4.5E-2, 5.9E-2, 7.5E-2, \
     9.5E-2, 1.15E-1, 1.197E-1, 1.214E-1, 1.218E-1, 1.195E-1, 1.125E-1, 9.75E-2,\
     8.71E-2, 7.91E-2, 7.35E-2, 6.88E-2, 6.5E-2, 6.1E-2, 5.71E-2, 5.4E-2,\
     5.15E-2, 4.88E-2, 4.59E-2, 4.31E-2, 4.2E-2, 0.0] 

xVec = [0] + xVec
fVec = [0] + fVec


def interpolate(xVec,yVec,x):
    if x < xVec[0] or x > xVec[-1]:
        return 0.0
    for i in range(len(xVec)-1):
        if x > xVec[i] and x < xVec[i+1]:
            b = yVec[i]
            m = (yVec[i+1]-yVec[i])/(xVec[i+1]-xVec[i])
            return m*(x-xVec[i])+b
    return 0.0

uniformX = np.linspace(0,0.17,17*4+1)
uniformY = [interpolate(xVec,fVec,x) for x in uniformX]

uniformX = [float("{0:.7f}".format(round(a,7))) for a in uniformX]
uniformY = [float("{0:.5f}".format(round(a,5))) for a in uniformY]
#print(len(uniformX))
#print(" "); print(uniformX)
#print(" "); print(uniformY); print(" ")

fig = plt.figure()
plt.plot(xVec,fVec)
plt.plot(uniformX,uniformY)

#print(""); print(uniformX[1:-1])
#print(""); print(uniformY[1:-1]); print("")

X_to_Print = uniformX[1:-1]
Y_to_Print = uniformY[1:-1]
#for i in range(len(X_to_Print)):
#    print("X("+str(i+1)+") = "+str(X_to_Print[i]))
#    print("Q("+str(i+1)+") = "+str(Y_to_Print[i]))

#plt.draw();plt.pause(1); raw_input("<Hit Enter To Close>"); plt.close(fig)



