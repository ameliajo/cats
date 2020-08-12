import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


def getInterpolatedDist(temp):
    invArray = np.transpose(array) 
    fullG_array = []
    for i in range(len(invArray)):
        g1 = interp1d(temperature,invArray[i],bounds_error=True,kind=len(temperature)-1)
        fullG_array.append(g1(temp))
    return fullG_array



if __name__=="__main__":
    import matplotlib.pyplot as plt
    from h2oDists import *

    temp = 410.0
    array = np.array(yVec)
    fullG_array = getInterpolatedDist(temp)

    for i,vec in enumerate(array):
        plt.plot(x,vec,label=str(temperature[i]))


    plt.plot(x,fullG_array,'ko',label='Interpolated at '+str(temp))
    plt.legend(loc='best')
    plt.show()




