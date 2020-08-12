import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from math import sin,cos

def f1(x): return sin(1.0*x)
def f2(x): return sin(1.02*x)+1
def f3(x): return sin(1.03*x)+1.3
def f4(x): return sin(1.08*x)+1.8

x1 = np.linspace(0,15,100)
y1 = [f1(x) for x in x1]
x2 = np.linspace(0.1,15.1,100)
y2 = [f2(x) for x in x2]
x3 = np.linspace(-0.1,14.9,100)
y3 = [f3(x) for x in x3]
y4 = [f4(x) for x in x3]

x = np.linspace(0.5,12.5,50)

temp = 15.0 
temperature = [10.0,20.0,30.0,40]

f1 = interp1d(x1,y1,bounds_error=False,fill_value=0.0)#'extrapolate')
f2 = interp1d(x2,y2,bounds_error=False,fill_value=0.0)#'extrapolate')
f3 = interp1d(x3,y3,bounds_error=False,fill_value=0.0)#'extrapolate')
f4 = interp1d(x3,y4,bounds_error=False,fill_value=0.0)#'extrapolate')
y1_new = f1(x)
y2_new = f2(x)
y3_new = f3(x)
y4_new = f4(x)
#y3_new[3]*=5

array = np.array([y1_new,y2_new,y3_new,y4_new])

plt.plot(x,array[0],'r-',label=str(temperature[0]))
plt.plot(x,array[1],'b-',label=str(temperature[1]))
plt.plot(x,array[2],'k-',label=str(temperature[2]))
plt.plot(x,array[3],'k-',label=str(temperature[3]))
#plt.legend(loc='best')
#plt.show()
#exit()

invArray = np.transpose(array) 

fullG_array = []
for i in range(len(invArray)):
    g1 = interp1d(temperature,invArray[i],bounds_error=True,kind=len(temperature)-1)
    fullG_array.append(g1(temp))
    #plt.plot(temperature,invArray[i],'bo')#,'r-')
    #plt.plot(temp,g1(temp),'ro')

plt.plot(x,fullG_array,'y-')
plt.show()




