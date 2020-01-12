import matplotlib.pyplot as plt

dwf_graphite = [2.1997,  2.7448,  3.2912,  3.8510,  4.4210,  4.9969, 6.1624,  7.3387,  9.6287, 11.992]
dwf_Be       = [3.16663, 3.88842, 4.62944, 5.40517, 6.19880, 7.0042, 8.63665, 10.2865, 0.0,    0.0   ]
dwf_BeO      = [2.153,   2.6374,  3.1348,  3.6513,  4.1798,  4.7164, 5.8052,  6.9068,  0.0,    0.0   ]

temps = [ 296.0, 400.0, 500.0, 600.0, 700.0, 800.0, 1000.0, 1200.0, 1600.0, 2000.0 ]

njoy_crystalline_graphite = [0.868056, 1.490619, 2.256850, 3.189282, 4.288908, 5.556286, 8.595482, 12.308253, 21.756507, 33.902898]
njoy_reactor_graphite_10p = [2.242738, 3.998034, 6.173456, 8.828603, 11.964347, 15.581185, 24.259217, 34.863973, 61.855523, 96.557605]
njoy_reactor_graphite_30p = [2.417081, 4.307196, 6.651144, 9.512948, 12.893411, 16.792958, 26.150167, 37.585555, 66.692219, 104.114189]
njoy_Be_metal = [0.603339, 0.982390, 1.449021, 2.017434, 2.688255, 3.461769, 5.317377, 7.584794, 13.355706, 20.775080]
njoy_Be_in_BeO = [0.451394, 0.726037, 1.054004, 1.451808, 1.920369, 2.460133, 3.754105, 5.334634, 9.356524, 14.526837]
crystalline_graphite_A = 11.898
reactor_graphite_10p_A = 11.898
reactor_graphite_30p_A = 11.898
Be_metal_A = 8.93478
Be_in_BeO_A = 8.93478

kb = 8.617333e-5

for t in range(len(temps)):
    njoy_crystalline_graphite[t] /= (crystalline_graphite_A*kb*temps[t])
    njoy_reactor_graphite_10p[t] /= (reactor_graphite_10p_A*kb*temps[t])
    njoy_reactor_graphite_30p[t] /= (reactor_graphite_30p_A*kb*temps[t])
    njoy_Be_metal[t] /= (Be_metal_A*kb*temps[t])
    njoy_Be_in_BeO[t] /= (Be_in_BeO_A*kb*temps[t])


#plt.plot(temps,dwf_graphite,label='hard coded')
#plt.plot(temps,njoy_crystalline_graphite,label='crystal')
#plt.plot(temps,njoy_reactor_graphite_10p,label='10p')
#plt.plot(temps,njoy_reactor_graphite_30p,label='30p')

#plt.plot(temps,dwf_Be,label='hard coded Be')
#plt.plot(temps,njoy_Be_metal,label='njoy generated Be')
plt.plot(temps,dwf_BeO,label='hard coded Be in BeO')
plt.plot(temps,njoy_Be_in_BeO,label='njoy Be in BeO')

plt.legend(loc='best')
plt.show()





