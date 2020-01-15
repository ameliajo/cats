import matplotlib.pyplot as plt

colors = ["#99B898", "#FECEAB", "#FF847C", "#E84A5F", "#2A363B"]

kb = 8.617333e-5


temps = [ 296.0, 400.0, 500.0, 600.0, 700.0, 800.0, 1000.0, 1200.0 ]
dwf_h_in_zrh  = [ 8.4795, 9.0854, 9.8196, 10.676, 11.625, 12.643, 14.822, 17.125 ]
dwf_zr_in_zrh = [ 1.9957, 2.6546, 3.2946, 3.9380, 4.5835, 5.2302, 6.5260, 7.8236 ]

c11a=162.88; c11b=296.;  c11c=34.957; c11d=350.; 
c11e=40.282; c12a=81.44; c13a=6.3366; 

njoy_h_in_zrh  = [0.216300, 0.313176, 0.423111, 0.552013, 0.701266, 0.871643, 1.277359, 1.770990]
njoy_zr_in_zrh = [4.572379, 8.347289, 12.950096, 18.575446, 25.223441, 32.894127, 51.303645, 73.804083]
temp_ch2 = [77.0, 196.0, 233.0, 293.6, 300.0, 303.0, 313.0, 323.0, 333.0, 343.0, 350.0]
njoy_h_in_ch2  = [0.082896, 0.343301, 0.460748, 0.691311, 0.718451, 0.731358, 0.775230, 0.820414, 0.866912, 0.914725, 0.948979]

h_in_zrh_A = 0.99917
zr_in_zrh_A = 90.436
h_in_ch2_A = 0.9991673


for t in range(len(temps)):
    njoy_h_in_zrh[t] /= (h_in_zrh_A*kb*temps[t])
    njoy_zr_in_zrh[t] /= (zr_in_zrh_A*kb*temps[t])

for t in range(len(temp_ch2)):
    njoy_h_in_ch2[t] /= (h_in_ch2_A*kb*temp_ch2[t])

dwf_polyethylene = [c11c+(temp-c11b)*(c11e-c11c)/(c11d-c11b) for temp in temp_ch2]
plt.plot(temps,dwf_h_in_zrh,label='hard coded H(ZrH)',color=colors[0],linestyle='dashed',linewidth=3)
plt.plot(temps,njoy_h_in_zrh,label='njoy H(ZrH)',color=colors[0],linewidth=3)
plt.plot(temps,dwf_zr_in_zrh,label='hard coded Zr(ZrH)',color=colors[1],linestyle='dashed',linewidth=3)
plt.plot(temps,njoy_zr_in_zrh,label='njoy Zr(ZrH)',color=colors[1],linewidth=3)
plt.plot(temp_ch2,dwf_polyethylene,label='hard coded polyethylene',color=colors[2],linestyle='dashed',linewidth=3)
plt.plot(temp_ch2,njoy_h_in_ch2,label='njoy polyethylene',color=colors[2],linewidth=3)

plt.legend(loc='best')
plt.show()




