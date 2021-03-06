import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

title = 'H in H2O'
 
X = [ 0.0020301,0.0037546,0.0043421, 0.0054803, 0.0071531, 0.0076952, \
      0.0087993,0.0093379,0.010435, 0.011548, 0.013213, 0.015464, 0.018276, \
      0.020525, 0.023901, 0.026159, 0.027850, 0.029541, 0.031229, 0.032357,   \
      0.034022, 0.035708, 0.038503, 0.039637, 0.041329, 0.042458, 0.044150,   \
      0.045291, 0.047007, 0.049860, 0.052143, 0.054446, 0.056727, 0.059582,   \
      0.061284, 0.062428, 0.064112, 0.065798, 0.067487, 0.068601, 0.069718,   \
      0.071398, 0.071946, 0.073058, 0.074730, 0.075836, 0.078060, 0.080831,   \
      0.083626, 0.086955, 0.089729, 0.093618, 0.096944, 0.10084, 0.10418,     \
      0.106960, 0.10970, 0.11308, 0.11640, 0.11912, 0.12140, 0.12366, 0.12584, \
      0.129285, 0.132592, 0.138713, 0.145400, 0.156194, 0.174679, 0.176866,   \
      0.181339, 0.186439, 0.187581, 0.190318, 0.19269, 0.19433, 0.19497,      \
      0.196080, 0.19729, 0.19779, 0.19894, 0.19963, 0.20079, 0.20192, 0.20254, \
      0.204380, 0.20430, 0.20481, 0.20544, 0.20603, 0.20764, 0.20823, 0.20933, \
      0.209840, 0.21035, 0.21149, 0.21190, 0.21305, 0.21468, 0.21513, 0.21735, \
      0.219437, 0.222847, 0.226708, 0.232313, 0.2424175, 0.2497561, 0.260300, \
      0.270453, 0.2979112, 0.3125884, 0.329365, 0.337281, 0.343313, 0.361861, \
      0.374212, 0.378169, 0.383293, 0.386092, 0.39005, 0.39231, 0.39409,      \
      0.394640, 0.39579, 0.39748, 0.39747, 0.39801, 0.39929, 0.40030, 0.40157, \
      0.403220, 0.40491, 0.40721, 0.40831, 0.41003, 0.41112, 0.41233, 0.41281, \
      0.414090, 0.41512, 0.41734, 0.41905, 0.41968, 0.42013, 0.42184, 0.42288, \
      0.424570, 0.42613, 0.42780, 0.42942, 0.43114, 0.43273, 0.43491, 0.43778, \
      0.439960, 0.44102, 0.44213, 0.44321, 0.44558, 0.448262, 0.449911,        \
      0.452137, 0.454425, 0.459415, 0.464554, 0.471793, 0.479019, 0.48588,    \
      0.494250, 0.496]

Q = [0.0023176, 0.0029208, 0.0033018, 0.0035400, 0.0034131, 0.0031432, 0.0028893, 0.0025719, 0.0022228, 0.0020958, 0.0019213, 0.0019214, 0.0020486, 0.0021915, 0.0023981, 0.0025569, 0.0026840, 0.0028111, 0.0028429, 0.0028747, 0.0028272, 0.0027638, 0.0027163, 0.0027957, 0.0029069, 0.0030498, 0.0032880, 0.0035421, 0.0039707, 0.0047010, 0.0053043, 0.0060346, 0.0066697, 0.0073682, 0.0077334, 0.0079398, 0.0080192, 0.0080828, 0.0081147, 0.0080671, 0.0079878, 0.0078768, 0.0077816, 0.0076388, 0.0074960, 0.0073056, 0.0070200, 0.0066551, 0.0063536, 0.0058935, 0.0054809, 0.0048938, 0.0044654, 0.0039736, 0.0035928, 0.0032596, 0.0029740, 0.0025298, 0.0020220, 0.0016412, 0.0013080, 0.0011176, 0.00091140, 0.00075285, 0.00064192, 0.00053113, 0.00046800, 0.00042094, 0.00043779, 0.00046966, 0.00053339, 0.00069239, 0.00077181, 0.00099418, 0.0012959, 0.0016293, 0.0019627, 0.0024548, 0.0028199, 0.0032009, 0.0041534, 0.0049153, 0.0058519, 0.0068202, 0.0074234, 0.0082648, 0.0083600, 0.0084711, 0.0085664, 0.0085823, 0.0084554, 0.0082173, 0.0078523, 0.0072174, 0.0065825, 0.0060429, 0.0054080, 0.0042335, 0.0030748, 0.0023923, 0.0015035, 0.00078937, 0.00042447, 0.00025008, 0.00015513, 0.000092175, 0.000092561, 0.00010900, 0.00010953, 0.000095112, 0.000095884, 0.00011265, 0.00014481, 0.00014513, 0.00020961, 0.00035312, 0.00043269, 0.00068693, 0.00095692, 0.0016555, 0.0022271, 0.0027669, 0.0029574, 0.0034019, 0.0039258, 0.0041797, 0.0043861, 0.0046719, 0.0051640, 0.0054657, 0.0060372, 0.0063071, 0.0067041, 0.0069581, 0.0071963, 0.0074580, 0.0077043, 0.0078631, 0.0080219, 0.0081172, 0.0081014, 0.0079745, 0.0078158, 0.0076413, 0.0071334, 0.0065938, 0.0061177, 0.0055463, 0.0049750, 0.0044672, 0.0040704, 0.0035943, 0.0032611, 0.0026898, 0.0023566, 0.0020392, 0.0017853, 0.0015631, 0.0011188, 0.00065862, 0.00048411, 0.00034137, 0.00024625, 0.00016715, 0.00015154, 0.00015193, 0.00016819, 0.00013680, 0.00010550, 0.0]


########
import numpy as np
import matplotlib.pyplot as plt
X_uniform = [0] + X
Q_uniform = [0] + Q

def interpolate(xVec,yVec,x):
    if x < 0 or x > xVec[-1]: return 0.0
    for i in range(len(xVec)-1):
        if x > xVec[i] and x < xVec[i+1]:
            b = yVec[i]
            m = (yVec[i+1]-yVec[i])/(xVec[i+1]-xVec[i])
            return m*(x-xVec[i])+b
    return 0.0

X_orig = X[:]
Q_orig = Q[:]
uniformX = np.linspace(0,X[-1],300)
uniformY = [interpolate([0.0]+X,[0.0]+Q,x) for x in uniformX]

uniformX = [float("{0:.7f}".format(round(a,7))) for a in uniformX]
uniformY = [float("{0:.5f}".format(round(a,5))) for a in uniformY]

X = uniformX[1:]
Q = uniformY[1:]
########


#print(uniformX)
#print(uniformY)
#exit()


#kbT = 0.0255

rho_f = interp1d([0.0]+X,[0.0]+Q,bounds_error=False,fill_value=0.0,kind='cubic')
#uniform_x = np.linspace(0,X[-1],3*len(X))
#uniform_y = rho_f(uniform_x)
#uniform_x = [x/kbT for x in uniform_x]






if __name__=="__main__":
    plot_just_damian = True
    plot_damian_beo = False

    from colors import colors
    plt.rcParams.update({'font.size': 14})

    if plot_damian_beo == True:
        from beoData import X as beoX
        from beoData import Q as beoQ
        watX = X[:]
        watQ = Q[:]

        invArea = 1.0/np.trapz(watQ,watX); watQ = [invArea * Qval for Qval in watQ]
        invArea = 1.0/np.trapz(beoQ,beoX); beoQ = [invArea * Qval for Qval in beoQ]

        watX = [1000*xval for xval in watX]
        beoX = [1000*xval for xval in beoX]

        plt.plot(beoX,beoQ,color='#f69b59',linewidth=1.5)
        plt.fill_between(beoX,beoQ,color=colors[1],alpha=1.0,label='Be in BeO')

        plt.plot(watX,watQ,color=colors[1],linewidth=1)
        plt.fill_between(watX,watQ,color=colors[3],alpha=0.4,label='H in H2O')

        plt.tick_params(axis='y',which='both',left=False,labelleft=False)
        plt.title('Normalized Phonon Distributions for BeO and H2O')
        plt.xlabel('Energy (meV)'); plt.ylabel('Normalized')
        plt.legend(loc='best');     plt.show()


    if plot_just_damian == True:
        #plt.plot(X,Q,color=colors[3],linewidth=1)
        #plt.fill_between(X,Q,color=colors[3],alpha=0.4,label='H in H2O')
        #plt.tick_params(axis='y',which='both',left=False,labelleft=False)
        #plt.title('Normalized Phonon Distribution for H in H2O')
        plt.xlabel('Energy (meV)'); plt.ylabel('Normalized')
        #plt.show()


    from math import sinh
    plt.rcParams.update({'font.size': 12})
    invTemp = 1.0/0.0255
    betas = [x*invTemp for x in X]
    plt.plot(betas,Q,color=colors[3],linewidth=1,label="Phonon DOS")
    T1 = [Q[b]/(2.0*betas[b]*sinh(betas[b]*0.5)) for b in range(len(betas))]
    plt.plot(betas,T1,color=colors[4],linewidth=1,label="T1")


    T2 = np.convolve(T1,T1)
    betas2 = [(betas[1]-betas[0])*i for i in range(len(T2))]
    plt.plot(betas2,T2,color=colors[0],linewidth=1,label="T2")



    plt.xlabel('beta (energy transfer)'); plt.ylabel('Normalized')
    plt.legend(loc='best')
    plt.yscale('log')
    plt.show()






