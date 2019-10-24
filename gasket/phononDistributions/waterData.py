title = "H in H2O"
X = [0.0025*i for i in range(69)]
Q = [0.0, 0.00049, 0.00098, 0.00191, 0.00338, 0.00485, 0.00721, 0.00966, 0.01245, 0.01588, 0.01931, 0.02353, 0.02794, 0.0326, 0.03799, 0.04338, 0.04884, 0.05433, 0.05994, 0.06622, 0.07249, 0.07971, 0.08755, 0.09539, 0.10324, 0.11108, 0.1173, 0.1205, 0.12158, 0.12081, 0.11662, 0.11015, 0.10426, 0.09838, 0.09403, 0.08995, 0.08616, 0.08302, 0.07988, 0.07745, 0.07526, 0.07313, 0.07129, 0.06945, 0.06783, 0.06634, 0.06484, 0.06327, 0.06171, 0.06016, 0.05863, 0.057, 0.05588, 0.05467, 0.05356, 0.05258, 0.0516, 0.05055, 0.04949, 0.0484, 0.04726, 0.04613, 0.04502, 0.04392, 0.04299, 0.04256, 0.04213, 0.02471, 0.0]
continWgt = 0.444444
freeGasWgt = 0.055556
oscWgt = 1.0-continWgt-freeGasWgt
oscEnergies = [0.205   , 0.48    ]
oscWeights  = [0.166667, 0.333333]

X = X[1:]
Q = Q[1:]
if __name__=="__main__":
    colors = ["#F8B195", "#F67280", "#C06C84", "#6C5B7B", "#355C7D"]
    import matplotlib.pyplot as plt
    import numpy as np
    from beoData import X as beoX
    from beoData import Q as beoQ

    invArea = 1.0/np.trapz(Q,X)
    Q = [invArea * Qval for Qval in Q]
    invArea = 1.0/np.trapz(beoQ,beoX)
    beoQ = [invArea * Qval for Qval in beoQ]

    X = [1000*xval for xval in X]
    beoX = [1000*xval for xval in beoX]
    plt.rcParams.update({'font.size': 14})
    plt.plot(X,Q,color=colors[0],linewidth=2)
    plt.fill_between(X,Q,color=colors[0],alpha=0.5,label='H in H2O')
    plt.plot(beoX,beoQ,color=colors[1],linewidth=2)
    plt.fill_between(beoX,beoQ,color=colors[1],alpha=0.5,label='Be in BeO')
    plt.tick_params(axis='y',which='both',left=False,labelleft=False)
    plt.legend(loc='best')
    plt.title('Continuous Components of Phonon Distributions')
    plt.xlabel('Energy (meV)')
    plt.ylabel('Normalized')

    plt.show()

