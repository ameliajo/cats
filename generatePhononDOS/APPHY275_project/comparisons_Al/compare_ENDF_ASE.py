import matplotlib.pyplot as plt
import numpy as np
from Al_PhononDists_EAM_LJ import *
from Al_PhononDists_ENDF import *



if __name__=="__main__":
    plt.title('Comparing FCC Aluminum Phonon Distribution (ENDF vs EAM vs. LJ)')
    invArea = 1.0/np.trapz(ENDF_y,ENDF_x)
    ENDF_y = [y*invArea for y in ENDF_y]
    
    i = 0
    for j in range(len(Al_PhononDOS_x)):
        if Al_PhononDOS_x[j] >= 0:
            i = j
            break

    plt.plot(Al_PhononDOS_x[j:],Al_PhononDOS_EAM[j:],label='EAM')
    plt.plot(ENDF_x,ENDF_y,label='ENDF')
    plt.plot(Al_PhononDOS_x,Al_PhononDOS_LJ,label='LJ')
    plt.xlabel('Energy (eV)')
    plt.ylabel('Phonon DOS (unitless)')
    plt.legend(loc='best')
    plt.show()
 






