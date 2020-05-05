import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import numpy as np
import sys
from ENDF_ZrH_ZrH2 import *
from EAM_Zr import *

options = [0,1,2]
""" 
0 - Plot the phonon distribution comparing DFT-generated ZrH
    with DFT-generated ZrH2. 
1 - Option 0, but with ENDF-generated ZrH (actually is ZrH2 though) alongside
2 - Plot EAM-generated Zr alongside DFT-generated ZrH and ZrH2
"""



def normalize(x,y,area=1.0):
    invArea = area/np.trapz(y,x)
    return [invArea*val for val in y]




def plotThisFile(fileName,label,color):
    THz_to_eV = 4.13567
    print("Reading data from "+fileName+"...")
    with open(fileName, "r") as f:
        E, dos   = [], []
        f.readline()
        for line in f:
            E.append(  float(line.split()[0])*THz_to_eV)
            dos.append(float(line.split()[1]))
        dos = normalize(E,dos)
        x = [e for e in E if e >= 0 ]
        y = [dos[i] for i in range(len(dos)) if E[i] >= 0 ]

        y2 = [y[i]*float(i)/30 if i < 30 else y[i]  for i in range(len(y)) ] 
        plt.plot(x,y2,linewidth=1.5,color=color)
        plt.fill(x,y2,label=label,linewidth=1.5,color=color,alpha=0.5)


def plot_ZrH_ZrH2():
    names = [ ['DFT_ZrH.dat' ,'ZrH (DFT)'],\
              ['DFT_ZrH2.dat' ,'ZrH2 (DFT)'],\
            ]
    colors = ['#F65058FF', '#FDD20EFF', '#28334AFF']
    for i,name in enumerate(names):
        plotThisFile(names[i][0],names[i][1],colors[i])
    plt.xlabel('Energy (meV)')
    plt.ylabel('Phonon Distribution (normalized)')


if 0 in options:
    plot_ZrH_ZrH2()
    plt.title('DFT-generated ZrHx Phonon Density of States')
    plt.legend(loc='best')
    plt.show()

if 1 in options:
    plot_ZrH_ZrH2()
    plt.title('DFT-generated ZrHx Phonon Density of States')
    endf_x,endf_y = get_ENDF_ZrH_dos()
    plt.plot(endf_x,endf_y,color='#28334AFF')
    plt.fill(endf_x,endf_y,color='#28334AFF',label='ZrH (ENDF)',alpha=0.3)
    plt.legend(loc='best')
    plt.show()
    
if 2 in options:
    plot_ZrH_ZrH2()
    plt.title('DFT-generated ZrHx Phonon Density of States')

    x = [x1*1000 for x1 in x_EAM if x1 >= 0]
    y = [y_EAM[i] for i in range(len(y_EAM)) if x_EAM[i] >= 0]
    y = normalize(x,y,0.4)

    plt.plot(x,y,color='#28334AFF')
    plt.fill(x,y,label='Zr (EAM)',color='#28334AFF',alpha=0.3)
    plt.xlabel('Energy (meV)')
    plt.ylabel('Phonon Density of States (arbitrary)')
    plt.title('EAM Zr and DFT-ZrHx Phonon Distributions')
    plt.legend(loc='best')
    plt.show()



