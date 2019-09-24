import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import numpy as np
import sys
sys.path.append("NJOY/")
sys.path.append("GASKET/")





def prepPlot(vector):
    cnorm = colors.Normalize(vmin=0,vmax=2*len(vector)+1)
    scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('tab20')) #hot autumn tab10
    mymap = colors.LinearSegmentedColormap.from_list('funTestColors',\
            [scalarMap.to_rgba(a) for a in range(len(vector))])
    colorBar = plt.contourf([[0,0],[0,0]], vector, cmap=mymap)
    plt.clf()
    return scalarMap, colorBar

def endPlot(colorBar):
    ax = plt.gca()
    plt.colorbar(colorBar).ax.set_ylabel('alpha values')
    plt.yscale('log')
    plt.show()

def checkThatAlphasBetasMatch(gasket_alphas,gasket_betas,njoy_alphas,njoy_betas):
    assert(len(njoy_alphas) == len(gasket_alphas))
    assert(len(njoy_betas) == len(gasket_betas))
    for a in range(len(njoy_alphas)): assert(abs(njoy_alphas[a]-gasket_alphas[a])<1e-5)
    for b in range(len(njoy_betas)):  assert(abs(njoy_betas[b]-gasket_betas[b])<1e-5)


def doThePlotting(alphas,betas,sab1,sab2):
    for a in range(len(alphas)):
        #njoy_sym = [njoy_sab[b+a*len(betas)] * \
        #            np.exp(-betas[b]*0.5) for b in range(len(betas))]
        #gasket_sym = [0.02550730568*gasket_sab[b+a*len(betas)] for b in range(len(betas))]
        plt.plot(betas,[sab1[b+a*len(betas)] for b in range(len(betas))],color=scalarMap.to_rgba(a),linewidth=3)
        plt.plot(betas,[sab1[b+a*len(betas)] for b in range(len(betas))],color=scalarMap.to_rgba(a),\
                 marker='o',markersize=3,linewidth=2,linestyle='dashed')
        #difference = [(gasket_sym[i] - njoy_sym[i])/(njoy_sym[i]) for i in range(len(njoy_sym))] 
        #plt.plot(betas[:-1],difference[:-1],color=scalarMap.to_rgba(a))



def getValuesFromInput(name):
    
    name = name.split('.') if '.' in name else \
           name.split('_') if '_' in name else name.split()

    if name[0] == 'g':
        fileName = 'gasket_output_'+name[1]
        alphas = importMod(fileName).alphas; betas  = importMod(fileName).betas
        temp   = importMod(fileName).temp;   sab    = importMod(fileName).sab
        sab    = [temp*val for val in sab]
        return temp,alphas,betas,sab

    elif name[0] == 'n':
        fileName = 'njoy_output_'+name[1]
        alphas = importMod(fileName).alphas; betas  = importMod(fileName).betas
        temp   = importMod(fileName).temp;   sab    = importMod(fileName).sab
        sab = [sab[b+a*len(betas)]*np.exp(-betas[b]*0.5) for b in range(len(betas)) for a in range(len(alphas))]
        return temp,alphas,betas,sab






if __name__=='__main__':
    from importlib import import_module as importMod
    if len(sys.argv) == 3:
        T1,alphas1,betas1,sab1 = getValuesFromInput(sys.argv[1])
        T2,alphas2,betas2,sab2 = getValuesFromInput(sys.argv[2])
        checkThatAlphasBetasMatch(alphas1,betas1,alphas2,betas2)
        alphas, betas = alphas1[:], betas1[:]
        fig = plt.figure(num=None,figsize=(17,12),dpi=80,facecolor='w',edgecolor='k') #fig = plt.figure()
        scalarMap, colorBar = prepPlot(alphas)
        doThePlotting(alphas,betas,sab1,sab2)
        endPlot(colorBar)


   

    


