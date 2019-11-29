import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import numpy as np
from numpy import sin, cos, exp, sinh, cosh
import sys
sys.path.append("NJOY/")
sys.path.append("GASKET/")
sys.path.append("GASKET_correction/")
sys.path.append("../../../phononDistributions"); 

def prepPlot(vector):
    cnorm = colors.Normalize(vmin=0,vmax=len(vector)+2)
    scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('tab10')) 
    mymap = colors.LinearSegmentedColormap.from_list('funTestColors',\
            [scalarMap.to_rgba(a) for a in range(len(vector))])
    colorBar = plt.contourf([[0,0],[0,0]], vector, cmap=mymap)
    return scalarMap, colorBar

def endPlot(colorBar):
    ax = plt.gca()
    plt.colorbar(colorBar).ax.set_ylabel('alpha values')
    plt.show()

def checkThatAlphasBetasMatch(alphas2,betas2,alphas1,betas1):
    assert(len(alphas1) == len(alphas2) and len(betas1) == len(betas2))
    for a in range(len(alphas1)): assert(abs(alphas1[a]-alphas2[a])<1e-4)
    for b in range(len(betas1 )): assert(abs(betas1[b] -betas2[b]) <1e-4)

def doThePlottingReverse(alphas,betas,sab,linestyle,scalarMap):
    b = [betas[i] for i in range(0,len(betas)-2,2)]
    for a in range(len(alphas)):
        chunk = [sab[b+(len(alphas)-a-1)*len(betas)] for b in range(len(betas))]
        for b in range(0,len(chunk)-2,10):
            for bp in range(10):
                if chunk[b+bp] < 1e-8:
                    for bp2 in range(10): chunk[b+bp2] = 1e-8
                    break

        plt.plot(betas[:-1],chunk[:-1],color=scalarMap.to_rgba(len(alphas)-a-1),\
                 linewidth=2,linestyle=linestyle)
        plt.fill_between(betas[:-1],chunk[:-1],linewidth=2,linestyle=linestyle,\
                 color=scalarMap.to_rgba(len(alphas)-a-1),alpha=0.6)

def doThePlottingSingleAlpha(alphas,betas,sab,linestyle):
    for a in range(len(alphas)):
        chunk = [sab[b+a*len(betas)] for b in range(len(betas))]
        plt.plot(betas,chunk,linewidth=2,linestyle=linestyle)
    #for b in range(0,len(chunk)-2,10):
    #    for bp in range(10):
    #        if chunk[b+bp] < 0.0:
    #            for bp2 in range(10): chunk[b+bp2] = 1e-8
    #            break

    #plt.plot(betas[1:],chunk[1:],linewidth=2,linestyle=linestyle)



def doThePlotting(alphas,betas,sab,linestyle,scalarMap):
    b = [betas[i] for i in range(0,len(betas)-2,2)]
    for a in range(len(alphas)):
        chunk = [sab[b+a*len(betas)] for b in range(len(betas))]
        for b in range(0,len(chunk)-2,10):
            for bp in range(10):
                if chunk[b+bp] < 0.0:
                    for bp2 in range(10): chunk[b+bp2] = 1e-8
                    break

        plt.plot(betas[:-1],chunk[:-1],color=scalarMap.to_rgba(a),linewidth=2,\
                 linestyle=linestyle)

def plot_H_or_F(time,func,linestyle,color):
    plt.plot(time,func,color=color,linewidth=2,linestyle=linestyle)
    plt.fill_between(time, func, color=color,alpha =0.2)

def getValuesFromInput(name):
    name = name.split('.') if '.' in name else name.split() 
    if name[0] == 'g':
        fName = 'gasket_'+name[1]
        alphas = importMod(fName).alphas; betas   = importMod(fName).betas
        temp   = importMod(fName).temp;   sab     = importMod(fName).sab
        num_t  = importMod(fName).num_t;  delta_t = importMod(fName).delta_t 
        H      = importMod(fName).H    ;  F       = importMod(fName).F
        name   = importMod(fName).name
        #sab    = [temp*val for val in sab] 
        H, F   = [x*temp for x in H], [x*temp for x in F]
        return temp,alphas,betas,sab,num_t,delta_t,H,F,name,None,None

    elif name[0] == 'gc':
        fName = 'gasket_corrected_'+name[1]
        alphas = importMod(fName).alphas; betas   = importMod(fName).betas
        temp   = importMod(fName).temp;   sab     = importMod(fName).sab
        num_t  = importMod(fName).num_t;  delta_t = importMod(fName).delta_t 
        H      = importMod(fName).H    ;  F       = importMod(fName).F
        oscBegin        = importMod(fName).oscBegin 
        correctionBegin = importMod(fName).correctionBegin 
        name   = importMod(fName).name
        #H, F   = [x*temp for x in H], [x*temp for x in F]
        sab = [sab[b+a*len(betas)]*np.exp(betas[b]*0.5) for a in range(len(alphas)) for b in range(len(betas))]
        return temp,alphas,betas,sab,num_t,delta_t,H,F,name,oscBegin,correctionBegin

    elif name[0] == 'n':
        fName = 'njoy_output_'+name[1]
        alphas = importMod(fName).alphas; betas  = importMod(fName).betas
        temp   = importMod(fName).temp;   sab    = importMod(fName).sab
        #name   = importMod(fName).name
        name = ''
        sab = [sab[b+a*len(betas)]*np.exp(-betas[b]*0.5) for a in range(len(alphas)) for b in range(len(betas))]
        return temp,alphas,betas,sab,None,None,None,None,name,None,None












if __name__=='__main__':
    from importlib import import_module as importMod
    fig = plt.figure(num=None,figsize=(10,6),dpi=120,facecolor='w',edgecolor='k') 
    plt.rcParams.update({'font.size': 17})
    if len(sys.argv) == 1: print('Tell me what to plot!'); exit()

    option = 0

    if option == 0:
        cnorm = colors.Normalize(vmin=0,vmax=len(sys.argv))
        scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('tab10')) 

        for i,name in enumerate(sys.argv[1:]):
            T,alphas,betas,sab,nt,dt,H,F,title,oscBegin,correctionBegin = getValuesFromInput(name)
            label = 'gasket'           if name.split('.')[0] == 'g'  else \
                    'corrected gasket' if name.split('.')[0] == 'gc' else \
                    'njoy'             if name.split('.')[0] == 'n'  else '?'
            for a in range(len(alphas)):
                chunk = [sab[b+a*len(betas)] for b in range(len(betas))]
                plt.plot(betas,chunk,linewidth=1.5,label=label,color=scalarMap.to_rgba(i))
            if oscBegin != 'None' and oscBegin != None:
                plt.plot([float(oscBegin),float(oscBegin)],[min(chunk),max(chunk)],'k',linewidth=1.5)
                plt.plot([float(correctionBegin),float(correctionBegin)],[min(chunk),max(chunk)],'k',linewidth=1.5)

        plt.xlabel('beta'); plt.ylabel('Symmetric S(a,b)'); 
        plt.title('S(a,b) for '+title+' at '+str(int(T/8.617e-5))+'K for alpha='+str(alphas[0]))
        plt.legend(loc='best')
        plt.yscale('log'); 
        plt.show()




