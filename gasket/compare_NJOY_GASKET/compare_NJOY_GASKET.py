import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import numpy as np
import sys
sys.path.append("NJOY/")
sys.path.append("GASKET/")





def prepPlot(vector):
    cnorm = colors.Normalize(vmin=0,vmax=len(vector)+1)
    scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('tab20')) #hot autumn tab10
    mymap = colors.LinearSegmentedColormap.from_list('funTestColors',\
            [scalarMap.to_rgba(a) for a in range(len(vector))])
    colorBar = plt.contourf([[0,0],[0,0]], vector, cmap=mymap)
    return scalarMap, colorBar

def endPlot(colorBar):
    ax = plt.gca()
    plt.colorbar(colorBar).ax.set_ylabel('alpha values')
    #plt.yscale('log')
    plt.show()

def checkThatAlphasBetasMatch(gasket_alphas,gasket_betas,njoy_alphas,njoy_betas):
    assert(len(njoy_alphas) == len(gasket_alphas))
    assert(len(njoy_betas) == len(gasket_betas))
    for a in range(len(njoy_alphas)): assert(abs(njoy_alphas[a]-gasket_alphas[a])<1e-5)
    for b in range(len(njoy_betas)):  assert(abs(njoy_betas[b]-gasket_betas[b])<1e-5)


def doThePlotting(alphas,betas,sab,linestyle):
    for a in range(1,len(alphas)):
        plt.plot(betas[:-2],[sab[b+a*len(betas)] for b in range(len(betas)-2)],\
                 color=scalarMap.to_rgba(a),linewidth=2,linestyle=linestyle)

def doThePlottingSingleAlpha(betas,sab,linestyle,color):
    plt.plot(betas[:-2],[sab[b] for b in range(len(betas)-2)],\
             color=color,linewidth=2,linestyle=linestyle)

def plot_H_or_F(time,func,linestyle,color):
    plt.plot(time,func,color=color,linewidth=2,linestyle=linestyle)
    plt.fill_between(time, func, color=color,alpha =0.2)


def getValuesFromInput(name):
    
    name = name.split('.') if '.' in name else \
           name.split('_') if '_' in name else name.split()

    if name[0] == 'g':
        fName = 'gasket_output_'+name[1]
        alphas = importMod(fName).alphas; betas   = importMod(fName).betas
        temp   = importMod(fName).temp;   sab     = importMod(fName).sab
        num_t  = importMod(fName).num_t;  delta_t = importMod(fName).delta_t 
        H      = importMod(fName).H    ;  F       = importMod(fName).F
        sab    = [temp*val for val in sab]
        H, F   = [x*temp for x in H], [x*temp for x in F]
        return temp,alphas,betas,sab,num_t,delta_t,H,F

    elif name[0] == 'n':
        fName = 'njoy_output_'+name[1]
        alphas = importMod(fName).alphas; betas  = importMod(fName).betas
        temp   = importMod(fName).temp;   sab    = importMod(fName).sab
        sab = [sab[b+a*len(betas)]*np.exp(-betas[b]*0.5) for b in range(len(betas)) for a in range(len(alphas))]
        return temp,alphas,betas,sab,None,None,None,None






if __name__=='__main__':
    from importlib import import_module as importMod
    fig = plt.figure(num=None,figsize=(10,6),dpi=120,facecolor='w',edgecolor='k') 
    if len(sys.argv) == 1: exit()

    option = 0 # Here, I'm just plotting a simple F and H across time 
    option = 1 # Here, we look at how time step things affect the H and F functions
               # from pg 37 in the GASKET manual
    option = 2 # Here, we look at a few different sab datasets and see how the 
               # timestep sizes affect them. These should all be GASKET inputs
    option = 3 # Here, we plot one or more S(a,b) against beta for various alphas

    if option == 0:
        T,alphas,betas,sab,nt,dt,H,F = getValuesFromInput(sys.argv[1])
        time = [i*dt for i in range(nt)]
        scalarMap, colorBar = prepPlot(list(range(len(sys.argv)))); plt.clf(); 
        plt.fill_between(time[1:],H[1:],linewidth=2,color=scalarMap.to_rgba(0),label='H(t)',alpha=0.5)
        plt.fill_between(time,F,linewidth=2,color=scalarMap.to_rgba(1),label='F(t)',alpha=0.5)
        plt.xlabel('t'); plt.ylabel('Normalized'); plt.legend(loc='best')
        plt.title('H(t) and F(t) for water at '+str(int(T/8.617e-5))+'K')
        plt.show()




    if option == 1:
        tSpacings = []; maxVal = 0.0
        scalarMap, colorBar = prepPlot(list(range(len(sys.argv)))); plt.clf();
        for i,name in enumerate(sys.argv[1:]):
            T,alphas,betas,sab,nt,dt,H,F = getValuesFromInput(name)
            if max(max(H),max(F)) > maxVal: maxVal = max(max(H),max(F))
            time = [i*dt for i in range(nt)]
            tSpacings.append(dt)
            plot_H_or_F(time,H,'solid',scalarMap.to_rgba(i))
            plot_H_or_F(time,F,'dashed',scalarMap.to_rgba(i))
        scalarMap, colorBar = prepPlot(tSpacings+[tSpacings[-1]+1])
        plt.ylim([-0.1*maxVal,1.1*maxVal])
        endPlot(colorBar)




    if option == 2:
        tSpacings = []
        scalarMap, colorBar = prepPlot(list(range(len(sys.argv)))); plt.clf()
        for i,name in enumerate(sys.argv[1:]):
            T,alphas,betas,sab,nt,dt,H,F = getValuesFromInput(name)
            tSpacings.append(dt)
            assert(len(betas)==len(sab))
            if i > 0: checkThatAlphasBetasMatch(alphasOld,betasOld,alphas,betas)
            alphasOld, betasOld = alphas[:], betas[:]
            plt.plot(betas[:-1],sab[:-1],linewidth=1.5,color=scalarMap.to_rgba(i),label='dt = '+str(dt))
        plt.xlabel('beta'); plt.ylabel('S(a,b)'); plt.legend(loc='best')
        plt.title('S(a,b) for water at '+str(int(T/8.617e-5))+'K for various dt values. alpha = '+str(alphas[0]))
        plt.show()

    if option == 3:
        for i,name in enumerate(sys.argv[1:]):
            T,alphas,betas,sab,nt,dt,H,F = getValuesFromInput(name)
            if i == 0:
                scalarMap, colorBar = prepPlot(alphas); plt.clf()
            else:
                checkThatAlphasBetasMatch(alphasOld,betasOld,alphas,betas)
            alphasOld, betasOld = alphas[:], betas[:]
            doThePlotting(alphas,betas,sab,'solid')
            #plt.plot(betas[:-1],sab[:-1],linewidth=1.5,color=scalarMap.to_rgba(i),label='dt = '+str(dt))
        plt.xlabel('beta'); plt.ylabel('S(a,b)'); 
        plt.title('S(a,b) for water at '+str(int(T/8.617e-5))+'K as prepared by GASKET (t = 0->600, dt = 0.1)')
        plt.yscale('log')
        endPlot(colorBar)

