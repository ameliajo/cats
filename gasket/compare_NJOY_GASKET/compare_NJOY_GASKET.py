import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import numpy as np
from numpy import sin, cos, exp, sinh, cosh
import sys
sys.path.append("NJOY/")
sys.path.append("GASKET/")
sys.path.append("SIMPLE_GASKET/")
sys.path.append("QUAD_GASKET/")

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
        #for b in range(0,len(chunk)-2,10):
        #    for bp in range(10):
        #        if chunk[b+bp] < 0.0:
        #            for bp2 in range(10): chunk[b+bp2] = 1e-8
        #            break

        plt.plot(betas[1:],chunk[1:],linewidth=2,\
                 linestyle=linestyle)



def doThePlotting(alphas,betas,sab,linestyle,scalarMap):
    b = [betas[i] for i in range(0,len(betas)-2,2)]
    for a in range(len(alphas)):
        chunk = [sab[b+a*len(betas)] for b in range(len(betas))]
        for b in range(0,len(chunk)-2,10):
            for bp in range(10):
                if chunk[b+bp] < 0.0:
                    for bp2 in range(10): chunk[b+bp2] = 1e-8
                    break

        plt.plot(betas[:-1],chunk[:-1],color=scalarMap.to_rgba(a),linewidth=0.5,\
                 linestyle=linestyle,markersize=5,marker='o')

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
        sab    = [temp*val for val in sab] 
        H, F   = [x*temp for x in H], [x*temp for x in F]
        return temp,alphas,betas,sab,num_t,delta_t,H,F,name

    elif name[0] == 's':
        fName = 'simple_'+name[1]
        alphas = importMod(fName).alphas; betas   = importMod(fName).betas
        temp   = importMod(fName).temp;   sab     = importMod(fName).sab
        num_t  = importMod(fName).num_t;  delta_t = importMod(fName).delta_t 
        H      = importMod(fName).H    ;  F       = importMod(fName).F
        #         color=scalarMap.to_rgba(a),linewidth=2,linestyle=linestyle)
        name   = importMod(fName).name
        H, F   = [x*temp for x in H], [x*temp for x in F]
        return temp,alphas,betas,sab,num_t,delta_t,H,F,name

    elif name[0] == 'q':
        fName = 'quad_'+name[1]
        alphas = importMod(fName).alphas; betas   = importMod(fName).betas
        temp   = importMod(fName).temp;   sab     = importMod(fName).sab
        num_t  = importMod(fName).num_t;  delta_t = importMod(fName).delta_t 
        H      = importMod(fName).H    ;  F       = importMod(fName).F
        #         color=scalarMap.to_rgba(a),linewidth=2,linestyle=linestyle)
        name   = importMod(fName).name
        H, F   = [x*temp for x in H], [x*temp for x in F]
        return temp,alphas,betas,sab,num_t,delta_t,H,F,name



    elif name[0] == 'n':
        fName = 'njoy_output_'+name[1]
        alphas = importMod(fName).alphas; betas  = importMod(fName).betas
        temp   = importMod(fName).temp;   sab    = importMod(fName).sab
        name   = importMod(fName).name
        sab = [sab[b+a*len(betas)]*np.exp(-betas[b]*0.5) for a in range(len(alphas)) for b in range(len(betas))]
        return temp,alphas,betas,sab,None,None,None,None,name




if __name__=='__main__':
    from importlib import import_module as importMod
    fig = plt.figure(num=None,figsize=(10,6),dpi=120,facecolor='w',edgecolor='k') 
    plt.rcParams.update({'font.size': 17})
    if len(sys.argv) == 1: exit()

    option = 0 # Here, I'm just plotting a simple F and H across time 
    option = 1 # Here, we look at how time step things affect the H and F functions
               # from pg 37 in the GASKET manual
    option = 2 # Here, we look at a few different sab datasets and see how the 
               # timestep sizes affect them. These should all be GASKET inputs
    option = 3 # Here, we plot one or more S(a,b) against beta for various alphas
    #option = 4 # % difference betwen two S(a,b) grids, plotted against beta for
    #           # various alphas
    #option = 5 # Plots F(t), H(t), Q and R
    #option = 6 # Single alpha sab values
    
    if option == 0:
        T,alphas,betas,sab,nt,dt,H,F,name = getValuesFromInput(sys.argv[1])
        time = [i*dt for i in range(nt)]
        invArea = 1.0/np.trapz(H,time); H = [val*invArea for val in H]
        invArea = 1.0/np.trapz(F,time); F = [val*invArea for val in F]
        scalarMap, colorBar = prepPlot(list(range(7))); plt.clf(); 
        plt.fill_between(time,H,linewidth=2,color=scalarMap.to_rgba(0),\
                         label='H(t)',alpha=0.5)
        plt.fill_between(time,F,linewidth=2,color=scalarMap.to_rgba(3),\
                         label='F(t)',alpha=0.5)
        plt.xlabel('t'); plt.ylabel('Normalized'); plt.legend(loc='best')
        plt.title('H(t) and F(t) for '+name+' at '+str(int(T/8.617e-5))+'K')
        plt.tick_params(axis='y',which='both',left=False,labelleft=False)
        plt.show()




    if option == 1:
        tNum = []; maxVal = 0.0
        scalarMap, colorBar = prepPlot(list(range(len(sys.argv)))); plt.clf();
        for i,name in enumerate(sys.argv[1:]):
            T,alphas,betas,sab,nt,dt,H,F,title = getValuesFromInput(name)
            if max(max(H),max(F)) > maxVal: maxVal = max(max(H),max(F))
            time = [i*dt for i in range(nt)]
            tNum.append(dt)
            plot_H_or_F(time,H,'solid',scalarMap.to_rgba(i))
            plot_H_or_F(time,F,'dashed',scalarMap.to_rgba(i))
        scalarMap, colorBar = prepPlot(tNum+[tNum[-1]+1])
        plt.show()




    if option == 2:
        tNum = []
        scalarMap, colorBar = prepPlot(list(range(len(sys.argv)))); plt.clf()
        for i,name in enumerate(sys.argv[1:]):
            T,alphas,betas,sab,nt,dt,H,F,name = getValuesFromInput(name)
            tNum.append(nt-1)

        scalarMap, colorBar = prepPlot(tNum); plt.clf()
        colorList = [colors.rgb2hex(scalarMap.to_rgba(i)) for i in range(len(tNum))]
        mymap = colors.ListedColormap(colorList)
        colorBar = plt.contourf([[0,0],[0,0]], list(range(len(tNum)+1)), cmap=mymap)
        plt.rcParams.update({'font.size': 14}); plt.clf()
        cbar = plt.colorbar(colorBar);     cbar.ax.get_yaxis().set_ticks([])
        cbar.ax.get_yaxis().labelpad = 50; cbar.ax.set_ylabel('tSpacing')
        for j, lab in enumerate(["    "+str(a) for a in tNum[:]]):
            cbar.ax.text(1.7,(2*j+1)/(2*len(tNum)),lab,ha='center',va='center')


        for i,name in enumerate(sys.argv[1:]):
            T,alphas,betas,sab,nt,dt,H,F,name = getValuesFromInput(name)
            chunk = [sab[b] for b in range(len(betas))]
            for b in range(0,len(chunk)-2,10):
                for bp in range(10):
                    if chunk[b+bp] < 0.0:
                        for bp2 in range(10): chunk[b+bp2] = 1e-8
                        break

            plt.fill_between(betas[:-1],chunk[:-1],color=scalarMap.to_rgba(i),\
                             linewidth=2-0.2*i,alpha=0.5)

        plt.xlabel('beta'); plt.ylabel('S(a,b)'); 
        plt.title('S(a,b) for '+name+' at '+str(int(T/8.617e-5))+\
                  'K for various dt values. alpha = '+str(alphas[0]))
        plt.yscale('log'); plt.show()

    if option == 3:
        patterns = ['solid','dashed','dotted']
        for i,name in enumerate(sys.argv[1:]):
            T,alphas,betas,sab,nt,dt,H,F,name = getValuesFromInput(name)

            if i == 0:
                scalarMap, colorBar = prepPlot(alphas); plt.clf()
                colorList = [colors.rgb2hex(scalarMap.to_rgba(i)) for i in range(len(alphas))]
                mymap = colors.ListedColormap(colorList)
                colorBar = plt.contourf([[0,0],[0,0]], list(range(len(alphas)+1)), cmap=mymap)
                plt.clf()
                cbar = plt.colorbar(colorBar);     cbar.ax.get_yaxis().set_ticks([])
                cbar.ax.get_yaxis().labelpad = 50; cbar.ax.set_ylabel('alpha')


                plt.rcParams.update({'font.size': 14})
                for j, lab in enumerate(["    "+str(a) for a in alphas[:]]):
                    cbar.ax.text(1.7,(2*j+1)/(2*len(alphas)),lab,ha='center',va='center')
            else:
                checkThatAlphasBetasMatch(alphasOld,betasOld,alphas,betas)
            alphasOld, betasOld = alphas[:], betas[:]
            doThePlotting(alphas,betas,sab,patterns[i],scalarMap)
        plt.xlabel('beta'); plt.ylabel('S(a,-b)'); 
        #plt.title('S(a,-b) for '+name+' at '+str(int(T/8.617e-5))+'K (NJOY vs. GASKET)')
        plt.title('Solid = non-quad gasket, dotted = quad gasket')
        plt.yscale('log'); plt.show()



    if option == 4 and len(sys.argv) == 3:

        T1,alphas1,betas1,sab1,nt1,dt1,H1,F1,name1 = getValuesFromInput(sys.argv[1])
        T2,alphas2,betas2,sab2,nt2,dt2,H2,F2,name2 = getValuesFromInput(sys.argv[2])
        assert(name1 == name2)
        checkThatAlphasBetasMatch(alphas1,betas1,alphas2,betas2)
        alphas = alphas1[:]; betas = betas1[:]

        scalarMap = cmx.ScalarMappable(norm=colors.Normalize(0,len(alphas)+2),\
                                       cmap=plt.get_cmap('tab20')) 

        fig, axs = plt.subplots(3,figsize=(6,8),dpi=100)
        for a in range(1,len(alphas)):
            chunk1 = [sab1[b+a*len(betas)] for b in range(len(betas)-2)]
            chunk2 = [sab2[b+a*len(betas)] for b in range(len(betas)-2)]
            pct_diff = [100.0*(chunk1[b]-chunk2[b])/chunk1[b] for b in \
                        range(len(betas)-2)]
            rel_diff = [chunk1[b]-chunk2[b] for b in range(len(betas)-2)]

            axs[0].plot(betas[:-2],chunk1,color=scalarMap.to_rgba(a),\
                     linewidth=2,label='alpha = '+str(alphas[a]),linestyle='solid')
            axs[0].plot(betas[:-2],chunk2,color=scalarMap.to_rgba(a),\
                     linewidth=2,linestyle='dashed')
            axs[1].plot(betas[:-2],rel_diff,color=scalarMap.to_rgba(a),\
                     linewidth=2,label='alpha = '+str(alphas[a]))
            axs[2].plot(betas[:-2],pct_diff,color=scalarMap.to_rgba(a),\
                     linewidth=2,label='alpha = '+str(alphas[a]))

        axs[0].legend(loc='best')
        axs[1].legend(loc='best')
        axs[2].legend(loc='best')
        axs[0].set(ylabel='S(a,b) Values'); 
        axs[1].set(ylabel='Relative Diff.'); 
        axs[2].set(xlabel='beta', ylabel='Percent Diff.'); 
        plt.title('S(a,b) for '+name1+' at '+str(int(T1/8.617e-5))+\
                  'K as prepared by GASKET (t = 0->600, dt = 0.1)')
        plt.show()

    if option == 5:
        T,alphas,betas,sab,nt,dt,H,F,name = getValuesFromInput(sys.argv[1])
        alphas[0] = 5.0
        Q = [cos(alphas[0]*F[i])*exp(alphas[0]*H[i])-1 for i in range(len(F))]
        R = [sin(alphas[0]*F[i])*exp(alphas[0]*H[i]) for i in range(len(F))]
        time = [i*dt for i in range(nt)]
        scalarMap, colorBar = prepPlot(list(range(5))); plt.clf(); 
        plt.fill_between(time[1:],H[1:],linewidth=2,color=scalarMap.to_rgba(0),\
                         label='H(t)',alpha=0.3)
        plt.fill_between(time,F,linewidth=2,color=scalarMap.to_rgba(1),\
                         label='F(t)',alpha=0.3)
        plt.fill_between(time,Q,linewidth=2,color=scalarMap.to_rgba(2),\
                         label='Q(t)',alpha=0.3)
        plt.fill_between(time,R,linewidth=2,color=scalarMap.to_rgba(3),\
                         label='R(t)',alpha=0.3)

        plt.xlabel('t'); plt.ylabel('Normalized'); plt.legend(loc='best')
        plt.title('H(t) and F(t) for '+name+'at '+str(int(T/8.617e-5))+'K')
        plt.show()

    if option == 6:
        patterns = ['solid','dashed','dotted']
        for i,name in enumerate(sys.argv[1:]):
            T,alphas,betas,sab,nt,dt,H,F,name = getValuesFromInput(name)
            doThePlottingSingleAlpha(alphas,betas,sab,patterns[i])
        plt.xlabel('beta'); plt.ylabel('S(a,-b)'); 
        #plt.title('S(a,-b) for '+name+' at '+str(int(T/8.617e-5))+'K (NJOY vs. GASKET)')
        plt.title('Solid = non-quad gasket, dotted = quad gasket')
        plt.yscale('log'); 
        plt.show()




