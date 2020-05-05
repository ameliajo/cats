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
        sab = [sab[b+a*len(betas)]*np.exp(betas[b]*0.5) \
               for a in range(len(alphas)) for b in range(len(betas))]
        return temp,alphas,betas,sab,num_t,delta_t,H,F,name,oscBegin,correctionBegin

    elif name[0] == 'n':
        fName = 'njoy_output_'+name[1]
        alphas = importMod(fName).alphas; betas  = importMod(fName).betas
        temp   = importMod(fName).temp;   sab    = importMod(fName).sab
        name = ''
        sab = [sab[b+a*len(betas)]*np.exp(-betas[b]*0.5) \
               for a in range(len(alphas)) for b in range(len(betas))]
        return temp,alphas,betas,sab,None,None,None,None,name,None,None












if __name__=='__main__':
    from importlib import import_module as importMod
    fig = plt.figure(num=None,figsize=(10,6),dpi=120,facecolor='w',edgecolor='k') 
    plt.rcParams.update({'font.size': 17})
    if len(sys.argv) == 1: print('Tell me what to plot!'); exit()

    num = sys.argv[1].split('.')[-1]
    arguments = ['n.'+str(num), 'g.'+str(num), 'gc.'+str(num)] if sys.argv[1][0] == 'a' else sys.argv[1:]


    option = 0
    #option = 1

    if option == 0:
        cnorm = colors.Normalize(vmin=0,vmax=len(arguments))
        #mapp = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('tab20')) 
        mapp = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('Set1')) 

        for i,name in enumerate(arguments):
            T,alphas,betas,sab,nt,dt,H,F,title,oscBegin,correctionBegin = \
              getValuesFromInput(name)
            label = 'gasket'           if name.split('.')[0] == 'g'  else \
                    'corrected gasket' if name.split('.')[0] == 'gc' else \
                    'njoy'             if name.split('.')[0] == 'n'  else '?'

            if label == 'corrected gasket':
                print(label)
                print('alpha * lambda = ',alphas[0]*H[0])
                print('exp(a*l)-1     = ',np.exp(alphas[0]*H[0])-1.0)
                print('')

            if label == 'gasket':
                print(label)
                print('alpha * lambda = ',alphas[0]*H[0]/T)
                print('exp(a*l)-1     = ',np.exp(alphas[0]*H[0]/T)-1.0)
                print('')
            if label == 'corrected gasket' and oscBegin == 'None':
                print('No correction needed!')
                continue
            for a in range(len(alphas)):
                chunk = [sab[b+a*len(betas)] for b in range(len(betas))]
                plt.plot(betas,chunk,linewidth=1.5,label=label,color=mapp.to_rgba(i))
            if oscBegin != 'None' and oscBegin != None:
                plt.plot([float(oscBegin),float(oscBegin)],\
                         [min(chunk),max(chunk)],'k',linewidth=1.5)
                plt.plot([float(correctionBegin),float(correctionBegin)],\
                         [min(chunk),max(chunk)],'k',linewidth=1.5)

        plt.xlabel('beta'); plt.ylabel('Symmetric S(a,b)'); 
        plt.title('S(a,b) for '+title+' at '+str(int(T/8.617e-5))+'K for alpha='+str(alphas[0]))
        plt.legend(loc='best')
        plt.yscale('log'); 
        plt.show()

    if option == 1:
        assert(len(arguments) == 2)
        cnorm = colors.Normalize(vmin=0,vmax=len(arguments))
        mapp = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('tab10')) 

        T,alphas,betas,sab,nt,dt,H,F,title,oscBegin,correctionBegin = getValuesFromInput(arguments[0])
        time = [i*dt for i in range(len(H))]
        invArea = 1.0/np.trapz(F,time); F = [invArea*x for x in F]
        invArea = 1.0/np.trapz(H,time); H = [invArea*x for x in H]
        plt.plot(time,H,label='H(t)')
        plt.plot(time,F,label='F(t)')

        plt.xlabel('time'); plt.ylabel('Key Gasket Functions'); 
        plt.title('F(t) and H(t) for '+title+' at '+str(int(T/8.617e-5))+'K (normalized)')
        plt.legend(loc='best')
        plt.show()




