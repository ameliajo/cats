import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import numpy as np
import sys
sys.path.append("NJOY/")
sys.path.append("GASKET/")

#from njoy_output_1 import *
#from gasket_output_1 import *
from njoy_output_2 import *
from gasket_output_2 import *
#from njoy_output_3 import *
#from gasket_output_3 import *





def prepPlot(vector):
    cnorm = colors.Normalize(vmin=0,vmax=2*len(vector)+1)
    scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('tab20')) #hot autumn tab10
    mymap = colors.LinearSegmentedColormap.from_list('funTestColors',\
            [scalarMap.to_rgba(a) for a in range(len(vector))])
    colorBar = plt.contourf([[0,0],[0,0]], vector, cmap=mymap)
    plt.clf()
    return scalarMap, colorBar


assert(len(alphas) == len(gasket_alphas))
assert(len(betas) == len(gasket_betas))
for a in range(len(alphas)): assert(abs(alphas[a]-gasket_alphas[a])<1e-5)
for b in range(len(betas)):  assert(abs(betas[b]-gasket_betas[b])<1e-5)



#fig = plt.figure(num=None, figsize=(17, 12), dpi=80, facecolor='w', edgecolor='k')
fig = plt.figure()
scalarMap, colorBar = prepPlot(alphas)


for a in range(0,len(alphas)):
    njoy_sym = [sab_nonsym_negative_side[b+a*len(betas)] * \
                      np.exp(-betas[b]*0.5) for b in range(len(betas))]
    gasket_sym = [0.02550730568*gasket_sab[b+a*len(betas)] for b in range(len(betas))]
    difference = [gasket_sym[i] - njoy_sym[i] for i in range(len(njoy_sym))] 
    plt.plot(betas,njoy_sym,color=scalarMap.to_rgba(a),linewidth=3)
    plt.plot(betas[:-1],gasket_sym[:-1],color=scalarMap.to_rgba(a),\
             marker='o',markersize=3,linewidth=2,linestyle='dashed')
    #plt.plot(betas[:-1],difference[:-1],color=scalarMap.to_rgba(a))



ax = plt.gca()
plt.colorbar(colorBar).ax.set_ylabel('alpha values')
#plt.yscale('log')
plt.show()
#plt.draw()
#plt.pause(1) 
#input("<Hit Enter To Close>")
#plt.close(fig)


