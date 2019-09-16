import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import numpy as np
import sys
sys.path.append("NJOY/")
sys.path.append("GASKET/")

from njoy_output_1 import *
from gasket_output_1 import *

def prepPlot(vector):
    cnorm = colors.Normalize(vmin=0,vmax=2*len(vector)+1)
    #scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('hot')) #hot autumn tab10
    scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('autumn')) #hot autumn tab10
    mymap = colors.LinearSegmentedColormap.from_list('funTestColors',\
            [scalarMap.to_rgba(a) for a in range(len(vector))])
    colorBar = plt.contourf([[0,0],[0,0]], vector, cmap=mymap)
    plt.clf()
    return scalarMap, colorBar







assert(len(njoy_alphas) == len(gasket_alphas))
assert(len(njoy_betas) == len(gasket_betas))
for a in range(len(njoy_alphas)): assert(abs(njoy_alphas[a]-gasket_alphas[a])<1e-5)
for b in range(len(njoy_betas)):  assert(abs(njoy_betas[b]-gasket_betas[b])<1e-5)

alphas = njoy_alphas[:]
betas  = njoy_betas[:]


fig = plt.figure(num=None, figsize=(17, 12), dpi=80, facecolor='w', edgecolor='k')
#fig = plt.figure()
scalarMap, colorBar = prepPlot(alphas)


for a in range(0,len(alphas),int(len(alphas)/20)):
    njoySAB_chunk   = [njoy_sab_output[b+a*len(betas)]*np.exp(-betas[b]*0.5) for b in range(len(betas))]
    gasketSAB_chunk = [0.0255*gasket_sab_output[b+a*len(betas)] for b in range(len(betas))]
    plt.plot(betas,njoySAB_chunk,color=scalarMap.to_rgba(a),linewidth=2)
    plt.plot(betas,gasketSAB_chunk,color=scalarMap.to_rgba(a),marker='o',markersize=3,linewidth=0)



ax = plt.gca()
plt.colorbar(colorBar).ax.set_ylabel('alpha values')
plt.yscale('log')
plt.draw()
plt.pause(1) 
raw_input("<Hit Enter To Close>")
plt.close(fig)


