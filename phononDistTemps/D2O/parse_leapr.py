




with open('tsl-DinD2O.leapr','r') as f:
    numTemps = int(f.readline().split()[0])
    
    temps = []
    bigDict = {}

    for itemp in range(numTemps):
        thisDict = {}
        temp = float(f.readline().split()[0])
        temps.append(temp)

        [dx,len_rho] = f.readline().split('/')[0].split()
        thisDict['dx'] = float(dx)
        thisDict['len_rho'] = int(len_rho)
    
        rhoVec = []
        while len(rhoVec) < int(len_rho):
            thisLine = f.readline().split()
            if len(thisLine)+ len(rhoVec) >= int(len_rho):
                thisLine = thisLine[:int(len_rho)-len(rhoVec)]
            rhoVec += [float(x) for x in thisLine]

        thisDict['rho'] = rhoVec
        thisDict['twt'], thisDict['c'], thisDict['tbeta'] = \
                [float(x) for x in f.readline().split()[:3]]
        nd = int(f.readline().split()[0])
        oscE = [float(x) for x in f.readline().split()[:nd]]
        oscW = [float(x) for x in f.readline().split()[:nd]]
    
        [nka, dka] = f.readline().split('/')[0].split()
        thisDict['nka'] = int(nka)
        thisDict['dka'] = float(dka)
    
        kappaVec = []
        while len(kappaVec) < int(nka):
            thisLine = f.readline().split()
            if len(thisLine)+ len(kappaVec) >= int(nka):
                thisLine = thisLine[:int(nka)-len(kappaVec)]
            kappaVec += [float(x) for x in thisLine]
    
        thisDict['kappa'] = kappaVec
        thisDict['cfrac'] = float(f.readline().split()[0])
    
        assert(len(thisDict['rho']) == thisDict['len_rho'])
        assert(len(thisDict['kappa']) == thisDict['nka'])
     
        bigDict[temp] = thisDict



if __name__=="__main__":
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import matplotlib.cm as cmx


    def prepPlot(vec):
        cnorm = colors.Normalize(vmin=vec[0],vmax=vec[-1]+100)
        #scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('hot')) #hot autumn tab10
        scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('tab20')) #hot autumn tab10
        scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('inferno')) #hot autumn tab10
        mymap = colors.LinearSegmentedColormap.from_list('test',\
                [scalarMap.to_rgba(vec[a]) for a in range(len(vec))])
        colorBar = plt.contourf([[0,0],[0,0]], vec, cmap=mymap)
        plt.clf()
        return scalarMap, colorBar
 

    def finishPlotting(colorBar):
        ax = plt.gca()
        plt.colorbar(colorBar).ax.set_ylabel('Temperature [K]')
        plt.show()

    
    #toPlot = ['tbeta','twt','c']
    toPlot = []
    toPlot = ['kappa']

    if 'DOS' in toPlot:
        scalarMap, colorBar = prepPlot(temps)
        for temp in temps:
            rho_y = bigDict[temp]['rho']
            rho_x = [i*bigDict[temp]['dx'] for i in range(len(rho_y))]
            plt.plot(rho_x,rho_y,color=scalarMap.to_rgba(temp))
        plt.title('Solid-type Contribution to the Phonon Distribution\nfor D in D$_2$O at Various Temperatures')
        plt.xlabel('Energy (eV)')
        plt.ylabel('Phonon DOS (unitless)')
        finishPlotting(colorBar)

    if 'kappa' in toPlot:
        scalarMap, colorBar = prepPlot(temps)
        for temp in temps:
            kappa_y = bigDict[temp]['kappa']
            kappa_x = [i*bigDict[temp]['dka'] for i in range(len(kappa_y))]
            plt.plot(kappa_x,kappa_y,color=scalarMap.to_rgba(temp))
        plt.title('Skold Correction Factors for D$_2$O')
        plt.xlabel('$\kappa$ [Angstrom]')
        plt.ylabel('Skold Correction Factors')
        plt.yscale('log')
        finishPlotting(colorBar)

    if 'tbeta' in toPlot:
        tbetaVec = [bigDict[temp]['tbeta'] for temp in temps]
        plt.plot(temps,tbetaVec)
        plt.title('tbeta Values for D$_2$O')
        plt.xlabel('Temperature (K)$')
        plt.ylabel('tbeta')
        plt.show()


    if 'twt' in toPlot:
        twtVec   = [bigDict[temp]['twt'] for temp in temps]
        plt.plot(temps,twtVec)
        plt.title('tbeta Values for D$_2$O')
        plt.xlabel('Temperature (K)$')
        plt.ylabel('tbeta')
        plt.show()



    if 'c' in toPlot:
        cVec     = [bigDict[temp]['c'] for temp in temps]
        plt.plot(temps,cVec)
        plt.title('tbeta Values for D$_2$O')
        plt.xlabel('Temperature (K)$')
        plt.ylabel('tbeta')
        plt.show()



exit()

scalarMap, colorBar = prepPlot(temps)
for temp in temps:
    rho_y = bigDict[temp]['rho']
    rho_x = [i*bigDict[temp]['dx'] for i in range(len(rho_y))]
    plt.plot(rho_x,rho_y,color=scalarMap.to_rgba(temp))
plt.title('Solid-type Contribution to the Phonon Distribution\nfor D in D$_2$O at Various Temperatures')
plt.xlabel('Energy (eV)')
plt.ylabel('Arbitrary')
finishPlotting(colorBar)



 

