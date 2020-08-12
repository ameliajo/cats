




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

import matplotlib.pyplot as plt
for temp in temps:
    rho_y = bigDict[temp]['rho']
    rho_x = [i*bigDict[temp]['dx'] for i in range(len(rho_y))]
    plt.plot(rho_x,rho_y)
plt.show()

for temp in temps:
    kappa_y = bigDict[temp]['kappa']
    kappa_x = [i*bigDict[temp]['dka'] for i in range(len(kappa_y))]
    plt.plot(kappa_x,kappa_y)
plt.show()



    
