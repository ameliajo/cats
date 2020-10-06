import matplotlib.pyplot as plt
from math import ceil 

def makeFloat(thisNum):
    if   thisNum[-3] == '-':
        thisNum = thisNum[:-3]+'e-'+thisNum[-2:]
    if   thisNum[-2] == '-':
        thisNum = thisNum[:-2]+'e-'+thisNum[-1:]
    elif thisNum[-3] == '+':
        thisNum = thisNum[:-3]+'e+'+thisNum[-2:]
    elif thisNum[-2] == '+':
        thisNum = thisNum[:-2]+'e+'+thisNum[-1:]
    return float(thisNum)


def getVal(line,index,numberType='float'):
    if numberType == 'float':
        return float(line.split()[index].replace('-','e-').replace('+','e+'))
    return int(line.split()[index].replace('-','e-').replace('+','e+'))



# info = {}
# with open('leaprRun/tape24','r') as f:
#     lines = [line for line in f.readlines() if line[71:75] == '7  4']

#     nbeta  = getVal(lines[4],5,'int')
#     nalpha = getVal(lines[6],5,'int')
#     numAlphaLines = ceil(nalpha/3)
#     lines = lines[6:]

#     for ibeta in range(nbeta):
#         offset = 2
#         beta = getVal(lines[0],1)
#         alphas, sab = [], []
#         for i in range(numAlphaLines):
#             for j in range(3):
#                 alphas.append(getVal(lines[offset+i],2*j  ))
#                 sab.append(   getVal(lines[offset+i],2*j+1))
#         info[beta] = [alphas,sab]
#         lines = lines[numAlphaLines+offset:]

#     betas = (list(info.keys()))

#     #for ibeta in [1,3,5,10,50,100]:
#     #    beta = betas[ibeta]
#     #    alphas = info[beta][0]
#     #    sab    = info[beta][1]
#     #    plt.plot(alphas,sab)

#     sab = [info[beta][1][10] for beta in betas]
#     plt.plot(betas,sab)
#     sab = [info[beta][1][50] for beta in betas]
#     plt.plot(betas,sab)
#     sab = [info[beta][1][100] for beta in betas]
#     plt.plot(betas,sab)



#     #plt.plot(alphas,sab)
#     plt.yscale('log')
#     plt.show()


"""
info = {}
with open('leaprRun/tape24','r') as f:
    counter = 0
    alphas = []
    sab    = []
    nalpha = None
    beta = None
    #lines = []
    beginningOffset = 6

    lines = []
    for line in f.readlines():
        if line[71:75] == '7  4':
            counter += 1
            if counter == 5:
                nalpha = int(line.split()[5])

                #print(ceil(nalpha/3))

            if counter > 8 and counter < 8+ceil(nalpha/3):
                for i in range(6):
                    lines.append(getVal(line,i))
    print(lines[:6])
    print(lines[6:12])

"""

#exit()
"""
    for line in f.readlines():
        if line[71:75] == '7  4':
            counter += 1
            if counter < 5: continue
            if counter == 5:
                nalpha = int(line.split()[5])
            if counter == 7:
                beta   = getVal(line,1)
                nbeta  = int(line.split()[5])
            if 3*counter > 8+nalpha/3.0:
                info[beta] = [alphas,sab]

                break
            if counter > 8:
                for i in range(3):
                    alphas.append(getVal(line,2*i  ))
                    sab.append(   getVal(line,2*i+1))


"""


    #from pprint import pprint
    #pprint(alphas)
    #import matplotlib.pyplot as plt
    #plt.plot(alphas,sab)
    #plt.yscale('log')
    #plt.show()

