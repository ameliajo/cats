import sys

def writeOutput(fileName,alphas,betas,sab,H,F,T,t,title,oscBegin=None,correctionBegin=None):
    f = open(fileName,'w')
    f.write("alphas = "+str(alphas)   +"\n")
    f.write("betas  = "+str(betas)    +"\n")
    f.write("sab    = "+str(sab)      +"\n")
    f.write("H      = "+str(H)        +"\n")
    f.write("F      = "+str(F)        +"\n")
    f.write("temp   = "+str(T)        +"\n")
    f.write("num_t  = "+str(len(t))   +"\n")
    f.write("delta_t= "+str(t[1]-t[0])+"\n")
    f.write("name   = '"+title+"'\n")
    f.write("oscBegin        = '"+oscBegin+"'\n")
    f.write("correctionBegin = '"+correctionBegin+"'\n")
    f.close()

