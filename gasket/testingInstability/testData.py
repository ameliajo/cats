import numpy as np

title = "H in H2O"
Q  = [\
     1e-5, 1e-5, 1e-5, 10.0, 1e-5, 1e-5, \
     1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5,  \
     1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5,  \
     1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5,  \
     1e-5, 1e-5, 1e-5, 0.0]
Q2 = [\
     1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5,  \
     1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5,  \
     1e-5, 1e-5, 1e-5, 10.0, 0.0 ]
Q3 = [\
     1e-5, 1e-5, 1e-5, 10.0, 1e-5, 1e-5, \
     1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5,  \
     1e-5, 1e-5, 1e-5, 10.0, 1e-5, 1e-5, \
     1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5,  \
     1e-5, 1e-5, 1e-5, 0.0]

Q4 = [\
     5e-6, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 5.000005, 10.0, 5.000005, 1e-5, 1e-5, 1e-5, \
     1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5,  \
     1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5,  \
     1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5,  \
     1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 5e-6, 0.0]

X  = [0.01*i for i in range(1,len(Q)+1)]
X2 = [0.005*i for i in range(1,len(Q2)+1)]
X3 = [0.01*i for i in range(1,len(Q3)+1)]
X4 = [0.005*i for i in range(1,len(Q4)+1)]

Q1 = Q[:]
X1 = X[:]

invArea = 1.0/np.trapz(Q1,X1); Q1 = [invArea * Qval for Qval in Q1]
invArea = 1.0/np.trapz(Q2,X2); Q2 = [invArea * Qval for Qval in Q2]
invArea = 1.0/np.trapz(Q3,X3); Q3 = [invArea * Qval for Qval in Q3]
invArea = 1.0/np.trapz(Q4,X4); Q4 = [invArea * Qval for Qval in Q4]


