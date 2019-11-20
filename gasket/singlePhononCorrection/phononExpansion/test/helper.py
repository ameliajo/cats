import unittest
import sys 
sys.path.append('../')

def approxEqual(x,y,tol):
    assert(abs((x-y)/y) < tol)
def approxEqualVec(a,b,tol):
    assert(len(a)==len(b))
    for i in range(len(a)):
        if abs(a[i]) < 1e-10 and abs(b[i]) < 1e-10:
            continue
        assert(abs((a[i]-b[i])/b[i]) < tol)


