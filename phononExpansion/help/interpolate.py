
def interpolate1(yVec, xVal, delta):
    if xVal >= len(yVec)*delta or xVal < 0: return 0.0
    for i in range(len(yVec)-1):
        if i*delta <= xVal and xVal < (i+1)*delta:
            return yVec[i] + (xVal-delta*i) * (yVec[i+1]-yVec[i]) / delta
    return 0.0

def interpolate(xVec, yVec, xVal):
    if xVal >= xVec[-1] or xVal < 0: return 0.0
    for i in range(len(yVec)-1):
        if xVec[i] <= xVal and xVal < xVec[i+1]:
            return yVec[i] + (xVal-xVec[i]) * (yVec[i+1]-yVec[i]) / (xVec[i+1]-xVec[i])
    return 0.0

