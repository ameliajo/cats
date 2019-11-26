
def interpolate(yVec, xVal, delta):
    if xVal >= len(yVec)*delta or xVal < 0: return 0.0
    for i in range(len(yVec)-1):
        if i*delta <= xVal and xVal < (i+1)*delta:
            return yVec[i] + (xVal-delta*i) * (yVec[i+1]-yVec[i]) / delta
    return 0.0

