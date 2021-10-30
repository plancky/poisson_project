dict1 = {"D":500,"C":100,"L":50,"X":10,"V":5,"I":1}

inp = "DXLLII"

def conv(x):
    l = []
    for i in x :
        l.append(dict1[i])
    return(l)

def max_(l):
    for j in range(len(l)) :
        if l[j]==max(l):
            return(j)
def recur(d):
    if d == [] or None :
        return(0)
    ini = max_(d)
    a,b = d[:ini],d[ini+1:]
    return(max(d)-recur(a)+recur(b))

def deci(x):
    d = conv(x)
    print(recur(d))

deci(inp)