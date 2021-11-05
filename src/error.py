import numpy as np
from lib.poisson import *
import matplotlib.pyplot as plt
n = np.arange(1,10000,1000)
lst=[]
rms=[]
'''f = lambda x : -1*x 
for i in n :
    h = 1/(i)
    m = Mesh((0,1),h=h,gtype="1D")
    mesh = m.get()
    mesh[0],mesh[-1] = 0,-1
    B = gauss1d(mesh,f(m.x_dom))
    er =max(abs(f(m.x_dom)-B))
    print(er)
    lst.append(er)
    rms.append(np.linalg.norm((f(m.x_dom)-B))/i)
    #for j in B[1]:
    #    plt.plot(m.x_dom,j)
'''
f = lambda x : np.sin(np.pi*x) 
for i in n :
    h = 1/(i)
    m = mesh((0,1),h=h,gtype="1D")
    m.grid[0],m.grid[-1] = 0,0
    B = model1d(m.grid,m.x_dom,rtol=1e-4)
    er =abs(f(m.x_dom)-B)
    lst.append(max(er))
    rms.append(np.linalg.norm((f(m.x_dom)-B))/len(m.grid))
    #plt.plot(m.x_dom,er)
    #for j in B[1]:
    #    plt.plot(m.x_dom,j)

print(rms)
plt.plot(n,lst)
plt.show()
plt.plot(n,rms)
plt.xscale("log")
plt.show()