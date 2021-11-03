import numpy as np
from Laplace import *


n = [10,20,30,50,80,100,150,200]
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
    m = Mesh((0,1),h=h,gtype="1D")
    mesh = m.get()
    mesh[0],mesh[-1] = 0,0
    B = ne1d(mesh,m.x_dom)
    er =abs(f(m.x_dom)-B)
    #print(er)
    lst.append(max(er))
    rms.append(np.linalg.norm((f(m.x_dom)-B))/len(mesh))
    #plt.plot(m.x_dom,er)
    #for j in B[1]:
    #    plt.plot(m.x_dom,j)

print(rms)
plt.plot(n,lst)
plt.show()
plt.plot(n,rms)
plt.xscale("log")
plt.show()