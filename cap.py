from Laplace import *
import numpy as np
from mpl_toolkits import mplot3d

ep = 8.854e-3
#us,xs=(8.85*40**2)*10**-12,1/40
us = 1e+40
xs=np.sqrt(us*ep)
mm = Mesh((0,4),(0,4.4),0.1)
h= np.diff(mm.x_dom)[0]/xs
#lower_bound_y = np.ones((mm.x_dim,))
#upper_bound_y = np.ones((mm.x_dim,))
lower_bound_x = np.ones((mm.y_dim,))*5/us
upper_bound_x = np.ones((mm.y_dim,))*5/us
mm.dirichlet([lower_bound_x,upper_bound_x,None,None])
ini= mm.get()
X,Y = np.meshgrid(mm.x_dom,mm.y_dom)
distri = np.zeros(ini.shape)
print(X)
bool_B = np.any([X==0.5,X==1.5,X==2.5,X==3.5],axis=0) & (Y >= 0.3)
bool_A = (X==2) & (Y <=4)
#distri[bool_B] = -5/4.1
#distri[bool_A] = 5/4
#print(distri[bool_B])
#print(bool_B[:,25])
#print(INPUT2D)
#print(mm.x_dom)
pfield = sor2dpoisson(ini,h,overcf=1.9,p=distri)*us
#vect= np.gradient(-1*pfield)
np.savetxt("potential.dat",pfield)

#pfield =np.loadtxt("potential.dat")
#plt.contourf(X,Y,pfield,50,cmap="coolwarm")
#print(pfield)
#plt.quiver(X,Y,vect[1],vect[0])
fig = plt.figure()
ax = plt.axes(projection="3d")
ax.plot_surface(X,Y,pfield,cmap = "coolwarm", edgecolor='none')
plt.savefig("plot.png")
plt.show()
#x=np.linspace(5,100,100)
#plt.plot(x,gauss1d(INPUT))
t1= t.time()
print(t1-t0)
