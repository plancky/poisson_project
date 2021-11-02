from Laplace import *
import numpy as np
from mpl_toolkits import mplot3d

ep = 8.854e-12
us = 1e+40
xs=np.sqrt(us*ep)
mm = Mesh((0,4),(0,4.4),0.1)
h= np.diff(mm.x_dom)[0]/xs

## Set Dirichlet 
'''
------------upper_y_bound-----------
lower                           upper
_x_                             _x_
bound                           bound
------------lower_y_bound-----------
'''
#lower_bound_y = np.ones((mm.x_dim,))
#upper_bound_y = np.ones((mm.x_dim,))
lower_bound_x = np.ones((mm.y_dim,))*5/us
upper_bound_x = np.ones((mm.y_dim,))*5/us
mm.dirichlet([lower_bound_x,upper_bound_x,None,None])


## get the mesh potential (initial) guess 
ini= mm.get()
X,Y = np.meshgrid(mm.x_dom,mm.y_dom)

## charge distributiion
distri = np.zeros(ini.shape)
xr,yr = np.round_(X,2),np.round_(Y,2)
bool_A = ((xr==2) | (xr==1) | (xr==3)) & (yr <=4) 
bool_B = np.any([xr==0.5,xr==1.5,xr==2.5,xr==3.5],axis=0) & (yr >= 0.3)
distri[bool_A] = 80/5*1e-6*8.854/us
distri[bool_B] = -20*1e-6*8.854/us

pfield = sor2dpoisson(ini,h,p=distri)*us

dat = np.array([X.flatten(),Y.flatten(),pfield.flatten()],dtype =float)
np.savetxt("potential_laplace.dat",dat.T)

#dat = np.loadtxt("potential_laplace.dat").T
#dat[2]*=1e+40


#plt.contourf(X,Y,pfield,50,cmap="coolwarm")
#print(pfield)
#plt.quiver(X,Y,vect[1],vect[0])
#fig = plt.figure()
#ax = plt.axes(projection="3d")
#ax.plot_surface(X,Y,pfield,cmap = "coolwarm", edgecolor='none')
#plt.savefig("plot.png")
#plt.show()
#x=np.linspace(5,100,100)
#plt.plot(x,gauss1d(INPUT))
#t1= t.time()
#print(t1-t0)
