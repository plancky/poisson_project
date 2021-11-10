#from mpl_toolkits import mplot3d
from lib.poisson import *
import numpy as np
import time as t

'''Set constants'''
ep = 8.854e-12
us = 1
xs = 1
'''Set Domain'''
mm = mesh((0,4),(0,4.4),0.1)
h = np.diff(mm.x_dom)[0]
'''Set Dirichlet Boundary value conditions''' 
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

'''Set charge distributiion'''
distri = np.zeros(mm.grid.shape)
xr,yr = np.round_(mm.X,2),np.round_(mm.Y,2)
bool_A = ((xr==2) | (xr==1) | (xr==3)) & (yr <=4) 
bool_B = np.any([xr==0.5,xr==1.5,xr==2.5,xr==3.5],axis=0) & (yr >= 0.4) 
distri[bool_A] = 2*1e+5*xs**2/us
distri[bool_B] = -2*1e+5*xs**2/us

'''Create a folder where data gets stored'''
mm.set_loc("t2")

'''
#Solve poisson
t1= t.time()
s1=mm.sor_poisson2d(distri,1.9,nu=us,rtol=1e-8)
print(s1[1])
t2= t.time()
s2=mm.sor_poisson2d(distri,1,nu=us,rtol=1e-8)
t3= t.time()
print(s2[1])
s3=mm.jacobi_poisson2d(distri,us,rtol=1e-8)
t4= t.time()
print(s3[1])

#mm.save_omega()
print(t2-t1)
print(t3-t2)
print(t4-t3)

mean = (s1[0]+s2[0]+s3[0])/3
var = 1/3*((s1[0]-mean)**2+(s2[0]-mean)**2+(s3[0]-mean)**2) 
np.savetxt("var2.dat",np.array([mm.X.flatten("F"),mm.Y.flatten("F"),var.flatten("F")]).T)
convert_gnuplotgrid3d("var2.dat")
'''

# Variation of iterations with omega
iter = []
for i in np.arange(1,2,0.01):
    iter.append(mm.sor_poisson2d(distri,i,nu=us,rtol=1e-8)[1])
np.savetxt("omega_variation.dat",np.array([np.arange(1,2,0.01),iter]).T)