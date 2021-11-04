#from numba import jit
import numpy as np
import time as t
import os
from .methods import poisson_module as pm

class mesh:
    def __init__(self,x=(0,1),y=(0,1),h=0.1,gtype="2D"):
        self.x_dim,self.y_dim=int((x[1]-x[0])/h +1),int((y[1]-y[0])/h +1) 
        self.y_dom=np.linspace(y[0],y[1],self.y_dim)
        self.x_dom=np.linspace(x[0],x[1],self.x_dim)
        self.h = h
        self.omega = dict()
        self.gtype= gtype 
        self.X,self.Y = np.meshgrid(self.x_dom,self.y_dom)
        if self.gtype == "1D":
            self.grid=np.ones((self.x_dim,),dtype="single")
        elif self.gtype == "2D":
            self.grid=np.ones((self.y_dim,self.x_dim),dtype="double")

    def dirichlet(self,x,excluded = "y"):
        if self.gtype=="2D" and len(x)==4:
            if excluded != "y":
                self.grid[0,:],self.grid[-1,:] = x[2] , x[3]
            if excluded !="x":
                self.grid[:,0],self.grid[:,-1] =  x[0] , x[1]
        elif self.gtype == "1D" and len(x)==2:
            self.grid[0],self.grid[-1] = x[0] , x[1]
        else:
            raise ValueError("Expected {0} but recieved {1} Boundary conditions.".format(self.gtype[0],len(x)))
    def get(self):
        return(self.grid)
    def set_loc(self,loc):
        self.loc = os.getcwd()+f"/dats/{loc}"
        os.mkdir(self.loc)

    def jacobi_poisson2d(self,p,nu=1,rtol=None):
        self.u = pm.jacobi2d(self.grid,self.h,p,rtol)[0]*nu
        self.dat = np.array([self.X.flatten(),self.Y.flatten(),self.u.flatten()],dtype =float)
        np.savetxt(f"{self.loc}/jacobi_poisson2d_{self.h}.dat",self.dat.T,fmt="%.20g")

    def sor_poisson2d(self,p,w,nu=1,rtol=None):
        solve = pm.sor2dpoisson(self.grid,self.h,w,p,rtol)
        self.u = solve[0]*nu
        self.omega[w] = solve[1]
        self.dat = np.array([self.X.flatten(),self.Y.flatten(),self.u.flatten()],dtype =float)
        np.savetxt(f"{self.loc}/sor_poisson2d_{self.h}_{w}_{nu}.dat",self.dat.T,fmt="%.20g")
    
    def save_omega(self):
        np.savetxt(f"/home/planck/Desktop/secc_project_proposal/dats/{self.loc}/omega_variation.dat",
        np.array([list(self.omega.keys()),list(self.omega.values())]).T)

def gauss1d(x,an):
    y=np.array((x,),dtype=float)
    for j in range(20000):
        y=np.vstack((y,y[-1]))
        for i in range(1,len(y[-1])-1):
            y[-1][i]=(y[-1][i-1]+y[-1][i+1])/2
        er = max(abs((y[-1]-y[-2])/y[-1])[1:-1])
        if er<=0.5e-14 and y.shape[0]>=2:
            print(j)
            break
    return(y[-1])

def ne1d(x,g):
    h= g[1] - g[0]
    y=np.array((x,),dtype=float)
    for j in range(int(5e+4)):
        y=np.vstack((y,y[-1]))
        for i in range(1,len(y[-1])-1):
            y[-1][i]=(y[-1][i-1]+y[-1][i+1]+2*(h**2)*(np.pi**2)*np.sin(np.pi*g[i]))/(2 + (np.pi*h)**2)
        er = max(abs((y[-1]-y[-2])/y[-1])[1:-1])
        if er<=0.5e-8 and y.shape[0]>=2:
            print(j)
            break
    return(y[-1])
'''
@jit("Tuple((f8[:,:],f8))(f8[:,:],f8,f8,f8[:,:],f8)",nopython=True)
def sor2dpoisson(x,h,overcf=1.9,p=None,rtol=None):
    k,m = x.shape[0],x.shape[1]
    if p is None :
        p = np.zeros(x.shape)
    if rtol is None :
        rtol = h**2
    for f in np.arange(1,15000,1):
        new_x= x.copy()
        for i in np.arange(0,k,1):
            for j in np.arange(1,m-1,1):
                left = new_x[i][j-1]
                right = new_x[i][j+1]
                ##Neumann 0 wrt y-axis at y_upper and y_lower bounds 
                if i == k-1 :
                    up = new_x[i-1][j]
                else :
                    up = new_x[i+1][j]
                if i == 0 :
                    down = new_x[i+1][j]
                else:
                    down = new_x[i-1][j]
                new_x[i][j] = new_x[i][j] + ((up+down+right+left + h**2*p[i][j])/4 - new_x[i][j]) * overcf
        er = np.abs((new_x - x )/ new_x).max()
        if er<=rtol: #Why relative error instead of abs ?# significant digits and Decimal places
            print(f)
            break
        else:
            x = new_x.copy()
    return(new_x,f)

@jit("Tuple((f8[:,:],f8))(f8[:,:],f8,f8[:,:],f8)",nopython=True)
def jacobi2d(x,h,p=None,rtol=None):
    k,m = x.shape[0],x.shape[1]
    if p is None :
        p = np.zeros(x.shape)
    if rtol is None :
        rtol = h**2
    for f in np.arange(1,15000,1):
        new_x= x.copy()
        for i in np.arange(0,k,1):
            for j in np.arange(1,m-1,1):
                left = x[i][j-1]
                right = x[i][j+1]
                ##Neumann 0 wrt y-axis at y_upper and y_lower bounds 
                if i == k-1 :
                    up = x[i-1][j]
                else :
                    up = x[i+1][j]
                if i == 0 :
                    down = x[i+1][j]
                else:
                    down = x[i-1][j]
                new_x[i][j] = ((up+down+right+left + h**2*p[i][j])/4) 
        er = abs((new_x - x) / new_x).max()
        if er<=rtol: #Why relative error instead of abs ?# significant digits and Decimal places
            break
        else:
            x = new_x.copy()
    return(new_x,[f,er])
'''
if __name__ == "__main__":
    pass
