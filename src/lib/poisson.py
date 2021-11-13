from numba import jit
import numpy as np
import time as t
import os
from .methods import poisson_module as pm

def convert_gnuplotgrid3d(FILE):
    f = open(FILE,'r')
    prev_x = f.readline().split()[0]
    lst1= ""
    for line in f :
        new_x = line.split(' ')[0]
        if new_x!=prev_x:
            lst1+="\n"+line
            prev_x = new_x
        else:
            lst1+=line
    f.close()
    f2 =open(FILE,"w")
    f2.write(lst1)
    f2.close()

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
            self.grid=np.ones((self.x_dim,),dtype="double")
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
        try:
            os.mkdir(self.loc)
        except:
            pass
    def jacobi_poisson2d(self,p,nu=1,rtol=None):
        solve = jacobi2d(self.grid,self.h,p,rtol)
        self.u = solve[0]*nu
        self.dat = np.array([self.X.flatten('F'),self.Y.flatten('F'),self.u.flatten('F')],dtype =float)
        np.savetxt(f"{self.loc}/jacobi_poisson2d_{self.h}.dat",self.dat.T,fmt="%.20g")
        convert_gnuplotgrid3d(f"{self.loc}/jacobi_poisson2d_{self.h}.dat")
        return(solve)
    def sor_poisson2d(self,p,w,nu=1,rtol=None):
        solve = sor2dpoisson(self.grid,self.h,w,p,rtol)
        self.u = solve[0]*nu
        self.omega[w] = solve[1]
        self.dat = np.array([self.X.flatten('F'),self.Y.flatten('F'),self.u.flatten('F')],dtype =float)
        np.savetxt(f"{self.loc}/sor_poisson2d_{self.h}_{w}_{nu}.dat",self.dat.T,fmt="%.20g")
        convert_gnuplotgrid3d(f"{self.loc}/sor_poisson2d_{self.h}_{w}_{nu}.dat")
        return(solve)
    
    def save_omega(self):
        np.savetxt(f"/home/planck/Desktop/secc_project_proposal/dats/{self.loc}/omega_variation.dat",
        np.array([list(self.omega.keys()),list(self.omega.values())]).T)


@jit("(f8[:],f8[:],f8,f8)",nopython=True)
def sor1dpoisson(x,an,overcf=1,rtol=1e-8):
    old_x = x.copy()
    for j in np.arange(1,1e+4):
        y= old_x.copy()
        for i in np.arange(1,len(y)-1):
            y[i]=y[i] + ((y[i-1]+y[i+1])/2 - y[i])*overcf
        er = np.abs((y-old_x)/y).max()
        if er<=rtol :
            print(j)
            break
        else:
            old_x= y.copy()
    return(y,j)

@jit("f8[:](f8[:],f8[:],f8)",nopython=True)
def model1d(x,g,rtol=1e-8):
    h= g[1] - g[0]
    old_x= x.copy()
    for j in np.arange(1,5e+4,1):
        y = old_x.copy()
        for i in np.arange(1,len(y)-1):
            y[i]=(y[i-1]+y[i+1]+2*(h**2)*(np.pi**2)*np.sin(np.pi*g[i]))/(2 + (np.pi*h)**2)
        er = np.abs((y-old_x)/y).max()
        if er<=rtol :
            print(j)
            break
        else:
            x_old = y.copy()
    return(y)

@jit("Tuple((f8[:,:],f8))(f8[:,:],f8,f8,f8[:,:],f8)",nopython=True,cache=True)
def sor2dpoisson(x,h,overcf=1.9,p=None,rtol=None):
    k,m = x.shape[0],x.shape[1]
    if p is None :
        p = np.zeros(x.shape)
    if rtol is None :
        rtol = h**2
    for f in np.arange(1,1e+6,1):
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
        if er<=rtol:
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
    for f in np.arange(1,1e+6,1):
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
        er = np.abs((new_x - x) / new_x).max()
        if er<=rtol:
            break
        else:
            x = new_x.copy()
    return(new_x,f)

if __name__ == "__main__":
    pass
