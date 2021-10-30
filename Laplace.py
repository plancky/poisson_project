#from mpl_toolkits import mplot3d
#from numba import jit
import numpy as np
import matplotlib.pyplot as plt
import time as t

class Mesh:
    def __init__(self,x=(0,1),y=(0,1),h=0.1,gtype="2D"):
        self.x_dim,self.y_dim=int((x[1]-x[0])/h +1),int((y[1]-y[0])/h +1) 
        self.y_dom=np.linspace(y[0],y[1],self.y_dim)
        self.x_dom=np.linspace(x[0],x[1],self.x_dim)
        self.gtype= gtype 
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

#@jit("f8[:,:](f8[:,:],f4,f8[:,:])",nopython=True,nogil=True)
def sor2dpoisson(x,h,overcf=1.9,p=None):
    k,m = x.shape[0],x.shape[1]
    if p is None :
        p = np.zeros(x.shape)
    for f in range(500):
        new_x= np.array(x)
        for i in range(0,k):
            for j in range(1,m-1):
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
        er = (abs(new_x - x) / new_x).max()
        print(er)
        if er<=0.5e-4: #Why relative error instead of abs ?# significant digits and Decimal places
            print(f)
            break
        else:
            x = new_x
    return(new_x)

if __name__ == "__main__":

    t0 = t.time()
    '''    
    m = Mesh((0,1),h=0.1,gtype="1D")
    mesh = m.get()
    mesh[0],mesh[-1] = 0,-1
    f = lambda x : -1*x 
    B = gauss1d(mesh,f(m.x_dom))
    print(max(abs(f(m.x_dom)-B)))
    plt.plot(m.x_dom,B)
    #for j in B[1]:
    #    plt.plot(m.x_dom,j)
    t1= t.time()
    print(t1-t0)
    plt.show()
    '''
    mm = Mesh((0,2),(0,2),h=0.1)
    bound_x = np.ones((mm.y_dim,))
    bound_y = np.ones((mm.x_dim,))
    mm.dirichlet([bound_x*2,bound_x*2.1,None,None])
    INPUT2D = mm.get()
    #print(INPUT2D)
    X, Y = np.meshgrid(mm.x_dom,mm.y_dom)
    #print(mm.x_dom)
    
    pfield = sor2dpoisson(INPUT2D,0.1)
    vect= np.gradient(-1*pfield)
    plt.contourf(X,Y,pfield,38,cmap="coolwarm")

    #plt.quiver(X,Y,vect[1],vect[0])
    #fig = plt.figure()
    #ax = plt.axes(projection="3d")
    #ax.plot_surface(X,Y,sor2dpoisson(INPUT2D,charge=None),cmap = "coolwarm", rstride=1, cstride=1, edgecolor='none')
    plt.savefig("plot.png")
    plt.show()
    #x=np.linspace(5,100,100)
    #plt.plot(x,gauss1d(INPUT))
    t1= t.time()
    print(t1-t0)
    

