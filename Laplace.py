from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
import time as t

gd =20
INPUT=np.ones((gd,),dtype="single")*5

class Mesh:
    def __init__(self,x,y):
        self.input2d=np.zeros((x,y),dtype="single")
        self.xdim,self.ydim=x,y
    def dirichlet(self,x):
        self.input2d[0],self.input2d[-1] = x[0] , x[1]
        self.input2d.T[0],self.input2d.T[-1] =  x[2] ,x[3]
    def get(self):
        return(self.input2d)


def gauss1d(x):
    y = np.array(x)
    for j in range(10000):
        for i in range(1,len(y)-1):
            y[i]=(y[i-1]+y[i+1])/2
    return(y)

def sor2dpoisson(x,overcf=1.9,charge=[[20,22,-2],[20,24,2]]):
    k,m = x.shape[0],x.shape[1]
    h = np.zeros(x.shape)
    if charge is not None :
        for jk in charge:
            h[jk[0]][jk[1]] = jk[2]

    for f in range(500):
        new_x= np.array(x)
        for i in range(1,k-1):
            for j in range(1,m-1):
                new_x[i][j] += ((new_x[i-1][j]+new_x[i+1][j]+new_x[i][j+1]+new_x[i][j-1] + h[i][j])/4 - new_x[i][j]) * overcf
    
        if np.allclose(x,new_x,atol=1e-5):
            print(f)
            break
        else:
            x = new_x
    return(new_x)

if __name__ == "__main__":
    mymesh = Mesh(40,40)
    bound = np.ones((mymesh.xdim,))
    mymesh.dirichlet([bound*0,bound*0,bound*2,bound*-2])
    INPUT2D = mymesh.get()
    pfield = sor2dpoisson(INPUT2D)
    vect= np.gradient(pfield)
    #fig = plt.figure()
    #ax = plt.axes(projection="3d")
    #x=np.linspace(5,100,100)
    #plt.plot(x,gauss1d(INPUT))
    x = np.linspace(0,mymesh.xdim,mymesh.xdim)
    y = np.linspace(0,mymesh.xdim,mymesh.xdim)
    #x,y = np.meshgrid(np.linspace(-5,5,10),np.linspace(-5,5,10))
    X, Y = np.meshgrid(x, y)
    #print(func2d(INPUT2D))
    t0 = t.time()
    plt.pcolormesh(pfield,cmap="coolwarm")
    #ax.plot_surface(X,Y,sor2dpoisson(INPUT2D,charge=None),cmap = "coolwarm", rstride=1, cstride=1, edgecolor='none')
    plt.quiver(X,Y,vect[1]*-1,vect[0]*-1)
    t1= t.time()
    print(t1-t0)

    plt.show()