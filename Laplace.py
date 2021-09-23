from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
from numpy.lib.function_base import average
import time as t


def gauss1d(x):
    y = np.array(x)
    for j in range(10000):
        for i in range(1,len(y)-1):
            y[i]=(y[i-1]+y[i+1])/2
    return(y)

def gauss2d(x):
    k,m = x.shape[0],x.shape[1]
    for f in range(500):
        for i in range(1,k-1):
            for j in range(1,m-1):
                x[i][j]= (x[i-1][j]+x[i+1][j]+x[i][j+1]+x[i][j-1])/4
    return(x)

def sor2d(x,overcf=1.93):
    k,m = x.shape[0],x.shape[1]
    for f in range(50):
        for i in range(1,k-1):
            for j in range(1,m-1):
                x[i][j] += ((x[i-1][j]+x[i+1][j]+x[i][j+1]+x[i][j-1])/4 - x[i][j]) * overcf
    return(x)

def sor2dpoisson(x,overcf=1.93,charge=[43,52,2]):
    k,m = x.shape[0],x.shape[1]
    h = np.zeros(x.shape)
    h[charge[0],charge[1]] = charge[2]
    for f in range(50):
        for i in range(1,k-1):
            for j in range(1,m-1):
                x[i][j] += ((x[i-1][j]+x[i+1][j]+x[i][j+1]+x[i][j-1] + h[i][j])/4 - x[i][j]) *overcf
    return(x)

if __name__ == "__main__":

    #starting potentials 
    INPUT=np.ones((100,),dtype="single")*5
    INPUT2D = np.ones((100,100),dtype="single") 
    ##Dirichlet Boundary values
    boundary_var = [1,0]
    boundary_var2d = np.ones(100,)*0
    INPUT2D[0],INPUT2D[-1] = boundary_var2d ,boundary_var2d /2
    INPUT[0],INPUT[-1] = boundary_var[0],boundary_var[-1]
    
    #fig = plt.figure()
    #ax = plt.axes(projection="3d")
    #x=np.linspace(5,100,100)
    #plt.plot(x,gauss1d(INPUT))
    x = np.linspace(-5,5,100)
    y = np.linspace(-5,5,100)
    X, Y = np.meshgrid(x, y)
    #print(func2d(INPUT2D))
    #ax.plot_surface(X,Y,sor2dpoisson(INPUT2D), rstride=1, cstride=1,cmap='winter', edgecolor='none')
    t0 = t.time()
    plt.pcolormesh(sor2dpoisson(INPUT2D))
    t1= t.time()
    print(t1-t0)
    
    plt.show()