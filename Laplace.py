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

def sor2dpoisson(x,overcf=1.9,charge=[[23,10,10],[35,23,-10]]):
    k,m = x.shape[0],x.shape[1]
    h = np.zeros(x.shape)
    if charge is not None :
        for jk in charge:
            h[jk[0]][jk[1]] = jk[2]

    for f in range(1000):
        new_x= np.array(x)    #IMP : duplicates array and avoids linking 
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
    gd = 50
    #starting potentials
    INPUT=np.ones((gd,),dtype="single")*5
    INPUT2D = np.ones((gd,gd),dtype="single") 
    #Dirichlet Boundary values
    boundary_var = [1,0]
    boundary_var2d = np.ones(gd,)*10
    INPUT2D[0],INPUT2D[-1] = boundary_var2d ,-1*boundary_var2d
    INPUT2D.T[0],INPUT2D.T[-1] =  boundary_var2d/2 ,-1*boundary_var2d/2
    INPUT[0],INPUT[-1] = boundary_var[0],boundary_var[-1]
    
    fig = plt.figure()
    ax = plt.axes(projection="3d")
    #x=np.linspace(5,100,100)
    #plt.plot(x,gauss1d(INPUT))
    x = np.linspace(-5,5,gd)
    y = np.linspace(-5,5,gd)
    X, Y = np.meshgrid(x, y)
    #print(func2d(INPUT2D))
    t0 = t.time()
    #plt.pcolormesh(sor2dpoisson(INPUT2D))
    ax.plot_surface(X,Y,sor2dpoisson(INPUT2D,charge=None), rstride=1, cstride=1, edgecolor='none')
    t1= t.time()
    print(t1-t0)

    plt.show()