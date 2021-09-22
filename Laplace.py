
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
from numpy.lib.function_base import average

boundary_var = [1,0]
boundary_var2d = np.ones(100,)*4
INPUT=np.ones((100,))*5
INPUT2D = np.ones((100,100)) 
INPUT2D[0],INPUT2D[-1] = boundary_var2d ,boundary_var2d /2
INPUT[0],INPUT[-1] = boundary_var[0],boundary_var[-1]

def func(x):
    y = np.array(x)
    for j in range(10000):
        for i in range(1,len(y)-1):
            y[i]=(y[i-1]+y[i+1])/2
    return(y)

def func2d(x):
    k,m = x.shape[0],x.shape[1]
    for f in range(500):
        for i in range(1,k-1):
            for j in range(1,m-1):
                x[i][j]= (x[i-1][j]+x[i+1][j]+x[i][j+1]+x[i][j-1])/4
    return(x)

def sor2d(x,overcf=1.93):
    k,m = x.shape[0],x.shape[1]
    for f in range(500):
        for i in range(1,k-1):
            for j in range(1,m-1):
                correction =  (x[i-1][j]+x[i+1][j]+x[i][j+1]+x[i][j-1])/4 - x[i][j] 
                x[i][j] += correction* overcf
    return(x)

fig = plt.figure()
ax = plt.axes(projection="3d")
#x=np.linspace(5,100,100)
#plt.plot(x,func(INPUT))
x = np.linspace(-5,5,100)
y = np.linspace(-5,5,100)
X, Y = np.meshgrid(x, y)
print(func2d(INPUT2D))
ax.plot_surface(X,Y,sor2d(INPUT2D))

plt.show()