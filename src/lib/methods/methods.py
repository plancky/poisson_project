from numba.pycc import CC
import numpy as np

cc = CC("poisson_module")
cc.verbose = True

@cc.export("sor2dpoisson","Tuple((f8[:,:],f8))(f8[:,:],f8,f8,f8[:,:],f8)")
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

@cc.export("jacobi2d","Tuple((f8[:,:],f8))(f8[:,:],f8,f8[:,:],f8)")
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
        er = np.abs((new_x - x) / new_x).max()
        if er<=rtol: #Why relative error instead of abs ?# significant digits and Decimal places
            break
        else:
            x = new_x.copy()
    return(new_x,f)

if __name__ == "__main__":
    cc.compile()
