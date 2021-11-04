import numpy as np
import matplotlib.pyplot as plt
from matplotlib import style
style.use("ggplot")
import Laplace as lp
import pandas as pd

def get_a_b(n):
    h = 1 / (n+1)
    N = n+2
    x = np.linspace(0, 1, N)
    # Get A
    A = np.zeros((N, N))
    A[0, 0] = 1
    A[n+1, n+1] = 1
    for i in range(1, n+1):
        A[i, i - 1] = 1
        A[i, i] = -(2 + (np.pi*h)**2)
        A[i, i + 1] = 1
    # Get b
    b = np.zeros((N,))
    for i in range(1, N-1):
        b[i] = -2 *(h**2)*(np.pi**2)*np.sin(np.pi*x[i])
    return x, A, b


def error(a, b):
    temp1, temp2 = [], 0
    for i in range(len(b) -1):
        temp3 = abs(a[i] - b[i])
        temp1.append(temp3)
    temp2 = max(temp1)
    rms = 0
    for j in temp1:
        rms = rms + j**2
    return [temp2, np.sqrt(rms/len(a))]

marker = [ "1", "2"]


def scatter_plot(ax,a, b,d, title, markers):
    ax.scatter(a, b, color = "red" , marker = markers[0], label = "Numerical")
    ax.plot(a, d,  marker=markers[1], label = "Exact")
    ax.set_title(title)
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.legend()
    ax.grid(color = "k")


n_s = [2, 4, 8, 16, 32,64,128,256,512]
# n_s = [10,100,1000,10000]

solution, num_solution = [], []
x_values = []
iterations = []
maxerr, rmserr = [], []


for n in n_s:
    x, A, b = get_a_b(n)
    jl = np.zeros((len(x),))
    f, iteration = lp.ne1d(jl, np.sin(np.pi * x), x)
    v = np.sin(np.pi * x)
    x_values.append(x)
    solution.append(v)
    num_solution.append(f)
    iterations.append(iteration)
    c, d = error(f, v)
    maxerr.append(c), rmserr.append(d)




print(maxerr,"/",
      rmserr)
plt.plot(x_values[2], solution[2],"g" , label = "Exact")
plt.scatter(x_values[2], num_solution[2], c = "red", label = "Numerical")
plt.legend()
plt.grid(c = "k")
plt.show()

dict1 = {"x" : x_values[2], "Y_numerical" : solution[2], "Y_exact" : num_solution[2], "Ei" : np.array(num_solution[2])- np.array(solution[2])}
data1 = pd.DataFrame(dict1)
print(data1)


emax_ratio = [ "none"]
erms_ratio = ["none"]
for t in range(len(maxerr)-1):
    t = t +1
    emax_ratio.append(maxerr[t-1]/maxerr[t])
    erms_ratio.append(rmserr[t-1]/rmserr[t])

dict2 = {"n" : n_s,"E_maximum"  : maxerr, "E_max_ratio" : emax_ratio, "E_rms" : rmserr, "E_rms_ratio" : erms_ratio}
data2 = pd.DataFrame({key: pd.Series(value) for key, value in dict2.items()})
print(data2)





fig1 = plt. figure()
(ax1, ax2), (ax3, ax4) = fig1.subplots(2,2)
scatter_plot(ax1,x_values[0], num_solution[0],solution[0],"n = 2",marker )
scatter_plot(ax2,x_values[1], num_solution[1],solution[1],"n = 4 ",marker )
scatter_plot(ax3,x_values[3], num_solution[3],solution[3],"n = 16",marker )
scatter_plot(ax4,x_values[4], num_solution[4],solution[4],"n = 32",marker )
plt.tight_layout()
plt.show()

fig2 = plt. figure()
(axx1, axx2), (axx3, axx4) = fig2.subplots(2,2)
scatter_plot(axx4,x_values[8], num_solution[8],solution[8],"n = 512",marker )
scatter_plot(axx1,x_values[5], num_solution[5],solution[5],"n = 64",marker )
scatter_plot(axx2,x_values[6], num_solution[6],solution[6],"n = 128",marker )
scatter_plot(axx3,x_values[7], num_solution[7],solution[7],"256",marker )
plt.tight_layout()
plt.show()



plt.plot(n_s, maxerr)
plt.yscale('log')
plt.xscale('log')
plt.grid(c = "k")
plt.show()
plt.plot(n_s, rmserr)
plt.yscale("log")
plt.xscale('log')
plt.grid(c = "k")
plt.show()