import numpy as np
import matplotlib.pyplot as plt
import Laplace as lp

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


n_s = [2, 4, 8, 16, 32, 64, 128, 256, 1000]
errors = []
max_error = []
y = []

x, A, b = get_a_b(5)
print(A)
for n in n_s:
    x, A, b = get_a_b(n)
    temp = np.linalg.inv(A)
    y1 = np.dot(temp, b)
    y.append(y1)
    x2 = np.linspace(0, 1, n+2)
    v = np.sin(np.pi * x2)
    for i in range(len(y1)):
        e = abs(v[i] - y1[i])
        errors.append(e)
    max_error.append(max(errors))

#x3 = np.linspace(0, 1, 9)
#v3 = np.sin(np.pi * x3)
#print(y[2],'/', v3 )

plt.figure(figsize=(10, 8))
#plt.plot(n_s, max_error)


x, A, b = get_a_b(100)
jl = np.zeros((len(x),))

f=lp.ne1d(jl,np.sin(np.pi*x),x)

plt.plot(x,f[0])
plt.plot(x,np.sin(np.pi*x),ls="--")

#plt.yscale('log')
plt.xlabel('n gird points')
plt.ylabel('errors at x = $\pi/2$')
plt.show()
