import numpy as np
import matplotlib.pyplot as plt


def get_a_b(n):
    h = 1 / n
    x = np.linspace(0, 1, n + 1)
    # Get A
    A = np.zeros((n + 1, n + 1))
    A[0, 0] = 1
    A[n, n] = 1
    A[n, n - 1] = 0
    for i in range(1, n):
        A[i, i - 1] = 1
        A[i, i] = -2 - (np.pi**2) * h ** 2
        A[i, i + 1] = 1

    # Get b
    b = np.zeros(n + 1)
    for i in range(1, n):
        b[i] = -2 * (h ** 2) *(np.pi**2)*np.sin(np.pi*x[i])

    return x, A, b


n_s = [2, 4, 8, 16, 32, 64, 128, 256, 1000]
errors = []
max_error = []
y = []
for n in n_s:
    x, A, b = get_a_b(n)
    temp = np.linalg.inv(A)
    y1 = np.dot(temp, b)
    y.append(y1)
    x2 = np.linspace(0, 1, n + 1)
    v = np.sin(np.pi * x2)
    for i in range(1,len(y1+1)):
        e = abs(v[i] - y1[i])
        errors.append(e)
    max_error.append(max(errors))
x3 = np.linspace(0, 1, 9)
v3 = np.sin(np.pi * x3)
print(y[2],'/', v3 )
plt.figure(figsize=(10, 8))
plt.plot(n_s, max_error)
#plt.yscale('log')
plt.xlabel('n gird points')
plt.ylabel('errors at x = $\pi/2$')
plt.show()
