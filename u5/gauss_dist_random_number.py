import random
import math
import matplotlib.pyplot as plt
import numpy as np

def gaussvar(mu, sigma):
    u1 = random.random()
    u2 = random.random()
    
    # Box Muller Method
    N = (math.sqrt(-2 * math.log(u1)) * math.cos(2 * math.pi * u2) * sigma) + mu
    return N

def histgauss(mu, sigma, k, m):
    N = np.zeros(k)
    for i in range(k):
        N[i] = gaussvar(mu, sigma)
    plt.hist(N, density=False, bins=m)
    plt.ylabel('occurances')
    plt.xlabel('Random Numbers')
    plt.show()

# T: Time of whole process
# n: Number of time steps
# N: Number of Wiener Processes
def wienerprocess(T, n, N):
    h = T / n
    mu = 0
    sigma = math.sqrt(h)

    W = np.zeros((n,N))
    print(W.shape)
    for i in range(N):
        W[i,0] = 0
        for j in range(1, n):
            d_W = gaussvar(mu, sigma)
            W[j,i] = W[j-1, i] + d_W

    t = list(range (n))

    fig, ax = plt.subplots()

    for w in range(N):
        ax.plot(t, W[:, w])
    # ax.plot(t, W[:][1])
    # ax.plot(t, W[:][2])
    plt.show()
    return W

# T: Time of whole process
# n: Number of steps
# N: Number of processes
# XO: The initial value of the process
# c: deterministic constant
# d: probabilistic constant
def eulermaruyama(T, n, N, X0, c, d):
    h = T/n
    X = np.zeros((n, N))
    X[0, :] = X0
    W = wienerprocess(T,n,N)

    for i in range(n-1):
        X[i+1, :] = np.multiply(X[i,:], (1 + c*h + np.multiply(d, (W[i+1, :] - W[i,:]))))
    
    fig, ax = plt.subplots()
    t = np.arange(0, T, h)
    print(t)
    print("size of numpy array: ", t.size)
    for j in range(N):
        ax.plot(t, X[:, j])
    plt.show()

# histgauss(0, 1, 10000, 30)
# wienerprocess(20, 100, 100)
eulermaruyama(1, 1000, 100, 1, 1, 0.3)