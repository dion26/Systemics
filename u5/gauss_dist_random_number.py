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
        W[0,i] = 0
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
    return X

# X0: Initial Value
# T0: Initial Time
# XT: Put value
# i: Deterministic SDE parameter
# r: Probabilistic SDE parameter
# T: terminal time

def put(X0, T0, XT, i, r, T):
    d1 = (math.log(X0/XT) + (i + 0.5 * r**2) * (T - T0)) / (r * math.sqrt(T - T0))
    d2 = d1 - r * math.sqrt(T - T0)
    n1 = 0.5 * (1 + math.erf(-d1 / math.sqrt(2)))
    n2 = 0.5 * (1 + math.erf(-d2 / math.sqrt(2)))
    result = XT * math.exp(-i * (T - T0)) * n2 - X0 * n1
    return result

def montecarlo():
    T=1
    n=100
    h=T/n
    N=1000
    X0=80
    XT=100
    i=0.08
    r=0.2

    Y = eulermaruyama(T, n, N, X0, i, r)

    ### Plot Euler-Maruyama solution paths
    time_grid = np.arange(0,T,h)

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)
    fig.suptitle('European Put Option')

    ax1.plot(time_grid, Y)
    ax1.set_title("Paths of Euler-Maruyama solutions")

    histogram = np.sort(Y[n-1, :])
    ax2.hist(histogram, bins=N)
    ax2.grid(True)

    # Evaluate payoff function 
    payoff = np.maximum(0, XT - Y[n-1, :])
    
    # Compute discounted expected value
    V0 = math.exp(-i*T) * (np.cumsum(payoff) / np.arange(0, N))

    # Plot discounted expected value
    ax3.plot(V0)
    ax3.set_title('Discounted expected value')
    ax3.grid(True)

    # Compute error in expected value
    Vexact = put(X0, 0, XT, i, r, T)
    ax4.plot(abs(V0 - Vexact * np.ones(np.size(V0))) / Vexact)
    ax4.set_title('Error in expected value')
    ax4.grid(True)

    plt.show()

# histgauss(0, 1, 10000, 30)
# wienerprocess(20, 100, 100)
# eulermaruyama(1, 1000, 100, 1, 1, 0.3)
montecarlo()