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

def wienerprocess(T, n, N):
    h = T / n
    mu = 0
    sigma = math.sqrt(h)

    W = np.zeros((n, N))
    for i in range(N):
        W[0, i] = 0
        for j in range(1, n):
            d_W = gaussvar(mu, sigma)
            W[j,i] = W[j-1, i] + d_W

    t = range (n)
    plt.plot(t, W[:1])
    plt.plot(t, W[:2])
    plt.plot(t, W[:3])
    plt.show()


histgauss(0, 1, 10000, 30)
wienerprocess(20, 10, 10)