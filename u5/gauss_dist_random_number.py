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

histgauss(0, 1, 10000, 30)