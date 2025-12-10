# source_code/firefly_simple.py
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def firefly_simple(instr=None):
    if instr is None:
        instr = [12, 50]
    
    n = int(instr[0])
    MaxGeneration = int(instr[1])
    
    np.random.seed(0)
    
    def f(x, y):
        return (np.exp(-(x-4)**2-(y-4)**2) + np.exp(-(x+4)**2-(y-4)**2) + 
                2*np.exp(-x**2-(y+4)**2) + 2*np.exp(-x**2-y**2))
    
    range_vals = [-5, 5, -5, 5]
    
    alpha = 0.2
    gamma = 1.0
    
    Ngrid = 100
    dx = (range_vals[1] - range_vals[0]) / Ngrid
    dy = (range_vals[3] - range_vals[2]) / Ngrid
    x = np.arange(range_vals[0], range_vals[1], dx)
    y = np.arange(range_vals[2], range_vals[3], dy)
    X, Y = np.meshgrid(x, y)
    Z = f(X, Y)
    
    fig1 = plt.figure(1)
    ax1 = fig1.add_subplot(111, projection='3d')
    ax1.plot_surface(X, Y, Z, cmap='viridis', alpha=0.8)
    plt.title('Objective Function')
    
    xn, yn, Lightn = init_ffa_simple(n, range_vals)
    
    fig2 = plt.figure(2)
    
    for i in range(MaxGeneration):
        plt.clf()
        plt.contour(X, Y, Z, 15)
        
        zn = f(xn, yn)
        
        sorted_idx = np.argsort(zn)
        Lightn = zn[sorted_idx]
        xn = xn[sorted_idx]
        yn = yn[sorted_idx]
        xo = xn.copy()
        yo = yn.copy()
        Lighto = Lightn.copy()
        
        plt.plot(xn, yn, 'g.', markersize=10)
        
        xn, yn = ffa_move_simple(xn, yn, Lightn, xo, yo, Lighto, alpha, gamma, range_vals)
        
        plt.pause(0.1)
    
    best = np.column_stack([xo, yo, Lighto])
    
    plt.show()
    
    return best

def init_ffa_simple(n, range_vals):
    xrange = range_vals[1] - range_vals[0]
    yrange = range_vals[3] - range_vals[2]
    xn = np.random.rand(n) * xrange + range_vals[0]
    yn = np.random.rand(n) * yrange + range_vals[2]
    Lightn = np.zeros(n)
    return xn, yn, Lightn

def ffa_move_simple(xn, yn, Lightn, xo, yo, Lighto, alpha, gamma, range_vals):
    ni = len(yn)
    nj = len(yo)
    
    for i in range(ni):
        for j in range(nj):
            r = np.sqrt((xn[i] - xo[j])**2 + (yn[i] - yo[j])**2)
            
            if Lightn[i] < Lighto[j]:
                beta0 = 1
                beta = beta0 * np.exp(-gamma * r**2)
                xn[i] = xn[i] * (1 - beta) + xo[j] * beta + alpha * (np.random.rand() - 0.5)
                yn[i] = yn[i] * (1 - beta) + yo[j] * beta + alpha * (np.random.rand() - 0.5)
    
    xn, yn = findrange_simple(xn, yn, range_vals)
    
    return xn, yn

def findrange_simple(xn, yn, range_vals):
    for i in range(len(yn)):
        if xn[i] <= range_vals[0]:
            xn[i] = range_vals[0]
        if xn[i] >= range_vals[1]:
            xn[i] = range_vals[1]
        if yn[i] <= range_vals[2]:
            yn[i] = range_vals[2]
        if yn[i] >= range_vals[3]:
            yn[i] = range_vals[3]
    
    return xn, yn

if __name__ == "__main__":
    best = firefly_simple([12, 50])
    print("Best solutions:")
    print(best)