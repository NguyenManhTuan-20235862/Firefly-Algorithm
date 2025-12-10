# ifa_wsn.py
import numpy as np
from init_ffa import init_ffa
from ifa_move import ifa_move
from coverage import coverage
from findlimits import findlimits

def ifa_wsn(u0, Lb, Ub, para, q):
    n = int(para[0])
    MaxGeneration = int(para[1])
    alpha0 = para[2]
    betamin = para[3]
    gamma = para[4]
    beta0 = 1
    b = beta0 - betamin
    NumEval = n * MaxGeneration
    
    alpha = alpha0
    theta = 0.97
    
    if len(Lb) != len(Ub):
        print('Simple bounds/limits are improper!')
        return None
    
    d = len(u0)
    zn = np.ones(n)
    
    nsx, nsy, Lightn = init_ffa(n, d, Lb, Ub, u0)
    
    maxzn = []
    
    for iter in range(MaxGeneration):
        print(f'IFA Iteration iter = {iter+1}')
        
        alpha = alpha0 * (theta ** (iter + 1))
        print(f'Current alpha = {alpha}')
        
        q = q + 1
        for i in range(n):
            Solution = np.vstack([nsx[i, :], nsy[i, :]])
            zn[i] = coverage(Solution, 100, 7)
            Lightn[i] = zn[i]
        
        maxzn.append(np.max(zn))
        print(f'IFA coverage of current solution: {maxzn[-1]}')
        
        nsxo = nsx.copy()
        nsyo = nsy.copy()
        Lighto = Lightn.copy()
        
        progress = (iter + 1) / MaxGeneration
        b_dynamic = b * (1 - 0.5 * progress)
        
        nsx, nsy = ifa_move(n, nsx, nsy, Lightn, nsxo, nsyo, Lighto, alpha, betamin, gamma, b_dynamic)
    
    Index = np.argsort(zn)
    nsx = nsx[Index, :]
    nsy = nsy[Index, :]
    Lightn = zn[Index]
    
    nxbest = nsx[-1, :]
    nybest = nsy[-1, :]
    Lightbest = Lightn[-1]
    fbest = Lightbest
    
    nsx, nsy = findlimits(n, nsx, nsy, Lb, Ub)
    
    return nxbest, nybest, fbest, NumEval, maxzn