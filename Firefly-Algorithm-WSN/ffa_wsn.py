# ffa_wsn.py
import numpy as np
from init_ffa import init_ffa
from ffa_move import ffa_move
from coverage import coverage
from findlimits import findlimits

def ffa_wsn(u0, Lb, Ub, para, q):
    n = int(para[0])
    MaxGeneration = int(para[1])
    alpha = para[2]
    betamin = para[3]
    gamma = para[4]
    beta0 = 1
    b = beta0 - betamin
    NumEval = n * MaxGeneration
    
    if len(Lb) != len(Ub):
        print('Simple bounds/limits are improper!')
        return None
    
    d = len(u0)
    zn = np.ones(n)
    
    nsx, nsy, Lightn = init_ffa(n, d, Lb, Ub, u0)
    
    maxzn = []
    
    for iter in range(MaxGeneration):
        print(f'Iteration iter = {iter+1}')
        
        q = q + 1
        for i in range(n):
            Solution = np.vstack([nsx[i, :], nsy[i, :]])
            zn[i] = coverage(Solution, 100, 7)
            Lightn[i] = zn[i]
        
        maxzn.append(np.max(zn))
        print(f'coverage of current solution: {maxzn[-1]}')
        
        nsxo = nsx.copy()
        nsyo = nsy.copy()
        Lighto = Lightn.copy()
        
        nsx, nsy = ffa_move(n, nsx, nsy, Lightn, nsxo, nsyo, Lighto, alpha, betamin, gamma, b)
    
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