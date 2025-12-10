# mofa_wsn.py
import numpy as np
from init_ffa import init_ffa
from mofa_move import mofa_move
from coverage import coverage
from findlimits import findlimits

def mofa_wsn(u0, Lb, Ub, para, q):
    n = int(para[0])
    MaxGeneration = int(para[1])
    alpha = para[2]
    betamin = para[3]
    gamma = para[4]
    beta0 = 1
    b = beta0 - betamin
    NumEval = n * MaxGeneration
    
    elite_size = max(3, int(np.floor(n * 0.1)))
    
    if len(Lb) != len(Ub):
        print('Simple bounds/limits are improper!')
        return None
    
    d = len(u0)
    zn = np.ones(n)
    
    nsx, nsy, Lightn = init_ffa(n, d, Lb, Ub, u0)
    
    best_solution_x = None
    best_solution_y = None
    best_fitness = 0
    
    maxzn = []
    
    for iter in range(MaxGeneration):
        print(f'MOFA Iteration iter = {iter+1}')
        
        gamma_adaptive = gamma * (1 + 0.5 * np.sin(np.pi * (iter + 1) / MaxGeneration))
        print(f'Current gamma = {gamma_adaptive}')
        
        q = q + 1
        for i in range(n):
            Solution = np.vstack([nsx[i, :], nsy[i, :]])
            zn[i] = coverage(Solution, 100, 7)
            Lightn[i] = zn[i]
        
        current_best_idx = np.argmax(zn)
        current_best = zn[current_best_idx]
        if current_best > best_fitness:
            best_fitness = current_best
            best_solution_x = nsx[current_best_idx, :].copy()
            best_solution_y = nsy[current_best_idx, :].copy()
        
        maxzn.append(np.max(zn))
        print(f'MOFA coverage of current solution: {maxzn[-1]}')
        
        nsxo = nsx.copy()
        nsyo = nsy.copy()
        Lighto = Lightn.copy()
        
        nsx, nsy = mofa_move(n, nsx, nsy, Lightn, nsxo, nsyo, Lighto, alpha, betamin, 
                             gamma_adaptive, b, best_solution_x, best_solution_y, best_fitness)
        
        sorted_idx = np.argsort(zn)[::-1]
        for e in range(min(elite_size, n)):
            idx = sorted_idx[e]
            current_fitness = coverage(np.vstack([nsx[idx, :], nsy[idx, :]]), 100, 7)
            if current_fitness < zn[idx]:
                nsx[idx, :] = nsxo[idx, :]
                nsy[idx, :] = nsyo[idx, :]
    
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