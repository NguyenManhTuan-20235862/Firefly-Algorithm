# init_ffa.py
import numpy as np

def init_ffa(n, d, Lb, Ub, u0):
    Range = Ub - Lb
    Lower = Lb
    
    if len(Lb) > 0:
        nsx = np.zeros((n, d))
        nsy = np.zeros((n, d))
        for i in range(n):
            nsx[i, :] = Lower + Range * np.random.rand(d)
            nsy[i, :] = Lower + Range * np.random.rand(d)
    else:
        nsx = np.zeros((n, d))
        nsy = np.zeros((n, d))
        for i in range(n):
            nsx[i, :] = u0 + np.random.randn(d)
            nsy[i, :] = u0 + np.random.randn(d)
    
    Lightn = np.ones(n)
    
    return nsx, nsy, Lightn