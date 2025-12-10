# findlimits.py
import numpy as np

def findlimits(n, nsx, nsy, Lb, Ub):
    for i in range(n):
        ns_tmp = nsx[i, :]
        I = ns_tmp < Lb
        ns_tmp[I] = Lb[I]
        
        J = ns_tmp > Ub
        ns_tmp[J] = Ub[J]
        nsx[i, :] = ns_tmp
        
        ns_tmp = nsy[i, :]
        I = ns_tmp < Lb
        ns_tmp[I] = Lb[I]
        
        J = ns_tmp > Ub
        ns_tmp[J] = Ub[J]
        nsy[i, :] = ns_tmp
    
    return nsx, nsy