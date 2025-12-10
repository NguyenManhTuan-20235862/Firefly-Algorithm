# ffa_move.py
import numpy as np
from coverage import coverage

def ffa_move(n, nsx, nsy, Lightn, nsxo, nsyo, Lighto, alpha, betamin, gamma, b):
    for i in range(n):
        for j in range(n):
            k = np.random.randint(0, 100)
            l = np.random.randint(0, 100)
            
            rx = abs(nsx[i, k] - nsx[j, l])
            ry = abs(nsy[i, k] - nsy[j, l])
            
            rx1 = abs(nsx[i, l] - nsx[j, k])
            ry1 = abs(nsy[i, l] - nsy[j, k])
            
            if Lightn[i] < Lighto[j]:
                betax = b * np.exp(-gamma * rx**2) + betamin
                betay = b * np.exp(-gamma * ry**2) + betamin
                
                betax1 = b * np.exp(-gamma * rx1**2) + betamin
                betay1 = b * np.exp(-gamma * ry1**2) + betamin
                
                tmpf = alpha * (np.random.rand() - 0.5)
                
                nsx[i, k] = nsx[i, k] * (1 - betax) + nsxo[j, l] * betax + tmpf
                nsy[i, k] = nsy[i, k] * (1 - betay) + nsyo[j, l] * betay + tmpf
                
                nsx[i, l] = nsx[i, l] * (1 - betax1) + nsxo[j, k] * betax + tmpf
                nsy[i, l] = nsy[i, l] * (1 - betay1) + nsyo[j, k] * betay + tmpf
                
                Solution_temp = np.vstack([nsx[i, :], nsy[i, :]])
                Lightn_temp = coverage(Solution_temp, 100, 7)
                
                if Lightn_temp < Lightn[i]:
                    nsx[i, :] = nsxo[i, :]
                    nsy[i, :] = nsyo[i, :]
                else:
                    Lightn[i] = Lightn_temp
                    Lighto[i] = Lightn[i]
    
    return nsx, nsy