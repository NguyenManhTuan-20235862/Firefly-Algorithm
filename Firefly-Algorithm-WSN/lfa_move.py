# lfa_move.py
import numpy as np
from coverage import coverage
from scipy.special import gamma

def levy_flight(beta):
    numerator = gamma(1 + beta) * np.sin(np.pi * beta / 2)
    denominator = gamma((1 + beta) / 2) * beta * 2**((beta - 1) / 2)
    sigma = (numerator / denominator)**(1 / beta)
    
    u = np.random.randn() * sigma
    v = np.random.randn()
    
    step = u / abs(v)**(1 / beta)
    
    MAX_STEP_SIZE = 100 * 0.1
    
    if abs(step) > MAX_STEP_SIZE:
        step = np.sign(step) * MAX_STEP_SIZE
    
    return step

def lfa_move(n, nsx, nsy, Lightn, nsxo, nsyo, Lighto, alpha, betamin, gamma_val, b, levy_beta):
    for i in range(n):
        for j in range(n):
            k = np.random.randint(0, 100)
            l = np.random.randint(0, 100)
            
            rx = abs(nsx[i, k] - nsx[j, l])
            ry = abs(nsy[i, k] - nsy[j, l])
            
            rx1 = abs(nsx[i, l] - nsx[j, k])
            ry1 = abs(nsy[i, l] - nsy[j, k])
            
            if Lightn[i] < Lighto[j]:
                betax = b * np.exp(-gamma_val * rx**2) + betamin
                betay = b * np.exp(-gamma_val * ry**2) + betamin
                
                betax1 = b * np.exp(-gamma_val * rx1**2) + betamin
                betay1 = b * np.exp(-gamma_val * ry1**2) + betamin
                
                levy_step_x = levy_flight(levy_beta)
                levy_step_y = levy_flight(levy_beta)
                
                tmpf_x = alpha * levy_step_x
                tmpf_y = alpha * levy_step_y
                
                nsx[i, k] = nsx[i, k] * (1 - betax) + nsxo[j, l] * betax + tmpf_x
                nsy[i, k] = nsy[i, k] * (1 - betay) + nsyo[j, l] * betay + tmpf_y
                
                levy_step_x2 = levy_flight(levy_beta)
                levy_step_y2 = levy_flight(levy_beta)
                tmpf_x2 = alpha * levy_step_x2
                tmpf_y2 = alpha * levy_step_y2
                
                nsx[i, l] = nsx[i, l] * (1 - betax1) + nsxo[j, k] * betax1 + tmpf_x2
                nsy[i, l] = nsy[i, l] * (1 - betay1) + nsyo[j, k] * betay1 + tmpf_y2
                
                Solution_temp = np.vstack([nsx[i, :], nsy[i, :]])
                Lightn_temp = coverage(Solution_temp, 100, 7)
                
                if Lightn_temp < Lightn[i]:
                    nsx[i, :] = nsxo[i, :]
                    nsy[i, :] = nsyo[i, :]
                else:
                    Lightn[i] = Lightn_temp
                    Lighto[i] = Lightn[i]
    
    return nsx, nsy