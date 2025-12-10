# source_code/fa_mincon.py
import numpy as np

def fa_mincon():
    para = [40, 150, 0.5, 0.2, 1]
    
    print('Solve the simple spring design problem ...')
    Lb = np.array([0.05, 0.25, 2.0])
    Ub = np.array([2.0, 1.3, 15.0])
    
    u0 = (Lb + Ub) / 2
    
    u, fval, NumEval = ffa_mincon(cost, constraint, u0, Lb, Ub, para)
    
    print(f"Best solution: {u}")
    print(f"Best objective: {fval}")
    print(f"Total evaluations: {NumEval}")
    
    return u, fval, NumEval

def cost(x):
    z = (2 + x[2]) * x[0]**2 * x[1]
    return z

def constraint(x):
    g = np.zeros(4)
    g[0] = 1 - x[1]**3 * x[2] / (71785 * x[0]**4)
    
    tmpf = (4 * x[1]**2 - x[0] * x[1]) / (12566 * (x[1] * x[0]**3 - x[0]**4))
    g[1] = tmpf + 1 / (5108 * x[0]**2) - 1
    g[2] = 1 - 140.45 * x[0] / (x[1]**2 * x[2])
    g[3] = x[0] + x[1] - 1.5
    
    geq = np.array([])
    
    return g, geq

def ffa_mincon(fhandle, nonhandle, u0, Lb, Ub, para):
    n = int(para[0])
    MaxGeneration = int(para[1])
    alpha = para[2]
    betamin = para[3]
    gamma = para[4]
    
    NumEval = n * MaxGeneration
    
    if len(Lb) != len(Ub):
        print('Simple bounds/limits are improper!')
        return None
    
    d = len(u0)
    
    zn = np.ones(n) * 1e100
    
    ns, Lightn = init_ffa_mincon(n, d, Lb, Ub, u0)
    
    for k in range(MaxGeneration):
        alpha = alpha_new_mincon(alpha, MaxGeneration)
        
        for i in range(n):
            zn[i] = Fun(fhandle, nonhandle, ns[i, :])
            Lightn[i] = zn[i]
        
        sorted_idx = np.argsort(zn)
        Lightn = zn[sorted_idx]
        ns = ns[sorted_idx, :]
        
        nso = ns.copy()
        Lighto = Lightn.copy()
        nbest = ns[0, :]
        Lightbest = Lightn[0]
        
        fbest = Lightbest
        
        ns = ffa_move_mincon(n, d, ns, Lightn, nso, Lighto, nbest, 
                             Lightbest, alpha, betamin, gamma, Lb, Ub)
    
    return nbest, fbest, NumEval

def init_ffa_mincon(n, d, Lb, Ub, u0):
    ns = np.zeros((n, d))
    
    if len(Lb) > 0:
        for i in range(n):
            ns[i, :] = Lb + (Ub - Lb) * np.random.rand(d)
    else:
        for i in range(n):
            ns[i, :] = u0 + np.random.randn(d)
    
    Lightn = np.ones(n) * 1e100
    
    return ns, Lightn

def ffa_move_mincon(n, d, ns, Lightn, nso, Lighto, nbest, Lightbest, 
                    alpha, betamin, gamma, Lb, Ub):
    scale = np.abs(Ub - Lb)
    
    for i in range(n):
        for j in range(n):
            r = np.sqrt(np.sum((ns[i, :] - ns[j, :])**2))
            
            if Lightn[i] > Lighto[j]:
                beta0 = 1
                beta = (beta0 - betamin) * np.exp(-gamma * r**2) + betamin
                tmpf = alpha * (np.random.rand(d) - 0.5) * scale
                ns[i, :] = ns[i, :] * (1 - beta) + nso[j, :] * beta + tmpf
    
    ns = findlimits_mincon(n, ns, Lb, Ub)
    
    return ns

def alpha_new_mincon(alpha, NGen):
    delta = 1 - (0.005 / 0.9)**(1 / NGen)
    alpha = (1 - delta) * alpha
    return alpha

def findlimits_mincon(n, ns, Lb, Ub):
    for i in range(n):
        ns_tmp = ns[i, :]
        I = ns_tmp < Lb
        ns_tmp[I] = Lb[I]
        
        J = ns_tmp > Ub
        ns_tmp[J] = Ub[J]
        
        ns[i, :] = ns_tmp
    
    return ns

def Fun(fhandle, nonhandle, u):
    z = fhandle(u)
    z = z + getnonlinear(nonhandle, u)
    return z

def getnonlinear(nonhandle, u):
    Z = 0
    lam = 1e15
    lameq = 1e15
    
    g, geq = nonhandle(u)
    
    for k in range(len(g)):
        Z = Z + lam * g[k]**2 * getH(g[k])
    
    for k in range(len(geq)):
        Z = Z + lameq * geq[k]**2 * geteqH(geq[k])
    
    return Z

def getH(g):
    if g <= 0:
        H = 0
    else:
        H = 1
    return H

def geteqH(g):
    if g == 0:
        H = 0
    else:
        H = 1
    return H

if __name__ == "__main__":
    fa_mincon()