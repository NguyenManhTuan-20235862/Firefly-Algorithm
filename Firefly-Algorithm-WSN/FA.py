# FA.py
import numpy as np
import matplotlib.pyplot as plt
from ffa_wsn import ffa_wsn
from draw import draw

def main():
    w = 100
    d = 100
    point = d
    r = 7
    q = 0
    para = [25, 25, 0.7, 0.2, 1]
    
    Ub = np.ones(d) * w
    Lb = np.zeros(d)
    
    u0 = (Lb + Ub) / 2
    
    ux, uy, fval, NumEval, maxzn = ffa_wsn(u0, Lb, Ub, para, q)
    
    draw(ux, uy, 100, 7)
    
    print(f"Best solution x: {ux}")
    print(f"Best solution y: {uy}")
    print(f"Best objective: {fval}")
    print(f"Total evaluations: {NumEval}")
    
    plt.show()

if __name__ == "__main__":
    main()