import numpy as np
import matplotlib.pyplot as plt
from lfa_wsn import lfa_wsn
from draw import draw

def main():
    w = 100
    d = 100 
    r = 7
    q = 0
    # N=25, MaxGen=25, alpha=0.7, betamin=0.2, gamma=1
    para = [25, 25, 0.7, 0.2, 1] 
    
    Ub = np.ones(d) * w
    Lb = np.zeros(d)
    u0 = (Lb + Ub) / 2
    
    print("--- Running LFA-WSN (Levy Flight Firefly Algorithm) ---")
    ux, uy, fval, NumEval, maxzn = lfa_wsn(u0, Lb, Ub, para, q)
    
    draw(ux, uy, 100, 7)
    
    print(f"\n--- LFA Results ---")
    print(f"LFA Best objective (fbest): {fval}")
    print(f"LFA Total evaluations: {NumEval}")
    print(f"LFA Convergence History (maxzn): {maxzn}") # Dùng để vẽ biểu đồ
    
    plt.show()

if __name__ == "__main__":
    main()