# source_code/alpha_new.py
def alpha_new(alpha, NGen):
    delta = 1 - (0.005 / 0.9)**(1 / NGen)
    alpha = (1 - delta) * alpha
    return alpha