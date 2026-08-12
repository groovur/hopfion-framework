import mpmath as mp

mp.mp.dps = 60

phi = (1 + mp.sqrt(5)) / 2
pi = mp.pi

D1 = 16*phi**8/5 - 360/phi**2
D1_alt = (2136*phi - 3392)/5
C = 4*pi + mp.mpf(15)/8 - phi
delta = D1 - C

def rel_err(x, target):
    return abs(x - target) / abs(target)
