import numpy as np


def calc_cubic(p,q):
    
    D = - p**3/27 - q**2/4
    
    if D<=0:
        x = np.sign(-q/2-(-D)**0.5)*(np.abs(-q/2-(-D)**0.5))**(1/3)\
          + np.sign(-q/2+(-D)**0.5)*(np.abs(-q/2+(-D)**0.5))**(1/3)
    else:
        phi = np.arccos( -q/2*(-3/p)**(1.5) )
        x = np.array([ 2*(-p/3)**0.5*np.cos((phi+2*np.pi*k)/3) for k in [0,1,2] ]).max()
        
    return x


def calc_ground_state(la, ga):
    
    if ga <= -la/(4*(2*3**0.5)**0.5):
        A_bar = 1
        L_bar = (24/3**0.5)**0.5
    else:
        r0 = calc_cubic(2/9*(12*ga-3**0.5), 2/9*la)
        A_bar = 3**0.5*3*0.5*r0**2
        L_bar = r0*6
        
    return A_bar, L_bar


def calc_param(n, params):
    
    la = params[0]
    ga = params[1]
    A_bar = params[2]
    L_bar = params[3]
    
    sin_n = np.sin(np.pi/n)
    cos_n = np.cos(np.pi/n)
    tan_n = np.tan(np.pi/n)
    
    p1 = (A_bar - 2*ga*n*tan_n)/3
    p2 = ( la*(2*sin_n-1) + 2*ga*L_bar*(sin_n-1) )/2/(4*sin_n*cos_n/n)**0.5
    
    return p1, p2


def calc_area_1(n, params):
    
    p1, p2 = calc_param(n,params)        
    sqrt_a = calc_cubic(-3*p1, 2*p2)
    
    if sqrt_a < 0:
        a = 0
    else:
        a = sqrt_a**2
    
    return a


def calc_shear(la,ga,L):
    
    if L>0:
        return 12*3**0.5*ga + 3**0.5*la/(L/6)
    else:
        return -1


def calc_bulk(la,ga,L):

    if L>0:        
        return 9*3**0.5*(L/6)**2 + 8*3**0.5*ga - 2
    else:
        return -1



