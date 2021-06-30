import numpy as np
import matplotlib.pyplot as plt
import calc_func as cf

area_mce_dict = {(0.01,0.025):0.0181, (0.01,0.047):0.0131, (0.01,0.069):0.00996,
                 (0.06,0.015):0.101, (0.06,0.035):0.0865, (0.06,0.055):0.0771,
                 (0.12,0.002):0.188, (0.12,0.021):0.168, (0.12,0.04):0.148}

A_N_dict = {(0.01,0.025):0.649, (0.01,0.047):0.41, (0.01,0.069):0.114,
            (0.06,0.015):0.668, (0.06,0.035):0.432, (0.06,0.055):0.156,
            (0.12,0.002):0.716, (0.12,0.021):0.458, (0.12,0.04):0.21}


def plot_analy_solution(la, ga, N, A_bar=0, L_bar=0, a_mce=0):
    
    n = (2*3**0.5*np.pi*N)**0.5*2

    n_in = 6
    k = 2*(np.tan(np.pi/n_in)/n_in)**0.5
    
    if A_bar > 0:
        if L_bar > 0:
            pass
        else:
            b = 3*ga/A_bar**0.5*(1-2*np.pi/n)
            c = -( A_bar + 6*ga*(1+2*np.pi/n) + 2*np.pi*la/n/k/A_bar**0.5 )
            A_bar = ( (b**2 - c)**0.5 - b)**2
            L_bar = 6*k*( (b**2 - c)**0.5 - b)
    else:
        A_bar, L_bar = cf.calc_ground_state(la, ga)
    
    a = np.linspace(0,1,1000)
    a_cri = A_bar/3 - n_in/3*ga - n_in*2/3*ga*np.pi/n
    mu = ( -k*a**1.5 + 3*k*a_cri*a**0.5 + ga*L_bar*(1-2*np.pi/n) ) * n /2/np.pi
    
    if a_mce > 0:
        mu_mce = ( -k*a_mce**1.5 + 3*k*a_cri*a_mce**0.5 + ga*L_bar*(1-2*np.pi/n) ) * n /2/np.pi
    
    if a_cri < 0:
        a_cri = 0
    
    fig = plt.figure(figsize=(5,3.5))
    ax = fig.add_subplot(111)

    ax.plot(mu, a, c='k')
    ax.plot(mu.max(), a_cri, 'ko')
    ax.plot(mu.max()*np.ones((100,)), np.linspace(0,a_cri,100), 'k--')
    ax.plot(np.linspace(0,mu.max(),100), a_cri*np.ones((100,)), 'k--')
    
    if a_mce > 0:
        ax.plot(mu_mce, a_mce, 'ro')
        ax.plot(np.linspace(0,mu_mce,100), a_mce*np.ones((100,)), 'r--')
        
    ax.set_xlim(la, mu.max()+0.2)
    ax.set_ylim(-0.02, A_bar+0.05)



# To compute analytical solutions (Fig 4) 

la = 0.12                         # Input Lambda
ga = 0.04                         # Input Gamma
N = 100                           # Input cluster size N 
A_N = 0; L_N = 0                  # Input area of normal cells: 0 to use ground state (for initial condition 1) 
# A_N = A_N_dict[(la,ga)]         # to use simulation result (for initial condition 2)
a_mce = area_mce_dict[(la,ga)]    # Input area of MCE

plot_analy_solution(la, ga, N, A_N, L_N, a_mce)

