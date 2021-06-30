import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import calc_func as cf

area_mce_dict = {(0.01,0.025):0.0181, (0.01,0.047):0.0131, (0.01,0.069):0.00996,
                 (0.06,0.015):0.101, (0.06,0.035):0.0865, (0.06,0.055):0.0771,
                 (0.12,0.002):0.188, (0.12,0.021):0.168, (0.12,0.04):0.148}

area_poly_exist_dict = {(0.01,0.025):(0.0862,0.00581,0.00061),
                        (0.01,0.047):(0.0244,0.00171,0.00021),
                        (0.01,0.069):(0.0114,0.00081,0.00011),
                        (0.06,0.015):(0.396,0.226,0.0524),
                        (0.06,0.035):(0.523,0.106,0.00971),
                        (0.06,0.055):(0.6,0.043,0.00391),
                        (0.12,0.002):(0.465,0.327,0.196),
                        (0.12,0.021):(0.609,0.357,0.107),
                        (0.12,0.04):(0.738,0.324,0.0295)}


def make_polygon_exist_area(la, ga, n):
    
    A_list = np.linspace(0.00001,1,10000)
    a_list = []
    for A_bar in A_list:
        L_bar = 2*(2*3**0.5*A_bar)**0.5
        params = [la, ga, A_bar, L_bar]
        a_list.append(cf.calc_area_1(n,params))
    a_list = np.array(a_list)    
    
    return A_list[np.where(a_list>0)[0][0]]


# To make area_poly_exist_dict dataset
la = 0.12
ga = 0.04
n = 3
make_polygon_exist_area(la, ga, n)



def plot_area_mce_poly_exist(area_mce_dict, area_poly_exist_dict):

    fig = plt.figure(figsize=(6,4))
    ax = fig.add_subplot(111)
    param_list = list(area_mce_dict.keys())
    
    for index, param in zip(range(len(param_list)), param_list):
        A, L = cf.calc_ground_state(param[0], param[1])
        x = np.linspace(0,0.8,100)
        y = np.ones_like(x)*index
        ax.plot(x,y,'k-',lw=1)
        ax.plot(area_mce_dict[param], index, 'r|')
        ax.plot(area_poly_exist_dict[param][0], index, 'b|')
        ax.plot(area_poly_exist_dict[param][1], index, 'g|')
        ax.plot(area_poly_exist_dict[param][2], index, 'k|')
    
    ax.set_yticks(range(len(param_list)))
    ax.set_yticklabels(param_list)

    
def plot_area_corr(area_mce_dict, area_poly_exist_dict, flag='reg'):
    
    fig = plt.figure(figsize=(4,4))
    ax = fig.add_subplot(111)
    
    param_list = list(area_poly_exist_dict.keys())
    
    X = []
    Y = []
    
    for index, param in zip(range(len(param_list)), param_list):
        if flag == '3':
            ax.plot(area_mce_dict[param], area_poly_exist_dict[param][0], 'bo', alpha=0.4)
        if flag == '4':
            ax.plot(area_mce_dict[param], area_poly_exist_dict[param][1], 'go', alpha=0.4)
        if flag == '5':
            ax.plot(area_mce_dict[param], area_poly_exist_dict[param][2], 'ko', alpha=0.4)
        if flag == 'reg':
            X.append([area_poly_exist_dict[param][0], area_poly_exist_dict[param][1], area_poly_exist_dict[param][2]])
            Y.append(area_mce_dict[param])

    if flag == 'reg':
        X = np.array(X)
        Y = np.array(Y)
        mlr = LinearRegression()
        mlr.fit(X, Y)
        Y_ = mlr.intercept_ + mlr.coef_[0]*X[:,0] + mlr.coef_[1]*X[:,1] + mlr.coef_[2]*X[:,2]
        ax.plot(Y,Y_,'ko', alpha=0.4)
        ax.plot(np.linspace(-0.01,0.22,10),np.linspace(-0.01,0.22,10),'k--')
        ax.set_xlim(-0.01,0.22)
        ax.set_ylim(-0.01,0.22)
        


# To plot about area of mce and existance condition of the polygon (Fig 3)
plot_area_mce_poly_exist(area_mce_dict, area_poly_exist_dict)
plot_area_corr(area_mce_dict, area_poly_exist_dict, flag='reg')
