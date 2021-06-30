import numpy as np
import matplotlib.pyplot as plt
import calc_func as cf

param_set = [(0.01,0.025),(0.01,0.047),(0.01,0.069),
             (0.06,0.015),(0.06,0.035),(0.06,0.055),
             (0.12,0.002),(0.12,0.021),(0.12,0.04)]
    
def plot_parameter_diagram(la_lim=(0.001,0.2), ga_lim=(0,0.1), flag='bs', param_set=[(0.12,0.04)]):
    
    fig = plt.figure(figsize=(3.5,3))
    ax = fig.add_subplot(111)
    ax.set_xlim(la_lim)
    ax.set_ylim(ga_lim)

    la_list = np.linspace(la_lim[0],la_lim[1],1000)
    ga_list = np.linspace(ga_lim[0],ga_lim[1],1000)
    la_mesh,ga_mesh = np.meshgrid(la_list,ga_list)
    L_mesh = np.array([cf.calc_ground_state(la,ga) for la,ga in np.c_[la_mesh.ravel(),ga_mesh.ravel()]])[:,1].reshape(la_mesh.shape)
    L_mesh = np.where(L_mesh<0,0,L_mesh)
    
    if 'b' in flag:
        LL_mesh = np.array([cf.calc_bulk(la,ga,L) for la,ga,L in np.c_[la_mesh.ravel(),ga_mesh.ravel(),L_mesh.ravel()]]).reshape(la_mesh.shape)    
        z = np.ma.masked_array(LL_mesh,mask=LL_mesh<-0.1)
        lv = np.arange(0.2,10,0.6)
        cs = ax.contour(la_list,ga_list,z,levels=lv, colors='k',linewidths=1)
        ax.clabel(cs, fmt='%1.1f')
        
    if 's' in flag:
        LL_mesh = np.array([cf.calc_shear(la,ga,L) for la,ga,L in np.c_[la_mesh.ravel(),ga_mesh.ravel(),L_mesh.ravel()]]).reshape(la_mesh.shape)
        z = np.ma.masked_array(LL_mesh,mask=LL_mesh<-0.1)
        lv = np.arange(0,5,0.5)
        cs = ax.contour(la_list,ga_list,z,levels=lv, colors='k',linewidths=2, alpha=0.2)
        ax.clabel(cs, fmt='%1.1f')
        
    ax.contour(la_list,ga_list,L_mesh>0, colors='k',linewidths=1)
    
    for param in param_set:
        ax.scatter(param[0], param[1], edgecolors='none', facecolors='r', s=50)


# To plot Lambda-Gamma parameter space (Fig 1)
plot_parameter_diagram(flag='bs', param_set=param_set)
