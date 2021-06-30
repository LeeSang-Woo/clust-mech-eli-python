import csv
import numpy as np
import matplotlib.pyplot as plt

data_path = "./init_data/"


def make_cell_stat(param_set, data_path=data_path, fig='off', cells_marker=-1,
                   cell_num_flag=False, save_flag=False):
    
    vertex_list, cell_list = read_init_cond(param_set, data_path)
    if fig == 'on':
        plot_init_cond(vertex_list, cell_list, cells_marker,
                          cell_num_flag, save_flag)


def read_init_cond(param_set, data_path):

    la = param_set[0]
    ga = param_set[1]    
    cell_num = param_set[2]
    
    str1 = 'in' + str(cell_num)
    str2 = 'L' + '{:04d}'.format(int(la*1000)) + 'G' + '{:04d}'.format(int(ga*1000))
    
    fn_c = data_path+'/initial2D_RD05_'+str1+'_c_'+str2+'.csv'
    fn_p = data_path+'/initial2D_RD05_'+str1+'_p_'+str2+'.csv'
    
    vertex_list = []
    with open(fn_c, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter = ',')
        for row in reader:
            vertex_list.append([float(x) for x in row])
    
    cell_list = []
    with open(fn_p, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter = ',')
        for row in reader:
            cell_list.append([int(x) for x in row])
            
    return vertex_list, cell_list


def plot_init_cond(vertex_list, cell_list, cells_mark=-1,
                   cell_num_flag=False, save_flag=False):
    
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)
    ax.axis('off')
    
    for cell in cell_list:
        polys = plt.Polygon([vertex_list[x-1] for x in cell[:-1]], fill=False, edgecolor='black')
        ax.add_patch(polys)
    
    if type(cells_mark) is int:
        if cells_mark > -1:
            target_polys = plt.Polygon([vertex_list[x-1] for x in cell_list[cells_mark][:-1]], color='r', alpha=0.5)
            ax.add_patch(target_polys)
    elif type(cells_mark) is list:
        if len(cells_mark)>0:
            for cell in cells_mark:
                target_polys = plt.Polygon([vertex_list[x-1] for x in cell_list[cell][:-1]], color='r', alpha=0.5)
                ax.add_patch(target_polys)
    
    plt.xlim((np.array(vertex_list)[:,0].min(),np.array(vertex_list)[:,0].max()))
    plt.ylim((np.array(vertex_list)[:,1].min(),np.array(vertex_list)[:,1].max()))
    ax.set_aspect('equal')
    
    if cell_num_flag:
        index = 0
        fs = 10-len(cell_list)//200
        for cell in cell_list:
            center = np.array([vertex_list[x-1] for x in cell[:-1]]).mean(axis=0)
            ax.text(center[0]-0.2, center[1], str(index), color='black', fontsize=fs)
            index += 1
    
    if save_flag:
        fig.savefig("./init_cond.pdf", bbox_inches='tight')



# To plot initial condition (Fig 1)
la = 0.12
ga = 0.04
cell_num = 1000
make_cell_stat((la,ga,cell_num), data_path, 'on', -1, False, False)


    

