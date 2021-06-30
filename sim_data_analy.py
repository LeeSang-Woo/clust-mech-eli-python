import numpy as np
import matplotlib.pyplot as plt
import csv
import os
import itertools
from scipy import optimize

data_path = "./sim_data/"

file_index = {"L001G0025_IC1":["InitCond1",0.01,0.025,np.linspace(120,400,29),[100,150,200,250]],
              "L001G0047_IC1":["InitCond1",0.01,0.047,np.linspace(120,230,12),[100]],
              "L001G0069_IC1":["InitCond1",0.01,0.069,np.linspace(80,220,15),[100]],
              "L006G0015_IC1":["InitCond1",0.06,0.015,np.linspace(24,44,11),[100]],
              "L006G0035_IC1":["InitCond1",0.06,0.035,np.linspace(16,40,13),[100]],
              "L006G0055_IC1":["InitCond1",0.06,0.055,np.linspace(8,34,14),[100]],
              "L012G0002_IC1":["InitCond1",0.12,0.002,np.linspace(13,23,11),[100]],
              "L012G0021_IC1":["InitCond1",0.12,0.021,np.linspace(8,20,13),[100]],
              "L012G0040_IC1":["InitCond1",0.12,0.04,np.linspace(2,20,19),[100,150,200,250]],
              "L001G0025_IC2":["InitCond2",0.01,0.025,np.linspace(150,210,13),[100]],
              "L001G0047_IC2":["InitCond2",0.01,0.047,np.linspace(100,150,11),[100]],
              "L001G0069_IC2":["InitCond2",0.01,0.069,np.linspace(3,69,23),[100]],
              "L006G0015_IC2":["InitCond2",0.06,0.015,np.linspace(18,40,12),[100]],
              "L006G0035_IC2":["InitCond2",0.06,0.035,np.linspace(8,30,12),[100]],
              "L006G0055_IC2":["InitCond2",0.06,0.055,np.linspace(2,15,14),[100]],
              "L012G0002_IC2":["InitCond2",0.12,0.002,np.linspace(12,18,13),[100]],
              "L012G0021_IC2":["InitCond2",0.12,0.021,np.linspace(5,11,13),[100]],
              "L012G0040_IC2":["InitCond2",0.12,0.04,np.linspace(1.5,7,12),[100]],
              "IC1_0010r":["IC1 10r",0.12,0.04,np.linspace(5,15,11),[100]],
              "IC1_0100r":["IC1 100r",0.12,0.04,np.linspace(5,15,11),[100]],
              "IC1_1000r":["IC1 1000r",0.12,0.04,np.linspace(5,15,11),[100]],
              "IC2_0010r":["IC2 10r",0.12,0.04,np.linspace(1.5,7,12),[100]],
              "IC2_0100r":["IC2 100r",0.12,0.04,np.linspace(1.5,7,12),[100]],
              "IC2_1000r":["IC2 1000r",0.12,0.04,np.linspace(1.5,7,12),[100]],
              "AB_G001":["Ab Gamma 0.01",0.12,0.04,np.linspace(1.5,13,24),[100]],
              "AB_G004":["Ab Gamma 0.04",0.12,0.01,np.linspace(2,8,7),[100]] }

ckey = {'timestep':0, 'recon_num':3, 'ab_cells':6, 'ab_out_cells':7, 'ab_in_cells':8, 'edge_num':10,
        'all_area':12, 'mean_all_area':13, 'mean_out_area':14, 'mean_in_area':15,
        'poly_freq_all':17, 'poly_area_all':24, 'poly_freq_out':31, 'poly_area_out':38,
        'poly_freq_in':45, 'poly_area_in':52, 'nei_cells':59, 'mean_nei_area':60,
        'mean_sq_nei_area':61, 'mean_nei_st':62, 'mean_all_st':64, 'mean_out_st':65,
        'mean_in_st':66, 'mean_sq_all_area':67, 'mean_sq_out_area':68, 'mean_sq_in_area':69,
        'mean_all_peri':71, 'mean_all_side':72, 'mean_out_peri':73, 'mean_out_side':74,
        'mean_in_peri':75, 'mean_in_side':76, 'mean_nei_peri':77, 'mean_nei_side':78}


def read_file(folder_name, dn, pop_size, onset_mu):
    
    file_name = folder_name + "%02d/data2D_N%03dmu%03d.csv" % (dn, pop_size, int(onset_mu*100))
    data_list = []
    if os.path.exists(file_name):
        with open(file_name, 'r') as csvfile:
            reader = csv.reader(csvfile, delimiter = ',')
            for row in reader:
                row_data = [int(row[ckey['timestep']]), #0
                            int(row[ckey['recon_num']]), #1
                            int(row[ckey['ab_cells']]), #2
                            int(row[ckey['ab_out_cells']]), #3
                            int(row[ckey['ab_in_cells']]), #4
                            float(row[ckey['all_area']]), #5
                            float(row[ckey['mean_all_area']]), #6
                            float(row[ckey['mean_out_area']]), #7
                            float(row[ckey['mean_in_area']]) #8
                            ]
                row_data += [float(x) for x in row[ckey['poly_freq_all']:ckey['poly_freq_all']+6]] #9-14
                row_data += [float(x) for x in row[ckey['poly_area_all']:ckey['poly_area_all']+6]] #15-20
                row_data += [float(x) for x in row[ckey['poly_freq_out']:ckey['poly_freq_out']+6]] #21-26
                row_data += [float(x) for x in row[ckey['poly_area_out']:ckey['poly_area_out']+6]] #27-32
                row_data += [float(x) for x in row[ckey['poly_freq_in']:ckey['poly_freq_in']+6]] #33-38
                row_data += [float(x) for x in row[ckey['poly_area_in']:ckey['poly_area_in']+6]] #39-44
                
                row_data += [float(row[ckey['mean_all_st']]), #45
                             float(row[ckey['mean_out_st']]), #46
                             float(row[ckey['mean_in_st']]), #47
                             int(row[ckey['nei_cells']]), #48
                             float(row[ckey['mean_nei_area']]), #49
                             float(row[ckey['mean_nei_st']]), #50
                             float(row[ckey['mean_sq_nei_area']]), #51
                             float(row[ckey['mean_sq_all_area']]), #52
                             float(row[ckey['mean_sq_out_area']]), #53
                             float(row[ckey['mean_sq_in_area']]) #54
                             ]
                row_data += [int(row[ckey['edge_num']])]    #55
                
                if len(row) > ckey['mean_nei_peri']:                    
                    row_data += [float(row[ckey['mean_nei_peri']])]     #56
                    row_data += [float(row[ckey['mean_nei_side']])]     #57
                
                data_list.append(row_data)
    
    return np.array(data_list)


def plot_series(delta_t, data_dict, dn, pop_size, onset_mu, ylim=0, flag=0):
    
    fig = plt.figure(figsize=(3,3))
    ax = fig.add_subplot(111)
    
    if flag == 0:
        y_list = data_dict[(dn, pop_size, onset_mu)][:,2]        
    elif flag == 1:
        y_list = data_dict[(dn, pop_size, onset_mu)][:,5]
    elif flag == 2:
        y_list = data_dict[(dn, pop_size, onset_mu)][:,6]
    x = data_dict[(dn, pop_size, onset_mu)][:,0]*delta_t*0.2
    ax.plot(x, y_list, 'k', lw=1)
    ax.set_xlim(0,x.max()+1)
    x = np.linspace(data_dict[(dn, pop_size, onset_mu)][:,0].min()*delta_t,
                    data_dict[(dn, pop_size, onset_mu)][:,0].max()*delta_t,
                    100)
    
    if ylim > 0:
        ax.set_ylim(-10, ylim)    


def plot_lewis_law_1_ab(delta_t, data_dict, dn, pop_size, onset_mu):
    
    n_range = np.array([4,5,6,7,8])    
    cell_freq_normal = np.array([0.04804805,0.28328328,0.38638639,0.21021021,0.06406406])
    cell_freq_fixed = np.array([0.009161891,0.389046012,0.366121099,0.178352162,0.04781197])
    data_list = data_dict[(dn, pop_size, onset_mu)]
    onset_time = np.where(data_list[:,2]>=pop_size)[0][0]
    time_data = data_list[onset_time:,0]
    cell_freq_abnormal = np.array([0,0,0,0,0], dtype=float)
    for n in n_range:
        cell_freq_abnormal[n-4] = calc_time_ave(time_data,data_list[onset_time:,n+29])
    
    fig = plt.figure(figsize=(4,3))
    ax = fig.add_subplot(111)
    bar_width = 0.3
    ax.bar(n_range, cell_freq_normal, bar_width, color='w', edgecolor='k')
    ax.bar(n_range+bar_width, cell_freq_abnormal, bar_width, color='k', edgecolor='k')
    ax.bar(n_range+bar_width*2, cell_freq_fixed, bar_width, color='gray', edgecolor='gray')
    
    ax.set_xticks(n_range+bar_width)
    ax.set_xticklabels(n_range)
    ax.set_ylim((0,0.45))
    ax.tick_params(direction='in', top=False, right=True, bottom=False)


def calc_time_ave(time_data, target_data):
    return (np.diff(time_data)*target_data[:-1]).sum() / (time_data[-1]-time_data[0])


def make_type_dict(delta_t, data_dict, dn_list, pop_size_list, onset_mu_list, sim_end_list):        
    
    type_dict = {}
    
    for sim_end in sim_end_list:
        
        mu_dict = {}
        
        for pop_size in pop_size_list:
            temp_dict = {}
            for onset_mu in onset_mu_list:
                temp_dict[onset_mu] = [0,0,0]
            mu_dict[pop_size] = temp_dict    
    
        for dn, pop_size, onset_mu in itertools.product(dn_list, pop_size_list, onset_mu_list):
            
            data_list = data_dict[(dn, pop_size, onset_mu)]
        
            onset_time = np.where(data_list[:,2]>=pop_size)[0][0]
            sim_end_arr = np.where((data_list[:,0]-data_list[onset_time,0]) > sim_end/delta_t)
        
            if len(sim_end_arr[0]) > 0:
                if data_list[-1,2] < 6 and data_list[-1,0]-data_list[-2,0] > 100/delta_t:
                    eli_type = 2
                else:
                    eli_type = 1
            else:                
                if data_list[-1,2] < 10:
                    eli_type = 2
                else:                
                    eli_type = 0
                
            mu_dict[pop_size][onset_mu][eli_type] += 1
            
        type_dict[sim_end] = mu_dict
  
    return type_dict


def plot_type_dict(type_dict, sim_end, pop_size, onset_mu_list, sample_n, la, flag='mf', fit_param=(2,2.5,75)):
    
    fig = plt.figure(figsize=(4,2.5))
    ax = fig.add_subplot(111)
    lw = 2
    
    mu_list = []
    for mu_key in onset_mu_list:
        mu_list.append([mu_key]+type_dict[sim_end][pop_size][mu_key])
    mu_list = np.array(mu_list)       
    
    if 'l' in flag:
        ax.plot(mu_list[:,0]*la,mu_list[:,2]/sample_n,'-',color='k',lw=lw)
        ax.plot(mu_list[:,0]*la,mu_list[:,1]/sample_n,'-',color='tab:red',lw=lw)
        ax.plot(mu_list[:,0]*la,mu_list[:,3]/sample_n,'-',color='tab:green',lw=lw)
        
    if 'm' in flag:
        ax.plot(mu_list[:,0]*la,mu_list[:,2]/sample_n, ls='', marker='^', ms=5, color='k')
        ax.plot(mu_list[:,0]*la,mu_list[:,1]/sample_n, ls='', marker='o', ms=5, color='tab:red')
        ax.plot(mu_list[:,0]*la,mu_list[:,3]/sample_n, ls='', marker='x', ms=5, color='tab:green')
    
    if 'f' in flag:
        if fit_param[0] == 1:
            x, t0 = make_eli_fitting(mu_list[:,0]*la, mu_list[:,1]/sample_n, 0, fit_param[1:])
            x, t1 = make_eli_fitting(mu_list[:,0]*la, mu_list[:,3]/sample_n, 1, fit_param[1:])
        if fit_param[0] == 2:
            x, t0 = make_eli_fitting_p2(mu_list[:,0]*la, mu_list[:,1]/sample_n, 0, fit_param[1:])
            x, t1 = make_eli_fitting_p2(mu_list[:,0]*la, mu_list[:,3]/sample_n, 1, fit_param[1:])
        t2 = 1-t0-t1

        ax.plot(x, t2, color='k', lw=lw)
        ax.plot(x, t0, color='tab:red', lw=lw)
        ax.plot(x, t1, color='tab:green', lw=lw)
        
    ax.set_ylim(-0.05,1.05)
    ax.set_yticks((0,0.5,1))


def plot_normalize_type_dict(type_dict, sim_end, pop_size_list, onset_mu_list, sample_n, la, fit_param=(0.25,74,74)):
    
    fig = plt.figure(figsize=(6,3))
    ax = fig.add_subplot(111)
    lw=0.5
    
    mu_plot_arr = []
    for pop_size in pop_size_list:
        mu_plot_list = [[1,sample_n,0,0]]
        for mu_key in onset_mu_list:
            mu_plot_list.append([(mu_key-1)/pop_size**0.5]\
                                 +type_dict[sim_end][pop_size][mu_key])
        
        mu_plot_arr.append(mu_plot_list)
        mu_plot_list = np.array(mu_plot_list)        
        
        x, t0 = make_eli_fitting_p2(mu_plot_list[:,0]*la, mu_plot_list[:,1]/sample_n, 0, (fit_param[0], fit_param[1]))
        x, t1 = make_eli_fitting_p2(mu_plot_list[:,0]*la, mu_plot_list[:,3]/sample_n, 1, (fit_param[0], fit_param[2]))
        t2 = 1-t0-t1
        
        ax.plot(x, t2, color='k', lw=lw)
        ax.plot(x, t0, color='tab:red', lw=lw)
        ax.plot(x, t1, color='tab:green', lw=lw)
                
    mu_plot_arr = np.array(mu_plot_arr)
    
    x, t0 = make_eli_fitting_p2(mu_plot_arr[:,1:,0].ravel()*la,
                                mu_plot_arr[:,1:,1].ravel()/sample_n, 0, (fit_param[0], fit_param[1]))
    x, t1 = make_eli_fitting_p2(mu_plot_arr[:,1:,0].ravel()*la,
                                mu_plot_arr[:,1:,3].ravel()/sample_n, 1, (fit_param[0], fit_param[2]))
    t2 = 1-t0-t1

    ax.plot(x, t2, color='k', lw=4)    
    ax.plot(x, t0, color='tab:red', lw=4)    
    ax.plot(x, t1, color='tab:green', lw=4)


def make_eli_fitting(x_list, y_list, data_type, p0):
    
    if data_type == 0:
        fit_func = lambda x,p1,p2: 1-1/(1+(p1/x)**p2)
    if data_type == 1:
        fit_func = lambda x,p1,p2: 1/(1+(p1/x)**p2)
    
    popt, pcov = optimize.curve_fit(fit_func, x_list, y_list, p0 = p0)
    
    x = np.linspace(x_list.min(), x_list.max(), 1000)
    
    return x, fit_func(x,popt[0],popt[1])


def make_eli_fitting_p2(x_list, y_list, data_type, p0):
    
    if data_type == 0:
        fit_func = lambda x,p1: 1-1/(1+(p1/x)**p0[1])
    if data_type == 1:
        fit_func = lambda x,p1: 1/(1+(p1/x)**p0[1])
   
    popt, pcov = optimize.curve_fit(fit_func, x_list, y_list, p0 = p0[0])
    
    x = np.linspace(x_list.min(), x_list.max(), 1000)
    
    return x, fit_func(x,popt[0])

    
def make_mean_area_list(delta_t, data_dict, sim_end, dn_list, pop_size_list, onset_mu_list):
    
    result_list = []    
    
    for dn, pop_size, onset_mu in itertools.product(dn_list, pop_size_list, onset_mu_list):
        
        data_list = data_dict[(dn, pop_size, onset_mu)]
        
        onset_time = np.where(data_list[:,2]>=pop_size)[0][0]
        sim_end_arr = np.where((data_list[:,0]-data_list[onset_time,0]) > sim_end/delta_t)
        eli_type = 0
        
        if len(sim_end_arr[0]) > 0:
            if data_list[-1,2] < 6 and data_list[-1,0]-data_list[-2,0] > 100/delta_t:
                eli_type = 2
            else:
                eli_type = 1
        else:                
            if data_list[-1,2] < 10:
                eli_type = 2
            else:                
                eli_type = 0
        
        if eli_type == 1 or eli_type == 0:
            data_list_end = data_list[np.where(data_list[:,0] > (data_list[-1,0]-5/delta_t))]
            time_data = data_list_end[:,0]
            ab_cells = calc_time_ave(time_data,data_list_end[:,2])
            all_area = calc_time_ave(time_data,data_list_end[:,5])            
            result_list.append([onset_mu, ab_cells, all_area, all_area/ab_cells])
        else:
            result_list.append([onset_mu, 0, 0, 0])
                
    return np.array(result_list)

            
def plot_mean_area_list(area_list, la, a_cri=0, a_mh=0):
    
    fig = plt.figure(figsize=(5,4))
    ax = fig.add_subplot(111)
    rng = np.random.default_rng()
    if la == 0.01:
        r = rng.random(len(area_list))*5
    elif la == 0.12:
        r = rng.random(len(area_list))*0.25
        
    ax.plot((area_list[:,0]+r)*la, area_list[:,3], 'k.', ms=2, alpha=0.2)  

    if a_cri > 0:
        ax.plot(area_list[:,0]*la, np.ones_like(area_list[:,0])*a_cri, 'r--')
    if a_mh > 0:
        ax.plot(area_list[:,0]*la, np.ones_like(area_list[:,0])*a_mh, 'b--')
 
    if la == 0.01:
        ax.set_xlim(1.,4.2)
        ax.set_ylim(-0.01, 0.5)
    elif la == 0.12:
        ax.set_xlim(0.15,2.5)
        ax.set_ylim(-0.01, 0.25)
        

# To read simulation result files
fn = "L001G0025_IC1"                            # Input file name
folder_name = data_path+fn+"/"
delta_t = 0.0001
sim_end = 490                                   # Input simulation finish time
sim_end_list = np.arange(10,500,10)
dn_list = list(range(1,21))                     # Input simulation runs
la, ga = file_index[fn][1], file_index[fn][2]   # Input Lambda, Gamma (Auto)
onset_mu_list = file_index[fn][3]               # Input onset contraction force (mu) list (Auto)
pop_size_list = file_index[fn][4]               # Input onset population size (N_theta) list (Auto)

data_dict = {}
for dn, pop_size, onset_mu in itertools.product(dn_list, pop_size_list, onset_mu_list):
    data_list = read_file(folder_name, dn, pop_size, onset_mu)
    data_dict[(dn, pop_size, onset_mu)] = data_list
    
# To plot time course of each simulation run (Fig 2)
# plot_series(delta_t, data_dict, 1, 100, 10, 0, 0)

# To plot distribution of the number of cell sides (Fig 3)
# plot_lewis_law_1_ab(delta_t, data_dict, 1, 100, 10)


# To plot the frequencies of elimination failure, growth suspension, and elimination success (Fig 2)
type_dict = make_type_dict(delta_t, data_dict, dn_list, pop_size_list, onset_mu_list, sim_end_list)
plot_type_dict(type_dict, sim_end, 100, onset_mu_list, len(dn_list), la, 'f')
plot_normalize_type_dict(type_dict, sim_end, pop_size_list, onset_mu_list, len(dn_list), la)


# To plot the time averages of cell area (Fig 2)
area_list = make_mean_area_list(delta_t, data_dict, sim_end, dn_list, [100], onset_mu_list)
plot_mean_area_list(area_list, la)

 
