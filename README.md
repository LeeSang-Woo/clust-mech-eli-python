# clust-mech-eli-python
 This python codes was used in the paper "Two types of critical cell density for mechanical elimination of abnormal cell clusters from epithelial tissue" for analytical solution and data analysis

#### Environment
 Python 3.8.5

#### Python code files
- sim_data_analy.py : Data analysis of simulation results, plotting (simulation result files required)
- analy_sol.py : Analytical solution for the paper, plotting
- param_diagram.py : Plotting the parameter space of vertex dynamics model
- init_cond.py : Plotting the initial condition of simulations
- area_mce.py : About mechanical homeostatic density, plotting
- calc_func.py : Common functions for computation

#### Data folders
- init_data/ : Initial conditions of simulations (for init_cond.py)
- sim_data/ : Simulation result files (for sim_data_analy.py)

#### Simulation
 To simulate, put the release file of C++ ([clust-mech-eli.exe](https://github.com/LeeSang-Woo/clust-mech-eli-cpp/releases/download/v1.0/clust-mech-eli.exe)) in **sim_data/(TargetFolder)/** and execute.
