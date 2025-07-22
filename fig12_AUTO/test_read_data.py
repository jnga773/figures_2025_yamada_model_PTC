#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 22 11:47:15 2025

@author: jnga773
"""

# Load extra functions
import auto
import data_functions as data_funcs
from os import listdir
import numpy as np
import matplotlib.pyplot as plt

# %%
# This run name
run_in = 'run08_PR_DTC_scan'

# Get list of directories inside run_in
dir_A = sorted(listdir('./data/{}/'.format(run_in)))

# Empty arrays
A_perturb_data = np.zeros((3, 2, 2))
theta_old_data = np.zeros((3, 2, 2))

# Bifurcation data
theta_perturb_A1 = []
theta_perturb_A2 = []
theta_perturb_A3 = []
theta_new_A1     = []
theta_new_A2     = []
theta_new_A3     = []

# Cycle through A_perturb directories
for idx_A, sub_run_A in enumerate(dir_A):
    # Get list of directories inside sub_run_A
    dir_P = sorted(listdir('./data/{}/{}/'.format(run_in, sub_run_A)))
    
    # Cycle through theta_perturb directories
    for idx_P, sub_run_P in enumerate(dir_P):
        # Get list of directories inside sub_run_P
        dir_O = sorted(listdir('./data/{}/{}/{}/'.format(run_in, sub_run_A, sub_run_P)))
        
        # Cycle through theta_old directories
        for idx_O, sub_run_O in enumerate(dir_O):
            # Create run name
            sub_run_name = '{}/{}/{}/{}'.format(run_in, sub_run_A, sub_run_P, sub_run_O)
        
            # Read bifurcation data
            bd_read = data_funcs.bd_read(sub_run_name)
            # Get first solution
            sol_read = bd_read(1)
            
            # Read data
            A_perturb_data[idx_A, idx_P, idx_O] = sol_read['A_perturb']
            theta_old_data[idx_A, idx_P, idx_O] = sol_read['theta_old']
            theta_perturb_read = bd_read['theta_perturb']
            theta_new_read     = bd_read['theta_new']
            
            # Append data
            if idx_A == 0:
                theta_perturb_A1.append(theta_perturb_read)
                theta_new_A1.append(theta_new_read)
            elif idx_A == 1:
                theta_perturb_A2.append(theta_perturb_read)
                theta_new_A2.append(theta_new_read)
            elif idx_A == 2:
                theta_perturb_A3.append(theta_perturb_read)
                theta_new_A3.append(theta_new_read)
                
# %%
plt.close('all')
fig = plt.figure(figsize=[6, 6])
ax = plt.gca()

theta_perturb_plot = theta_perturb_A3
theta_new_plot     = theta_new_A3

theta_old_plot = theta_old_data[0]
A_perturb_plot = A_perturb_data[0, 0, 0]

# Plot
# ax.plot(theta_perturb_plot[0], theta_new_plot[0], color='C0', ls='solid')
# ax.plot(theta_perturb_plot[1], theta_new_plot[1], color='C1', ls='dashed')
# ax.plot(theta_perturb_plot[2], theta_new_plot[2], color='C2', ls='dotted')
ax.plot(theta_perturb_plot[3], theta_new_plot[3], color='C3', ls='dashdot')

ax.set_xlim(-1, 1)
ax.set_ylim(-2, 2)

fig.tight_layout()
fig.show()
