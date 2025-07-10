#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 10 14:13:12 2025

@author: jnga773
"""

import numpy as np
import matplotlib.pyplot as plt
# plt.

plt.close('all')

#-----------------------------------------------------------------------------#
#                                  FUNCTIONS                                  #
#-----------------------------------------------------------------------------#
def sort_data_folders(run_str_in):
    from os import listdir
    from numpy import argsort
    
    # List all directories in ./data/run08_isochron_scan/
    dir_main = './data/{}/'.format(run_str_in)
    dir_list = listdir(dir_main)

    dir_list = [f for f in dir_list if not f.startswith('.')]
    
    # Cyycle through and get the integer number
    int_list = []
    for f in dir_list:
        int_read = f[4:]
        int_list.append(int(int_read))
        
    # Sort int_list and get indices
    idx_sorted = argsort(int_list)
    
    dir_list_sorted = [dir_list[i] for i in idx_sorted]
    
    return dir_list_sorted

def read_scan_bd(run_str_in):
    """
    Reads all of the b.dat files from the scan run and outputs a list of
    bd files.
    """
    from continuation_scripts.data_functions import bd_read

    # Get sorted data folders
    data_dir = sort_data_folders(run_str_in)
    
    #-------------------#
    #     Read Data     #
    #-------------------#
    bd_out = []
    
    # Read data
    for idx, run in enumerate(data_dir):
        # Read bifurcation data
        bd = bd_read('{}/{}'.format(run_str_in, run))
        bd_out.append(bd)
        
    #----------------#
    #     Output     #
    #----------------#
    return bd_out

#-----------------------------------------------------------------------------#
#                                  DATA THINGS                                #
#-----------------------------------------------------------------------------#
# Run name
run_name = 'run10_DTC_scan'
# Read bifurcation data
bd_read  = read_scan_bd(run_name)

#-----------------------#
#     Read DTC Data     #
#-----------------------#
idx = 3

sol = bd_read[idx](1)

# Read A_perturb and theta_old
A_perturb = sol['A_perturb']
theta_old = sol['theta_old']

import auto
bd1 = auto.merge(bd_read[idx])

# Read DTC data
theta_perturb = bd1['theta_perturb']
theta_new     = bd1['theta_new'] - 1.0

#-----------------------------------------------------------------------------#
#                                 PLOT DTC                                    #
#-----------------------------------------------------------------------------#
# Figure size in centimetres
figsize = np.array([7.4, 9.1])
figsize *= 1 / 2.54

# Create figure
fig = plt.figure(num='DTC', figsize=figsize.tolist())
ax = plt.gca()

# Plot
ax.plot(theta_perturb, theta_new, color='C0', ls='solid')

# Fundamental domain
from matplotlib.patches import Rectangle
fund = Rectangle((-3, 0), 6, 1, linewidth=0, edgecolor='none', facecolor='C2', alpha=0.25)
ax.add_patch(fund)

# Labels
ax.set_xlabel(r'$\varphi_{d}$')
ax.set_ylabel(r'$\vartheta_{\mathrm{n}}$')
title_str = '$A = {:.3f}, \\vartheta_{{\mathrm{{o}}}} = {:.3f}$'.format(A_perturb, theta_old)
ax.set_title(title_str)

# Limits
ax.set_xlim(-1.0, 1.0)

# Figure stuff
fig.tight_layout()
fig.show()

