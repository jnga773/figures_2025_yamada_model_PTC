#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 21 12:16:27 2025

@author: jnga773
"""

# # Append Python path to add AUTO functions
# import sys
# sys.path.append('$HOME/auto/07p/python')
# sys.path.append('$HOME/auto/07p/python/auto')

# Load extra functions
import auto
import data_functions as data_funcs

# %%
#==============================================================================#
##                            INITIAL CONTINUATION                            ##
#==============================================================================#
# We calculate the initial periodic orbit of the Yamada model.

# Set parameters
gamma_PO = 3.5e-2
A_PO     = 7.4
B        = 5.8
a        = 1.8

# Parameter vector
p0 = {1: gamma_PO, 2: A_PO, 3: B, 4: a}

# Initial solution is the 'off' state
x0 = [10, 10, 10]

# State-space variable names
unames = {1: 'x1', 2: 'x2', 3: 'x3'}

# %%
#------------------------------------------------------------------------------#
#                    Confirm ODE45 Periodic Orbit Solution                     #
#------------------------------------------------------------------------------#
# Calculate the periodic orbit using MATLAB's ode45 function.

#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'run01_initial_PO_ode45'

#--------------------------#
#     Print to Console     #
#--------------------------#
print('=====================================================================')
print('Initial Periodic Orbit: First Run')
print('Find new periodic orbit')
print('---------------------------------------------------------------------')
print('This run name           : {}'.format(run_new_str))
print('Continuation parameters : {}'.format('gamma, T, A'))
print('=====================================================================')

#----------------------------#
#     Calculate Solution     #
#----------------------------#
# Solve using scipy's solve_ivp
x_init_solve_ivp = data_funcs.calc_initial_solution_solve_ivp(x0, p0)

# Parameter names
pnames_PO = {1: 'gamma', 2: 'A', 3: 'B', 4: 'a', 11: 'T'}

#-------------------------------#
#     Run AUTO Continuation     #
#-------------------------------#
# Saved points
UZR = {'A': A_PO, 'gamma': gamma_PO}
# Continuation parameters
pcont = ['gamma', 'T', 'A']
# Parameter range
prange = {'A': [5.0, 20.0], 'gamma': [0.0, 0.4]}

# Run continuation
run_new = auto.run(x_init_solve_ivp, e='./functions/yamada', IPS=2, IRS=0,
                   NPAR=len(pnames_PO), PAR=p0, parnames=pnames_PO, NDIM=len(unames), unames=unames,
                   ICP=pcont, UZSTOP=prange, UZR=UZR,
                   JAC=1, NBC=0, NINT=0,
                   NTST=50, DSMIN=1e-1, DS=1e-1, DSMAX=1e-1,
                   NMX=500, NPR=50)

#-------------------#
#     Save Data     #
#-------------------#
# Save data
data_funcs.save_move_data(run_new, run_new_str)

# Print clear line
print('\n')

# %%
#------------------------------------------------------------------------------#
#            Shifted Periodic Orbit with Zero-Phase Point Condition            #
#------------------------------------------------------------------------------#
# Continuing from a specified point in the previous continuation, we shift the
# origin of the periodic orbit to match the 'zero-phase point' condition, and
# continue the periodic orbit as a new BVP.

#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'run02_initial_PO'
# Previous run name
run_old_str = 'run01_initial_PO_ode45'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
label_old = run_old('UZ1')
label_old = label_old['LAB']

#--------------------------#
#     Print to Console     #
#--------------------------#
print('=====================================================================')
print('Initial Periodic Orbit: Second Run')
print('Continue periodic orbit with shifted phase condition')
print('---------------------------------------------------------------------')
print('This run name           : {}'.format(run_new_str))
print('Previous run name       : {}'.format(run_old_str))
print('Previous label_solution : {}'.format(label_old))
print('Continuation parameters : {}'.format('A, gamma'))
print('=====================================================================')

#-------------------#
#     Read Data     #
#-------------------#
# Calculate initial solution from previous run
x_init_PO, p_PO, pnames_PO = data_funcs.calc_initial_solution_PO(run_old(label_old))

#-------------------------------#
#     Run AUTO Continuation     #
#-------------------------------#
# Saved points
UZR = {'A': A_PO, 'gamma': gamma_PO}
# Continuation parameters
pcont = ['gamma', 'T']
# Parameter range
prange = {'gamma': [0.0, 0.4]}

#-------------------------------#
#     Run AUTO Continuation     #
#-------------------------------#
# Run continuation
run_new = auto.run(x_init_PO, e='./functions/yamada_VAR', IPS=4, IRS=0, LAB=1,
                   NPAR=len(pnames_PO), PAR=p_PO, parnames=pnames_PO,
                   NDIM=len(unames), unames=unames,
                   ICP=pcont, UZSTOP=prange, UZR=UZR,
                   JAC=0, NBC=4, NINT=0,
                   NMX=10, NPR=1,
                   DSMIN=1e-3, DS=-1e-3, DSMAX=1e-2,
                   NCOL=4, IAD=1, NTST=50)

#-------------------#
#     Save Data     #
#-------------------#
# Save data
data_funcs.save_move_data(run_new, run_new_str)

# Print clear line
print('\n')

# %%
#==============================================================================#
##                 COMPUTE FLOQUET BUNDLE AT ZERO PHASE POINT                 ##
#==============================================================================#
# Here we compute the stable Floquet bundle of the periodic orbit, as well
# as the perpendicular vector, w.

#------------------------------------------------------------------------------#
#               Continue the Floquet multiplier until mu_s = 1.0               #
#------------------------------------------------------------------------------#
# We now add the adjoint function and Floquet boundary conditions to
# compute the adjoint (left or right idk) eigenvectors and eigenvalues.
# This will give us the perpendicular vector to the tangent of the periodic
# orbit. However, this will only be for the eigenvector corresponding to
# the eigenvalue \mu = 1. Hence, here we continue in \mu (mu_f) until
# mu_f = 1.

#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'run03_floquet_mu'
# Previous run name
run_old_str = 'run02_initial_PO'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
label_old = 1

#--------------------------#
#     Print to Console     #
#--------------------------#
print('=====================================================================')
print('Floquet Bundle: First Run')
print('Calculate stable Floquet bundle eigenvalue')
print('---------------------------------------------------------------------')
print('This run name           : {}'.format(run_new_str))
print('Previous run name       : {}'.format(run_old_str))
print('Previous label_solution : {}'.format(label_old))
print('Continuation parameters : {}'.format('mu_s, w_norm'))
print('=====================================================================')

#-------------------#
#     Read Data     #
#-------------------#
# Calculate initial solution from previous run
x_init_VAR, p_VAR, unames_VAR, pnames_VAR = data_funcs.calc_initial_solution_VAR(run_old(label_old))

# Parameter names
pnames_PO = {1: 'gamma', 2: 'A', 3: 'B', 4: 'a', 11: 'T'}

#-------------------------------#
#     Run AUTO Continuation     #
#-------------------------------#
# Saved points
UZR = {'mu_s': 1.0}
# Continuation parameters
pcont = ['mu_s', 'w_norm', 'gamma']
# Parameter range
prange = {'mu_s': [0.0, 1.1]}

#-------------------------------#
#     Run AUTO Continuation     #
#-------------------------------#
# Run continuation
run_new = auto.run(x_init_VAR, e='./functions/yamada_VAR', IRS=0, IPS=4, ISW=1,
                   NPAR=len(p_VAR), PAR=p_VAR, parnames=pnames_VAR,
                   NDIM=len(unames_VAR), unames=unames_VAR,
                   ICP=pcont, UZSTOP=prange, UZR=UZR,
                   JAC=0, NBC=8, NINT=0,
                   NTST=50, NCOL=4, IAD=1,
                   DSMIN=5e-4, DS=1e-3, DSMAX=1e-2,
                   NMX=300, NPR=10)

#-------------------#
#     Save Data     #
#-------------------#
# Save data
data_funcs.save_move_data(run_new, run_new_str)

# Print clear line
print('\n')

# %%
#------------------------------------------------------------------------------#
#                   Continue the w-vector until w_norm = 1.0                   #
#------------------------------------------------------------------------------#
# Having found the solution (branching point 'BP') corresponding to
# \mu = 1, we can continue in the norm of the vector w (w_norm), until the
# norm is equal to zero. Then we will have the correct perpendicular
# vector.

#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'run04_floquet_wnorm'
# Previous run name
run_old_str = 'run03_floquet_mu'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
label_old = run_old('BP1')
label_old = label_old['LAB']

#--------------------------#
#     Print to Console     #
#--------------------------#
print('=====================================================================')
print('Floquet Bundle: Second Run')
print('Grow norm of stable Floquet bundle vector')
print('---------------------------------------------------------------------')
print('This run name           : {}'.format(run_new_str))
print('Previous run name       : {}'.format(run_old_str))
print('Previous label_solution : {}'.format(label_old))
print('Continuation parameters : {}'.format('mu_s, w_norm'))
print('=====================================================================')

#-------------------------------#
#     Run AUTO Continuation     #
#-------------------------------#
# Run continuation
run_new = auto.run(run_old(label_old), LAB=1,
                   ISW=-1, ILP=0,
                   DSMIN=5e-2, DS=1e-1, DSMAX=5e-1,
                   NMX=500, NPR=100,
                   UZSTOP={'w_norm': [0.0, 1.1]}, UZR={'w_norm': 1.0})

#-------------------#
#     Save Data     #
#-------------------#
# Save data
data_funcs.save_move_data(run_new, run_new_str)

# Print clear line
print('\n')

# %%
#==============================================================================#
##               PHASE RESPONSE - DIRECTIONAL TRANSITION CURVES               ##
#==============================================================================#
# We set up the phase resetting problem by creating four segments of the
# periodic orbit, with boundary conditions described in a paper somewhere.

#------------------------------------------------------------------------------#
##                    Change theta_old Along Periodic Orbit                   ##
#------------------------------------------------------------------------------#
#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'run05_PR_move_theta_old'
# Previous run name
run_old_str = 'run04_floquet_wnorm'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
label_old = run_old('UZ1')
label_old = label_old['LAB']

#--------------------------#
#     Print to Console     #
#--------------------------#
print('=====================================================================')
print('Directional Transition Curve: First Run')
print('Move along periodic orbit')
print('---------------------------------------------------------------------')
print('This run name           : {}'.format(run_new_str))
print('Previous run name       : {}'.format(run_old_str))
print('Previous label_solution : {}'.format(label_old))
print('Continuation parameters : {}'.format('theta_old, theta_new, eta, mu_s, T'))
print('=====================================================================')

#-------------------#
#     Read Data     #
#-------------------#
# Set initial phase resetting parameters
# Periodicity
k = 60

# Perturbation direction (in units of 2 \pi)
theta_perturb = 0.0

# Calculate initial solution
x_init_PR, p_PR, unames, pnames_PR = \
    data_funcs.calc_initial_solution_PR(run_old(label_old), k, theta_perturb)

#-------------------------------#
#     Run AUTO Continuation     #
#-------------------------------#
# Saved points perturbation amplitudes
SP_points = [0.339386, 1.339386]

# Set saved points
UZR = {'theta_old': SP_points}

# Set continuation parameters
pcont = ['theta_old', 'theta_new', 'eta', 'mu_s', 'T']
# Set continuation stop points
prange = {'theta_old': [0.0, 1.0]}

# Run continuation
run_new = auto.run(dat='./initial_solution_PR.dat', PAR=p_PR, parnames=pnames_PR,
                   NPAR=len(p_PR), NDIM=len(unames), unames=unames,
                   e='./functions/yamada_PR', IRS=0, IPS=4, ISW=1, ILP=0, ISP=0,
                   JAC=1, NBC=22, NINT=0,
                   ICP=pcont, UZSTOP=prange, UZR=UZR,
                   NTST=k*50, NCOL=4, IAD=10,
                   DSMIN=1e-2, DS=-5e-1, DSMAX=1e0,
                   NMX=1000, NPR=25)

#-------------------#
#     Save Data     #
#-------------------#
# Save data
data_funcs.save_move_data(run_new, run_new_str)

# Print clear line
print('\n')

# %%
#------------------------------------------------------------------------------#
##                         Change Perturbation Angle                          ##
#------------------------------------------------------------------------------#
#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'run06_PR_change_angle'
# Previous run name
run_old_str = 'run05_PR_move_theta_old'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
label_old = run_old('UZ')
label_old = [sol['LAB'] for sol in label_old[:-1]]
label_old = label_old[0]

#--------------------------#
#     Print to Console     #
#--------------------------#
print('=====================================================================')
print('Directional Transition Curve: Second Run')
print('Change perturbation angle')
print('---------------------------------------------------------------------')
print('This run name           : {}'.format(run_new_str))
print('Previous run name       : {}'.format(run_old_str))
print('Previous label_solution : {}'.format(label_old))
print('Continuation parameters : {}'.format('theta_perturb, theta_new, eta, mu_s, T'))
print('=====================================================================')

#-------------------------------#
#     Run AUTO Continuation     #
#-------------------------------#
# Saved points perturbation amplitudes
# SP_points = [0.0, 0.25, 0.5, 0.75]

from numpy import arange
SP_points = arange(0.0, 1.0, 0.125)

# Set saved points
UZR = {'theta_perturb': SP_points}

# Set continuation parameters
pcont = ['theta_perturb', 'theta_new', 'eta', 'mu_s', 'T']
# Set continuation stop points
# prange = {'theta_perturb': [-0.1, 0.5]}
prange = {'theta_perturb': [0.0, 1.0]}

# Run continution
data_funcs.run_PR_continuation(run_new_str, run_old, label_old,
                               pcont, prange,
                               UZR=UZR, NMX=4000, NPR=100,
                               DSMIN=1e-3, DS=1e-1, DSMAX=1e0,
                               reverse=False)

# %%
#------------------------------------------------------------------------------#
##                      Increase Perturbation Amplitude                       ##
#------------------------------------------------------------------------------#
#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'run07_PR_increase_perturbation'
# Previous run name
run_old_str = 'run06_PR_change_angle'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
label_old = run_old('UZ')
label_old = [sol['LAB'] for sol in label_old]

#-------------------------------------#
#     Cycle Through Previous Runs     #
#-------------------------------------#
for run in range(len(label_old)):
    # This label
    this_run_label = label_old[run]

    # Data directory for this run
    this_run_name = '{}/theta_perturb_{}'.format(run_new_str, str(run+1).zfill(2))

    #--------------------------#
    #     Print to Console     #
    #--------------------------#
    print('=====================================================================')
    print('Directional Transition Curve: Third Run')
    print('Increase perturbation amplitude')
    print('---------------------------------------------------------------------')
    print('This run name           : {}'.format(this_run_name))
    print('Previous run name       : {}'.format(run_old_str))
    print('Previous label_solution : {}'.format(this_run_label))
    print('Continuation parameters : {}'.format('A_perturb, theta_new, eta, mu_s, T'))
    print('=====================================================================')
    
    #-------------------------------#
    #     Run AUTO Continuation     #
    #-------------------------------#
    # Set saved points for perturbation
    SP_points = [0.1, 0.724236, 10.0]
    # UZR = {'theta_old': SP_points}
    UZR = {'A_perturb': SP_points}
    
    # Set continuation parameters
    pcont = ['A_perturb', 'theta_new', 'eta', 'mu_s', 'T']
    # Set continuation stop points
    prange = {'A_perturb': [0.0, max(SP_points)+0.1],
              'eta': [-1e-4, 1e-2]}
    True
    # Run continution
    data_funcs.run_PR_continuation(this_run_name, run_old, this_run_label,
                                    pcont, prange, UZR=UZR,
                                    NMX=10000, NPR=100,
                                    DSMIN=1e-3, DS=1e-1, DSMAX=1e0,
                                    reverse=False)
    
# %%
#------------------------------------------------------------------------------#
##                               Calculate DTCs                               ##
#------------------------------------------------------------------------------#
#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'run08_PR_DTC_scan'
# Previous run name
run_old_str = 'run07_PR_increase_perturbation'

# Read directories inside run_old
from os import listdir

#-------------------------------------#
#     Cycle Through Previous Runs     #
#-------------------------------------#
# Cycle through theta_perturb directories
for idx_theta in range(len(listdir('./data/{}/'.format(run_old_str)))):
    # Set run string for this run to read
    theta_run_name = '{}/theta_perturb_{}'.format(run_old_str, str(idx_theta+1).zfill(2))

    # Read bifurcation data
    run_old = data_funcs.bd_read(theta_run_name)
    
    # Read previous labels
    label_old = run_old('UZ')
    label_old = [sol['LAB'] for sol in label_old]
    
    for run in range(len(label_old)):
        # This label
        this_run_label = label_old[run]

        # Data directory for this run
        this_run_name = ('{}/theta_perturb_{}/DTC_{}'
                            ).format(run_new_str,
                                    str(idx_theta+1).zfill(2),
                                    str(run+1).zfill(2))

        #--------------------------#
        #     Print to Console     #
        #--------------------------#
        print('=====================================================================')
        print('Directional Transition Curve: Fourth Run')
        print('Calculate DTC')
        print('---------------------------------------------------------------------')
        print('This run name           : {}'.format(this_run_name))
        print('Previous run name       : {}'.format(theta_run_name))
        print('Previous label_solution : {}'.format(this_run_label))
        print('Continuation parameters : {}'.format('theta_perturb, theta_new, eta, mu_s, T'))
        print('=====================================================================')
        
        #-------------------------------#
        #     Run AUTO Continuation     #
        #-------------------------------#    
        # Set continuation parameters
        pcont = ['theta_perturb', 'theta_new', 'eta', 'mu_s', 'T']
        # Set continuation stop points
        prange = {'theta_perturb': [0.0, 1.0]}
        
        # Run continution
        data_funcs.run_PR_continuation(this_run_name, run_old, this_run_label,
                                        pcont, prange,
                                        NMX=2000, NPR=100,
                                        DSMIN=1e-3, DS=1e-2, DSMAX=1e0,
                                        reverse=True)

# %%
#==============================================================================#
#                                 END OF FILE                                  #
#==============================================================================#
data_funcs.clean_directories()
