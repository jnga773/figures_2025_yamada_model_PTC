#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 16:42:11 2025

@author: jnga773
"""

# # Append Python path to add AUTO functions
# import sys
# sys.path.append('$HOME/auto/07p/python')
# sys.path.append('$HOME/auto/07p/python/auto')

# Load extra functions
import auto
import continuation_scripts.data_functions as data_funcs

# %%
#==============================================================================#
##                            INITIAL CONTINUATION                            ##
#==============================================================================#
# We calculate the initial periodic orbit of the Yamada model.

# Set parameters
gamma = 0.1
A     = 6.6
B     = 5.8
a     = 1.8
p0 = {1: gamma, 2: A, 3: B, 4: a}

# Initial solution is the 'off' state
x0 = [A, B, 0]

#------------------------------------------------------------------------------#
#                           Compute Equilibrium Point                          #
#------------------------------------------------------------------------------#
# We compute and continue the equilibrium point of the model.

#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'run01_initial_EP'

# Print to console
print('~~~ Initial Periodic Orbit: First Run (c.initial_EP) ~~~')
print('Initial continuation from some point x0')
print('Run name: {}'.format(run_new_str))

#-------------------------------#
#     Run AUTO Continuation     #
#-------------------------------#
# # Copy continuation script
auto.copy('./continuation_scripts/', 'initial_EP')

# Run the first continuation from the initial equilibrium point
run_new = auto.run(x0, PAR=p0, c='initial_EP')

#-------------------#
#     Save Data     #
#-------------------#
# Save data
data_funcs.save_move_data(run_new, run_new_str)

# Print clear line
print('\n')

# %%
#------------------------------------------------------------------------------#
#                         Continue from Branching Point                        #
#------------------------------------------------------------------------------#
# Continue from the branching point until we detect a Hopf bifurcation.

#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'run02_branching_point'
# Previous run name
run_old_str = 'run01_initial_EP'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
label_old = run_old('BP1')
label_old = label_old['LAB']

# Print to console
print('~~~ Initial Periodic Orbit: Second Run ~~~')
print('Continue bifurcations from the branching point')
print('Run name: {}'.format(run_new_str))
print('Continuing from point {} in run: {}'.format(label_old, run_old_str))

#-------------------------------#
#     Run AUTO Continuation     #
#-------------------------------#
# Continue from branching point in run01_initial_EP
run_new = auto.run(run_old('BP1'), ISW=-1)

#-------------------#
#     Save Data     #
#-------------------#
# Save data
data_funcs.save_move_data(run_new, run_new_str)

# Print clear line
print('\n')

# %%
#------------------------------------------------------------------------------#
#                           Continue Hopf To A = 7.3757                        #
#------------------------------------------------------------------------------#
# We continue the Hopf bifurcation, varying the 'A' parameter until A = 7.3757

#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'run03_hopf'
# Previous run name
run_old_str = 'run02_branching_point'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
label_old = run_old('HB1')
label_old = label_old['LAB']

# Print to console
print('~~~ Initial Periodic Orbit: Third Run ~~~')
print('Follow Hopf birfucation until z=-0.8')
print('Run name: {}'.format(run_new_str))
print('Continuing from point {} in run: {}'.format(label_old, run_old_str))

#-------------------------------#
#     Run AUTO Continuation     #
#-------------------------------#
# Continue from Hop bifurcation
run_new = auto.run(run_old(label_old), ISW=2,
                   ICP=['A', 'gamma'],
                   UZSTOP={'gamma': [0.0, 0.4], 'A': [5.0, 20.0]},
                   DSMIN=5e-3, DS=5e-3, DSMAX=5e-3,
                   NMX=500, NPR=50,
                   UZR={'A': 7.37570000})

#-------------------#
#     Save Data     #
#-------------------#
# Save data
data_funcs.save_move_data(run_new, run_new_str)

# Print clear line
print('\n')

# %%
#------------------------------------------------------------------------------#
#                            Hopf to Periodic Orbit                            #
#------------------------------------------------------------------------------#
# We compute a family of periodic orbits originating off the Hopf bifurcation.

#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'run04_hopf_to_PO'
# Previous run name
run_old_str = 'run03_hopf'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
label_old = run_old('GH1')
label_old = label_old['LAB']

# Print to console
print('~~~ Initial Periodic Orbit: Fourth Run ~~~')
print('Continue periodic orbits from the Hopf bifurcation')
print('Run name: {}'.format(run_new_str))
print('Continuing from point {} in run: {}'.format(label_old, run_old_str))

#-------------------------------#
#     Run AUTO Continuation     #
#-------------------------------#
# Gamma value for saved point
SP = 3.54e-2

# Follow periodic orbits
run_new = auto.run(run_old(label_old), ISW=-1, IPS=2, LAB=1,
                   ICP=['gamma'],
                   NTST=50, DSMIN=1e-1, DS=1e-1, DSMAX=1e-1,
                   NMX=500, NPR=50,
                   UZR={'gamma': SP})

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
run_new_str = 'run05_initial_PO'
# Previous run name
run_old_str = 'run04_hopf_to_PO'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
label_old = run_old('UZ1')
label_old = label_old['LAB']

# Print to console
print('~~~ Initial Periodic Orbit: Fifth Run (c.initial_PO) ~~~')
print('Continue periodic orbit with shifted phase condition')
print('Run name: {}'.format(run_new_str))
print('Continuing from point {} in run: {}'.format(label_old, run_old_str))

#-------------------#
#     Read Data     #
#-------------------#
# Calculate initial solution from previous run
x_init_PO, p_PO, pnames_PO = data_funcs.calc_initial_solution_PO(run_old(label_old))

#-------------------------------#
#     Run AUTO Continuation     #
#-------------------------------#
# Copy continuation script
auto.copy('./continuation_scripts/', 'initial_PO')

# Run continuation
run_new = auto.run(x_init_PO, PAR=p_PO, parnames=pnames_PO,
                   c='initial_PO', LAB=1)

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
run_new_str = 'run06_floquet_mu'
# Previous run name
run_old_str = 'run05_initial_PO'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
label_old = 1

# Print to console
print('~~~ Floquet Bundle: First Run (c.floquet_variational) ~~~')
print('Calculate Floquet bundle (mu) ')
print('Run name: {}'.format(run_new_str))
print('Continuing from point {} in run: {}'.format(label_old, run_old_str))

#-------------------#
#     Read Data     #
#-------------------#
# Calculate initial solution from previous run
x_init_VAR, p_VAR, pnames_VAR = data_funcs.calc_initial_solution_VAR(run_old(label_old))

#-------------------------------#
#     Run AUTO Continuation     #
#-------------------------------#
# Copy continuation script
auto.copy('./continuation_scripts/', 'floquet_variational')

# Run continuation
run_new = auto.run(x_init_VAR, PAR=p_VAR, parnames=pnames_VAR,
                   c='floquet_variational', NPR=50,
                   DSMIN=1e-4, DS=1e-4, DSMAX=1e-3)

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
run_new_str = 'run07_floquet_wnorm'
# Previous run name
run_old_str = 'run06_floquet_mu'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
label_old = run_old('BP1')
label_old = label_old['LAB']

# Print to console
print('~~~ Floquet Bundle: Second Run (c.floquet_variational) ~~~')
print('Calculate Floquet bundle (w_norm) ')
print('Run name: {}'.format(run_new_str))
print('Continuing from point {} in run: {}'.format(label_old, run_old_str))

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
# Save solution to MATLAB .mat file
# old_functions.save_floquet_data_matlab(run_new('UZ1'))

# Save data
data_funcs.save_move_data(run_new, run_new_str)

# Print clear line
print('\n')

# %%
#==============================================================================#
##                  PHASE RESPONSE - PHASE TRANSITION CURVES                  ##
#==============================================================================#
# We set up the phase resetting problem by creating four segments of the
# periodic orbit, with boundary conditions described in a paper somewhere.

#------------------------------------------------------------------------------#
##                 First Continuation: Perturbation Amplitude                 ##
#------------------------------------------------------------------------------#
# We first compute the phase response to a perturbation in a fixed direction
# applied at the zero-phase point, theta_old = 1.0 (or 0.0). We free  A_perturb
# and theta_old.

#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'run08_phase_reset_perturbation'
# Previous run name
run_old_str = 'run07_floquet_wnorm'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
label_old = run_old('UZ1')
label_old = label_old['LAB']

# Print to console
print('~~~ Phase Reset: First Run ~~~')
print('Continue in the perturbation ampltiude A_perturb')
print('Run name: {}'.format(run_new_str))

#-------------------#
#     Read Data     #
#-------------------#
# Set initial phase resetting parameters
# Periodicity
k = 55

# Perturbation direction
from numpy import pi
theta_perturb = 0.0
# theta_perturb = 0.5 * pi

# Calculate initial solution
x_init_PR, p_PR, pnames_PR = \
    data_funcs.calc_initial_solution_PR(run_old(label_old), k, theta_perturb)

#-------------------------------#
#     Run AUTO Continuation     #
#-------------------------------#
# Set saved points
from numpy import linspace, concatenate, unique

# Saved points for large scan of G perturbation
SP_points = concatenate((linspace(0.0, 0.25, 25),
                         linspace(0.30, 1.0, 25),
                         linspace(1.0, 1.3, 25),
                         linspace(1.3, 2.0, 25)))
SP_points = concatenate((SP_points, [0.05, 0.1, 0.15, 0.55, 1.0, 1.5, 2.0]))
SP_points = unique(SP_points)

# Copy continuation script
auto.copy('./continuation_scripts/', 'PTC_initial')

# Try set up phase reset calculation lol
run_new = auto.run(dat='./initial_solution_PR.dat', PAR=p_PR, parnames=pnames_PR,
                   c='PTC_initial',
                   NMX=2000, NTST=k * 50,
                   UZR={'A_perturb': SP_points},
                   UZSTOP={'A_perturb': max(SP_points) + 0.1})

#-------------------#
#     Save Data     #
#-------------------#
# Save data
data_funcs.save_move_data(run_new, run_new_str)

# Print clear line
print('\n')


# %%
#------------------------------------------------------------------------------#
##           Second Continuation: Phase Transition Curve (Multiple)           ##
#------------------------------------------------------------------------------#
# We then fix the perturbation amplitude A_perturb and free theta_old and
# theta_new to compute the phase transition curve (PTC).

#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'run09_PTC_scan'
# Previous run name
run_old_str = 'run08_phase_reset_perturbation'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
label_old = run_old('UZ')
label_old = [sol['LAB'] for sol in label_old[:-1]]

# Print to console
print('~~~ Phase Reset: Second Run ~~~')
print('Compute phase transition curve (PTC)')
print('Run name: {}'.format(run_new_str))
print('Continuing from (Saved Points) in run: {}'.format(run_old_str))

#--------------------------------------#
#     Define Continuation Function     #
#--------------------------------------#
# Define function for parallelising
def calculate_PTC(i):
    """
    Run PTC continuation run for label 'i' in run_old.
    """
    from numpy import arange
    
    # This label
    this_label = label_old[i]

    # Run string identifier
    this_run = 'sol_' + str(i+1).zfill(3)

    # Print run information
    print('Continuing from point {} in run: {}'.format(this_label, run_old_str))
    
    # Set saved solutions for theta_old
    SP_points = arange(0.1, 2.0, 0.1)
    theta_old_stop = [0.0, 2.0]
    theta_new_stop = [-1.0, 3.0]

    # Run continuation
    run_scan = auto.run(run_old(this_label), LAB=1,
                        ICP=['theta_old', 'theta_new', 'eta', 'mu_s', 'T'],
                        UZSTOP={'theta_old': theta_old_stop, 'theta_new': theta_new_stop},
                        UZR={'theta_old': SP_points},
                        DSMIN=1e-2, DS=1e-1, DSMAX=1e0,
                        EPSL=1e-7, EPSU=1e-7, EPSS=1e-4,
                        NMX=8000, NPR=100)
    run_scan += auto.run(DS='-')
    
    #-------------------#
    #     Save Data     #
    #-------------------#
    # Append runs and save data
    run_scan = auto.relabel(run_scan)
    # run_scan = auto.merge(run_scan)
    
    # Save data
    data_funcs.save_move_data(run_scan, '{}/{}'.format(run_new_str, this_run))

    # Print new line
    print('\n')

#--------------------------------------------#
#     Run Continuation: Regular For Loop     #
#--------------------------------------------#
# Regular for loop run
for i in range(len(label_old)):
    # Run continuation
    calculate_PTC(i)

#--------------#
#     Plot     #
#--------------#
# Save data
import save_fig7_data as data_PTC
data_PTC.save_PTC_scan(run_new_str, '../data_files/fig7_data.mat')

# %%
#==============================================================================#
#                                 END OF FILE                                  #
#==============================================================================#
data_funcs.clean_directories()