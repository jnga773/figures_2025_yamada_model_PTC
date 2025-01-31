# # Append Python path to add AUTO functions
# import sys
# sys.path.append('/Users/jnga773/auto/07p/python')
# sys.path.append('/Users/jnga773/auto/07p/python/auto')

# Load extra functions
import auto
import data_functions as data_funcs
import continuation_scripts.old_functions as old_functions

#==============================================================================#
##                  PHASE RESPONSE - PHASE TRANSITION CURVES                  ##
#==============================================================================#
# We set up the phase resetting problem by creating four segments of the
# periodic orbit, with boundary conditions described in a paper somewhere.

# %%
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
run_new_str = 'run01_phase_reset_perturbation'
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
from numpy import pi
k = 25
theta_perturb = 0.5 * pi


# Calculate initial solution
x_init_PR, p_PR, pnames_PR = \
    data_funcs.calc_initial_solution_PR(run_old(label_old), k, theta_perturb)

#-------------------------------#
#     Run AUTO Continuation     #
#-------------------------------#
# Set saved points
from numpy import linspace, concatenate, arange

# Saved points for large scan of G perturbation
SP_points = concatenate((linspace(0.0, 0.25, 25), linspace(0.30, 2.0, 25)))

# # Saved points for large scan of I perturbation
# SP_points = concatenate((linspace(0.0, 1.0, 35), linspace(0.30, 25.0, 35)))


# SP_points = [0.0, 0.025, 0.05, 0.075, 0.1, 0.2, 0.3, 0.4, 0.5,
#              1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
# SP_points = [0.0, 0.025, 0.05, 0.075, 0.1, 0.2, 0.3, 0.4, 0.5,
#              1.0, 1.2, 1.4, 1.6, 1.8, 2.0,
#              3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
#              12.0, 14.0, 16.0, 18.0, 20.0]

# Copy continuation script
auto.copy('./constant_files/', 'PTC_initial')

# Try set up phase reset calculation lol
# run_new = auto.run(x_init_PR, PAR=p_PR, parnames=pnames_PR,
#                    c='PTC_initial',
#                    NMX=2000, NTST=p_PR[6] * 50,
#                    UZR={'A_perturb': SP_points},
#                    UZSTOP={'A_perturb': max(SP_points) + 0.1})
run_new = auto.run(dat='./data_mat/initial_solution_PR.dat', PAR=p_PR, parnames=pnames_PR,
                   c='PTC_initial',
                   NMX=2000, NTST=p_PR[6] * 50,
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
#==============================================================================#
#                                 END OF FILE                                  #
#==============================================================================#
data_funcs.clean_directories()
