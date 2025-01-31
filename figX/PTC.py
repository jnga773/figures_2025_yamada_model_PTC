# # Append Python path to add AUTO functions
# import sys
# sys.path.append('/Users/jnga773/auto/07p/python')
# sys.path.append('/Users/jnga773/auto/07p/python/auto')

# Load extra functions
import auto
import continuation_scripts.data_functions as data_funcs
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

# Print to console
print('~~~ Phase Reset: First Run ~~~')
print('Continue in the perturbation ampltiude A_perturb')
print('Run name: {}'.format(run_new_str))

#-------------------#
#     Read Data     #
#-------------------#
# Define function for reading previous parameters and setting phase resetting
# parameters
def set_parameters_PR():
    from scipy.io import loadmat
    # Read Floquet data
    data_floquet = loadmat('./data_mat/floquet_solution.mat')
    
    # Read parameters
    gamma = data_floquet['gamma'].item()
    A     = data_floquet['A'].item()
    B     = data_floquet['B'].item()
    a     = data_floquet['a'].item()
    T     = data_floquet['T'].item()
    mu_s  = data_floquet['mu_s'].item()

    #----------------------------#
    #     Initial Parameters     #
    #----------------------------#
    from numpy import pi

    # Integer for period
    k             = 50
    # \theta_old (where perturbation starts)
    theta_old     = 1.0
    # \theta_new (where segment comes back to \Gamma)
    theta_new     = 1.0
    # Distance from perturbed segment to \Gamma
    eta           = 0.0
    # Size of perturbation
    A_perturb     = 0.0
    # Angle at which perturbation is applied?
    # theta_perturb = 0.0
    theta_perturb = 0.5 * pi
    # Azimuthal angle at which perturbation is applied?
    phi_perturb   = 0.0

    #----------------#
    #     Output     #
    #----------------#
    # Parameter vector
    p_out = { 1: gamma, 2: A, 3: B, 4: a,
              5: T, 6: k, 7: theta_old, 8: theta_new,
              9: mu_s, 10: eta,
             11: A_perturb, 12: theta_perturb, 13: phi_perturb}
    
    # Parameter names
    pnames_out = { 1: 'gamma', 2: 'A', 3: 'B', 4: 'a',
                   5: 'T', 6: 'k', 7: 'theta_old', 8: 'theta_new',
                   9: 'mu_s', 10: 'eta',
                  11: 'A_perturb', 12: 'theta_perturb', 13: 'phi_perturb'}
    
    return p_out, pnames_out

# Read parameters from previous run
par_PR, pnames_PR = set_parameters_PR()

# Calculate and write initial solution from previous run
old_functions.write_initial_solution_phase_reset(par_PR[6])

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
run_new = auto.run(c='PTC_initial', PAR=par_PR, parnames=pnames_PR,
                   NMX=2000,
                   NTST=par_PR[6] * 60,
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
