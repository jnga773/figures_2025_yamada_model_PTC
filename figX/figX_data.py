# # Append Python path to add AUTO functions
# import sys
# sys.path.append('$HOME/auto/07p/python')
# sys.path.append('$HOME/auto/07p/python/auto')

# Load extra functions
import auto
import continuation_scripts.data_functions as data_funcs
import continuation_scripts.initial_PO_functions as initial_PO
import continuation_scripts.PTC_functions as phase_reset

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

# %%
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
# Copy continuation script
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
# SP = 2.5e-2

# Follow periodic orbits
print('Continuing periodic orbits ...')
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

# Write periodic orbit data to run06_initial_solution.dat
initial_PO.write_initial_solution_PO(run_new('UZ1'))

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
def read_parameters_PO(sol_in):
    # Read parameters
    gamma = sol_in['gamma']
    A     = sol_in['A']
    B     = sol_in['B']
    a     = sol_in['a']
    T     = sol_in['PAR(11)']


    # Parameter values
    p_out = {1: gamma, 2: A, 3: B, 4: a,
             5: T}
    # Parameter names
    pnames_out = {1: 'gamma', 2: 'A', 3: 'B', 4: 'a',
                  5: 'T'}
    
    return p_out, pnames_out

# Read parameters from previous run
par_PO, pnames_PO = read_parameters_PO(run_old(label_old))

#-------------------------------#
#     Run AUTO Continuation     #
#-------------------------------#
# Copy continuation script
auto.copy('./continuation_scripts/', 'initial_PO')

# Run continuation
run_new = auto.run(c='initial_PO', PAR=par_PO, parnames=pnames_PO)

# Save data
data_funcs.save_move_data(run_new, run_new_str)

# Write shifted periodic data and 'zero-ed' variational data as
# initial solution to the variational problem to ./data/run07_initial_solution.dat'
# initial_PO.write_initial_solution_floquet(run_new('EP1'))

#--------------#
#     Plot     #
#--------------#
# Save solution to MATLAB .mat file
# initial_PO.save_PO_data_matlab(run_new('EP1'))

# Print clear line
print('\n')

#==============================================================================#
##                 COMPUTE FLOQUET BUNDLE AT ZERO PHASE POINT                 ##
#==============================================================================#
# Here we compute the stable Floquet bundle of the periodic orbit, as well
# as the perpendicular vector, w.

# %%
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
label_old = run_old('EP1')
label_old = label_old['LAB']

# Print to console
print('~~~ Floquet Bundle: First Run (c.floquet_variational) ~~~')
print('Calculate Floquet bundle (mu) ')
print('Run name: {}'.format(run_new_str))
print('Continuing from point {} in run: {}'.format(label_old, run_old_str))

#-------------------#
#     Read Data     #
#-------------------#
# Define function for reading previous parameters
def read_parameters_VAR(sol_in):
    # Read parameters
    gamma = sol_in['gamma']
    A     = sol_in['A']
    B     = sol_in['B']
    a     = sol_in['a']
    T     = sol_in['T']

    # Initial parameter values
    mu_s  = 0.8
    wnorm = 0.0

    # Parameter values
    p_out = {1: gamma, 2: A, 3: B, 4: a,
             5: T,
             6: mu_s, 7: wnorm}
    # Parameter names
    pnames_out = {1: 'gamma', 2: 'A', 3: 'B', 4: 'a',
                  5: 'T',
                  6: 'mu_s', 7: 'w_norm'}
    
    return p_out, pnames_out

# Read parameters from previous run
par_VAR, pnames_VAR = read_parameters_VAR(run_old(label_old))

#-------------------------------#
#     Run AUTO Continuation     #
#-------------------------------#
# Copy continuation script
auto.copy('./continuation_scripts/', 'floquet_variational')

# Run continuation
run_new = auto.run(c='floquet_variational', PAR=par_VAR, parnames=pnames_VAR,
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
# Try run this bitch
run_new = auto.run(run_old(label_old), LAB=1,
                   ISW=-1, ILP=0,
                   DSMIN=5e-2, DS=1e-1, DSMAX=5e-1,
                   NMX=500, NPR=100,
                   UZSTOP={'w_norm': [0.0, 1.1]}, UZR={'w_norm': 1.0})

#-------------------#
#     Save Data     #
#-------------------#
# Save solution to MATLAB .mat file
initial_PO.save_floquet_data_matlab(run_new('UZ1'))

# Save data
data_funcs.save_move_data(run_new, run_new_str)

# Print clear line
print('\n')

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
run_new_str = 'run08_phase_reset_perturbation'

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
phase_reset.write_initial_solution_phase_reset(par_PR[6])

#-------------------------------#
#     Run AUTO Continuation     #
#-------------------------------#
# Set saved points
from numpy import linspace, concatenate, unique

# Saved points for large scan of I perturbation
SP_points = concatenate((linspace(0.0, 1.0, 20), 
                         linspace(1.0, 10.0, 15),
                         linspace(10.0, 13.0, 15),
                         linspace(13.0, 25.0, 20)))
SP_points = unique(SP_points)

# Copy continuation script
auto.copy('./continuation_scripts/', 'PTC_initial')

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
import continuation_scripts.phase_reset.save_PTC_scan_data as data_PTC
data_PTC.save_PTC_scan(run_new_str)

# %%
#==============================================================================#
#                                 END OF FILE                                  #
#==============================================================================#
data_funcs.clean_directories()
