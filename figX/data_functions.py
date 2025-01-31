#==============================================================================#
#                    SELF-DEFINED FUNCTIONS FOR AUTO THINGS                    #
#==============================================================================#
# This Python module contains the following functions:
#
# - calc_stationary_points            : Calculates the three stationary
#                                       points of the Yamada model.
#
# - clean_directories                 : Deletes all of the AUTO generated
#                                       files that are yuck as.
#
# - save_move_data                    : Saves the data from run [RUN_IN] to
#                                       ./data/RUN_NAME_IN/.
#
# - bd_read                           : Loads the bifucation and solution data
#                                       in ./data_run_name_in/XXX.run_name_in.
#
# - calc_initial_solution_PO          : Reads the periodic orbit data from
#                                       the AUTO solution, rotates it, and sets
#                                       it as the initial solution.
#
# - calc_initial_solution_VAR         : Reads the data from the previous AUTO
#                                       solution, and sets the initial
#                                       conditions.
#
# - calc_initial_solution_PR          : Calculates and sets the initial
#                                       solution to solve for the phase
#                                       resetting problems, and write it
#                                       to ./data/[XXX].dat file.

#------------------------------------------------------------------------------#
def clean_directories():
    """
    Cleans up fort files, .o and .exe files in ./functions/, and the
    GLOBAL_VARIABLE .o file.
    """
    from os import remove, listdir

    # Cycle through files in current directory
    function_files = listdir('./')
    for file in function_files:
        # Remove c.XXX continuation scripts
        if file.startswith('c.'):
            remove('./{}'.format(file))
        
        # Remove Fortran .mod file
        if file.endswith('.mod'):
            remove('./{}'.format(file))
        
        # Remove fort.* data files
        if file.startswith('fort.'):
            remove('./{}'.format(file))


    # Remove things inside ./functions/ folder
    function_files = listdir('./functions/')
    for file in function_files:
        if file.endswith('.o') or file.endswith('.exe'):
            remove('./functions/{}'.format(file))

    print("Removed fort.*, *.o, *.exe, c.*, and *.mod files")
    
#------------------------------------------------------------------------------#
def save_move_data(run_in, run_name_in):
    """
    Saves the data from run [RUN_IN] to './data/RUN_NAME_IN/'.
    To make things easy, make sure that run_in and run_name_in are
    the same thing, except run_name_in will have '' around it :)
    """
    from os.path import isdir
    from os import makedirs
    import auto

    # Check for './data/' directory
    if not isdir('./data'):
        makedirs('./data')

    # Check for './data/run_name_in/' directory
    if not isdir('./data/{}'.format(run_name_in)):
        makedirs('./data/{}'.format(run_name_in))

    # Save the data
    auto.save(run_in, 'dat')

    # Move the data to the './data/run_name_in/' directory
    auto.move('dat', './data/{}'.format(run_name_in))

#------------------------------------------------------------------------------#
def bd_read(run_name_in):
    """
    Loads the bifurcation and solution data in ./data/run_name_in/XXX.run_name_in.
    """
    from os import remove, chdir
    from os.path import abspath
    import auto

    # Get current directory complete path
    main_dir = abspath('./')

    # Move into data directory
    chdir('./data/{}/'.format(run_name_in))

    # Load it
    bd_out = auto.loadbd('dat')

    # Change back to working directory
    chdir(main_dir)

    return bd_out

#------------------------------------------------------------------------------#
def calc_initial_solution_PO(sol_in):
    """
    Reads the periodic orbit data from [sol_in], rotates the periodic orbit,
    and writes the initial time and state space solution to
    './data_mat/initial_solution_PO.dat'.
    """
    from numpy import argmax, concatenate, array
    
    #-------------------#
    #     Read Data     #
    #-------------------#
    # Read time data
    t_read = sol_in['t']

    # Read data
    x1_read = sol_in['x1']
    x2_read = sol_in['x2']
    x3_read = sol_in['x3']

    #-------------------------#
    #     Read Parameters     #
    #-------------------------#
    # Read parameters
    gamma = sol_in['gamma']
    A     = sol_in['A']
    B     = sol_in['B']
    a     = sol_in['a']
    T     = sol_in['PAR(11)']

    #--------------------#
    #     Shift Data     #
    #--------------------#
    # Find the point where the first component is the maximum
    max_idx = argmax(x1_read)

    # Shift data
    if max_idx > 0:
        # Shift around arrays
        x1 = [x1_read[max_idx:], x1_read[1:max_idx+1]]
        x2 = [x2_read[max_idx:], x2_read[1:max_idx+1]]
        x3 = [x3_read[max_idx:], x3_read[1:max_idx+1]]

        # Shift time array
        t = [t_read[max_idx:] - t_read[max_idx],
             t_read[1:max_idx+1] + (t_read[-1] - t_read[max_idx])]

        # Concatenate arrays
        x1 = concatenate((x1[0], x1[1]))
        x2 = concatenate((x2[0], x2[1]))
        x3 = concatenate((x3[0], x3[1]))
        t = concatenate((t[0], t[1]))

    #----------------#
    #     Output     #
    #----------------#
    # State space solution
    x_init_out = [t, x1, x2, x3]
    x_init_out = array(x_init_out)

    # Parameter values
    p_out = {1: gamma, 2: A, 3: B, 4: a,
             5: T}
    # Parameter names
    pnames_out = {1: 'gamma', 2: 'A', 3: 'B', 4: 'a',
                  5: 'T'}

    return x_init_out, p_out, pnames_out


#------------------------------------------------------------------------------#
# Calculate the initial solution to the adjoint BVP from
# the previous BVP run.
def calc_initial_solution_VAR(sol_in):
    """
    Calculates and sets the initial solution to solve for the
    adjoint problem and write it to a .dat file
    """
    from numpy import zeros, array
    #-------------------#
    #     Read Data     #
    #-------------------#
    # Read time data
    t = sol_in['t']

    # Read data
    x1_read = sol_in['x1']
    x2_read = sol_in['x2']
    x3_read = sol_in['x3']

    # Zeros for perpindicular solution
    w1 = zeros(shape=(len(x1_read)))
    w2 = zeros(shape=(len(x1_read)))
    w3 = zeros(shape=(len(x1_read)))

    #-------------------------#
    #     Read Parameters     #
    #-------------------------#
    # Read parameters
    gamma = sol_in['gamma']
    A     = sol_in['A']
    B     = sol_in['B']
    a     = sol_in['a']
    T     = sol_in['T']

    # Initial parameter values
    mu_s  = 0.8
    wnorm = 0.0

    #----------------#
    #     Output     #
    #----------------#
    # State space solution
    x_init_out = [t, x1_read, x2_read, x3_read, w1, w2, w3]
    x_init_out = array(x_init_out, dtype='float')

    # Parameter values
    p_out = {1: gamma, 2: A, 3: B, 4: a,
             5: T,
             6: mu_s, 7: wnorm}
    # Parameter names
    pnames_out = {1: 'gamma', 2: 'A', 3: 'B', 4: 'a',
                  5: 'T',
                  6: 'mu_s', 7: 'w_norm'}

    return x_init_out, p_out, pnames_out

#------------------------------------------------------------------------------#
# Calculate the initial solutions to the phase resetting problem.
def calc_initial_solution_PR(sol_in, k_in, theta_in):
    """
    Reads solution from sol_in, and calculates initial solution to the phase
    resetting problem.
    """
    from numpy import pi, concatenate, interp, ones, array
    #-------------------#
    #     Read Data     #
    #-------------------#
    # Read time data
    t_read = sol_in['t']

    # Periodic orbit solution
    x1_read = sol_in['x1']
    x2_read = sol_in['x2']
    x3_read = sol_in['x3']

    # Perpindicular vector solution
    wn1_read = sol_in['wn_1']
    wn2_read = sol_in['wn_2']
    wn3_read = sol_in['wn_3']

    # Initial points
    x1_0 = x1_read[0]
    x2_0 = x2_read[0]
    x3_0 = x3_read[0]
    wn1_0 = wn1_read[0]
    wn2_0 = wn2_read[0]
    wn3_0 = wn3_read[0]

    #-------------------------#
    #     Read Parameters     #
    #-------------------------#
    # Read parameters
    gamma = sol_in['gamma']
    A     = sol_in['A']
    B     = sol_in['B']
    a     = sol_in['a']
    T     = sol_in['T']
    mu_s  = sol_in['mu_s']

    #----------------------------#
    #     Initial Parameters     #
    #----------------------------#
    # Integer for period
    k             = k_in
    # \theta_old (where perturbation starts)
    theta_old     = 1.0
    # \theta_new (where segment comes back to \Gamma)
    theta_new     = 1.0
    # Distance from perturbed segment to \Gamma
    eta           = 0.0
    # Size of perturbation
    A_perturb     = 0.0
    # Angle at which perturbation is applied?
    theta_perturb = theta_in
    # Azimuthal angle at which perturbation is applied?
    phi_perturb   = 0.0

    #--------------------------#
    #     "Periodise" Data     #
    #--------------------------#
    # Set initial time array
    t_period  = t_read

    # Set inintial state date
    x1_period = x1_read
    x2_period = x2_read
    x3_period = x3_read

    # If k > 1, cycle through and keep appending solutions
    if k_in > 1:
        # Cycle through number of periods
        for i in range(k_in-1):
            # Append time
            t_period = concatenate((t_read, 1 + t_period[1:]))
            
            # Segment 4: State space data
            x1_period = concatenate((x1_read, x1_period[1:]))
            x2_period = concatenate((x2_read, x2_period[1:]))
            x3_period = concatenate((x3_read, x3_period[1:]))
            
        # Normalise time data
        t_period *= 1 / k_in

    #--------------------------#
    #     Interpolate Data     #
    #--------------------------#
    # Interpolate: State space
    x1_interp  = interp(t_period, t_read, x1_read)
    x2_interp  = interp(t_period, t_read, x2_read)
    x3_interp  = interp(t_period, t_read, x3_read)

    # Interpolate: Perpindicular space
    wn1_interp = interp(t_period, t_read, wn1_read)
    wn2_interp = interp(t_period, t_read, wn2_read)
    wn3_interp = interp(t_period, t_read, wn3_read)

    # New ones matrix
    ones_mat = ones(shape=(len(t_period)))

    # Multiply through single values
    x1_extended  = x1_0 * ones_mat
    x2_extended  = x2_0 * ones_mat
    x3_extended  = x3_0 * ones_mat

    wn1_extended = wn1_0 * ones_mat
    wn2_extended = wn2_0 * ones_mat
    wn3_extended = wn3_0 * ones_mat
    
    #--------------------#
    #     Write Data     #
    #--------------------#
    # Time
    t_seg = t_period

    # Segment 1
    x_seg1 = [x1_interp, x2_interp, x3_interp,
              wn1_interp, wn2_interp, wn3_interp]
    x_seg1 = array(x_seg1).T

    # Segment 2
    x_seg2 = [x1_extended, x2_extended, x3_extended,
              wn1_extended, wn2_extended, wn3_extended]
    x_seg2 = array(x_seg2).T

    # Segment 3
    x_seg3 = [x1_extended, x2_extended, x3_extended]
    x_seg3 = array(x_seg3).T

    # Segment 4
    x_seg4 = [x1_period, x2_period, x3_period]
    x_seg4 = array(x_seg4).T

    # filename_data = './data_mat/initial_solution_PR.dat'
    # with open(filename_data, 'w') as file_data:
    #     for i in range(len(t_seg)):
    #         # Break up each column section into separate pieces

    #         # Time
    #         t_write    = "{:>15.10f}".format(t_seg[i])
    #         # Segment 1
    #         seg1_write = ("{:>15.10f} \t {:>15.10f} \t {:>15.10f} \t"
    #                       "{:>15.10f} \t {:>15.10f} \t {:>15.10f}"
    #                       ).format(x_seg1[i, 0], x_seg1[i, 1], x_seg1[i, 2],
    #                                x_seg1[i, 3], x_seg1[i, 4], x_seg1[i, 5])
    #         # Segment 2
    #         seg2_write = ("{:>15.10f} \t {:>15.10f} \t {:>15.10f} \t"
    #                       "{:>15.10f} \t {:>15.10f} \t {:>15.10f}"
    #                       ).format(x_seg2[i, 0], x_seg2[i, 1], x_seg2[i, 2],
    #                                x_seg2[i, 3], x_seg2[i, 4], x_seg2[i, 5])
    #         # Segment 3
    #         seg3_write = ("{:>15.10f} \t {:>15.10f} \t {:>15.10f}"
    #                       ).format(x_seg3[i, 0], x_seg3[i, 1], x_seg3[i, 2])
    #         # Segment 4
    #         seg4_write = ("{:>15.10f} \t {:>15.10f} \t {:>15.10f}"
    #                       ).format(x_seg4[i, 0], x_seg4[i, 1], x_seg4[i, 2])

    #         # Total string
    #         str_write = ("{} \t {} \t {} \t {} \t {} \n"
    #                      ).format(t_write, seg1_write, seg2_write, seg3_write, seg4_write)
            
    #         # Write to file
    #         file_data.write(str_write)

    # # Close file
    # file_data.close()

    #----------------#
    #     Output     #
    #----------------#
    # Initial vector
    x_init_out = [t_period,
                  # Segment '1'
                  x1_interp, x2_interp, x3_interp,
                  wn1_interp, wn2_interp, wn3_interp,
                  # Segment '2'
                  x1_extended, x2_extended, x3_extended, 
                  wn1_extended, wn2_extended, wn3_extended,
                  # Segment '3'
                  x1_extended, x2_extended, x3_extended,
                  # Segment '4'
                  x1_period, x2_period, x3_period]
    x_init_out = array(x_init_out, dtype='float')

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
    
    return x_init_out, p_out, pnames_out