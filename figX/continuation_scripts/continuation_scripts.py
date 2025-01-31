#==============================================================================#
#                FUNCTIONS USED IN THE PHASE RESET CALCULATIONS                #
#==============================================================================#
# This Python module contains the following functions:
#
# - calc_initial_solution_PO          : Reads the periodic orbit data from
#                                       the AUTO solution and writes it to
#                                       ./data/[XXX].dat.
#
# - calc_initial_solution_VAR         : Reads the data from the previous AUTO
#                                       solution, sets the initial conditions
#                                       and writes it to ./data/[XXX].dat.
#
# - save_floquet_data_matlab          : Saves the Floquet variational problem
#                                       as a MATLAB .mat file.
#
# - calc_initial_solution_PR          : Calculates and sets the initial
#                                       solution to solve for the phase
#                                       resetting problems, and write it
#                                       to ./data/[XXX].dat file.

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

    #-----------------------#
    #     Write to file     #
    #-----------------------#
    filename_data = './data_mat/initial_solution_PO.dat'
    with open(filename_data, 'w') as file_data:
        for i in range(len(t)):
            str_write = ("{:>20.15f} \t "
                         "{:>20.15f} \t {:>20.15f} \t {:>20.15f} \n"
                         ).format(t[i],
                                  x1[i], x2[i], x3[i])
            
            # Write to file
            file_data.write(str_write)

    # Close file
    file_data.close()

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

    #-----------------------#
    #     Write to file     #
    #-----------------------#
    filename_data = './data_mat/initial_solution_VAR.dat'
    with open(filename_data, 'w') as file_data:
        for i in range(len(t)):
            str_write = ("{:>20.15f} \t "
                         "{:>20.15f} \t {:>20.15f} \t {:>20.15f} \t"
                         "{:>20.15f} \t {:>20.15f} \t {:>20.15f} \n"
                         ).format(t[i],
                                  x1_read[i], x2_read[i], x3_read[i],
                                  w1[i], w2[i], w3[i])
            
            # Write to file
            file_data.write(str_write)

    # Close file
    file_data.close()

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