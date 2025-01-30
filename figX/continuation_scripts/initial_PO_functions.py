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

#------------------------------------------------------------------------------#
# Calculate new initial solution to rotated periodic orbit solution
def calc_initial_solution_PO(bd_data_in):
    """
    Reads the periodic orbit data from [bd_data_in], and rotates it to the
    new head point.
    """
    from numpy import argmax, concatenate, array
    
    #-------------------#
    #     Read Data     #
    #-------------------#
    # Grab the periodic orbit data file
    sol = bd_data_in

    # Read time data
    t_read = sol['t']

    # Read data
    x1_read = sol['x1']
    x2_read = sol['x2']
    x3_read = sol['x3']

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
    x_init_out = [t, x1, x2, x3]
    x_init_out = array(x_init_out, dtype='float')

    return x_init_out

#------------------------------------------------------------------------------#
# Calculate the initial solution to the adjoint BVP from
# the previous BVP run.
def calc_initial_solution_VAR(bd_data_in):
    """
    Calculates and sets the initial solution to solve for the
    adjoint problem.
    """
    from numpy import zeros, array

    #-------------------#
    #     Read Data     #
    #-------------------#
    # Grab the periodic orbit data file
    sol = bd_data_in

    # Read time data
    t = sol['t']

    # Read data
    x1_read = sol['x1']
    x2_read = sol['x2']
    x3_read = sol['x3']

    # Zeros for perpindicular solution
    w1 = zeros(len(x1_read))
    w2 = zeros(len(x1_read))
    w3 = zeros(len(x1_read))

    #----------------#
    #     Output     #
    #----------------#
    x_init_out = [t, x1_read, x2_read, x3_read, w1, w2, w3]
    x_init_out = array(x_init_out, dtype='float')

    return x_init_out

#------------------------------------------------------------------------------#
# Write Floquet data to .dat file to test in Spyder
def write_floquet_data(bd_data_in):
    """
    Read the final solution from run07_floquet_wnorm.
    """
    #-------------------#
    #     Read Data     #
    #-------------------#
    # Grab the periodic orbit data file
    sol = bd_data_in

    # Read time data
    t_read = sol['t']

    # Periodic orbit solution
    x1_read = sol['x1']
    x2_read = sol['x2']
    x3_read = sol['x3']

    # Perpindicular vector solution
    wn_1_read = sol['wn_1']
    wn_2_read = sol['wn_2']
    wn_3_read = sol['wn_3']

    #-----------------------#
    #     Write to file     #
    #-----------------------#
    filename_data = './data/VAR_data_test.dat'
    with open(filename_data, 'w') as file_data:
        for i in range(len(t_read)):
            str_write = ("{:>20.15f} \t "
                         "{:>20.15f} \t {:>20.15f} \t {:>20.15f} \t"
                         "{:>20.15f} \t {:>20.15f} \t {:>20.15f} \n"
                         ).format(t_read[i],
                                  x1_read[i], x2_read[i], x3_read[i],
                                  wn_1_read[i], wn_2_read[i], wn_3_read[i])
            
            # Write to file
            file_data.write(str_write)

    # Close file
    file_data.close()
