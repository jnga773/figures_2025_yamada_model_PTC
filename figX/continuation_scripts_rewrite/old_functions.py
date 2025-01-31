#==============================================================================#
#                FUNCTIONS USED IN THE PHASE RESET CALCULATIONS                #
#==============================================================================#
# This Python module contains the following functions:
#
# - calc_initial_solution_phase_reset : Calculates and sets the initial
#                                       solution to solve for the phase
#                                       resetting problems, and write it
#                                       to ./data/[XXX].dat file.

#------------------------------------------------------------------------------#
# Write initial periodic orbit data to [FILE].dat to then be loaded in to
# the following continuation.
def write_initial_solution_PO(bd_data_in):
    """
    Reads the periodic orbit data from [bd_data_in] and writes to to
    [FILE].dat.
    """
    from numpy import argmax, concatenate
    
    #-------------------#
    #     Read Data     #
    #-------------------#
    # Grab the periodic orbit data file
    sol = bd_data_in

    # Read parameters
    gamma = sol['gamma']
    A     = sol['A']
    B     = sol['B']
    a     = sol['a']

    # Read time data
    t_read = sol['t']

    # Read data
    x1_read = sol['x1']
    x2_read = sol['x2']
    x3_read = sol['x3']

    # Calculate stationary points
    x0, xpos, xneg = calc_stationary_points([gamma, A, B, a])

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
                         "{:>20.15f} \t {:>20.15f} \t {:>20.15f} \t"
                         "{:>20.15f} \t {:>20.15f} \t {:>20.15f} \t"
                         "{:>20.15f} \t {:>20.15f} \t {:>20.15f} \t"
                         "{:>20.15f} \t {:>20.15f} \t {:>20.15f} \n"
                         ).format(t[i],
                                  x1[i], x2[i], x3[i],
                                  x0[0], x0[1], x0[2],
                                  xneg[0], xneg[1], xneg[2],
                                  xpos[0], xpos[1], xpos[2])
            
            # Write to file
            file_data.write(str_write)

    # Close file
    file_data.close()

#------------------------------------------------------------------------------#
# Calculate the initial solutions to the phase resetting problem.
def write_initial_solution_phase_reset(k_in):
    """
    Calculates and sets the initial solution to solve for the
    phase resetting problem and write it to a .dat file
    """
    from scipy.io import loadmat

    #---------------------------------------------------------------------------#
    #                                 Read Data                                 #
    #---------------------------------------------------------------------------#
    # Load data from .mat file
    data_floquet = loadmat('./data_mat/floquet_solution.mat')

    # Read time data
    t_read = data_floquet['t_read'].flatten()

    #--------------------#
    #     State data     #
    #--------------------#
    # Periodic orbit solution
    x1_read = data_floquet['x1_read'].flatten()
    x2_read = data_floquet['x2_read'].flatten()
    x3_read = data_floquet['x3_read'].flatten()

    # Initial point
    x1_0 = x1_read[0]
    x2_0 = x2_read[0]
    x3_0 = x3_read[0]

    #------------------------------#
    #     Perpindicular vector     #
    #------------------------------#
    wn1_read = data_floquet['wn1_read'].flatten()
    wn2_read = data_floquet['wn2_read'].flatten()
    wn3_read = data_floquet['wn3_read'].flatten()

    # Initial point
    wn1_0 = wn1_read[0]
    wn2_0 = wn2_read[0]
    wn3_0 = wn3_read[0]

    #---------------------------------------------------------------------------#
    #                              "Periodise" Data                             #
    #---------------------------------------------------------------------------#
    from numpy import concatenate

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

    #---------------------------------------------------------------------------#
    #                              Interpolate Data                             #
    #---------------------------------------------------------------------------#
    from numpy import interp, ones

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
            
    #---------------------------------------------------------------------------#
    #                                Write Data                                 #
    #---------------------------------------------------------------------------#
    from numpy import array
    
    # Time
    t_seg = t_period

    # Segment 1
    x_seg1 = array([x1_interp, x2_interp, x3_interp,
                    wn1_interp, wn2_interp, wn3_interp])
    x_seg1 = x_seg1.T

    # Segment 2
    x_seg2 = array([x1_extended, x2_extended, x3_extended,
                    wn1_extended, wn2_extended, wn3_extended])
    x_seg2 = x_seg2.T

    # Segment 3
    x_seg3 = array([x1_extended, x2_extended, x3_extended])
    x_seg3 = x_seg3.T

    # Segment 4
    x_seg4 = array([x1_period, x2_period, x3_period])
    x_seg4 = x_seg4.T

    filename_data = './data_mat/initial_solution_PR.dat'
    with open(filename_data, 'w') as file_data:
        for i in range(len(t_seg)):
            # Break up each column section into separate pieces

            # Time
            t_write    = "{:>15.10f}".format(t_seg[i])
            # Segment 1
            seg1_write = ("{:>15.10f} \t {:>15.10f} \t {:>15.10f} \t"
                          "{:>15.10f} \t {:>15.10f} \t {:>15.10f}"
                          ).format(x_seg1[i, 0], x_seg1[i, 1], x_seg1[i, 2],
                                  x_seg1[i, 3], x_seg1[i, 4], x_seg1[i, 5])
            # Segment 2
            seg2_write = ("{:>15.10f} \t {:>15.10f} \t {:>15.10f} \t"
                          "{:>15.10f} \t {:>15.10f} \t {:>15.10f}"
                          ).format(x_seg2[i, 0], x_seg2[i, 1], x_seg2[i, 2],
                                  x_seg2[i, 3], x_seg2[i, 4], x_seg2[i, 5])
            # Segment 3
            seg3_write = ("{:>15.10f} \t {:>15.10f} \t {:>15.10f}"
                          ).format(x_seg3[i, 0], x_seg3[i, 1], x_seg3[i, 2])
            # Segment 4
            seg4_write = ("{:>15.10f} \t {:>15.10f} \t {:>15.10f}"
                          ).format(x_seg4[i, 0], x_seg4[i, 1], x_seg4[i, 2])

            # Total string
            str_write = ("{} \t {} \t {} \t {} \t {} \n"
                         ).format(t_write, seg1_write, seg2_write, seg3_write, seg4_write)
            
            # Write to file
            file_data.write(str_write)

    # Close file
    file_data.close()

#------------------------------------------------------------------------------#
# Calculate the three stationary points
def calc_stationary_points(p_in):
    """
    Calculates the three stationary points of the Yamada model from analytic
    expression.
    

    Input
    ----------
    p_in : float, array
        Equation parameters.
    
    Output
    ----------
    x0_out : float, array
        The 'off' state
    xpos_out : float, array
        The higher-amplitude intensity state.
    xneg_out : float, array
        The lower-amplitude intensity state.  
    """
    from numpy import sqrt

    #---------------#
    #     Input     #
    #---------------#
    # Read the parameters
    gamma = p_in[0]
    A     = p_in[1]
    B     = p_in[2]
    a     = p_in[3]

    #--------------------------#
    #     Calculate Things     #
    #--------------------------#
    # Components for stationary points
    temp1 = -B - 1 - a + (a * A);
    temp1 = temp1 / (2 * a);

    temp2 = (B + 1 + a - (a * A)) ** 2
    temp2 = temp2 - 4 * (1 + B - A) * a
    temp2 = sqrt(temp2)
    temp2 = temp2 / (2 * a)

    # The two non-trivial equilibrium values of intensity
    I_pos = temp1 + temp2
    I_neg = temp1 - temp2

    #----------------#
    #     Output     #
    #----------------#
    x0_out   = [A, B, 0]
    xpos_out = [A / (1 + I_pos), B / (1 + (a * I_pos)), I_pos]
    xneg_out = [A / (1 + I_neg), B / (1 + (a * I_neg)), I_neg]

    return x0_out, xpos_out, xneg_out


#------------------------------------------------------------------------------#
# Calculate the initial solution to the adjoint BVP from
# the previous BVP run.
def write_initial_solution_floquet(bd_data_in):
    """
    Calculates and sets the initial solution to solve for the
    adjoint problem and write it to a .dat file
    """
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
                                  0.0, 0.0, 0.0)
            
            # Write to file
            file_data.write(str_write)

    # Close file
    file_data.close()


#------------------------------------------------------------------------------#
# Write Floquet data to .MAT file to read for phase reset program
def save_floquet_data_matlab(bd_data_in):
    """
    Writes the periodic orbit and Floquet variation problem data to a MATLAB
    .mat file/

    Input
    -------
    sol_in : AUTO-continuation solution
        AUTO generate solution for the periodic orbit (like run('UZ')).
        Probably run05_initial_PO('UZ1'), along with the stationary points.
    """
    from scipy.io import savemat
    from numpy import array

    #-------------------#
    #     Read Data     #
    #-------------------#
    
    # Grab the periodic orbit data file
    sol = bd_data_in

    # Read time data
    t_read = sol['t']

    #--------------------#
    #     State data     #
    #--------------------#
    # Periodic orbit solution
    x1_read = sol['x1']
    x2_read = sol['x2']
    x3_read = sol['x3']

    #------------------------------#
    #     Perpindicular vector     #
    #------------------------------#
    wn1_read = sol['wn_1']
    wn2_read = sol['wn_2']
    wn3_read = sol['wn_3']

    #--------------------#
    #     Parameters     #
    #--------------------#
    # Read parameters
    gamma = sol['gamma']
    A     = sol['A']
    B     = sol['B']
    a     = sol['a']
    T     = sol['T']
    mu_s  = sol['mu_s']

    #-------------------#
    #     Save Data     #
    #-------------------#
    # Dictionary for data
    mat_out = {'t_read': t_read,
               'x1_read': x1_read, 'x2_read': x2_read, 'x3_read': x3_read,
               'wn1_read': wn1_read, 'wn2_read': wn2_read, 'wn3_read': wn3_read,
               'gamma': gamma, 'A': A, 'B': B, 'a': a, 'T': T, 'mu_s': mu_s}
    
    # Save data
    savemat('./data_mat/floquet_solution.mat', mat_out, oned_as='column') 
