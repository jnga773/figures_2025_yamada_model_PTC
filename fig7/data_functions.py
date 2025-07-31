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
# - yamada                            : Set of ODEs for the Yamada model.
#
# - calc_initial_solution_solve_ivp   : Calculates the periodic orbit using
#                                       ode45.
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
#
# - run_PR_continuation               : Runs the PTC continuation
#                                       for the Yamada model.

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
        
        # Remove initial solution .dat files
        if file.endswith('.dat'):
            remove('./{}'.format(file))


    # Remove things inside ./functions/ folder
    function_files = listdir('./functions/')
    for file in function_files:
        if file.endswith('.o') or file.endswith('.exe'):
            remove('./functions/{}'.format(file))

    print("Removed fort.*, *.o, *.exe, c.*, *.dat, and *.mod files")
    
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
def yamada(t_in, x_in, p_in):
  """
  Set of ODEs for the Yamada model.
  """
  # State space variables
  G = x_in[0]
  Q = x_in[1]
  I = x_in[2]

  # Parameters
  gamma, A, B, a = p_in

  # ODEs
  dG = gamma * (A - G - (G * I))
  dQ = gamma * (B - Q - (a * Q * I))
  dI = (G - Q - 1) * I

  return [dG, dQ, dI]

#------------------------------------------------------------------------------#
def calc_initial_solution_solve_ivp(x0_in, p0_in):
  """
  Calculate the periodic orbit using ode45.
  """
  from scipy.integrate import solve_ivp
  from numpy import linspace, array
  
  # Time span for the integration
  t_span = (0, 10000)  # Adjust as needed for the period

  # Solve the ODEs
  sol = solve_ivp(yamada, t_span, x0_in, args=(list(p0_in.values()),),
                  method='RK45', t_eval=None, dense_output=True)
  
  # Get final state
  x_long = sol.y[:, -1]
  
  # Guess period
  T_PO = 41
  
  # Solve again for close to a full period
  sol = solve_ivp(yamada, (0, T_PO), x_long, args=(list(p0_in.values()),),
                  method='RK45', dense_output=True)

  # Return the solution
  t_out = linspace(0, T_PO, 1000)
  x_out = sol.sol(t_out)
  
  sol_out = [t_out, x_out[0, :], x_out[1, :], x_out[2, :]]
  
  return sol_out

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
    T     = sol_in['T']

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

    # State-space variable names
    unames_out = {1: 'x1', 2: 'x2', 3: 'x3',
                  4: 'wn_1', 5: 'wn_2', 6: 'wn_3'}

    # Parameter values
    p_out = {1: gamma, 2: A, 3: B, 4: a,
             5: T,
             6: mu_s, 7: wnorm}
    
    # Parameter names
    pnames_out = {1: 'gamma', 2: 'A', 3: 'B', 4: 'a',
                  5: 'T',
                  6: 'mu_s', 7: 'w_norm'}

    return x_init_out, p_out, unames_out, pnames_out

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

    filename_data = './initial_solution_PR.dat'
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

    # State space variable names
    unames_out = { 1: 'seg1_x1',  2: 'seg1_x2',  3: 'seg1_x3',
                   4: 'seg1_w1',  5: 'seg1_w2',  6: 'seg1_w3',
                   7: 'seg2_x1',  8: 'seg2_x2',  9: 'seg2_x3',
                  10: 'seg2_w1', 11: 'seg2_w2', 12: 'seg2_w3',
                  13: 'seg3_x1', 14: 'seg3_x2', 15: 'seg3_x3',
                  16: 'seg4_x1', 17: 'seg4_x2', 18: 'seg4_x3'}

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
    
    return x_init_out, p_out, unames_out, pnames_out

#------------------------------------------------------------------------------#
def run_PR_continuation(run_new_str, run_old, label_old, pcont, prange, **kwargs):
    """
    Scan through SP labels from previous run (different values of A_perturb)
    and continue in \theta_{old} and \theta_{new}. Each run will save
    to the 'data/run10_phase_reset_PTC/' directory.

    Parameters
    ----------
    run_new_str : string
        The new run identifier for the main continuation problem.
    run_old : string
        The old run identifier for the sub continuation problem.
    label_old : integer
        The label identifier for the previous continuation problem.
    data_PR : struct
        Data structure containing the initial conditions for the trajectory segments.
    bcs_funcs : list of functions
        Structure containing boundary condition functions.
    pcont : cell
        Cell array containing additional parameters for the continuation.
    prange : cell
        Cell array containing the ranges for the continuation parameters.

    Optional Arguments
    ------------------
    UZ : dict
        Dictionary containing saved point parameters for the continuation.
        Default is an empty dictionary.
    DSMIN : float
        Minimum step size for the continuation. Default is 1e-4.
    DS : float
        Step size for the continuation. Default is 1e-3.
    DSMAX : float
        Maximum step size for the continuation. Default is 1e0.
    NMX : int
        Maximum number of continuation steps. Default is 2000.
    NPR : int
        Number of points to print during the continuation. Default is 100.
    IAD : int
        Interval for the continuation. Default is 10.
    reverse: logical
        Check if you want to run the continuation in the opposite direction
        with DS='-'.
    """
    from auto import run, merge, relabel
    # from data_functions import save_move_data

    #-------------------#
    #     Arguments     #
    #-------------------#
    # Default options
    options = {
        'UZR'    : {},
        'DSMIN'  : 1e-4,
        'DS'     : 1e-2,
        'DSMAX'  : 1e1,
        'NMX'    : 2000,
        'NPR'    : 100,
        'IAD'    : 10,
        'reverse': False
    }
    # Override defaults with any provided kwargs
    options.update(kwargs)

    #------------------------#
    #     Run Contination    #
    #------------------------#

    # Run continuation
    run_scan = run(run_old, LAB=1, IRS=label_old,
                   ICP=pcont, UZSTOP=prange, UZR=options['UZR'],
                   DSMIN=options['DSMIN'], DS=options['DS'], DSMAX=options['DSMAX'],
                   NMX=options['NMX'], NPR=options['NPR'],
                   IAD=options['IAD'])
    
    if options['reverse']:
        run_scan += run(DS='-')
    
    #-------------------#
    #     Save Data     #
    #-------------------#
    # Append runs and save data
    run_scan = merge(run_scan)
    run_scan = relabel(run_scan)
    
    # Save data
    save_move_data(run_scan, '{}'.format(run_new_str))

    # Print new line
    print('\n')

#------------------------------------------------------------------------------#