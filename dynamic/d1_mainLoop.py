import numpy as np
from dynamic.sistem_enacb import sistem_enacb
#from scipy.integrate import solve_ivp
from SaveAndLoad.SaveLoad import save_dynamic_process
from i_magnetic.i_le_fullList import i_le_fullList
from i_magnetic.i_Ip import i_Ip_all

#Using Assimulo's IDA for DAE Problems
from assimulo.solvers import IDA
from assimulo.problem import Implicit_Problem



# Define the ODE-ADE system as a function compatible with Assimulo
def make_system_dae(U, V, hh, ne_c, tipverige, e_zp, Cdi, Qe_gg, UV_indices,
                    qc, qv, Ce, ni, nd, M, odvisneU, odvisneUdot,
                    cab, le, le_full, Ip):
    """
    Returns a residual function res(t, y, ydot) -> np.ndarray
    with all context captured via closure.
    """
    def system_dae(t, y, ydot):
        dz = sistem_enacb(
            t, y, U, V, hh, ne_c, tipverige, e_zp, Cdi, Qe_gg, UV_indices,
            qc, qv, Ce, ni, nd, M, odvisneU, odvisneUdot,
            cab, le, le_full, Ip
        )

        return ydot - dz
    return system_dae


    #Solving system - equations of motion

def d1_mainLoop(zacpog, U,V, hh, ne_c, tipverige, Cdi,Qe_g, Ce_active, M,
                    cab, le,filename):


    # Create index arrays
    ni = V.size  # number independant  variable
    nd = U.size  # number of dependant variables
    #i = np.arange(ni)  # indexing for coordinates

    # Initial conditions for coordinates
    neodvisneV = zacpog[V.flatten()]
    odvisneU = zacpog[U.flatten()]

    # Initial conditions for velocities
    neodvisneVdot = np.zeros(ni)  # everything is in static (velocity is zero)
    odvisneUdot = np.zeros(nd)  # same
    e_zp = zacpog.copy()  # this is used for constant determination, new value e_zp just for better clarity

    # Create initial ODE vector
    zacpog_ode = np.zeros(2 * ni)
    zacpog_ode[0::2] = neodvisneV  # inital condition, in every second row/column
    zacpog_ode[1::2] = neodvisneVdot


    # Reorder matrices
    UV_indices = np.concatenate([U.flatten(), V.flatten()])
    Ce = Ce_active[:, UV_indices]

    # Reorder load vectors
    Qe_gg = Qe_g[UV_indices]

    # Reorder mass matrix
    M = M[:, UV_indices]
    M = M[UV_indices, :]

    # Other variables, derivation of constraints
    qc = np.zeros(nd)  # this is zero because constraints are independent of time
    qv = 0  # in ANCF Cariolius force is zero

    # prepare (pre-process) array le_full which are used in i_Qem
    le_full = i_le_fullList(ne_c, le)
    Ip=i_Ip_all(le_full,cab)

    #  Compute initial derivatives for IDA
    #For a DAE solver like IDA, you need BOTH the state AND its derivative
    #IDA(..., y0=zacpog_ode, ydot0=zacpog_ode_dot, ...)  # Both needed!
    # These should satisfy ydot = f(t, y) at t=0
    zacpog_ode_dot = sistem_enacb(
        0,  # t=0
        zacpog_ode,  # Current state [positions; velocities]
        U, V, hh, ne_c, tipverige,
        e_zp, Cdi, Qe_gg,
        UV_indices,
        qc, qv, Ce, ni, nd, M,
        odvisneU, odvisneUdot,
        cab, le, le_full, Ip
    )

    res = make_system_dae(
        U, V, hh, ne_c, tipverige, e_zp, Cdi, Qe_gg, UV_indices,
        qc, qv, Ce, ni, nd, M, odvisneU, odvisneUdot,
        cab, le, le_full, Ip
    )



    # Returns: [velocities; accelerations]  which is exactly ydot!


    # ========================
    # MAIN SIMULATION LOOP
    # ========================
    print("Starting simulation...")
    t_o = 0.0

    # ---- build once ----
    # Create the DAE problem
    # #Assimulo calls  residual as res(t, y, ydot)
    # Note: IDA requires consistent initial conditions
    prob = Implicit_Problem(res,
                            zacpog_ode.flatten(),  # Initial state
                            zacpog_ode_dot.flatten(),
                            t0=t_o  # Assimulo's simulate always starts from solver.t0
                            )

    solver = IDA(prob)

    # Set solver options (matching your original settings)
    solver.rtol = 1e-5  # Relative tolerance
    solver.atol = 1e-5  # Absolute tolerance
    # solver.max_step = 0.005  # Uncomment if you need to limit step size
    # solver.use_jacobian = True  # Uncomment if you want to provide a Jacobian
    #solver.initialize() # do the initialization once, before any simulate call


    for uu in range(cab.ste_delitev):
        print(f"Interval {uu} of {cab.ste_delitev - 1}")

        # Create time points for output
        #num_steps = int(round(cab.interval / cab.note)) + 1
        #t_final   = (uu + 1) * cab.interval          # avoid FP drift

        # Simulate
        try:
            # t, Z = solver.simulate(t_o + cab.interval)
            #or if control over output number oint is desired
            #t, Z, Zd = solver.simulate(t_final, ncp=num_steps)  #t_final is the final time, initialization is the key that keep solver same over whole time, after end of loop keep settings and take them as input into new loop
            t, Z, Zd = solver.simulate(cab.tfinal,
                                       ncp=0)  # t_final is the final time, ncp=0 mean internal number of time steps as output!

            # t: time points array
            # Z = state values (positions + velocities)
            # Zd = state derivatives


        except Exception as e:
            print(f"Warning: ODE solver failed at interval {uu}")
            print(f"Error: {e}")
            break




        # Save results for each interval separately
        filename1 = f'{uu}e_process_{filename}.h5'

        save_dynamic_process(filename1, Z, t, cab.ste_delitev)  # save calculation results

        # No need to rebuild solver; it continues from its internal state
        #If you need final state after interval
        #y_final = solver.y  # numpy array, internal state
        #yd_final = solver.yd  # numpy array, internal ydot



        # removed in this version
        # Update for next interval,  Set new initial conditions = last row of Z, transposed to column vector
        #zacpog_ode = Z[-1, :].copy()  # Last state
        #zacpog_ode_dot = Zd[-1, :].copy()   # use the returned derivative directly
        #t_o += cab.interval

        #print(f"  Saved {filename1}")

    #print("Simulation completed!")
    return cab.ste_delitev  # Or return whatever you need



