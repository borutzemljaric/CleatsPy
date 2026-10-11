import numpy as np



def d_open_menuDynamic(cab):
    """
    Add data about time and short circuit currents
    """
    cab.tfinal = 1  # (s) whole length of simulation
    cab.tsim = 0.05  #  (s) delay of simulation from start for better plots
    cab.trelay=0.5 # (s) length of short circuit - time required to  switch off
    cab.interval = cab.tfinal # (s) time length of sub-intervals, for easier memory management, must be smaller than tfinal, because ASSIMULO it is changed into one
    cab.ste_delitev = int(np.round(cab.tfinal / cab.interval))
    cab.note = 0.02  # (s ) interval of time for saved results

    #cab.Ik=9400  # (A) electric current Ik''
    cab.Ik = 11700  # (A) electric current Ik''
    cab.phi=89 # (degrees)surge impedance angle
    cab.Nseg=2 # number of sub-segment, minimum 2
    cab.Ncon=2 # number of contact points (checks) on element

    cab.a=30# # damping coefficient alpha
    cab.b=0.002# # damping coefficient beta
    return