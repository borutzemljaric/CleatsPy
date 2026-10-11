import numpy as np
from rutines.resfor_Re import Ke_vect
from dynamic.sistem_enacb import pdotdot_lagrange
from rutines.odvisneNeodvisne import odvi_U
from i_magnetic.i_Qem import i_Qem
from i_magnetic.i_le_fullList import i_le_fullList
from i_damping.Qd_vec import Qd_vect
from j_contact.j_Qcon import j_Qcon
from i_magnetic.i_Ip import i_Ip_all

def systemEquationPostProcess(t, Z, zacpog, e_a, ne_c, V, U, le,
                                 Qeg_0, M, tipverige, hh, Cdi, Ce_active,cab):
    """
    to produce all coordinates, dependant and independent in original order and for each time step
    """
    # Reconstruct
    le_full = i_le_fullList(ne_c, le)
    Ip = i_Ip_all(le_full, cab)

    # Extract independent coordinates and velocities
    neodvisneV = Z[0::2 ]  #
    neodvisneVdot = Z[1::2]  #

    # Dependent coordinates from previous step
    odvisneU=e_a[U] # prepišem odvisne koordinate, vzamem zadnje izračunane iz predhodnega koraka

    # Update dependent coordinates
    odvisneU = odvi_U(odvisneU, hh, ne_c, neodvisneV, V, U, tipverige, zacpog)
    # e_a[U.flatten()] = odvisneU

    # NAZAJ SESTAVIM VEKTOR e_a po ponovnem izračunu odvisnih koordinat
    #odvisne na novo prirejene tem neodvisnim
    e_a_new=np.zeros_like(e_a)
    e_a_new[U] = odvisneU
    e_a_new[V] = neodvisneV
    #%vrne izračunane koordinate v osnovnem vrstnem redu

    # %% sedaj izracunam se HITROSTI
    #derivation Ct and Ctt are equal zero

    # Update dependent velocities
    odvisneUdot = Cdi @ neodvisneVdot


    # velocity vector- back to original index order
    e_dot = np.zeros_like(e_a)
    e_dot[U]=odvisneUdot
    e_dot[V]=neodvisneVdot

    # Calculate qc and qv
    #nd = U.size
    qc = np.zeros(len(U))
    qv = 0  # in ANCF, Coriolis and centrifugal forces are zero

    # Determine load during integration

    # Calculate forces
    #Qeg_0 = Qe_g.copy()  # just to define
    # Qeg_0 = load_over_time(t, timeload, Qe, Qe_g, Qe_breme, Qe_gg)

    # note- because there is no need to include velocity( there in no damp or relative wind speed) the same function as in static can be used for qkqt forces
    qkqt = Ke_vect(hh, ne_c, e_a_new, cab, le)
    Qd = Qd_vect(hh, ne_c, cab, le, e_a, e_dot)
    Qem = i_Qem(t, e_a_new, cab, le_full,Ip)
    Qcon = j_Qcon(t, e_a, e_dot, cab, le_full)
    qkqt = qkqt + Qd - Qem - Qcon

    # Calculate accelerations
    e_dotdot, lambda_vec = pdotdot_lagrange(M, Qeg_0, qkqt, qv, qc, Ce_active)

    Qc_izhod = Ce_active.T @ lambda_vec

    return e_a_new, e_dot, e_dotdot, Qc_izhod