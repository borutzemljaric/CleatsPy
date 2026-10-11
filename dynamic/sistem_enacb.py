import numpy as np
from rutines.odvisneNeodvisne import odvi_U
from rutines.resfor_Re import Ke_vect
from dynamic.pdotdot_lagrange import pdotdot_lagrange
from i_magnetic.i_Qem import i_Qem
from i_damping.Qd_vec import Qd_vect
from j_contact.j_Qcon import j_Qcon



def sistem_enacb(t, z,U,V, hh, ne_c, tipverige, e_zp, Cdi,Qe_g, UV_indices, qc,qv, Ce, ni, nd,  M,
                   odvisneU,  odvisneUdot,  cab, le, le_full,Ip):
    """
    solver for time step
    """


    # Reconstruct coordinates
    e_a = np.zeros(len(U) + len(V))  # predefine vector
    e_a[U] = odvisneU
    e_a[V] = z[0::2]

    neodvisneV = z[0::2]
    neodvisneVdot = z[1::2]

    # assemble velocity vector - into original index order
    qdot = np.zeros(e_a.shape)  # predefine vector
    qdot[U] = odvisneUdot
    qdot[V] = neodvisneVdot

    # Update dependent coordinates
    odvisneU = odvi_U(odvisneU, hh, ne_c, neodvisneV, V, U, tipverige, e_zp)

    # e_a[U.flatten()] = odvisneU
    e_a[U] = odvisneU

    # some notes
    # % e_zp are all coordinates, ordered as cable
    # e_a are the current values, so they are temporarily transformed/rotated within the program

    # % since the constraints are time-invariant, I read the pre-computed 'Cdi','Ced' from the static analysis

    # the derivatives of Ct and Ctt are both zero and I do not consider them further
    # Ct=0, Ctt=0


    # Update dependent velocities
    odvisneUdot = Cdi @ neodvisneVdot

    # Calculate forces
    Qeg_0 = Qe_g.copy() # just to define
    #Qeg_0 = load_over_time(t, timeload, Qe, Qe_g, Qe_breme, Qe_gg)

    #note-
    qkqt = Ke_vect(hh, ne_c, e_a, cab, le)
    Qd = Qd_vect(hh, ne_c, cab, le, e_a, qdot)
    Qem=i_Qem(t,e_a, cab,le_full,Ip)
    Qcon=j_Qcon(t, e_a, qdot, cab, le_full)

    qkqt=qkqt + Qd - Qem -Qcon


    #E[V] = neodvisneV
    #E[U] = odvisneU
    # print(qkqt)
    # print(t)
    # breakpoint()

    qkqt = qkqt[UV_indices]

    # Calculate accelerations
    pdotdot, lambda_vec = pdotdot_lagrange(M, Qeg_0, qkqt, qv, qc, Ce)

    # Create derivative vector
    dz = np.zeros(2 * ni)
    #dz[2 * i] = z[2 * i + 1]
    #dz[2 * i + 1] = pdotdot[nd:nd + ni]
    dz[0::2] = z[1::2]  # position derivatives = velocities
    dz[1::2] = pdotdot[nd:nd + ni]  # velocity derivatives = accelerations
    #print(t)
    return dz