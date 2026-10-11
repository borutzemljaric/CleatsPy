import numpy as np
from scipy.optimize import brentq

from i_magnetic.i_S_Sk_matrike import i_S_Sk_matrike
from j_contact.j_overlap import j_overlap


def S_matrix(k, L):
    """
    Returns the 3x12 Hermite cubic shape function matrix evaluated at parameter k
    for an element of length L.  Used in sub-function called 'fun'
    """
    k2 = k * k
    k3 = k2 * k
    # Hermite basis functions
    N1 = 2 * k3 - 3 * k2 + 1
    N2 = L * (k3 - 2 * k2 + k)
    N3 = 3 * k2 - 2 * k3
    N4 = -L * (k2 - k3)  # = L*(k3 - k2)

    # Build 3x12 matrix (each row corresponds to x, y, z)
    # S = np.array([
    #     [N1, 0, 0, N2, 0, 0, N3, 0, 0, N4, 0, 0],
    #     [0, N1, 0, 0, N2, 0, 0, N3, 0, 0, N4, 0],
    #     [0, 0, N1, 0, 0, N2, 0, 0, N3, 0, 0, N4]
    # ])

    S = np.zeros((3, 12))
    S[0, 0::3] = N1, N2, N3, N4
    S[1, 1::3] = N1, N2, N3, N4
    S[2, 2::3] = N1, N2, N3, N4

    return S

def funPjPi(k, L, ei, pt, ei_dot, pt_dot):
    """
    Returns (PjPi, PjPi_dot) where:
        PjPi = S(k)*ei - pt
        PjPi_dot = S(k)*ei_dot - pt_dot
    """
    S = S_matrix(k, L)
    PjPi = S @ ei - pt
    PjPi_dot = S @ ei_dot - pt_dot
    return PjPi, PjPi_dot


def j_izracun_Qcon(ei, ej, Li, Lj, N_con, pK_con, cab, ei_dot, ej_dot):
    """
    Returns generalized contact forces (12x1 vector) for element i.
    """
    # convert d_c from mm to m
    d_c = cab.d * 1e-3

    # Obtain shape matrices at integration points
    Si, Ski, _ = i_S_Sk_matrike(pK_con, Li)  # Ski not used
    Sj, Skj, _ = i_S_Sk_matrike(pK_con, Lj)  # Skj needed for tangent

    # Preallocate arrays
    PTi = np.zeros((3, N_con))
    PTj = np.zeros((3, N_con))
    t_vec = np.zeros((3, N_con))
    PTj_dot = np.zeros((3, N_con))

    for k in range(N_con):
        PTj[:, k] = Sj[:, :, k] @ ej  # global vector of contact point on j element
        PTi[:, k] = Si[:, :, k] @ ei
        rdot = Skj[:, :, k] @ ej  # first derivation
        norm_rdot = np.linalg.norm(rdot)  # absolute value of rdot
        if norm_rdot > 1e-15:
            t_vec[:, k] = rdot / norm_rdot
        else:
            t_vec[:, k] = np.zeros(3)
        PTj_dot[:, k] = Sj[:, :, k] @ ej_dot

    # Check if pair of the elements are in the 'zone' of possible contact
    # Check overlap in y-direction (as in original code)
    c_y = j_overlap(PTj[1, 0], PTj[1, -1], PTi[1, 0], PTi[1, -1])


    # Inner function to compute contact force
    def silaKontakta():
        Qcontact = np.zeros(12)
        for kj in range(N_con):
            # Define the function whose root we seek (projection condition)
            def fun(k):
                S = S_matrix(k, Li)
                return np.dot((S @ ei - PTj[:, kj]), t_vec[:, kj])

            test1 = fun(0.0)  #left hand side point
            test2 = fun(1.0)   #right hand side point

            # If sign change, there is a root
            if np.sign(test1) != np.sign(test2):
                # Find root in [0, 1]
                k1 = brentq(fun, 0.0, 1.0)  #if the sign change the fzero has a solution

                # Compute vector from point on j to point on i and its derivative
                PjPi, delta_dot = funPjPi(k1, Li, ei, PTj[:, kj], ei_dot, PTj_dot[:, kj]) # vector of distance between points
                PjPi_abs = np.linalg.norm(PjPi)
                delta_dot_abs = np.linalg.norm(delta_dot)

                if PjPi_abs > 1e-15:  # normalisation of the vector
                    n_vec = PjPi / PjPi_abs
                else:
                    n_vec = np.zeros(3)

                if PjPi_abs < d_c:
                    delta = d_c - PjPi_abs
                    force_mag = cab.kn * (delta ** 1.1) + cab.mi * delta_dot_abs
                    force_vec = force_mag * n_vec
                    # Add contribution to generalized forces (using Sj for element j)
                    Qcontact += Sj[:, :, kj].T @ force_vec

        return Qcontact

    if c_y == 0:
        # no contact
        Qcon_i = np.zeros(12)
    else:
        # contact detected
        Qcon_i = silaKontakta()

    return Qcon_i.reshape(12, 1)