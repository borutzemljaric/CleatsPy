import numpy as np
from scipy.linalg import lstsq

from rutines.resfor_Re import Ke_vect

# ============================================================================
# Main function: reakcijeVerige
# ============================================================================
def reakctionNodes(zacpog, e_zp,Qe_g, hh, ne_c, le, tipverige,cab):
    """Main function to calculate reactions."""

  # Call for calculation of external forces (Ke-Qg)
    K_elementa = Ke_vect(hh, ne_c, zacpog, cab, le)
    Re = Qe_g - K_elementa
    #Qc = Re  # just rename so we don't modify the program further
    Flevo, Fdesno = sileMomenti(Re)

    return Flevo, Fdesno


# ============================================================================
# Helper function: sileMomenti
# ============================================================================
def sileMomenti(Qc):
    """Calculate internal forces in nodes in 3D."""
    nc = Qc.size  # number of generalized forces
    I = np.eye(3)
    O = np.zeros((3, 3))

    # For element A
    Sk1 = np.hstack([O, O, I, O])  # valid for r in point b of body A - ksi equals one
    Sk2 = np.hstack([I, O, O, O])  # valid for r in point a of body A - ksi equals zero

    # First three outside loop
    ii = 1
    # Extract coordinates for current element
    Qc_element = Qc[((ii - 1) * 12):((ii - 1) * 12 + 12)].reshape(-1, 1)

    # Solve linear systems: Sk2.T * x = Qc_element
    # Using least squares solution for minimum norm
    Flevo = lstsq(Sk2.T, Qc_element, cond=None)[0]  # dimensions: n × number of time records
    Fdesno = lstsq(Sk1.T, Qc_element, cond=None)[0]  # dimensions: n × number of time records

    for ii in range(2, nc // 12 + 1):
        Qc_element = Qc[((ii - 1) * 12):((ii - 1) * 12 + 12)].reshape(-1, 1)

        Flevo_n = lstsq(Sk2.T, Qc_element, cond=None)[0]
        Flevo = np.vstack([Flevo, Flevo_n])

        Fdesno_n = lstsq(Sk1.T, Qc_element, cond=None)[0]
        Fdesno = np.vstack([Fdesno, Fdesno_n])

    return Flevo, Fdesno


