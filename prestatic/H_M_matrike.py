import numpy as np
#import scipy.sparse as sp
from scipy.linalg import block_diag
from scipy.linalg import inv



def M_element(m_c, L):
    """Calculate the 12x12 mass matrix for a beam element."""
    M_ij = m_c * np.array([
        [(13 * L) / 35, 0, 0, (11 * L ** 2) / 210, 0, 0, (9 * L) / 70, 0, 0, -(13 * L ** 2) / 420, 0, 0],
        [0, (13 * L) / 35, 0, 0, (11 * L ** 2) / 210, 0, 0, (9 * L) / 70, 0, 0, -(13 * L ** 2) / 420, 0],
        [0, 0, (13 * L) / 35, 0, 0, (11 * L ** 2) / 210, 0, 0, (9 * L) / 70, 0, 0, -(13 * L ** 2) / 420],
        [(11 * L ** 2) / 210, 0, 0, L ** 3 / 105, 0, 0, (13 * L ** 2) / 420, 0, 0, -L ** 3 / 140, 0, 0],
        [0, (11 * L ** 2) / 210, 0, 0, L ** 3 / 105, 0, 0, (13 * L ** 2) / 420, 0, 0, -L ** 3 / 140, 0],
        [0, 0, (11 * L ** 2) / 210, 0, 0, L ** 3 / 105, 0, 0, (13 * L ** 2) / 420, 0, 0, -L ** 3 / 140],
        [(9 * L) / 70, 0, 0, (13 * L ** 2) / 420, 0, 0, (13 * L) / 35, 0, 0, -(11 * L ** 2) / 210, 0, 0],
        [0, (9 * L) / 70, 0, 0, (13 * L ** 2) / 420, 0, 0, (13 * L) / 35, 0, 0, -(11 * L ** 2) / 210, 0],
        [0, 0, (9 * L) / 70, 0, 0, (13 * L ** 2) / 420, 0, 0, (13 * L) / 35, 0, 0, -(11 * L ** 2) / 210],
        [-(13 * L ** 2) / 420, 0, 0, -L ** 3 / 140, 0, 0, -(11 * L ** 2) / 210, 0, 0, L ** 3 / 105, 0, 0],
        [0, -(13 * L ** 2) / 420, 0, 0, -L ** 3 / 140, 0, 0, -(11 * L ** 2) / 210, 0, 0, L ** 3 / 105, 0],
        [0, 0, -(13 * L ** 2) / 420, 0, 0, -L ** 3 / 140, 0, 0, -(11 * L ** 2) / 210, 0, 0, L ** 3 / 105]
    ])
    return M_ij


def H_M_matrike(hh, ne_c, cab, le ,  Ce_active):
    """Main function to calculate mass and H matrices."""

    # Calculate element mass matrices, for first conductor
    M_elementa_list = []

    for k in range(hh - 1):
        for i in range(ne_c[k]):
            M_elementa_list.append(M_element(cab.m, le[k, i]))

    # Create block diagonal matrix for first conductor
    M = M_elementa_list[0]
    for i in range(1, len(M_elementa_list)):
        M = block_diag(M, M_elementa_list[i])

    #print(type(M), M.dtype)
    # Convert to sparse and compute inverse
    #M_sparse = sp.csr_matrix(M)
    #if M_sparse.shape[0] > 0:

    # there are three (same) cables , so full block_diagonal is created
    M = block_diag(M, M, M)

    # inverse of mass matrix
    Minv = np.linalg.inv(M)

    Ce_1_active = Ce_active.T  # rezultat transponirana matrika Cq
    # Compute H matrices

    H_ll = inv(Ce_active @ Minv @ Ce_1_active)
    H_ee = Minv - (Minv @ Ce_1_active @ H_ll @ Ce_active @ Minv)
    H_el = Minv @ Ce_1_active @ H_ll
    H_le = H_el.T

    return H_ll, H_ee, H_le, H_el, M

