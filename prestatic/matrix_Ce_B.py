import numpy as np
from scipy.linalg import block_diag
#from typing import Tuple, List, Any
#from scipy.sparse import csr_matrix
from scipy.linalg import lstsq, solve, LinAlgError

from prestatic.constraintMatrices import constraint_C,constraint_Cgl,constraint_Cgr,constraint_CCL

def matrix_Ce(hh,ne_c,tipverige):
    """
    Matrix Ce_B function for constraint matrix assembly and Gaussian elimination
    """
    #print('tipVerige',tipverige)

    # Initialize variables
    Cq_elementa = np.empty((0, 0))

    konstante = []
    neod_fix = []

    # assembly constraints derivative matrix
    for k in range(0, hh):  #
        for i in range(0, ne_c[k]):  #

            if k == 0 and i == 0:  # A. first one, left hand side is given
                Cq_elementa = constraint_Cgl()
                #print(Cq_elementa)
                # [r,c]=size(Cq_elementa)
                konstante = [0, 1, 2]  #
                neod_fix = [3,4,5]  # list for fixed independent coordinates


            elif i >= 1 and i != ne_c[k]-1:  # B. cable-cable joint
                r, c = Cq_elementa.shape

                Cq_elementa = np.block([
                    [Cq_elementa, np.zeros((r, 12))],
                    [np.zeros((6, c - 12)), constraint_C()]
                ])
                # display('b')

            elif i >= 1 and i == ne_c[k]-1 and k != hh - 2:  # C. joint cable-cable-cleave
                # tipberige= SP or EP - no cleave-
                # tipverige=CC- - cable cleave
                # tipverige=IR- intermediate restraint

                if tipverige[k + 1] == 'CC':  # cable cleave

                    r, c = Cq_elementa.shape
                    konstante.extend([c - 6, c - 5, c - 4, c+0, c+1, c+2 ])  # end of element a
                    #neod_fix.extend([c + 9, c + 10, c + 11])  # inclinations
                    Cq_elementa = np.block([
                        [Cq_elementa, np.zeros((r, 12))],
                        [np.zeros((9, c - 12)), constraint_CCL()]
                    ])
                    # constrain from cable right-hand . connection to first proceeding cable part
                    r, c = Cq_elementa.shape
                    Cq_elementa = np.block([
                        [Cq_elementa, np.zeros((r, 12))],
                        [np.zeros((6, c - 12)), constraint_C()]
                    ])
                    # display('e2')

                elif tipverige[k + 1] == 'IR':  # intermediate restraint
                    r, c = Cq_elementa.shape
                    konstante.extend([c - 6, c - 5, c - 4, c+0, c+1, c+2])  # end of element a
                    # neod_fix.extend([c - 6, c - 5, c - 4, c - 3, c - 2, c - 1, c + 15, c + 16,
                    #                c + 17])  # subtract 1 from all indices
                    Cq_elementa = np.block([
                        [Cq_elementa, np.zeros((r, 12))],
                        [np.zeros((9, c - 12)), constraint_CCL()]
                    ])
                    # constrain from cable right-hand . connection to first proceeding cable part
                    r, c = Cq_elementa.shape
                    Cq_elementa = np.block([
                        [Cq_elementa, np.zeros((r, 12))],
                        [np.zeros((6, c - 12)), constraint_C()]
                    ])



            elif i >= 1 and i == ne_c[k]-1 and k == hh - 2:  # C. conductor fixed right-hand end
                r, c = Cq_elementa.shape
                # Append new rows instead of adding zeros matrix
                Cq_elementa = np.vstack([
                    Cq_elementa,
                    np.hstack([np.zeros((3, c - 12)), constraint_Cgr()])
                ])
                konstante.extend([c - 6, c - 5, c - 4])  #
                # display('z');


    Ce_active = block_diag(Cq_elementa,Cq_elementa,Cq_elementa)  # multiplied matrix because there are 3 cables
    i,j=Cq_elementa.shape
    del Cq_elementa # release memory

    # also for (konstante,neod_fix) extended list because there are 3 cables
    konstante=konstante + [x + j for x in konstante] + [x + 2*j for x in konstante]
    neod_fix=neod_fix + [x + j for x in neod_fix] + [x + 2*j for x in neod_fix]

    return Ce_active, konstante, neod_fix



def B_matrix(Ce_active,konstante,neod_fix):
    """
    Function that determines dependent U and independent V variables
    Python version converted from MATLAB
    """

    # Call UV_clenitev function (implement separately)
    U, V, Ce = UV_clenitev(Ce_active, konstante, neod_fix)

    # Get matrix dimensions
    nc, n = Ce.shape

    # Split Jacobian matrix into dependent and independent parts
    Ced = Ce[:, :nc]
    Cei = Ce[:, nc:]

    # Calculate Cdi = -inv(Ced) * Cei
    # Calculate Cdi = -inv(Ced) * Cei using least squares for robustness
    try:
        # Try direct solve first (faster if matrix is well-conditioned)
        Cdi = -solve(Ced, Cei)
        #print('izstop')
    except (LinAlgError, np.linalg.LinAlgError):
        print("Direct solve failed, using least squares")
        # Use least squares if direct solve fails
        Cdi = -lstsq(Ced, Cei, cond=None)[0]

    # Create identity matrix and construct Bdi
    I = np.eye(n - nc)
    Bdi = np.vstack([I, Cdi])

    # Find pairs (this part is complex )
    c = np.array([])

    # For first dependent variable after identity matrix
    row_idx = n - nc
    if row_idx < Bdi.shape[0]:
        # Find non-zero elements in this row
        non_zero_indices = np.flatnonzero(Bdi[row_idx, :])
        if non_zero_indices.size == 0:
            c = np.array([-1])
        else:
            c = non_zero_indices.copy()
    else:
        c = np.array([0])

    # For remaining dependent variables
    for i in range(n - nc + 1, n):
        if i < Bdi.shape[0]:
            non_zero_indices = np.flatnonzero(Bdi[i, :])
            if non_zero_indices.size == 0:
                c = np.append(c, -1)
            else:
                c = np.append(c, non_zero_indices)
        else:
            c = np.append(c, 0)

    # Create pairs array
    pari = np.column_stack([c, U])
    #prvi stolpec index koordinat odvisne in drugi stolpec index pripadajoče neodvisne koordinate

    return Bdi,U,V,pari,Cdi,Ced,Ce






def UV_clenitev(Ce_active, konstante, neod_fix):
    """
    Function that determines dependent U and independent V variables
    """

    # Get matrix dimensions
    nc, n = Ce_active.shape

    # Create index array
    index = np.arange(0, n )  # 0-based indexing

    # Remove fixed independent coordinates
    index_1 = np.setdiff1d(index, neod_fix)

    # Combine indices: first the non-fixed, then the fixed ones
    index_1 = np.concatenate([index_1, neod_fix])

    # Reorder Ce_active columns according to index_1
    # Note: MATLAB uses 1-based indexing, Python uses 0-based
    Ce_active = Ce_active[:, index_1]

    # Apply Gaussian elimination
    active_Ce, index_result = gausovaeliminacija_vectorized(Ce_active, nc, index_1)

    # Extract dependent and independent variables
    U_active = index_result[0:nc]
    V_active = index_result[nc:n]

    return U_active, V_active, active_Ce


#def gausovaeliminacija(Ce, nc, index): to je moja prvotna koda zamenjana spodaj s predlogom AI
  #  """
  #  Gaussian elimination algorithm with partial pivoting
   # """
    # n = Ce.shape[1]  # number of columns
    # Ce = Ce.copy()  # work with a copy
    # index = index.copy()  # work with a copy
    #
    # for k in range(nc):
    #     # Exchange step - find pivot in current row
    #     # Look only in columns from k to end for pivoting
    #     row_to_search = Ce[k, k:]
    #     if row_to_search.size > 0:
    #         max_idx = np.argmax(np.abs(row_to_search))
    #         i = max_idx + k  # adjust index to full column range
    #     else:
    #         i = k
    #
    #     # Swap columns if necessary
    #     if i > k:
    #         # Swap columns in Ce
    #         pom = Ce[:, k].copy()
    #         Ce[:, k] = Ce[:, i]
    #         Ce[:, i] = pom
    #
    #         # Swap indices
    #         pom_i = index[k]
    #         index[k] = index[i]
    #         index[i] = pom_i
    #
    #     # Gaussian elimination
    #     if k < nc - 1:  # elimination up to second last row
    #         for v in range(k + 1, nc):  # rows
    #             if Ce[k, k] == 0:
    #                 print('singularno')
    #                 m = 0
    #             else:
    #                 m = Ce[v, k] / Ce[k, k]
    #
    #             # Eliminate elements below pivot
    #             for s in range(k, n):  # columns
    #                 Ce[v, s] = Ce[v, s] - m * Ce[k, s]
    #
    # return Ce, index


# Alternative implementation using numpy operations for better performance
def gausovaeliminacija_vectorized(Ce, nc, index):
    """
    Vectorized version of Gaussian elimination for better performance
    """
    n = Ce.shape[1]
    Ce = Ce.copy().astype(float)  # ensure float type for division
    index = index.copy()

    for k in range(nc):
        # Partial pivoting
        if k < n:  # only if we have columns left
            max_idx = np.argmax(np.abs(Ce[k, k:])) + k
            if max_idx > k:
                # Swap columns
                Ce[:, [k, max_idx]] = Ce[:, [max_idx, k]]
                index[k], index[max_idx] = index[max_idx], index[k]

        # Gaussian elimination
        if k < nc - 1 and k < n and Ce[k, k] != 0:
            pivot = Ce[k, k]
            for v in range(k + 1, nc):
                factor = Ce[v, k] / pivot
                Ce[v, k:] = Ce[v, k:] - factor * Ce[k, k:]

    return Ce, index
