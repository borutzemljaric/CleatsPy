import numpy as np
from scipy.integrate import quad_vec

def i_Ip_all(le,cab):

    Nsub=cab.Nseg # number of subdivision of each element
    Nel=len(le)  # nuber of elements, whole system


    #For first element in span the Ip is calculated, the return is (12x36xNsum) array
    # the loop goes through all spans, and at final assemble final Ip_all is ((12*hh)x36xxNsub) array


    # List to collect all vertically stacked parts
    all_parts = []

    for i in range(Nel):
        # Call the function once for this row (passing only the first column of le_scaled)
        arr_3d = i_Ip(le[i], Nsub)  # returns a 3D array, shape (12, 36, Nsub)

        all_parts.append(arr_3d) #

    # Vertically concatenate all blocks along axis 0
    Ip_all = np.concatenate(all_parts, axis=0)  # final 3D array/ for first conductor

    #result is array in 3D, vertically stacked where each 12 rows goes to one element

    return Ip_all



# External function that calculates Ip integrals for individual segments
def i_Ip(L, N):
    """
    Calculate the integral of the shape function matrix over sub‑divisions.
    """

    # Integration limits over sub‑segments (0 to 1 divided into N parts)
    ss = np.linspace(0, 1, N + 1)

    # Pre‑allocate result: 12 rows, 36 columns, N layers
    Ip = np.zeros((12, 36, N))

    # Flatten the 12x36 output for quad_vec
    def fun_flat(k,L):
        return fun_Sb3(k,L).flatten()

    # Integrate over each sub‑segment
    for idx in range(N):
        a = ss[idx]
        b = ss[idx + 1]
        result, _ = quad_vec(fun_flat, a, b,args=(L,))
        Ip[:, :, idx] = result.reshape(12, 36)

    return Ip





def fun_Sb3(k, L):
    """
    #S'*b3
    Returns a 12x36 matrix with entries computed from k and L.
    """
    k2 = k * k
    k3 = k2 * k

    # ---- Block 1 (rows 0‑2) ----
    term1 = (-6 * k2 + 6 * k) * (2 * k3 - 3 * k2 + 1)
    term2 = (3 * k2 - 4 * k + 1) * (2 * k3 - 3 * k2 + 1)
    term3 = (-3 * k2 + 2 * k) * (2 * k3 - 3 * k2 + 1)

    # ---- Block 2 (rows 3‑5) ----
    term4 = L * (-6 * k2 + 6 * k) * (k3 - 2 * k2 + k)
    term5 = L**2 * (3 * k2 - 4 * k + 1) * (k3 - 2 * k2 + k)
    term6 = L**2 * (-3 * k2 + 2 * k) * (k3 - 2 * k2 + k)

    # ---- Block 3 (rows 6‑8) ----
    term8 = (-6 * k2 + 6 * k) * (-2 * k3 + 3 * k2)
    term9 = (-2 * k3 + 3 * k2) * (3 * k2 - 4 * k + 1)
    term10 = (-3 * k2 + 2 * k) * (-2 * k3 + 3 * k2)

    # ---- Block 4 (rows 9‑11) ----
    termA = (-6 * k2 + 6 * k) * (-k3 + k2)
    termB = (-k3 + k2) * (3 * k2 - 4 * k + 1)
    termC = (-3 * k2 + 2 * k) * (-k3 + k2)

    # Initialize matrix (12 rows, 36 columns)
    M = np.zeros((12, 36))

    # ---------- Row 0 ----------
    M[0, 14] = term1
    M[0, 17] = -L * term2
    M[0, 20] = -term1
    M[0, 23] = L * term3
    M[0, 25] = -term1
    M[0, 28] = L * term2
    M[0, 31] = term1
    M[0, 34] = -L * term3

    # ---------- Row 1 ----------
    M[1, 2] = -term1
    M[1, 5] = L * term2
    M[1, 8] = term1
    M[1, 11] = -L * term3
    M[1, 24] = term1
    M[1, 27] = -L * term2
    M[1, 30] = -term1
    M[1, 33] = L * term3

    # ---------- Row 2 ----------
    M[2, 1] = term1
    M[2, 4] = -L * term2
    M[2, 7] = -term1
    M[2, 10] = L * term3
    M[2, 12] = -term1
    M[2, 15] = L * term2
    M[2, 18] = term1
    M[2, 21] = -L * term3

    # ---------- Row 3 ----------
    M[3, 14] = term4
    M[3, 17] = -term5
    M[3, 20] = -term4
    M[3, 23] = term6
    M[3, 25] = -term4
    M[3, 28] = term5
    M[3, 31] = term4
    M[3, 34] = -term6

    # ---------- Row 4 ----------
    M[4, 2] = -term4
    M[4, 5] = term5
    M[4, 8] = term4
    M[4, 11] = -term6
    M[4, 24] = term4
    M[4, 27] = -term5
    M[4, 30] = -term4
    M[4, 33] = term6

    # ---------- Row 5 ----------
    M[5, 1] = term4
    M[5, 4] = -term5
    M[5, 7] = -term4
    M[5, 10] = term6
    M[5, 12] = -term4
    M[5, 15] = term5
    M[5, 18] = term4
    M[5, 21] = -term6

    # ---------- Row 6 ----------
    M[6, 14] = term8
    M[6, 17] = -L * term9
    M[6, 20] = -term8
    M[6, 23] = L * term10
    M[6, 25] = -term8
    M[6, 28] = L * term9
    M[6, 31] = term8
    M[6, 34] = -L * term10

    # ---------- Row 7 ----------
    M[7, 2] = -term8
    M[7, 5] = L * term9
    M[7, 8] = term8
    M[7, 11] = -L * term10
    M[7, 24] = term8
    M[7, 27] = -L * term9
    M[7, 30] = -term8
    M[7, 33] = L * term10

    # ---------- Row 8 ----------
    M[8, 1] = term8
    M[8, 4] = -L * term9
    M[8, 7] = -term8
    M[8, 10] = L * term10
    M[8, 12] = -term8
    M[8, 15] = L * term9
    M[8, 18] = term8
    M[8, 21] = -L * term10

    # ---------- Row 9 ----------
    M[9, 14] = -L * termA
    M[9, 17] = L**2 * termB
    M[9, 20] = L * termA
    M[9, 23] = -L**2 * termC
    M[9, 25] = L * termA
    M[9, 28] = -L**2 * termB
    M[9, 31] = -L * termA
    M[9, 34] = L**2 * termC

    # ---------- Row 10 ----------
    M[10, 2] = L * termA
    M[10, 5] = -L**2 * termB
    M[10, 8] = -L * termA
    M[10, 11] = L**2 * termC
    M[10, 24] = -L * termA
    M[10, 27] = L**2 * termB
    M[10, 30] = L * termA
    M[10, 33] = -L**2 * termC

    # ---------- Row 11 ----------
    M[11, 1] = -L * termA
    M[11, 4] = L**2 * termB
    M[11, 7] = L * termA
    M[11, 10] = -L**2 * termC
    M[11, 12] = L * termA
    M[11, 15] = -L**2 * termB
    M[11, 18] = -L * termA
    M[11, 21] = L**2 * termC

    return M

