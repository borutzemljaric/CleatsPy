import numpy as np

def h_Sk_matrike(k, L):
    """
    Skk matrices for bending. Returns array of shape (3, 12, n).
    """
    k = np.asarray(k, dtype=float)
    n = k.size

    # vectorized scalars, shape (n,)
    p = 12.0*k - 6.0          # = -(6 - 12k)
    q = L * (6.0*k - 4.0)
    s = L * (6.0*k - 2.0)

    Skk = np.zeros((3, 12, n))

    # row 0 block
    Skk[0, 0,  :] =  p
    Skk[0, 3,  :] =  q
    Skk[0, 6,  :] = -p
    Skk[0, 9,  :] =  s

    # row 1 is row 0 shifted right by 3 columns
    Skk[1, 1,  :] =  p
    Skk[1, 4,  :] =  q
    Skk[1, 7,  :] = -p
    Skk[1, 10, :] =  s

    # row 2 is row 0 shifted right by 6 columns
    Skk[2, 2,  :] =  p
    Skk[2, 5,  :] =  q
    Skk[2, 8,  :] = -p
    Skk[2, 11, :] =  s

    return Skk


# def _orig_Sk(k, L):
#     """
#     Calculate Skk matrices for bending calculation
#     Returns 3D array where each slice Skk[:, :, i] is a 3x12 matrix for the i-th integration point
#
#
#     Returns:
#     Skk: 3D array of shape (3, 12, n) where n is number of integration points
#     izracun matrik(i) B2-jev za izračun krivljenja, b-ji so tridimenzionalne matrike s po numel(k(i)(i)) globino
#     k je stolpni vektor
#     """
#     n = len(k)  #
#     Skk = np.zeros((3, 12, n))
#
#     for i in range(n):
#         ki = k[i]  # Get scalar value from array
#
#         # Build the 3x12 matrix for this integration point
#         Skk[:, :, i] = np.array([
#             [12 * ki - 6, 0, 0, L * (6 * ki - 4), 0, 0, 6 - 12 * ki, 0, 0, L * (6 * ki - 2), 0, 0],
#             [0, 12 * ki - 6, 0, 0, L * (6 * ki - 4), 0, 0, 6 - 12 * ki, 0, 0, L * (6 * ki - 2), 0],
#             [0, 0, 12 * ki - 6, 0, 0, L * (6 * ki - 4), 0, 0, 6 - 12 * ki, 0, 0, L * (6 * ki - 2)]
#         ])
#
#     return Skk
#
#
#
#
#
# rng = np.random.default_rng(0)
# for trial in range(5):
#     L = rng.uniform(0.5, 5.0)
#     k = rng.uniform(0.05, 0.95, size=4)
#     assert np.allclose(_orig_Sk(k, L), h_Sk_matrike(k, L), rtol=1e-12, atol=1e-14), \
#         f"Mismatch at L={L}, k={k}"
# print("Sk template verified ✔")


