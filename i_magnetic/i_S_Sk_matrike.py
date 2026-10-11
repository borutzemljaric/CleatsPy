import numpy as np

import numpy as np


def i_S_Sk_matrike(k, L):
    """
    Shape function matrix S and its derivatives Sk, Skk.
    Returns three (3, 12, n) arrays.
    """
    k = np.asarray(k, dtype=float)
    n = k.size

    k2 = k * k
    k3 = k2 * k

    S   = np.zeros((3, 12, n))
    Sk  = np.zeros((3, 12, n))
    Skk = np.zeros((3, 12, n))

    # ---- S ----
    s0 =  2.0*k3 - 3.0*k2 + 1.0
    s1 = L*(k3 - 2.0*k2 + k)
    s2 = -2.0*k3 + 3.0*k2
    s3 = L*(k3 - k2)

    # ---- Sk ----
    d0 =  6.0*k2 - 6.0*k
    d1 = L*(3.0*k2 - 4.0*k + 1.0)
    d2 =  6.0*k - 6.0*k2
    d3 = L*(3.0*k2 - 2.0*k)

    # ---- Skk ----
    dd0 =  12.0*k - 6.0
    dd1 = L*(6.0*k - 4.0)
    dd2 =  6.0 - 12.0*k
    dd3 = L*(6.0*k - 2.0)

    # Fill each row with the appropriate 4 slots
    # Row r has slots at columns: r, r+3, r+6, r+9
    for r in range(3):
        S  [r, r,      :] = s0
        S  [r, r + 3,  :] = s1
        S  [r, r + 6,  :] = s2
        S  [r, r + 9,  :] = s3

        Sk [r, r,      :] = d0
        Sk [r, r + 3,  :] = d1
        Sk [r, r + 6,  :] = d2
        Sk [r, r + 9,  :] = d3

        Skk[r, r,      :] = dd0
        Skk[r, r + 3,  :] = dd1
        Skk[r, r + 6,  :] = dd2
        Skk[r, r + 9,  :] = dd3

    return S, Sk, Skk

# def _orig_i_S_Sk(k, L):
#     """
#     Calculate shape function matrix S and its derivatives Sk, Skk.
#     Sub-function for electromagnetic force and contact force calculations.
#     Based on code from supplementary material to manuscript JSV,
#     author Borut Zemljaric, 2021.
#     Originally in Slovenian, translated to English.
#
#     Parameters:
#         k : 1D array of parameter values (integration points)
#         L : element length
#
#     Returns:
#         S  : shape function matrix, size (3, 12, n)
#         Sk : first derivative of S, size (3, 12, n)
#         Skk: second derivative of S, size (3, 12, n)
#     """
#     n = len(k)
#     # Pre-allocate 3D arrays
#     S = np.zeros((3, 12, n))
#     Sk = np.zeros((3, 12, n))
#     Skk = np.zeros((3, 12, n))
#
#     for i in range(n):
#         ki = k[i]
#         k2 = ki ** 2
#         k3 = ki ** 3
#
#         # Shape function matrix S (3x12) – exactly as in MATLAB
#         S[:, :, i] = np.array([
#             [2 * k3 - 3 * k2 + 1, 0, 0, L * (k3 - 2 * k2 + ki), 0, 0, -2 * k3 + 3 * k2, 0, 0, -L * (-k3 + k2), 0, 0],
#             [0, 2 * k3 - 3 * k2 + 1, 0, 0, L * (k3 - 2 * k2 + ki), 0, 0, -2 * k3 + 3 * k2, 0, 0, -L * (-k3 + k2), 0],
#             [0, 0, 2 * k3 - 3 * k2 + 1, 0, 0, L * (k3 - 2 * k2 + ki), 0, 0, -2 * k3 + 3 * k2, 0, 0, -L * (-k3 + k2)]
#         ])
#
#         # Derivative of shape function matrix Sk (3x12)
#         Sk[:, :, i] = np.array([
#             [6 * k2 - 6 * ki, 0, 0, L * (3 * k2 - 4 * ki + 1), 0, 0, 6 * ki - 6 * k2, 0, 0, -L * (2 * ki - 3 * k2), 0,
#              0],
#             [0, 6 * k2 - 6 * ki, 0, 0, L * (3 * k2 - 4 * ki + 1), 0, 0, 6 * ki - 6 * k2, 0, 0, -L * (2 * ki - 3 * k2),
#              0],
#             [0, 0, 6 * k2 - 6 * ki, 0, 0, L * (3 * k2 - 4 * ki + 1), 0, 0, 6 * ki - 6 * k2, 0, 0,
#              -L * (2 * ki - 3 * k2)]
#         ])
#
#         # Second derivative of shape function matrix Skk (3x12)
#         Skk[:, :, i] = np.array([
#             [12 * ki - 6, 0, 0, L * (6 * ki - 4), 0, 0, 6 - 12 * ki, 0, 0, L * (6 * ki - 2), 0, 0],
#             [0, 12 * ki - 6, 0, 0, L * (6 * ki - 4), 0, 0, 6 - 12 * ki, 0, 0, L * (6 * ki - 2), 0],
#             [0, 0, 12 * ki - 6, 0, 0, L * (6 * ki - 4), 0, 0, 6 - 12 * ki, 0, 0, L * (6 * ki - 2)]
#         ])
#
#     return S, Sk, Skk
#
# rng = np.random.default_rng(0)
# for trial in range(5):
#     L = rng.uniform(0.5, 5.0)
#     k = rng.uniform(0.05, 0.95, size=4)
#     S_old, Sk_old, Skk_old = _orig_i_S_Sk(k, L)
#     S_new, Sk_new, Skk_new = i_S_Sk_matrike(k, L)
#     assert np.allclose(S_old,   S_new,   rtol=1e-12, atol=1e-14), "S mismatch"
#     assert np.allclose(Sk_old,  Sk_new,  rtol=1e-12, atol=1e-14), "Sk mismatch"
#     assert np.allclose(Skk_old, Skk_new, rtol=1e-12, atol=1e-14), "Skk mismatch"
# print("i_S_Sk verified ✔fgggg")