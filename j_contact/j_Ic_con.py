import numpy as np
from scipy.integrate import quad_vec

from j_izracun_Qcon import S_matrix

def j_Ic_con(L, N):
    """
    Compute the integral of the shape function matrix over sub‑divisions.

    Ip : ndarray, shape (3, 12, N)
        Integral of the shape function matrix over each sub‑interval.
        Ip[:, :, k] = ∫_{a_k}^{b_k} L * S(k) dk, where a_k and b_k are the
        boundaries of the k‑th sub‑division.
    """

    # Integrand: L * S(k)
    def fun_flat(k):
        return (L * S_matrix(k,L)).flatten()   # flatten to 1D array of length 36

    # Sub‑division points (including boundaries)
    ss = np.linspace(0, 1, N + 1)

    # Preallocate output array (3x12xN)
    Ip = np.zeros((3, 12, N))

    for idx in range(N):
        a = ss[idx]
        b = ss[idx + 1]
        # Numerically integrate the flattened integrand over [a, b]
        res= quad_vec(fun_flat, a, b)[0]
        # Reshape back to (3, 12) and store
        Ip[:, :, idx] = res.reshape((3, 12))

    return Ip