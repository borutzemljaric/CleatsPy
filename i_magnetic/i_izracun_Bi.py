import numpy as np
from i_magnetic.i_S_Sk_matrike import i_S_Sk_matrike

def i_izracun_Bi(ei, ej, Li, Lj, Ij, N, pK):
    """
    Function to calculate the shape function matrix and derivative of shape function matrix.
    """


    # Get shape function matrices at integration points (pK)
    # i_S_Sk_matrike returns three arrays: S, Sk, Skk each of shape (3, 12, N+1)
    Si, Ski, Skki = i_S_Sk_matrike(pK, Li)  # for element i
    Sj, Skj, Skkj = i_S_Sk_matrike(pK, Lj)  # for element j

    # Coordinates of sub‑segment points for element i and j (global),there is N+1 point omn each element
    RTi = np.zeros((3, N + 1))
    RTj = np.zeros((3, N + 1))
    for k in range(N + 1):
        RTi[:, k] = Si[:, :, k] @ ei  # shape (3,12) * (12,) -> (3,)
        RTj[:, k] = Sj[:, :, k] @ ej

    # Pre‑allocate arrays for segment vectors and lengths
    rk = np.zeros((3, N))  # vector between endpoints of sub‑segments of j
    abs_Ljs = np.zeros(N)  # absolute length of each sub‑segment

    for kj in range(N):
        rk[:, kj] = RTj[:, kj + 1] - RTj[:, kj]
        abs_Ljs[kj] = np.sqrt(rk[:, kj].T @ rk[:, kj])

    # Transformation matrices for each sub‑segment (global → local)
    A = np.zeros((3, 3, N))
    for kj in range(N):
        # Check for special case: rk is purely in y direction
        if (rk[0, kj] == 0 and rk[2, kj] == 0):
            # Transformation matrix from global to local for a line along y‑axis
            A[:, :, kj] = np.array([[0, 1, 0],
                                    [0, 0, 1],
                                    [1, 0, 0]])
        else:
            vec_t = rk[:, kj] / np.sqrt(rk[:, kj].T @ rk[:, kj]) # this gives vector in i,j,k order, the components are tt,tn,tb,   shape (3,)

            # Local n vector (perpendicular to t in x‑z plane)
            denom = np.sqrt(vec_t[0] ** 2 + vec_t[2] ** 2)
            vec_n = np.array([-vec_t[2] / denom,  0,  vec_t[0] / denom])

            vec_b = np.cross(vec_t, vec_n)
            A[:, :, kj] = np.vstack([vec_t, vec_n, vec_b])  # rows: t, n, b

    # Compute magnetic field at each point on element i (if N=2 the return will be B in left, middle and right point )
    Bki = np.zeros((3, N + 1))
    for ki in range(N + 1):  # ki = 0..N (N+1 points)
        Bs = np.zeros(3)  # accumulated field from all sub‑segments of j
        for kj in range(N):
            # Vector from start of sub‑segment j to point ki on i, expressed in local coordinates
            RA = A[:, :, kj] @ (RTi[:, ki] - RTj[:, kj])
            # Add contribution of this sub‑segment, transformed back to global coordinates
            Bs += A[:, :, kj].T @ i_B_polje(RA, abs_Ljs[kj], Ij)
        Bki[:, ki] = Bs

    return Bki


def i_B_polje(r1, Ljs_a, Ij):
    """
    Internal function: magnetic field due to a straight filament segment.
    r1 : vector from segment start to field point (local coordinates)
    Ljs_a : length of the segment
    Ij : current in the segment
    Returns field vector [0; Hn; Hb] in local coordinates.
    """
    alfa2 = r1[1] ** 2 + r1[2] ** 2 + 1e-9  # stabilisation term to avoid singularity
    abs_r1 = np.sqrt(alfa2 + r1[0] ** 2)
    abs_r2 = np.sqrt(alfa2 + (r1[0] - Ljs_a) ** 2)

    Hn = -r1[0] * r1[2] / (abs_r1 * alfa2) + (r1[0] - Ljs_a) * r1[2] / (abs_r2 * alfa2)
    Hb = r1[0] * r1[1] / (abs_r1 * alfa2) - (r1[0] - Ljs_a) * r1[1] / (abs_r2 * alfa2)

    # Multiply by permeability factor (μ₀/(4π) gives 1e-7
    B = Ij * 1e-7 * np.array([0.0, Hn, Hb])
    return B

