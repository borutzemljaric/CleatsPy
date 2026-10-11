import numpy as np

def i_izracun_Bi_to_Qi(ei, Bki, Ii, Ip, N):
    """
    Compute generalized electromagnetic force Qi for element i.
    """
    # Arithmetic mean of component contributions to get piecewise constant functions
    # Bki is (3, N+1), we compute average over each segment
    B_s_ki = np.zeros((3, N))
    for ki in range(N):                    # MATLAB: ki=1:N  (0-based: 0 to N-1)
        B_s_ki[:, ki] = Bki[:, ki] / 2 + Bki[:, ki+1] / 2

    # Compute generalised force
    o = np.zeros((12, 1))
    # Ensure ei is a column vector (12,1)
    if ei.ndim == 1:
        ei = ei.reshape(-1, 1)

    # Build block matrix eb = [ei o o; o ei o; o o ei]  (36 x 3)
    eb = np.vstack([
        np.hstack([ei, o, o]),
        np.hstack([o, ei, o]),
        np.hstack([o, o, ei])
    ])

    # Pre-allocate Q as (N, 12) – each row will hold a 12‑component force
    Q = np.zeros((N, 12))

    for ki in range(N):
        # Ip[:,:,ki] is (12,36), eb is (36,3), B_s_ki[:,ki] is (3,)
        # Compute Ii * Ip * eb * B_s_ki  →  (12,)
        temp = Ii * (Ip[:, :, ki] @ eb) @ B_s_ki[:, ki]
        Q[ki, :] = temp

    # Sum over all segments and transpose to column vector (12,1)
    Qi = np.sum(Q, axis=0)#.reshape(-1, 1)

    return Qi