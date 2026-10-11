import numpy as np

def pdotdot_lagrange(M, Qeg_0, qkqt, qv, qc, Ce):
    """Calculate accelerations and Lagrange multipliers"""
    n_constraints = len(qc)
    n_M = M.shape[0]
    n_total = n_M + n_constraints

    # print('n_constraints', n_constraints)
    # print('n_M', n_M)

    # If M is sparse, convert to dense for np.block
    if hasattr(M, 'toarray'):
        M_dense = M.toarray()
    else:
        M_dense = M

    aug_matrix = np.block([
        [M_dense, Ce.T],
        [Ce, np.zeros((n_constraints, n_constraints))]
    ])
    # print('aug_matrix',aug_matrix.shape)

    rhs = np.concatenate([
        Qeg_0 - qkqt.flatten(),
        qc
    ])

    solution = np.linalg.solve(aug_matrix, rhs)
    pdotdot = solution[:n_M]
    lambda_vec = solution[n_M:]

    return pdotdot, lambda_vec