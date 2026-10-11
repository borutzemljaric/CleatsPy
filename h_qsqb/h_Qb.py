import numpy as np

from h_qsqb.h_Sk_matrike import  h_Sk_matrike
from h_qsqb.h_b1_matrike import h_b1_matrike
from h_qsqb.h_b2_matrike import h_b2_matrike
from h_qsqb.h_dcde import h_dcde


def h_Qb(e, EI, L):
    """
    Function to calculate strain energy due to bending
    Based on Shabana-Gertmaister-Zemljaric (T5) formulation

    Returns:
    Qb: bending energy vector
    """
    e = e.reshape(-1, 1)  # into column vector

    # Gauss-Legendre quadrature integration
    k = np.array([-0.774597 / 2 + 0.5, 0.5, 0.774597 / 2 + 0.5])  # # 3 quadrature points
    w = np.array([0.555556, 0.888889, 0.555556])  # corresponding quadrature coefficients

    # Moment of inertia with unit conversion (mm to m)
    # Ix = (np.pi * d**4 * 1e-6) / 64

    # Calculate matrices needed for continuation
    o = np.zeros(12)

    eb = np.block([[e.T, o, o], [o, e.T, o], [o, o, e.T]])

    # Get matrices for each integration point
    Skk = h_Sk_matrike(k, L)  # returns 3D array for integration points
    b1 = h_b1_matrike(k, L)  # returns 3D array for integration points
    B2 = h_b2_matrike(k, L)  # returns 3D array for integration points
    dcde = h_dcde()  # returns derivatives wrt coordinates 1-12

    # Initialize integral calculation
    del_I = np.zeros((12, 3))

    # Calculate integral at three points
    for i in range(3):
        g = (e.T @ b1[:, :, i] @ e) ** (3 / 2)
        dgde = 3 * (e.T @ b1[:, :, i] @ e) ** (1 / 2) * e.T @ b1[:, :, i].T

        f_del = e.T @ Skk[:, :, i].T @ eb @ B2[:, :, i] @ eb.T @ Skk[:, :, i] @ e / L ** 4

        f = f_del ** (1 / 2)

        dfde_del = h_dfde(e, Skk[:, :, i], eb, B2[:, :, i], dcde)
        dfde = (f_del ** (-1 / 2)) * dfde_del / (2 * L ** 4) if f_del != 0 else np.zeros(12)

        # Handle singularity in numerical calculation
        if f == 0:
            dfde = np.zeros(12)

        del_I[:, i] = (f * (dfde * g - f * dgde) * w[i] * EI * L / (2 * g ** 3)).flatten()

    Qb = np.sum(del_I, axis=1)  # sum along rows to get 'column' vector

    return Qb  # .reshape(-1, 1)


def h_dfde(e, Skk, eb, B2, dcde):
    """
    Helper function to calculate derivative for h_Qb
    """
    e = np.asarray(e).flatten()  # make it 1‑D
    dfde_del = np.zeros(12)

    for ii in range(12):
        term1 = Skk[:, ii].T @ eb @ B2 @ eb.T @ Skk @ e
        term2 = e.T @ Skk.T @ dcde[:, :, ii].T @ B2 @ eb.T @ Skk @ e
        term3 = e.T @ Skk.T @ eb @ B2 @ dcde[:, :, ii] @ Skk @ e
        term4 = e.T @ Skk.T @ eb @ B2 @ eb.T @ Skk[:, ii]

        dfde_del[ii] = term1 + term2 + term3 + term4

    return dfde_del


