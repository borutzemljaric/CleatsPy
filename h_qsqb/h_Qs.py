import numpy as np

from h_qsqb.h_b1_matrike import h_b1_matrike

def h_Qs(e, EA, L):
    """
    Function to calculate strain energy due to elongation
    Based on Shabana-Gertmaister-Zemljaric (L5 model) formulation - Strain energy


    Returns:
    Qs: strain energy vector
    """
    e = e.reshape(-1, 1)  # reshape into column vector
    # Gauss-Legendre quadrature integration
    k = np.array([-0.90618 / 2 + 0.5, -0.538469 / 2 + 0.5, 0.5,
                  0.538469 / 2 + 0.5, 0.90618 / 2 + 0.5])  # 5 quadrature points
    w = np.array([0.236927, 0.478629, 0.568889, 0.478629, 0.236927])  # quadrature coefficients

    # Calculate b matrices
    b1 = h_b1_matrike(k, L)

    skupno = EA * L / 2
    del_I1 = np.zeros((12, 5))

    for i in range(5):
        b1_i = b1[:, :, i]
        term = e.T @ b1_i.T * ((e.T @ b1_i @ e) - 1.)
        del_I1[:, i] = (term * w[i] / 2).flatten()

    Qsk1 = np.sum(del_I1, axis=1)  # sum along rows
    Qs = Qsk1 * skupno  # .reshape(-1, 1)

    return Qs  # shape (12,)


