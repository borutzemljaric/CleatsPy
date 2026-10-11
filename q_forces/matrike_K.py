import numpy as np
from scipy import sparse

def Kl_L1(e, EA, L):
    """
    Linear stiffness matrix for beam element
    """
    D = np.sqrt((e[6] - e[0])**2 + (e[7] - e[1])**2 + (e[8] - e[2])**2)
    Kl = EA * (D - L) * Kl2(L) @ e
    return Kl

def Kl2(L):
    """
    Helper function for Kl_L1 - stiffness matrix template
    """
    return np.array([
        [6/(5*L**2), 0, 0, 1/(10*L), 0, 0, -6/(5*L**2), 0, 0, 1/(10*L), 0, 0],
        [0, 6/(5*L**2), 0, 0, 1/(10*L), 0, 0, -6/(5*L**2), 0, 0, 1/(10*L), 0],
        [0, 0, 6/(5*L**2), 0, 0, 1/(10*L), 0, 0, -6/(5*L**2), 0, 0, 1/(10*L)],
        [1/(10*L), 0, 0, 2/15, 0, 0, -1/(10*L), 0, 0, -1/30, 0, 0],
        [0, 1/(10*L), 0, 0, 2/15, 0, 0, -1/(10*L), 0, 0, -1/30, 0],
        [0, 0, 1/(10*L), 0, 0, 2/15, 0, 0, -1/(10*L), 0, 0, -1/30],
        [-6/(5*L**2), 0, 0, -1/(10*L), 0, 0, 6/(5*L**2), 0, 0, -1/(10*L), 0, 0],
        [0, -6/(5*L**2), 0, 0, -1/(10*L), 0, 0, 6/(5*L**2), 0, 0, -1/(10*L), 0],
        [0, 0, -6/(5*L**2), 0, 0, -1/(10*L), 0, 0, 6/(5*L**2), 0, 0, -1/(10*L)],
        [1/(10*L), 0, 0, -1/30, 0, 0, -1/(10*L), 0, 0, 2/15, 0, 0],
        [0, 1/(10*L), 0, 0, -1/30, 0, 0, -1/(10*L), 0, 0, 2/15, 0],
        [0, 0, 1/(10*L), 0, 0, -1/30, 0, 0, -1/(10*L), 0, 0, 2/15]
    ])

# def Kl_L2a(e, E, A, L):
#     """
#     Nonlinear stiffness matrix variant A for beam element
#     """
#     ea1, ea2, ea3, ea4, ea5, ea6, ea7, ea8, ea9, ea10, ea11, ea12 = e.flatten()
#
#     Kl = E * A * Kl1(L, ea1, ea2, ea3, ea4, ea5, ea6, ea7, ea8, ea9, ea10, ea11, ea12) @ e
#     return Kl
#
# def Kl1(L, ea1, ea2, ea3, ea4, ea5, ea6, ea7, ea8, ea9, ea10, ea11, ea12):
#     """
#     Helper function for Kl_L2a - complex nonlinear stiffness matrix
#     """
#     # This is an extremely complex matrix - in practice, you might want to
#     # generate this symbolically or use a more efficient approach
#     # For now, I'll provide a placeholder structure
#     Kl1_matrix = np.zeros((12, 12))
#
#     # The actual implementation would contain the massive matrix from the MATLAB code
#     # This is simplified - you should implement the full matrix from the MATLAB code
#     for i in range(12):
#         for j in range(12):
#             # Placeholder - implement the actual complex expressions from MATLAB
#             Kl1_matrix[i, j] = 0  # Replace with actual expressions
#
#     return Kl1_matrix
#
# def Kl_L2b(e, E, A, L):
#     """
#     Nonlinear stiffness matrix variant B for beam element
#     """
#     ea1, ea2, ea3, ea4, ea5, ea6, ea7, ea8, ea9, ea10, ea11, ea12 = e.flatten()
#
#     Kl_2b = 0.5 * E * A * L * (Kl1(L, ea1, ea2, ea3, ea4, ea5, ea6, ea7, ea8, ea9, ea10, ea11, ea12) - Kl2(L)) @ e
#     return Kl_2b
#
# def Kl_L13(e, E, A, L):
#     """
#     Integral-based stiffness matrix for axial deformation
#     """
#     # Shape function first derivative
#     def Sk(k, L):
#         return np.array([
#             [6*k**2 - 6*k, 0, 0, L*(3*k**2 - 4*k + 1), 0, 0, 6*k - 6*k**2, 0, 0, -L*(2*k - 3*k**2), 0, 0],
#             [0, 6*k**2 - 6*k, 0, 0, L*(3*k**2 - 4*k + 1), 0, 0, 6*k - 6*k**2, 0, 0, -L*(2*k - 3*k**2), 0],
#             [0, 0, 6*k**2 - 6*k, 0, 0, L*(3*k**2 - 4*k + 1), 0, 0, 6*k - 6*k**2, 0, 0, -L*(2*k - 3*k**2)]
#         ])
#
#     def Depsde(k, L, e):
#         Sk_val = Sk(k, L)
#         return (e.T @ Sk_val.T @ Sk_val) / (L**2)
#
#     def epsi(k, L, e):
#         Sk_val = Sk(k, L)
#         return ((e.T @ Sk_val.T @ Sk_val @ e) / (L**2) - 1) / 2
#
#     # Integration points
#     N = 5
#     s = np.linspace(1/(N+1), 1-1/(N+1), N)
#
#     # Numerical integration
#     IA = np.zeros(12)
#     for si in s:
#         IA += epsi(si, L, e) * Depsde(si, L, e)
#
#     IA = (E * A * L / (2 * N)) * IA
#     return IA.T  # Return as column vector

def Kt_T1(e, EI, L):
    """
    Linear torsional stiffness matrix
    """
    Kt = EI  * L * I1(L) @ e
    return Kt

def I1(L):
    """
    Helper function for Kt_T1 - torsional stiffness matrix template
    """
    return np.array([
        [12/L**4, 0, 0, 6/L**3, 0, 0, -12/L**4, 0, 0, 6/L**3, 0, 0],
        [0, 12/L**4, 0, 0, 6/L**3, 0, 0, -12/L**4, 0, 0, 6/L**3, 0],
        [0, 0, 12/L**4, 0, 0, 6/L**3, 0, 0, -12/L**4, 0, 0, 6/L**3],
        [6/L**3, 0, 0, 4/L**2, 0, 0, -6/L**3, 0, 0, 2/L**2, 0, 0],
        [0, 6/L**3, 0, 0, 4/L**2, 0, 0, -6/L**3, 0, 0, 2/L**2, 0],
        [0, 0, 6/L**3, 0, 0, 4/L**2, 0, 0, -6/L**3, 0, 0, 2/L**2],
        [-12/L**4, 0, 0, -6/L**3, 0, 0, 12/L**4, 0, 0, -6/L**3, 0, 0],
        [0, -12/L**4, 0, 0, -6/L**3, 0, 0, 12/L**4, 0, 0, -6/L**3, 0],
        [0, 0, -12/L**4, 0, 0, -6/L**3, 0, 0, 12/L**4, 0, 0, -6/L**3],
        [6/L**3, 0, 0, 2/L**2, 0, 0, -6/L**3, 0, 0, 4/L**2, 0, 0],
        [0, 6/L**3, 0, 0, 2/L**2, 0, 0, -6/L**3, 0, 0, 4/L**2, 0],
        [0, 0, 6/L**3, 0, 0, 2/L**2, 0, 0, -6/L**3, 0, 0, 4/L**2]
    ])

# def Kt_T2(e, E, d, L):
#     """
#     Nonlinear torsional stiffness matrix
#     """
#     ea1, ea2, ea3, ea4, ea5, ea6, ea7, ea8, ea9, ea10, ea11, ea12 = e.flatten()
#
#     I = (np.pi * d**4 * 1e-12) / 64  # Moment of inertia
#
#     # Complex expressions for f_crt and kapa_crt
#     # These would need to be implemented based on the MATLAB expressions
#     f_crt = 1.0  # Placeholder - implement actual expression
#     kapa_crt = 1.0  # Placeholder - implement actual expression
#
#     Kt_T2 = E * I * L * ((I1(L) / f_crt**4) - 2 * (kapa_crt**2 / f_crt**2) * I2(L)) @ e
#     return Kt_T2
#
# def I2(L):
#     """
#     Helper function for Kt_T2 - additional torsional matrix
#     """
#     return np.array([
#         [6/(5*L**2), 0, 0, 1/(10*L), 0, 0, -6/(5*L**2), 0, 0, 1/(10*L), 0, 0],
#         [0, 6/(5*L**2), 0, 0, 1/(10*L), 0, 0, -6/(5*L**2), 0, 0, 1/(10*L), 0],
#         [0, 0, 6/(5*L**2), 0, 0, 1/(10*L), 0, 0, -6/(5*L**2), 0, 0, 1/(10*L)],
#         [1/(10*L), 0, 0, 2/15, 0, 0, -1/(10*L), 0, 0, -1/30, 0, 0],
#         [0, 1/(10*L), 0, 0, 2/15, 0, 0, -1/(10*L), 0, 0, -1/30, 0],
#         [0, 0, 1/(10*L), 0, 0, 2/15, 0, 0, -1/(10*L), 0, 0, -1/30],
#         [-6/(5*L**2), 0, 0, -1/(10*L), 0, 0, 6/(5*L**2), 0, 0, -1/(10*L), 0, 0],
#         [0, -6/(5*L**2), 0, 0, -1/(10*L), 0, 0, 6/(5*L**2), 0, 0, -1/(10*L), 0],
#         [0, 0, -6/(5*L**2), 0, 0, -1/(10*L), 0, 0, 6/(5*L**2), 0, 0, -1/(10*L)],
#         [1/(10*L), 0, 0, -1/30, 0, 0, -1/(10*L), 0, 0, 2/15, 0, 0],
#         [0, 1/(10*L), 0, 0, -1/30, 0, 0, -1/(10*L), 0, 0, 2/15, 0],
#         [0, 0, 1/(10*L), 0, 0, -1/30, 0, 0, -1/(10*L), 0, 0, 2/15]
#     ])
#
# def Kt_T13(e, E, d, L):
#     """
#     Integral-based torsional stiffness matrix
#     """
#     # Shape functions
#     def Sk(k, L):
#         return np.array([
#             [6*k**2 - 6*k, 0, 0, L*(3*k**2 - 4*k + 1), 0, 0, 6*k - 6*k**2, 0, 0, -L*(2*k - 3*k**2), 0, 0],
#             [0, 6*k**2 - 6*k, 0, 0, L*(3*k**2 - 4*k + 1), 0, 0, 6*k - 6*k**2, 0, 0, -L*(2*k - 3*k**2), 0],
#             [0, 0, 6*k**2 - 6*k, 0, 0, L*(3*k**2 - 4*k + 1), 0, 0, 6*k - 6*k**2, 0, 0, -L*(2*k - 3*k**2)]
#         ])
#
#     def Skk(k, L):
#         return np.array([
#             [12*k - 6, 0, 0, L*(6*k - 4), 0, 0, 6 - 12*k, 0, 0, L*(6*k - 2), 0, 0],
#             [0, 12*k - 6, 0, 0, L*(6*k - 4), 0, 0, 6 - 12*k, 0, 0, L*(6*k - 2), 0],
#             [0, 0, 12*k - 6, 0, 0, L*(6*k - 4), 0, 0, 6 - 12*k, 0, 0, L*(6*k - 2)]
#         ])
#
#     I_val = (np.pi * d**4 * 1e-12) / 64
#
#     def f1(k, L, e):
#         Sk_val = Sk(k, L)
#         Skk_val = Skk(k, L)
#         cross_prod = np.cross(Sk_val @ e, Skk_val @ e)
#         return np.sqrt(np.sum(cross_prod**2)) / L**3
#
#     def g1(k, L, e):
#         Sk_val = Sk(k, L)
#         return (np.sqrt(np.sum((Sk_val @ e)**2)) / L)**3
#
#     def odvod(fun, x):
#         dx = np.finfo(float).eps**(1/3)
#         nx = len(x)
#         nf = len(fun(x))
#         J = np.zeros((nf, nx))
#
#         for n in range(nx):
#             delta = np.zeros(nx)
#             delta[n] = dx
#             dF = fun(x + delta) - fun(x - delta)
#             J[:, n] = dF.flatten() / dx / 2
#
#         return J
#
#     def dkdk(k, L, e):
#         f1_val = f1(k, L, e)
#         g1_val = g1(k, L, e)
#
#         df1_de = odvod(lambda e_vec: f1(k, L, e_vec), e)
#         dg1_de = odvod(lambda e_vec: g1(k, L, e_vec), e)
#
#         return (1 / g1_val**2) * (g1_val * df1_de - f1_val * dg1_de)
#
#     # Integration points
#     N = 3
#     s = np.linspace(1/(N+1), 1-1/(N+1), N)
#
#     # Numerical integration
#     IB = np.zeros(12)
#     for si in s:
#         f1_val = f1(si, L, e)
#         g1_val = g1(si, L, e)
#         IB += (f1_val / g1_val) * dkdk(si, L, e).flatten()
#
#     IB = (E * I_val * L / N) * IB
#     return IB

# Export all functions
#__all__ = [
#    'Kl_L1', 'Kl2', 'Kl_L2a', 'Kl1', 'Kl_L2b', 'Kl_L13',
#    'Kt_T1', 'I1', 'Kt_T2', 'I2', 'Kt_T13'
#]