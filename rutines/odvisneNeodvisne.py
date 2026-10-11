import numpy as np

from q_forces.constrain_C_C import constrain_C_C
from rutines.resfor_Re import resfor_Re

def neod_U(neodvisneV, konstanteU, hh, ne_c, V, U, konstante, cab, le, Bdi, pari, Qe):
    """
    Function to calculate independent coordinates
    
    Parameters:
    neodvisneV: independent coordinates
    konstanteU: constant U coordinates
    hh, ne_c, V, U, konstante, E_c, A_c, le, d_c, Bdi, pari, Qe, ne_i, tipverige, izo, DK: various parameters
    
    Returns:
    neodvisneV: updated independent coordinates
    exitflag_Niter: exit flag indicating convergence status
    """

    Niter = 0
    eps_e = 1  # epsilon for step checking, dummy to enter the loop
    eps_R = 1  # epsilon for function residual checking, dummy to enter the loop


    #test=resfor_Re(neodvisneV, konstanteU, hh, ne_c, V, U, konstante, E_c, A_c, le, d_c, Bdi, pari, Qe, ne_i, tipverige, izo,
     #        DK)

    while (eps_e > 1e-5 or eps_R > 1e-4) and Niter < 200:
        # Calculate Jacobian using finite differences
        J = jacobian(lambda ei: resfor_Re(ei, konstanteU, hh, ne_c, V, U, konstante, cab, le, Bdi, pari, Qe), neodvisneV)
        #print("Rank:", np.linalg.matrix_rank(J))
        # Calculate function residual
        F = resfor_Re(neodvisneV, konstanteU, hh, ne_c, V, U, konstante, cab, le, Bdi, pari, Qe)
        #print("Rank:", np.linalg.matrix_rank(F))
        #print("Test:", F)

        try:
            delta_e = -np.linalg.solve(J, F)
        except np.linalg.LinAlgError:
            print("Warning: J is singular, using least-squares fallback.")
            delta_e, *_ = np.linalg.lstsq(J, F, rcond=None)
            delta_e = -delta_e

        # Calculate step
        #delta_e = -np.linalg.solve(J, F)
        
        # Update independent coordinates
        neodvisneV = neodvisneV + delta_e
        
        Niter += 1  # iteration count
        
        # Calculate convergence criteria
        eps_e = np.linalg.norm(delta_e)  # epsilon for step checking
        eps_R = np.linalg.norm(F)  # epsilon for function residual checking
        
        if Niter == 200:
            exitflag_Niter = 3  # if iteration count is exceeded, go to beginning of loop
        else:
            exitflag_Niter = 0
    
    return neodvisneV, exitflag_Niter


def odvi_U(odvisneU, hh, ne_c, neodvisneV, V, U, tipverige, e_zp):
    """
    Function to calculate dependent coordinates
    
    Parameters:
    odvisneU: dependent coordinates
    hh, ne_c, neodvisneV, V, U, tipverige, e_zp: various parameters
    
    Returns:
    odvisneU: updated dependent coordinates
    """
    Niter = 0
    eps_e = 1  # epsilon for step checking, dummy to enter the loop
    eps_R = 1  # epsilon for function residual checking, dummy to enter the loop


    while (eps_e > 1e-5 or eps_R > 1e-4) and Niter < 200:
        # Calculate Jacobian using finite differences
        J = jacobian(lambda edd: constrain_C_C(edd, hh, ne_c, neodvisneV, V, U, tipverige, e_zp), odvisneU)
        #print('J',J)
        #print("J shape:", J.shape)
        #print("Rank:", np.linalg.matrix_rank(J))
        # Calculate function residual
        F = constrain_C_C(odvisneU, hh, ne_c, neodvisneV, V, U, tipverige, e_zp)
        #print(F.size)
        # Calculate step
        delta_e = -np.linalg.solve(J, F)
        
        # Update dependent coordinates
        odvisneU = odvisneU + delta_e
        
        Niter += 1  # iteration count
        
        # Calculate convergence criteria
        eps_e = np.linalg.norm(delta_e)  # epsilon for step checking
        eps_R = np.linalg.norm(F)  # epsilon for function residual checking
    
    return odvisneU


def jacobian(func, x, epsilon=1e-8):
    """
    Calculate Jacobian matrix using finite differences
    
    Parameters:
    func: function that returns a vector
    x: point at which to calculate Jacobian
    epsilon: step size for finite differences
    
    Returns:
    J: Jacobian matrix
    """
    n = len(x)
    f0 = func(x)
    m = len(f0)
    J = np.zeros((m, n))
    
    for i in range(n):
        x_perturbed = x.copy()
        x_perturbed[i] += epsilon
        f_perturbed = func(x_perturbed)
        J[:, i] = (f_perturbed - f0) / epsilon
    
    return J


