import numpy as np

def ureditevparov(pari):
    """
    Function to arrange pairs
    
    Parameters:
    pari: numpy array with pairs information
    
    Returns:
    par_K: array for pairs where pari(k,0) == -1
    par_U: array for pairs where pari(k,0) != -1
    """

    N_parov = (pari.shape[0])

    par_K_list = []
    par_U_list = []
    
    for k in range(N_parov):
        if pari[k, 0] == -1:
            # First column - which independent belongs to dependent
            # Second column - dependent coordinates  
            # Third column - index where dependent is located in original order
            par_K_list.append([pari[k, 0], pari[k, 1], k])
        else:
            par_U_list.append([pari[k, 0], pari[k, 1], k])
    
    # Convert lists to numpy arrays
    if par_K_list:
        par_K = np.array(par_K_list)
    else:
        par_K = np.array([])
        
    if par_U_list:
        par_U = np.array(par_U_list)
    else:
        par_U = np.array([])
    
    return par_K, par_U