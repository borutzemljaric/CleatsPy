import numpy as np
import matplotlib.pyplot as plt

from rutines.ureditevparov import ureditevparov
from rutines.odvisneNeodvisne import odvi_U, neod_U
from q_forces.gravitation_Qe import gravitation_Qe


def calculateStaticCoordinates( cab, hh, ne_c, le,tipverige, e_zp, Bdi,U,V,konstante,pari):
    """
    Calculate static ANCF coordinates
    """
    g = 9.81  # m/s²

    # Arrange pairs
    par_K, par_U = ureditevparov(pari)
    #print(pari)
    #print(pari.shape)

    # Initialize coordinate arrays
    zacpog = np.zeros_like(e_zp)

    # Store independent coordinates
    neodvisneV = e_zp[V]

    # Store dependent coordinates with constants
    odvisneU = np.hstack([
        e_zp[par_U[:, 1].astype(int)],
        e_zp[np.array(konstante).astype(int)]
    ])


    # Reorder to original order
    red = np.hstack([par_U[:, 2], par_K[:, 2]]).astype(int)
    #odvisneU_temp = odvisneU.copy()
    odvisneU[red] = odvisneU#odvisneU_temp
    
    konstanteU = e_zp[np.array(konstante).astype(int) ].copy()
    #print('odvisneU',odvisneU)
    #print('neodvisneV',neodvisneV)

    # Prepare pairs dictionary
    pari_dict = {}
    pari_dict['V'] = par_U[:, 0]
    pari_dict['U'] = par_U[:, 1]
    #print('e_zp',e_zp)

    # Calculate gravitational forces
    Qe_g = gravitation_Qe(hh, ne_c, le, g, cab)
    #print(Qe_g)

    # STEP 1:
    print('Calculation of zero state started')
    # STATIC CALCULATION ALGORITHM

    # Calculate dependent coordinates
    odvisneU = odvi_U(odvisneU, hh, ne_c, neodvisneV, V, U, tipverige, e_zp)


    # Calculate independent coordinates
    neodvisneV, exitflag_Niter = neod_U(
            neodvisneV, konstanteU, hh, ne_c, V, U, konstante,
            cab, le, Bdi, pari_dict, Qe_g  )


    if exitflag_Niter==0:
        print('Solution found')

    # Prepare final results

    odvisneU = odvi_U(odvisneU, hh, ne_c, neodvisneV, V, U, tipverige, e_zp)
    zacpog[U ] = odvisneU.flatten()
    zacpog[V ] = neodvisneV.flatten()
    zacpog[konstante] = konstanteU.flatten()


    # Plot results
    #plot_results(zacpog, e_zp)
    
    return zacpog, Qe_g













