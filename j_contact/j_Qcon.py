import numpy as np

from j_contact.j_izracun_Qcon import j_izracun_Qcon



def j_Qcon(t,E, E_dot, cab, le_full):
    """
    Function to calculate the generalised contact force.
    """

    # Number of elements, for all three cables (each element has 12 DOFs, E is coordinates)
    Nel = len(E) // 12

    # Number of contact points, each segment is subdivided
    N_con = cab.Ncon  #

    # This creates N_con evenly spaced points starting at 1/(2N) and ending at 1 - 1/(2N)
    pK_con = np.linspace(1 / (2 * N_con), 1 - 1 / (2 * N_con), N_con)
    #


    # Extract element coordinates  for each element
    ee= np.zeros((12, Nel))  # pre‑allocate for convenience
    ee_dot = np.zeros((12, Nel))               # pre‑allocate for convenience
    #ee_dot = np.zeros((12, Nel))
    for i in range(Nel):
        start = i * 12
        end = start + 12
        ee[:, i] = E[start:end]
        ee_dot[:, i] = E_dot[start:end]            #


    # Sub‑segment approach
    Qcon_list = []    # will hold each element's contact force contribution

    for i in range(Nel):
        # Extract the 12-row from current Ip
        start_row = i * 12
        end_row = (i + 1) * 12
        #subIp = Ip[start_row:end_row, :, :]  # shape (12, 36, 2)

        Qcon_i = np.zeros((12, 1))          # reset to zero

        for j in range(Nel):
            if i == j:
                # No self‑interaction
                Qcon_i += np.zeros((12, 1))

            else:
                Qcon_i += j_izracun_Qcon(ee[:, i], ee[:, j], le_full[i], le_full[j],
                                   N_con, pK_con, cab ,ee_dot[:, i], ee_dot[:, j])

        # Accumulate  vectors
        Qcon_list.append(Qcon_i)

    # Vertically stack all element contributions
    Qcon = np.concatenate(Qcon_list).flatten()  # final array

    return Qcon

