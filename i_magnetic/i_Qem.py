import numpy as np
from i_magnetic.i_izracun_Bi import i_izracun_Bi
from i_magnetic.i_izracun_Bi_to_Qi import i_izracun_Bi_to_Qi
from i_magnetic.i_currents import i_currents




def i_Qem(t,E, cab,le_full,Ip):
    """
    Function to calculate the generalised electromagnetic force.
    """

    # Number of elements, for all three cables (each element has 12 DOFs, E is coordinates)
    Nel = len(E) // 12


    # Number of sub‑segments (minimum 2), each segment is subdivided
    # for better accuracy,
    N = cab.Nseg

    # External function that calculates Ip integrals for individual segments
    currents=i_currents(t, cab, Nel) # each element has current array(Nel,)

    # Parameter for segment endpoints
    pK = np.linspace(0, 1, N + 1)          #


    # Extract element coordinates  for each element
    ee = np.zeros((12, Nel))               # pre‑allocate for convenience
    #ee_dot = np.zeros((12, Nel))
    for i in range(Nel):
        start = i * 12
        end = start + 12
        ee[:, i] = E[start:end]            # the result is matrix 12xNel
        #ee_dot[:, i] = E_dot[start:end]

    # Sub‑segment approach
    Qem_list = []    # will hold each element's EM force contribution

    for i in range(Nel):

        # Extract the 12-row from current Ip
        start_row = i * 12
        end_row = (i + 1) * 12
        subIp = Ip[start_row:end_row, :, :]  # shape (12, 36, 2)

        Bi = np.zeros((3, N + 1))          # reset Bi to zero

        for j in range(Nel):
            if i == j:
                # No self‑interaction
                Bi += np.zeros((3, N + 1))

            else:
                # Calculate total B‑field (contribution from all elements to the current i-th element)
                # in the global coordinate system
                Bi += i_izracun_Bi(ee[:, i], ee[:, j], le_full[i], le_full[j],
                                   currents[j], N, pK)

        # Compute electromagnetic force for element i
        Qem_i = i_izracun_Bi_to_Qi(ee[:, i], Bi, currents[i], subIp, N)

        # Accumulate  vectors
        Qem_list.append(Qem_i)

    # Vertically stack all element contributions and add EM
    Qem = np.concatenate(Qem_list)  # final array


    return Qem

