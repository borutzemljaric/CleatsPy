import numpy as np


def constrain_C_C(edd, hh, ne_c, neodvisneV, V, U, tipverige, e_zp):
    """
    Constraints function -

    Parameters:
    edd : array
        Dependent coordinates
    hh : int
        Number of segments
    ne_c : array
        Number of elements per segment
    neodvisneV : array
        Independent variables
    V : array
        Independent variable indices
    U : array
        Dependent variable indices
    tipverige : array
        Chain type indicators
    e_zp : array
        Fixed coordinates

    Returns:
    numpy.ndarray
        Constraint vector Cc
    """

    # Initialize ed array (size needs to be determined based on your data)
    # You may need to adjust the initial size based on your maximum expected index
    max_size = len(V) + len(U)
    ed = np.zeros(max_size)

    # Set independent variables
    if len(V) > 0:
        ed[V] = neodvisneV.flatten()

    # Set dependent variables
    if len(U) > 0:
        ed[U] = edd.flatten()

    # Initialize constraint vector
    temp_list = []  # Use a Python list for accumulation

    a = 0

    for k in range(hh - 1):
        for i in range(ne_c[k]):

            if k == 0 and i == 0:  # A. column-conductor connection left
                constraints = np.array([
                    ed[a] - e_zp[a],# left point segment ground-conductor
                    ed[a + 1] - e_zp[a + 1],
                    ed[a + 2] - e_zp[a + 2],
                    ed[a + 6] - ed[a + 12], # right point segment first
                    ed[a + 7] - ed[a + 13],
                    ed[a + 8] - ed[a + 14],
                    ed[a + 9] - ed[a + 15],
                    ed[a + 10] - ed[a + 16],
                    ed[a + 11] - ed[a + 17]
                ])
                temp_list.append(constraints)  #
                a += 12

            elif i >= 1 and i != ne_c[k] - 1:  # B. conductor-conductor
                constraints = np.array([
                    ed[a + 6] - ed[a + 12],
                    ed[a + 7] - ed[a + 13],
                    ed[a + 8] - ed[a + 14],
                    ed[a + 9] - ed[a + 15],
                    ed[a + 10] - ed[a + 16],
                    ed[a + 11] - ed[a + 17]
                ])
                temp_list.append(constraints)
                a += 12

            elif i >= 1 and i == ne_c[k] - 1 and k != hh - 2:  # C. conductor-conductor-cleaves connection
                current_tip = tipverige[k + 1]

                if current_tip == 'CC':  # cable cleave
                    constraints = np.array([
                        ed[a + 6] - e_zp[a + 6],  # conductor-fixed point -conductor
                        ed[a + 7] - e_zp[a + 7],
                        ed[a + 8] - e_zp[a + 8],
                        ed[a + 6] - ed[a + 12],  # conductor-conductor continues
                        ed[a + 7] - ed[a + 13],
                        ed[a + 8] - ed[a + 14],
                        ed[a + 9] - ed[a + 15],  # conductor-conductor inclination
                        ed[a + 10] - ed[a + 16],
                        ed[a + 11] - ed[a + 17],
                        ed[a + 18] - ed[a + 24],  # joint to the next right element
                        ed[a + 19] - ed[a + 25],
                        ed[a + 20] - ed[a + 26],
                        ed[a + 21] - ed[a + 27],
                        ed[a + 22] - ed[a + 28],
                        ed[a + 23] - ed[a + 29],
                    ])
                    temp_list.append(constraints)
                    a += 24


                elif current_tip == 'IR':  #
                    constraints = np.array([
                        ed[a + 6] - e_zp[a + 6],  # conductor-fixed point -conductor
                        ed[a + 7] - e_zp[a + 7],
                        ed[a + 8] - e_zp[a + 8],
                        ed[a + 6] - ed[a + 12],  # conductor-conductor continuous
                        ed[a + 7] - ed[a + 13],
                        ed[a + 8] - ed[a + 14],
                        ed[a + 9] - ed[a + 15],  # conductor-conductor inclination
                        ed[a + 10] - ed[a + 16],
                        ed[a + 11] - ed[a + 17],
                        ed[a + 18] - ed[a + 24],  # joint to the next right element
                        ed[a + 19] - ed[a + 25],
                        ed[a + 20] - ed[a + 26],
                        ed[a + 21] - ed[a + 27],
                        ed[a + 22] - ed[a + 28],
                        ed[a + 23] - ed[a + 29],
                    ])
                    temp_list.append(constraints)
                    a += 24

            elif i >= 1 and i == ne_c[k] - 1 and k == hh - 2:  # D. conductor-right hand connection
                constraints = np.array([
                    ed[a + 6] - e_zp[a + 6],
                    ed[a + 7] - e_zp[a + 7],
                    ed[a + 8] - e_zp[a + 8]
                ])
                temp_list.append(constraints)

    #
    if temp_list:
        Cq_elementa = np.concatenate(temp_list)
    else:
        Cq_elementa = np.array([])

    # Now  triple concatenation as we have 3 cables
    Cq_elementa = np.concatenate([Cq_elementa, Cq_elementa, Cq_elementa])  # Tuple or list, both fine
    #Cq_elementa = Cq_elementa.reshape(-1, 1)  #


    return Cq_elementa #