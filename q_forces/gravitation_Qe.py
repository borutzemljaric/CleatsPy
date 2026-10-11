import numpy as np

def gravitation_Qe(hh, ne_c, le, g, cab):
    """
    Function Qe for calculating external gravitational forces
    """

    temp_list = []
    for k in range(hh - 1):
        for i in range(ne_c[k]):
            temp_list.append(Qg(g, cab.m, le[k, i]))

    # One single concatenation at the end (fast!)
    Qe_elementa = np.concatenate(temp_list) if temp_list else np.array([])

    Qe_elementa=np.concatenate([Qe_elementa, Qe_elementa, Qe_elementa])  # expand to 3 cables

    return Qe_elementa


# Qg function :
def Qg(g, m, L):
    """
    Returns gravitational force vector as column vector
    """
    Qg = -g * m * np.array([
        0, 0, L / 2,
        0, 0, L ** 2 / 12,
        0, 0, L / 2,
        0, 0, -L ** 2 / 12
    ])# in row vector ////  .reshape(-1, 1)  # Convert to column vector

    return Qg