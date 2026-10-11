import numpy as np
from i_damping.Qd_element import Qd_element

# assemble vector Qd, generalised  damping forces
def Qd_vect(hh, ne_c, cab, le, E, Edot):

    """
    Compute internal damping force for ALL elements.
    """
    # three conductors, so three times and then creation of the full vector
    offset = 0
    Qd_elementa0, vars0 =Qd_single(hh, ne_c, E, Edot, cab, le,offset)
    offset = vars0
    Qd_elementa1, vars1 =Qd_single(hh, ne_c, E, Edot, cab, le,offset)
    offset = vars1
    Qd_elementa2, vars2 =Qd_single( hh, ne_c, E, Edot, cab, le,offset)

    Qd_full = np.concatenate([Qd_elementa0, Qd_elementa1, Qd_elementa2])  # construct vector

    return Qd_full



# helper function for Ke, this goes through only one conductor
def Qd_single(hh, ne_c, E, Edot, cab, le,start_offset):

    ii = start_offset  # local offset for this conductor
    Qd_elementa = []  # predefined empty matrix

    for k in range(hh - 1):
        for i in range(ne_c[k]):
            # Get coordinates for current element, from-to by 12 coordinates
            e = E[(ii  * 12):(ii  * 12 + 12)]
            edot = Edot[(ii * 12):(ii * 12 + 12)]
            Qd_elementa.append(Qd_element(cab, le[k, i], e, edot))
            ii += 1  # counting conductor elements

    # Convert lists to arrays if needed
    Qd_elementa = np.hstack(Qd_elementa)

    return Qd_elementa,ii