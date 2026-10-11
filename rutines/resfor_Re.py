from q_forces.matrike_K import Kl_L1,Kt_T1
from q_forces.matrike_h import *
from h_qsqb.h_Qs import h_Qs
from h_qsqb.h_Qb import h_Qb



def resfor_Re(ei, konstanteU, hh, ne_c, V, U, konstante, cab, le,
               Bdi, pari_dict, Qe):
    """
    Function to calculate residual forces
    """
    E = np.hstack([ei, ei[pari_dict['V']], konstanteU])  # composition of vector with order: independent, pari(dependent), constants
    index = np.hstack([V, pari_dict['U'], konstante])#.flatten()
    E[index] = E#.flatten()  # placement in basic order

    K_elementa=Ke_vect(hh, ne_c, E, cab, le)
    Re = Qe - K_elementa
    # Re is still in correct order 1-n, so rearrange rows according to order V U
    reordered_indices = np.hstack([V, U]).flatten()
    #Re = sparse.csr_matrix(Re[reordered_indices, :])
    Re = Re[reordered_indices]

    #print(Qw)
    #print('Re',Re)
    #print('Re', Re.shape)
    #breakpoint()
    RTcrt = -Bdi.T @ Re

    return RTcrt


# assemble vector Ke, generalised  elastic forces
def Ke_vect(hh, ne_c, E, cab, le):
    # three conductors, so three times and then vector
    offset = 0
    K_elementa0, vars0 =Qe_Ke_vectors(hh, ne_c, E, cab, le,offset)
    offset = vars0
    K_elementa1, vars1 =Qe_Ke_vectors(hh, ne_c, E, cab, le,offset)
    offset = vars1
    K_elementa2, vars2 =Qe_Ke_vectors( hh, ne_c, E, cab, le,offset)

    K_elementa = np.concatenate([K_elementa0, K_elementa1, K_elementa2])  # construct vector

    return K_elementa



# helper function for Ke, this goes through only one conductor
def Qe_Ke_vectors(hh, ne_c, E, cab, le,start_offset):

    ii = start_offset  # local offset for this conductor
    K_elementa = []  # predefined empty matrix

    for k in range(hh - 1):
        for i in range(ne_c[k]):
            # Get coordinates for current element, from-to by 12 coordinates
            e = E[(ii  * 12):(ii  * 12 + 12)]
            #K_elementa.append(Kl_L1(e,cab.EA, le[k, i]) + Kt_T1(e, cab.EI, le[k, i]))
            K_elementa.append(h_Qs(e, cab.EA, le[k, i]) + h_Qb(e, cab.EI, le[k, i]))
            ii += 1  # counting conductor elements

    # Convert lists to arrays if needed
    K_elementa = np.hstack(K_elementa)

    return K_elementa,ii





