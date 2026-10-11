import numpy as np

def i_le_fullList(ne_c, le):
    # List to collect all vertically stacked parts

    hh_iter=len(ne_c)-1  # minus 1 because the ce_c has empty last entry

    all_parts = []
    for i in range(hh_iter):
        # Call the function once for this row
        arr_1d = np.full(ne_c[i], le[i,0])  #
        all_parts.append(arr_1d)  #

    # Vertically concatenate all blocks along axis 0
    le_final = np.concatenate(all_parts)  # final array

    le_full = np.concatenate([le_final, le_final, le_final])  # for all three conductors
    return le_full