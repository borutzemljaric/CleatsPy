import numpy as np
import os

from SaveAndLoad.SaveLoad import load_for_plot
from plots.findCoordIndex import findCoordIndex

def readForces(filename):
    """
    Read force data from file.
    """

    # Check if file exists
    if not os.path.exists(filename):
        raise FileNotFoundError(f"File {filename} not found")

    # Load  data
    hh, ne_c, tipverige, zacpog, t, _, _, izhod_Flevo, izhod_Fdesno= load_for_plot(filename)
        # Call function to get indices, the force results have different array as coordinates
    _, cleaveX1 = findCoordIndex(tipverige, hh, ne_c, 3)


    # Reshape into exactly 3 rows.
    # NumPy automatically calculates the number of columns: len(arr) // 3
    reshaped = cleaveX1.reshape(3, -1)
    # Subtract 3 from the last column (column index -1)
    #reshaped[:, -1] -= 3  # this is cleaveX1 just transformed in 3xn array

    split1L, split2L, split3L = np.split(izhod_Flevo, 3, axis=1)  # split matrices by cables 1,2,3
    split1D, split2D, split3D = np.split(izhod_Fdesno, 3, axis=1)
    # Each has shape (n, m//3)

    # Get the number of rows (all splits have the same shape)
    rows = split1L.shape[0]

    # Create a matrix of zeros with 3 columns
    zeros = np.zeros((rows, 3))

    # 1. Add 3 zero columns at the END of split L
    split1L_padded = np.hstack((split1L, zeros))
    split2L_padded = np.hstack((split2L, zeros))
    split3L_padded = np.hstack((split3L, zeros))

    # 2. Add 3 zero columns at the BEGINNING of split D
    split1D_padded = np.hstack((zeros, split1D))
    split2D_padded = np.hstack((zeros, split2D))
    split3D_padded = np.hstack((zeros, split3D))

    # 3. Sum both padded matrices together
    split1_sum = split1L_padded + split1D_padded  # cable 1 - the results are sum of forces by xyz in each node
    split2_sum = split2L_padded + split2D_padded  # cable 2
    split3_sum = split3L_padded + split3D_padded  # cable 3

    # Construct reaction forces in nodes-cleats

    idx = reshaped[0, :]  # 1D array of integer indices

    # Directly stack the columns horizontally
    FXcleave = np.hstack([split1_sum[:, idx], split2_sum[:, idx], split3_sum[:, idx]])
    FYcleave = np.hstack([split1_sum[:, idx + 1], split2_sum[:, idx + 1], split3_sum[:, idx + 1]])
    FZcleave = np.hstack([split1_sum[:, idx + 2], split2_sum[:, idx + 2], split3_sum[:, idx + 2]])


    #the output columns are forces by components that cable put on cleats, (rows are time)


    return t, FXcleave,FYcleave,FZcleave

