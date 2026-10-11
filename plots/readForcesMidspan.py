import numpy as np
import os

from SaveAndLoad.SaveLoad import load_for_plot
from plots.findCoordIndex import findCoordIndex

def readForcesMidspan(filename):
    """
    Read force data from file.
    """

    # Check if file exists
    if not os.path.exists(filename):
        raise FileNotFoundError(f"File {filename} not found")

    # Load  data
    hh, ne_c, tipverige, zacpog, t, _, _, izhod_Flevo, izhod_Fdesno= load_for_plot(filename)
        # Call function to get indices, the force results have different array as coordinates
    X1mid, _ = findCoordIndex(tipverige, hh, ne_c, 12)  # read coordinates in midspan
    X1mid = (X1mid / 4).astype(int)  # divided by 4 to get force index

    # Reshape into exactly 3 rows.
    # NumPy automatically calculates the number of columns: len(arr) // 3
    reshaped = X1mid.reshape(3, -1)



    # Assemble forces by components
    cols = []
    cols_y = []
    cols_z = []
    rows, m = reshaped.shape


    for i in range(rows):
        for j in range(m):
            idx = reshaped[i, j]
            cols.append(izhod_Flevo[:, idx])
            cols_y.append(izhod_Flevo[:, idx+1])
            cols_z.append(izhod_Flevo[:, idx+2])

    # Stack columns horizontally -> shape (t, rows * m)
    FXmid = np.column_stack(cols)
    FYmid = np.column_stack(cols_y)
    FZmid= np.column_stack(cols_z)


    #the output columns are matrices of forces by components in elements ends, (rows are time)


    return t, FXmid,FYmid,FZmid

