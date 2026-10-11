import numpy as np
import os

from SaveAndLoad.SaveLoad import load_for_plot
from plots.findCoordIndex import findCoordIndex


def readCoordinateALL(filename):
    """
    Read all coordinate data from file and create x,y,z.
    """

    # Check if file exists
    if not os.path.exists(filename):
        raise FileNotFoundError(f"File {filename} not found")

    # Load  data

    hh, ne_c, tipverige, zacpog,  t, izhod_e_a, _, _, _= load_for_plot(filename)

    # Call function to get indices
    X1, X2 = findCoordIndex(tipverige, hh, ne_c, 12)

    # here I use trick tro extract all coordinates along line
    a1 = np.array(X1) # this is mid span coordinates
    a2 = np.array(X2) # this is cleats

    # Reshape to exactly 3 rows (columns are automatic)
    m1 = a1.reshape(3, -1)  # shape (3, len(X1)//3)
    m2 = a2.reshape(3, -1)  # shape (3, len(X2)//3)

    # Join HORIZONTALLY (side by side) -> keeps 3 rows
    combined = np.hstack((m1, m2))  # shape (3, cols1 + cols2)

    # Sort each row independently in ascending order
    sorted_rows = np.sort(combined, axis=1)

    # Subtract 6 from the last column (column index -1)
    sorted_rows[:, -1] -= 6  # this is right-hand last coordinate


    # Flatten back to a list (row‑major order)
    X1=sorted_rows.flatten()#.tolist()

    X1 = X1.astype(int)

    # Extract  coordinates in one go (broadcasting)
    X = izhod_e_a[:, X1] #- zacpog[X1]  # shape (len(t), n_points)
    Y = izhod_e_a[:, X1 + 1] #- zacpog[X1 + 1]
    Z = izhod_e_a[:, X1 + 2] #- zacpog[X1 + 2]

    return t, X, Y, Z