import numpy as np
import os

from SaveAndLoad.SaveLoad import load_for_plot
from plots.findCoordIndex import findCoordIndex


def readCoordinate(filename):
    """
    Read displacement data from file.
    """

    # Check if file exists
    if not os.path.exists(filename):
        raise FileNotFoundError(f"File {filename} not found")

    # Load  data

    hh, ne_c, tipverige, zacpog,  t, izhod_e_a, _, _, _= load_for_plot(filename)


    # Call function to get indices
    X1, _= findCoordIndex(tipverige, hh, ne_c, 12)
    X1 = X1.astype(int)
    #cleaveX1= cleaveX1.astype(int)

    # Extract insulator coordinates
    # cleaveX = izhod_e_a[:, cleaveX1]
    # cleaveY = izhod_e_a[:, cleaveX1 + 1]
    # cleaveZ = izhod_e_a[:, cleaveX1 + 2]

    # Calculate relative displacements!!!

    # Extract relative displacements in one go (broadcasting)
    X = izhod_e_a[:, X1] #- zacpog[X1]  # shape (len(t), n_points)
    Y = izhod_e_a[:, X1 + 1] #- zacpog[X1 + 1]
    Z = izhod_e_a[:, X1 + 2] #- zacpog[X1 + 2]

    # If you also need cleave coordinates (insulators), do the same:
    # cleaveX = izhod_e_a[:, cleaveX1] - zacpog[cleaveX1]
    # cleaveY = izhod_e_a[:, cleaveX1 + 1] - zacpog[cleaveX1 + 1]
    # cleaveZ = izhod_e_a[:, cleaveX1 + 2] - zacpog[cleaveX1 + 2]

    return t, X, Y, Z