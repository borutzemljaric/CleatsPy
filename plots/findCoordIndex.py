import numpy as np

def findCoordIndex(tipverige: np.ndarray, hh: int, ne_c: np.ndarray, tip_K_F: int):
    """
    Function that determines indices for coordinates/forces in output files for plotting graphs.
indeksiKoordinat
    Parameters:
    -----------
    tipverige : np.ndarray
    hh : int
        Number of spans
    ne_c : np.ndarray
        Array of number of elements per span
    tip_K_F : int
        Type of data: 12 for coordinates, 3 for forces

    Returns:
    --------
    X1 : np.ndarray
        Indices for mid-span coordinates/forces
    insX1 : np.ndarray
        Indices for cleave coordinates/forces
    """

    # Find midpoint of each span of first conductor
    s = np.round(ne_c / 2).astype(int)  # Find midpoints of individual spans

    # Print coordinates of insulators/midpoints of spans for graph plotting
    # tip_K_F can be 12 (for coordinates) or 3 (for forces)
    qI = ne_c  * tip_K_F  # Number of  conductor coordinates, per spans

    # Total entries per conductor
    total_per_conductor = np.sum(qI)

    # predefine array
    X1 = np.array([])  #
    cleaveX1 = np.array([0])  # the first coordinate is starting cleave

    # Predefinitions of even matrices and indices
    previousSpan = 0

    # Count number of spans minus one (that's how many cleaves inside end points are)
    # Loop goes from i=1 to (hh-1)-1 (so hh-2 iterations)
    for i in range(0, hh-1 ):  # hh in number of cleaves
        X1 = np.append(X1, s[i] * tip_K_F + previousSpan)  # Valid for midpoint of span
        previousSpan += qI[i]
        cleaveX1 = np.append(cleaveX1,  previousSpan)  # Valid for midpoint of span

    # Append the final cleave (end of the last span)

    # extending coordinates on all three conductors
    offset2 = total_per_conductor
    offset3 = 2 * total_per_conductor

    X1_all = np.concatenate([X1, X1 + offset2, X1 + offset3])
    cleaveX1_all = np.concatenate([cleaveX1, cleaveX1 + offset2, cleaveX1 + offset3])

    return X1_all, cleaveX1_all