import numpy as np

def i_currents(t, cab, Nel):
    """
    Calculate currents in three conductors over integration time
    """

    # Current functions (transient signal for each phase)
    it1 = lambda t, ksi, fi, w: np.sin(w * t + ksi - fi) - np.exp(-t * w / np.tan(fi)) * np.sin(ksi - fi)
    it2 = lambda t, ksi, fi, w: np.sin(w * t + np.pi * 240 / 180 + ksi - fi) - np.exp(-t * w / np.tan(fi)) * np.sin(
        np.pi * 240 / 180 + ksi - fi)
    it3 = lambda t, ksi, fi, w: np.sin(w * t + np.pi * 120 / 180 + ksi - fi) - np.exp(-t * w / np.tan(fi)) * np.sin(
        np.pi * 120 / 180 + ksi - fi)

    # Parameters
    f = 50
    w = 2 * np.pi * f
    # t-time
    # ksi- switching angle
    # fi- impedance angle
    # w- angular frequency
    # Angle ranges

    fi = cab.phi * np.pi / 180  # radians
    ksi=0 # switching angle is taken to be 0°
    # ksi = ksi_angles * np.pi / 180  # radians
    mask=(t < cab.tsim) | (t > cab.tsim + cab.trelay)
    if np.all(mask):
        val1 = 0
        val2 = 0
        val3 = 0
    else:
        val1=it1(t-cab.tsim, ksi, fi, w)  # t-cab.tsim is used that time for short circuit starts from 0
        val2=it2(t-cab.tsim, ksi, fi, w)
        val3=it3(t-cab.tsim, ksi, fi, w)

        # Determine segment number with same current, because there are three equal cables
    n = Nel // 3
    #n2 = Nel // 3
    #n3 = Nel - n1 - n2  # remainder goes to the last segment

    # Build the array
    it = np.concatenate([
        np.full(n, val1),
        np.full(n, val2),
        np.full(n, val3)
    ])*np.sqrt(2)*cab.Ik


    return it  # array where each element has corresponding electric current

