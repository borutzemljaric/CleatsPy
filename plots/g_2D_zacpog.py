import numpy as np
import matplotlib.pyplot as plt
"""
2D plot, x,y or Z  along two axes
-particular for draw initial conditions
"""

def g_2D_zacpog(zacpog, e_zp):
    """Plot results"""
    # Example: assume zacpog is a flat list/array
    # zacpog = ...  (your data)

    # Convert to numpy array for easier slicing
    data = np.array(zacpog)

    # Total number of points per block
    n = len(data) // (3 * 12)  # because 3 blocks × n points × 12 values

    # Reshape into (3, n, 12)  → 3 blocks, each with n points of 12 coordinates
    blocks = data.reshape(3, n, 12)

    # Plot y vs z for each block
    plt.figure(figsize=(8, 6))

    colors = ['blue', 'green', 'red']
    labels = ['Block 1', 'Block 2', 'Block 3']
    markers = ['o', 'x', '*']
    markersizes=[3,5,3]

    for i in range(3):
        y = blocks[i, :, 1]  # index 1 = y coordinate
        z = blocks[i, :, 2]  # index 2 = z coordinate
        plt.plot(y, z, color=colors[i], label=labels[i], marker=markers[i], markersize=markersizes[i])

    plt.xlabel('y')
    plt.ylabel('z')
    # Add 5% padding to y‑axis (and optionally x‑axis)
    plt.margins(y=0.05, x=0.05)  # 5% on both axes
    plt.title('y‑z plot of the three blocks')
    plt.legend()
    plt.grid(True)
    plt.axis('equal')  # optional: keep aspect ratio equal
    plt.show()