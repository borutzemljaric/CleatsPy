import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def g_cForce_t(point_idx,force):
    # ----------------------------
    # Plot reaction forces in cleave joints vs time for all three cable, for a given cleave
    # ----------------------------

    point_idx = point_idx - 1  # to move on with python indexing

    #  X1, Y1, Z1 are 2D arrays with same shape, for all three cables, so we split them
    parts_F = np.array_split(force, 3, axis=1)  # list of 3 arrays
    # each array has number of rows equal to time records, and number of columns equal to number of left side elements



    x1 = parts_F[0][:, point_idx:point_idx + 1]
    x2 = parts_F[1][:, point_idx:point_idx + 1]
    x3 = parts_F[2][:, point_idx:point_idx + 1]


    # Group them into lists for easy plotting
    x_curves = [(x1, 'r-', 'C1'), (x2, 'g--', 'C2'), (x3, 'b:', 'C3')]

    # ----------------------------
    # Plot 1: Time history (3 subplots, each with three curves)
    # ----------------------------
    #fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 10), sharex=True)
    fig, ax= plt.subplots( figsize=(8, 4), sharex=True)

    # Subplot 1: X coordinate
    for data, style, label in x_curves:
        ax.plot(t, data, style, linewidth=2, label=label)
    ax.set_ylabel('F (N)')
    ax.grid(True, linestyle='--', alpha=0.6)  # make grid more visible
    ax.set_title('Force vs time')
    ax.legend()
    ax.set_xlabel('Time (s)')


    # ---- Set custom time ticks and labels on ALL subplots ----
    span = t[-1] - t[0]
    num_ticks = 10
    tick_positions = np.linspace(t[0], t[-1], num_ticks)

    # Choose format
    if span < 1:
        fmt = '{:.3f} s'
    elif span < 10:
        fmt = '{:.2f} s'
    else:
        fmt = '{:.0f} s'

    labels = [fmt.format(x) for x in tick_positions]


    ax.set_xticks(tick_positions)
    ax.set_xticklabels(labels)
    ax.tick_params(labelbottom=True)



    plt.tight_layout()
    plt.show()

