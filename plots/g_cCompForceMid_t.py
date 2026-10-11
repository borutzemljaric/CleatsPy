import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def g_cCompForceMid_t(point_idx,data,relative=False):
    # ----------------------------
    # Plot component forces in cleave vs time for all three cable
    # ----------------------------
    point_idx=point_idx-1  # to move on with python indexing
    t = data['t']
    # All keys in one list (corrected the typo)
    keys = ['FX_mid_1', 'FY_mid_1', 'FZ_mid_1']

    # Read and split all at once into a dictionary
    split_parts = {k: np.array_split(data.get(k), 3, axis=1) for k in keys}

    # Unpack to your original variable names (optional, keeps downstream code unchanged)
    parts_X, parts_Y, parts_Z = [split_parts[k] for k in keys]
    #parts_X or X or Z have the  depth 3, each layer present conductor 1,2, 3

    # resultant forces in node (summarize the end of two elements)


    # Now fx_node contains the sum -resultant force in cleaves

    # Y components
    x1 = parts_X[0][:, point_idx:point_idx + 1]
    x2 = parts_X[1][:, point_idx:point_idx + 1]
    x3 = parts_X[2][:, point_idx:point_idx + 1]

    # Y components
    y1 = parts_Y[0][:, point_idx:point_idx + 1]
    y2 = parts_Y[1][:, point_idx:point_idx + 1]
    y3 = parts_Y[2][:, point_idx:point_idx + 1]

    # Z components
    z1 = parts_Z[0][:, point_idx:point_idx + 1]
    z2 = parts_Z[1][:, point_idx:point_idx + 1]
    z3 = parts_Z[2][:, point_idx:point_idx + 1]


    # Group them into lists for easy plotting
    x_curves = [(x1, 'r-', 'C1'), (x2, 'g--', 'C2'), (x3, 'b:', 'C3')]
    y_curves = [(y1, 'r-', 'C1'), (y2, 'g--', 'C2'), (y3, 'b:', 'C3')]
    z_curves = [(z1, 'r-', 'C1'), (z2, 'g--', 'C2'), (z3, 'b:', 'C3')]

    # read maximum value from the curve, the list is return from this function
    yielfForce=120000 #N from 120 N/mm2 *1000 nn2
    max_values = [y1.max()/yielfForce, y2.max()/yielfForce, y3.max()/yielfForce]
    print('max_values of horizontal Y component, normalised by Aluminium yield strength')
    print('C1,C2,C3', max_values)
    # ----------------------------
    # Plot 1: Time history (3 subplots, each with three curves)
    # ----------------------------
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 10), sharex=True)

    # Subplot 1: X coordinate
    for data, style, label in x_curves:
        ax1.plot(t, data, style, linewidth=2, label=label)
    ax1.set_ylabel(f'$F_x$ (N)')
    ax1.grid(True, linestyle='--', alpha=0.6)  # make grid more visible
    ax1.set_title('Force midspan X vs time')
    ax1.legend()
    ax1.set_xlabel('Time (s)')

    # Subplot 2: Y coordinate
    for data, style, label in y_curves:
        ax2.plot(t, data, style, linewidth=2, label=label)
    ax2.set_ylabel(f'$F_y$ (N)')
    ax2.grid(True, linestyle='--', alpha=0.6)
    ax2.set_title('Force midspan Y vs time')
    ax2.legend()
    #ax2.set_ylim(-0.01, 0.01)

    # Subplot 3: Z coordinate
    for data, style, label in z_curves:
        ax3.plot(t, data, style, linewidth=2, label=label)
    ax3.set_ylabel(f'$F_z$ (N)')
    ax3.grid(True, linestyle='--', alpha=0.6)
    ax3.set_title('Force midspan Z  vs time')
    ax3.legend()

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

    for ax in (ax1, ax2, ax3):
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(labels)
        ax.tick_params(labelbottom=True)

    # Single x-axis label on the bottom subplot
    ax3.set_xlabel('Time (s)')

    plt.tight_layout()
    plt.show()



    # ----------------------------
    # Plot 2: 3D trajectory (optional)
    # Here we plot only the first set (Case A) as an example.
    # ----------------------------
    # fig = plt.figure(figsize=(8, 8))
    # ax = fig.add_subplot(111, projection='3d')
    # ax.plot(x1, y1, z1, 'k-', linewidth=2, label='Trajectory (Case A)')
    # ax.scatter(x1[0], y1[0], z1[0], color='green', s=100, label='Start')
    # ax.scatter(x1[-1], y1[-1], z1[-1], color='red', s=100, label='End')
    # ax.set_xlabel('X (m)')
    # ax.set_ylabel('Y (m)')
    # ax.set_zlabel('Z (m)')
    # ax.set_title('3D trajectory of the point')
    # ax.legend()
    # plt.show()