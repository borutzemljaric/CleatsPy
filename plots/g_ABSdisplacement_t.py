import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def g_ABSdisplacement_t(point_idx,data,relative, drawFigure):
    # ----------------------------
    # Plot xyz vs time for all three cable, for a given point
    # ----------------------------
    point_idx=point_idx-1  # to move on with python indexing

    t = data['t']
    cab_d= data['cab_d']*1e-3 # conversion to meters
    X_1 = data.get('X_1')
    Y_1 = data.get('Y_1')
    Z_1 = data.get('Z_1')

    #  X1, Y1, Z1 are 2D arrays with same shape, for alll three cables, so we split them
    parts_X = np.array_split(X_1, 3, axis=1)  # list of 3 arrays
    parts_Y = np.array_split(Y_1, 3, axis=1)
    parts_Z = np.array_split(Z_1, 3, axis=1)

    if relative==True:
        x1 = parts_X[0][:, point_idx:point_idx + 1]-parts_X[0][0, point_idx]
        x2 = parts_X[1][:, point_idx:point_idx + 1]-parts_X[1][0, point_idx]
        x3 = parts_X[2][:, point_idx:point_idx + 1]-parts_X[2][0, point_idx]

        # Y components
        y1 = parts_Y[0][:, point_idx:point_idx + 1]-parts_Y[0][0, point_idx]
        y2 = parts_Y[1][:, point_idx:point_idx + 1]-parts_Y[1][0, point_idx]
        y3 = parts_Y[2][:, point_idx:point_idx + 1]-parts_Y[2][0, point_idx]

        # Z components
        z1 = parts_Z[0][:, point_idx:point_idx + 1]-parts_Z[0][0, point_idx]
        z2 = parts_Z[1][:, point_idx:point_idx + 1]-parts_Z[1][0, point_idx]
        z3 = parts_Z[2][:, point_idx:point_idx + 1]-parts_Z[2][0, point_idx]
    else:
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

    r1 = np.sqrt(np.array(x1)**2 + np.array(y1)**2 + np.array(z1)**2) / cab_d # absolute distance
    r2 = np.sqrt(np.array(x2)**2 + np.array(y2)**2 + np.array(z2)**2) / cab_d
    r3 = np.sqrt(np.array(x3)**2 + np.array(y3)**2 + np.array(z3)**2) / cab_d
    # Group them into lists for easy plotting
    x_curves = [(r1, 'r-', 'C1'), (r2, 'g--', 'C2'), (r3, 'b:', 'C3')]

    #read maximum value from the curve, the list is return from this function
    max_values = [r1.max(), r2.max(), r3.max()]
    print(max_values)

    if drawFigure==True:

        # ----------------------------
        # Plot 1: Time history (3 subplots, each with three curves)
        # ----------------------------
        fig, ax = plt.subplots(figsize=(8, 6))

        # Subplot 1: X coordinate
        for data, style, label in x_curves:
            ax.plot(t, data, style, linewidth=2, label=label)
        ax.set_ylabel('Dimensionless displacement $h*=h/d_c$')
        ax.grid(True, linestyle='--', alpha=0.6)  # make grid more visible
        #ax.set_title('X coordinate vs time')
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

        plt.tight_layout()
        plt.show()

    return max_values



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