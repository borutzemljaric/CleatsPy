import numpy as np
import matplotlib.pyplot as plt

def g_xy_at_given_t(data):
    """
    Animate the cable displacement profiles (x vs y and z vs y) over time.
    """
    t = data['t']
    X_1 = data['allX_1']
    Y_1 = data['allY_1']
    Z_1 = data['allZ_1']

    # Split into 3 cables (assumed concatenated horizontally)
    parts_X = np.array_split(X_1, 3, axis=1)
    parts_Y = np.array_split(Y_1, 3, axis=1)
    parts_Z = np.array_split(Z_1, 3, axis=1)

    # Relative displacement (subtract initial time step)
    x1 = parts_X[0] - parts_X[0][0, :]
    x2 = parts_X[1] - parts_X[1][0, :]
    x3 = parts_X[2] - parts_X[2][0, :]

    z1 = parts_Z[0] - parts_Z[0][0, :]
    z2 = parts_Z[1] - parts_Z[1][0, :]
    z3 = parts_Z[2] - parts_Z[2][0, :]

    # Spatial y-coordinate for each cable (constant over time)
    y_spatial_0 = parts_Y[0][0, :]   # shape (points_per_cable,)
    y_spatial_1 = parts_Y[1][0, :]
    y_spatial_2 = parts_Y[2][0, :]
    # If the Y coordinates are identical for all cables, you can use any one.
    # For safety, we use a list.
    y_spatial_list = [y_spatial_0, y_spatial_1, y_spatial_2]

    # Group data with their corresponding y_spatial, style, and label
    # Each tuple: (y_spatial, displacement_array_2D, style, label)
    x_curves = [
        (y_spatial_0, x1, 'r-', 'Cable 1'),
        (y_spatial_1, x2, 'g--', 'Cable 2'),
        (y_spatial_2, x3, 'b:', 'Cable 3')
    ]
    z_curves = [
        (y_spatial_0, z1, 'r-', 'Cable 1'),
        (y_spatial_1, z2, 'g--', 'Cable 2'),
        (y_spatial_2, z3, 'b:', 'Cable 3')
    ]

    # Setup interactive plot
    plt.ion()
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10), sharex=True)

    # Pre-configure axes (static parts)
    ax1.set_ylabel('Displacement x (m)')
    ax1.set_title('X displacement vs Y coordinate')
    ax1.grid(True, linestyle='--', alpha=0.6)

    ax2.set_ylabel('Displacement z (m)')
    ax2.set_title('Z displacement vs Y coordinate')
    ax2.set_xlabel('Y coordinate along cable (m)')
    ax2.grid(True, linestyle='--', alpha=0.6)

    # Set x‑ticks based on the first cable's Y (they are the same)
    y_first = y_spatial_0
    num_ticks = 8
    tick_positions = np.linspace(y_first.min(), y_first.max(), num_ticks)
    for ax in (ax1, ax2):
        ax.set_xticks(tick_positions)
        ax.set_xticklabels([f'{val:.2f}' for val in tick_positions])

    # Animation loop
    for tt in range(len(t)):
        # Clear previous lines
        ax1.clear()
        ax2.clear()

        # Re‑apply static configurations
        ax1.set_ylabel('Displacement x (m)')
        ax1.set_title('X displacement vs Y coordinate')
        ax1.grid(True, linestyle='--', alpha=0.6)
        ax2.set_ylabel('Displacement z (m)')
        ax2.set_title('Z displacement vs Y coordinate')
        ax2.set_xlabel('Y coordinate along cable (m)')
        ax2.grid(True, linestyle='--', alpha=0.6)
        ax1.set_xticks(tick_positions)
        ax1.set_xticklabels([f'{val:.2f}' for val in tick_positions])
        ax2.set_xticks(tick_positions)
        ax2.set_xticklabels([f'{val:.2f}' for val in tick_positions])

        # Plot X profiles for each cable
        for y_spatial, data_arr, style, label in x_curves:
            ax1.plot(y_spatial, data_arr[tt, :], style, linewidth=2, label=label)

        # Plot Z profiles
        for y_spatial, data_arr, style, label in z_curves:
            ax2.plot(y_spatial, data_arr[tt, :], style, linewidth=2, label=label)

        # Add legends and current time
        ax1.legend(loc='best')
        ax2.legend(loc='best')
        fig.suptitle(f'Time = {t[tt]:.3f} s', fontsize=12)

        # Refresh GUI
        fig.canvas.draw()
        fig.canvas.flush_events()

        # Wait for keypress to advance
        plt.waitforbuttonpress()

    plt.ioff()
    plt.show()