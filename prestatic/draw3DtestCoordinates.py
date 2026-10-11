import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def draw3DtestCoordinates(e_zp):
    # Assuming e_zp is your 1D numpy array
    # e_zp = np.array(...)   # shape (n*12,)

    # Reshape into blocks of 12
    n = len(e_zp) // 12

    blocks = e_zp.reshape(n, 12)


    # Extract start coordinates (first three)
    start = blocks[:, 0:3]          # shape (n, 3)

    # Extract end coordinates (last three – indices 9,10,11)
    # If yours are at indices 7,8,9, change to: end = blocks[:, 7:10]
    end = blocks[:, 6:9]           # shape (n, 3)

    # Create 3D plot
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Plot each line segment
    for i in range(n):
        ax.plot([start[i, 0], end[i, 0]],
                [start[i, 1], end[i, 1]],
                [start[i, 2], end[i, 2]],
                color='b', alpha=0.7, linewidth=1)

    # Optional: mark the start points (optional)
    ax.scatter(start[:, 0], start[:, 1], start[:, 2], c='r', s=10, label='Start')

    # Labels and title
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('3D Lines from e_zp')
    ax.legend()

    plt.tight_layout()
    plt.show()