import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib as mpl
import os


def FD33_plotDynamic(FileToOpen):
    """
    Function for plotting dynamic simulation results.

    """




    # Set up plotting style
    plt.style.use('seaborn-v0_8-darkgrid')
    mpl.rcParams['figure.figsize'] = [10, 7]
    mpl.rcParams['font.size'] = 11

    # User parameters (from MATLAB code)
    k = 1  # span index to plot  (ničla pythona je upoštevana pri grafu)
    n_krivulj = 1  # number of curves to plot
    tip_grafa = 2  # type of graph
    smerGrafa = 3  # direction of graph
    # 1- z displacement
    # 2 x displacement
    # 3 zx displacement (vertical and horizontal)
    # 4 insulator displacement Y
    # 5 insulator displacement x
    # 6 time display (pressing adds graph one second at a time)
    # 7 time display for first ten seconds (simultaneous plot)
    # 8 forces at supports Left Right
    # 9 pole tip displacement
    # 10 forces in insulators
    # 11 insulator displacement in Y direction / for all insulators in span
    # choice between one or two curves


    # File names (hardcoded as in MATLAB)
    ime_datoteke_1 = 'I_40sek'
    ime_datoteke_2 = 'V_40sek'
    ime_datoteke_3 = 'BPA_M2_2_40sek'
    ime_datoteke_4 = 'A_40sek'
    ime_datoteke_5 = 'IB_40sek'

    file_names = ['yDR_fekr_post.h5', 'ime_datoteke_2', 'ime_datoteke_3',
                  'ime_datoteke_4', 'ime_datoteke_5']

    for i in range(n_krivulj):
        try:
            t_curr, X_i, Y_i, Z_i = readCoordinate(file_names[i])


    filename = Path(FileToOpen).stem  # only file name without extension
    FileToRead = f'ySR_{filename}.h5'  # Should be file name+extension, to read coordinate and forces

    # Usage
    data = load_curve_data(n_krivulj, file_names)
    t = data['t']
    X_1 = data.get('X_1')
    Y_1 = data.get('Y_1')
    Z_1 = data.get('Z_1')

    #     t, X_1, Y_1, Z_1, from first file , note that are coordinates for all three conductors
    #     t, X_2, Y_2,
    #     etc
    #     t, Flevo1, Fdesno1
    #     t, Flevo2, Fdesno2,
    #     etc.



    # Set graph limits
    mejagrafa = np.max(t)
    #print(mejagrafa)
    #breakpoint()
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 7))

    # Switch case for different plot types
    if tip_grafa == 1:
        # Z-displacement plot
        print('test')
        t = data['t']
        Z_1 = data.get('Z_1')
        ax.plot(t, Z_1[:, k - 1], '-k', linewidth=1.5, label='sredina 2 razpetine')

        if n_krivulj == 2:
            ax.plot(t, Z_2[:, k - 1], '-m', linewidth=1.5, label='sredina 2 razpetine-M2')
        elif n_krivulj == 5:
            ax.plot(t, Z_2[:, k - 1], 'xr', linewidth=1.5, label='I')
            ax.plot(t, Z_3[:, k - 1], '-.b', linewidth=1.5, label='V')
            ax.plot(t, Z_4[:, k - 1], '--m', linewidth=1.5, label='BPA')
            ax.plot(t, Z_5[:, k - 1], ':g', linewidth=1.5, label='IB')

        ax.set_ylabel('Displacement Z (m)', fontsize=12)
        ax.set_xlabel('Time (s)', fontsize=12)
        ax.set_xlim([0, 20])
        ax.set_xticks(np.arange(0, mejagrafa + 2, 2))
        ax.grid(True, which='both', linestyle='--', alpha=0.7)
        ax.legend(loc='best')

    elif tip_grafa == 2:
        # X-displacement plot
        t = data['t']
        X_1 = data.get('X_1')
        ax.plot(t, X_1[:, k - 1], '-k', linewidth=1.5, label='I-insu.')

        if n_krivulj == 2:
            ax.plot(t, X_2[:, k - 1], '--m', linewidth=1.5, label='I-insu.')
        elif n_krivulj == 5:
            ax.plot(t, X_2[:, k - 1], 'xr', linewidth=1.5, label='I')
            ax.plot(t, X_3[:, k - 1], '-.b', linewidth=1.5, label='V')
            ax.plot(t, X_4[:, k - 1], '--m', linewidth=1.5, label='BPA')
            ax.plot(t, X_5[:, k - 1], ':g', linewidth=1.5, label='IB')

        ax.set_ylabel('Pomik X (m)', fontsize=12)
        ax.set_xlabel('Čas (s)', fontsize=12)
        ax.set_xlim([0, 20])
        ax.set_xticks(np.arange(0, mejagrafa + 2, 2))
        ax.grid(True, which='both', linestyle='--', alpha=0.7)
        ax.legend(loc='best')

    elif tip_grafa == 3:
        # X-Z displacement plot with color coding
        # For color-coded line based on time
        points = np.array([X_1[:, k - 1], Z_1[:, k - 1]]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)

        # Create a color map based on time
        norm = plt.Normalize(t.min(), t.max())
        lc = mpl.collections.LineCollection(segments, cmap='viridis', norm=norm)
        lc.set_array(t)
        lc.set_linewidth(2)
        line = ax.add_collection(lc)

        ax.set_xlabel('Pomik X (m)', fontsize=12)
        ax.set_ylabel('Pomik Z (m)', fontsize=12)
        ax.grid(True, which='both', linestyle='--', alpha=0.7)

        # Add colorbar
        cbar = fig.colorbar(line, ax=ax)
        cbar.set_label('Čas (s)', fontsize=10)

    elif tip_grafa == 4:
        # Insulator displacement in Z direction
        ax.plot(t, insZ_1[:, k - 1], '-b', linewidth=1.5, label='I-insu.')

        if n_krivulj == 2:
            ax.plot(t, insZ_2[:, k - 1], '-m', linewidth=1.5, label='V-insu.')

        ax.set_ylabel('Pomik Z (m)', fontsize=12)
        ax.set_xlabel('Čas (s)', fontsize=12)
        ax.grid(True, which='both', linestyle='--', alpha=0.7)
        ax.legend(loc='best')

    elif tip_grafa == 5:
        # Insulator displacement in X direction
        ax.plot(t, insX_1[:, k - 1], '-b', linewidth=1.5, label='I-insu.')

        if n_krivulj == 2:
            ax.plot(t, insX_2[:, k - 1], '--.m', label='V-insu.')

        ax.set_ylabel('Pomik X (m)', fontsize=10)
        ax.set_xlabel('Čas (s)', fontsize=10)
        ax.set_xticks(np.arange(0, mejagrafa + 1, 1))
        ax.grid(True, which='both', linestyle='--', alpha=0.7)
        ax.legend(loc='best')

    elif tip_grafa == 8:
        # Forces at conductor supports
        ax.plot(t, Fdesno1, '-b', linewidth=1.5, label='Desni steber-TT2')
        ax.plot(t, Flevo1, '--om', markersize=4, linewidth=1, label='Levi steber-TT1')

        ax.set_ylabel('Sila absolutna (N)', fontsize=12)
        ax.set_xlabel('Čas (s)', fontsize=12)
        ax.set_xticks(np.arange(0, mejagrafa + 1, 1))
        ax.grid(True, which='both', linestyle='--', alpha=0.7)
        ax.legend(loc='best')

    elif tip_grafa == 9:
        # Pole tip displacements
        ax.plot(t, spX1[:, k - 1], '-b', linewidth=1.5, label='smer X')
        ax.plot(t, spY1[:, k - 1], '-r', linewidth=1.5, label='smer Y')
        ax.plot(t, spZ1[:, k - 1], '-g', linewidth=1.5, label='smer Z')

        ax.set_ylabel('Pomik (m)', fontsize=12)
        ax.set_xlabel('Čas (s)', fontsize=12)
        ax.set_xticks(np.arange(0, mejagrafa + 2, 2))
        ax.grid(True, which='both', linestyle='--', alpha=0.7)
        ax.legend(loc='best')

    elif tip_grafa == 10:
        # Forces at insulator supports
        ax.plot(t, Flevoins1[:, k - 1], '-b', linewidth=1.5, label='-P-M1')
        ax.plot(t, Fdesnoins1[:, k - 1], '-ob', linewidth=1.5, label='-B-M1')

        if n_krivulj == 2:
            ax.plot(t, Flevoins2[:, k - 1], '-m', linewidth=1.5, label='-P-M2')
            ax.plot(t, Fdesnoins2[:, k - 1], '-om', linewidth=1.5, label='-B-M2')

        ax.set_ylabel('Absolutna sila v obesišču na stebru (N)', fontsize=12)
        ax.set_xlabel('Čas (s)', fontsize=12)
        ax.set_xlim([0, 20])
        ax.set_xticks(np.arange(0, mejagrafa + 1, 1))
        ax.grid(True, which='both', linestyle='--', alpha=0.7)
        ax.legend(loc='best')

    elif tip_grafa == 11:
        # Insulator displacements for all insulators in span
        if n_krivulj == 1:
            ax.plot(t, insY_1, label=['insulator I1', 'insulator I2', 'insulator I3'])
        elif n_krivulj == 2:
            ax.plot(t, insY_1, label=['insulator I1-M1', 'insulator I2-M1', 'insulator I3-M1'])
            ax.plot(t, insY_2, label=['insulator I1-M2', 'insulator I2-M2', 'insulator I3-M2'])

        ax.set_ylabel('Displacement y (m)', fontsize=10)
        ax.set_xlabel('Time (s)', fontsize=10)
        ax.set_xticks(np.arange(0, mejagrafa + 1, 1))
        ax.grid(True, which='both', linestyle='--', alpha=0.7)
        ax.legend(loc='best')

    # Adjust layout and show plot
    plt.tight_layout()
    plt.show()

    return fig






# ============================================================================
# USAGE
# ============================================================================
if __name__ == "__main__":

    # Defining and Calling external files with  data to be presented
    # it is possible to define up to 5 files
    FileToOpen = 'TestFileCleavesX.txt'




    # File names
    file_names = ['TestFileCleavesX.txt', 'TestFileCleavesX.txt', None,
                  None, None]

    # User parameters
    k = 1  # span index to plot  (ničla pythona je upoštevana pri grafu)
    n_krivulj = 1  # number of curves to plot
    tip_grafa = 2  # type of graph
    smerGrafa = 3  # direction of graph
    # 1- xyz displacement
    # 2 x displacement
    # 3 zx displacement (vertical and horizontal)
    # 4 insulator displacement Y
    # 5 insulator displacement x
    # 6 time display (pressing adds graph one second at a time)
    # 7 time display for first ten seconds (simultaneous plot)
    # 8 forces at supports Left Right
    # 9 pole tip displacement
    # 10 forces in insulators
    # 11 insulator displacement in Y direction / for all insulators in span
    # choice between one or two curves

    FD33_plotDynamic(file_names)
