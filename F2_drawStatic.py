import numpy as np
import matplotlib.pyplot as plt
from SaveAndLoad.SaveLoad import load_StaticResults
from plots.g_2D_zacpog import g_2D_zacpog


def F2_drawStatic(ime_datoteke2):
    """
    Plots static  analysis results and calculates forces
    """

    # Extract variables
    cab, e_zp, e_zp_s, hh, ne_c, tipverige, le,zacpog=load_StaticResults(File_ToRead)

    #g_2D_zacpog(zacpog,e_zp)  # testing static results
    #breakpoint()

    # Indices for conductors and insulators
    ec = []      # indices for conductors
    ci = []      # indices for cleaves

    Fec = []     # force indices for conductors
    Fci = []     # force indices for insulators

    
    ii = 1  # auxiliary index
    # Convert ne_c to a 1D array for easier indexing
    # ne_c is in this shape [[3, 3, 0]] -> we want [3, 3, 0]
    ne_c_flat = ne_c.flatten()

    for k in range(hh - 1):
        # Get number of elements in this segment
        num_elements_in_segment = int(ne_c_flat[k])
        for i in range(num_elements_in_segment):
            # Conductor indices (12 DOFs per element)
            ec.extend(range(((ii-1)*12 + 0) , ((ii-1)*12 + 12)))  #
            Fec.extend(range(((ii-1)*3 + 0) , ((ii-1)*3 + 3)))    #

            if k != hh - 2 and i == num_elements_in_segment - 1:  # conductor-cleave element
                # Additional indices for insulators
                ei1 = list(range(((ii-1)*12 + 1 + 11) , ((ii-1)*12 + 12 + 12) ))
                Fei1 = list(range(((ii-1)*3 + 1 + 2) , ((ii-1)*3 + 3 + 3) ))
                ii += 1  # count conductor elements


                ci.extend(ei1)
                Fci.extend(Fei1)

            else:
                ii += 1  # count conductor elements
    
    # Convert to numpy arrays
    ec = np.array(ec) #index verižnice
    ci = np.array(ci) #%index izolacije

    Fec = np.array(Fec)
    Fci = np.array(Fci)

    # print('ec',ec)
    # print('ci',ci)
    # print('Fci', Fci)
    # breakpoint()
    # Get coordinates for plotting
    # Conductor points (every 12th element starting at 0, 1, 2 for x, y, z)
    tocke_x = ec[0::12]
    tocke_y = ec[1::12]
    tocke_z = ec[2::12]
    
    # Add the last points
    tocke_x = np.append(tocke_x, ec[-6])
    tocke_y = np.append(tocke_y, ec[-5])
    tocke_z = np.append(tocke_z, ec[-4])

    tocke_xi = ci[0::6]
    tocke_yi = ci[1::6]
    tocke_zi = ci[2::6]

    # print(tocke_zi)
    # breakpoint()

    n=int(len(e_zp)/3) # we have 3 cables, so every cable start at :0*n, 1*n, 2*n

    tocke_x_2d = np.zeros((3, len(tocke_x)))
    tocke_y_2d = np.zeros((3, len(tocke_x)))
    tocke_z_2d = np.zeros((3, len(tocke_x)))

    for i in range(3):
        # Compute the row
        tocke_x_2d[i, :] = tocke_x+i*n  # assign to the i‑th row
        tocke_y_2d[i, :] = tocke_y + i * n  # assign to the i‑th row
        tocke_z_2d[i, :] = tocke_z + i * n  # assign to the i‑th row


    #print(np.column_stack((zacpog[tocke_y_2d[0, :].astype(int) ], e_zp[tocke_y_2d[0, :].astype(int) ])))
    #breakpoint()

    # print(f"zacpog shape: {zacpog.shape}")
    # Plotting
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')

    ax.plot3D(zacpog[tocke_x_2d[0, :].astype(int) ], zacpog[tocke_y_2d[0, :].astype(int) ], zacpog[tocke_z_2d[0, :].astype(int) ],
              '-.r*', label='Zero state', markersize=6, linewidth=2)
    ax.plot3D(zacpog[tocke_x_2d[1, :].astype(int) ], zacpog[tocke_y_2d[1, :].astype(int) ], zacpog[tocke_z_2d[1, :].astype(int) ],
              '-.g*', label='Zero state', markersize=6, linewidth=2)
    ax.plot3D(zacpog[tocke_x_2d[2, :].astype(int) ], zacpog[tocke_y_2d[2, :].astype(int) ], zacpog[tocke_z_2d[2, :].astype(int) ],
              '-.b*', label='Zero state', markersize=6, linewidth=2)

    ax.scatter3D(
        zacpog[tocke_x_2d[0, :].astype(int)],
        zacpog[tocke_y_2d[0, :].astype(int)],
        zacpog[tocke_z_2d[0, :].astype(int)],
        marker='o',  # circular markers
        s=20,  # size of markers (points)
        color='r',  # or any color
        alpha=0.7  # transparency (optional)
    )


    # Configure plot
    ax.set_title('Cable Analysis', fontsize=14)
    ax.set_xlabel('x (m)', fontsize=12)
    ax.set_ylabel('y (m)', fontsize=12)
    ax.set_zlabel('z (m)', fontsize=12)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    # Save figure
    save_figure=False
    if save_figure is True:
        output_filename = f"{ime_datoteke2}_Razpetina.png"
        plt.savefig(output_filename, dpi=150)
        print(f"\nPlot saved as: {output_filename}")
    else:
        pass

    plt.show()

    # Force calculations
    print("\n" + "=" * 60)
    print("FORCE CALCULATIONS")
    print("=" * 60)

    if len(Fec) > 0:
        # Force indices (every 3rd element starting at 0, 1, 2 for x, y, z forces)
        indF_x = Fec[0::3]
        indF_y = Fec[1::3]
        indF_z = Fec[2::3]

        print(f"\nForce index counts:")
        print(f"  indF_x: {len(indF_x)} indices")
        print(f"  indF_y: {len(indF_y)} indices")
        print(f"  indF_z: {len(indF_z)} indices")

        # Check if Flevo is 1D or 2D
        if len(Flevo.shape) == 2 and Flevo.shape[1] > 1:
            # Use first column as in MATLAB code
            force_column = 0
            print(f"Using column {force_column} of Flevo (shape: {Flevo.shape})")
        else:
            force_column = 0

        # Left side conductor forces
        print("\nLeft side conductor forces:")
        Fprecno = Flevo[indF_x, force_column] if len(Flevo.shape) == 2 else Flevo[indF_x]
        Fvodoravno = Flevo[indF_y, force_column] if len(Flevo.shape) == 2 else Flevo[indF_y]
        Fnavpicno = Flevo[indF_z, force_column] if len(Flevo.shape) == 2 else Flevo[indF_z]
        Fskupna = np.sqrt(Fprecno ** 2 + Fvodoravno ** 2 + Fnavpicno ** 2)

        print(f"Transverse (x) forces: {Fprecno.flatten()}")
        print(f"Horizontal (y) forces: {Fvodoravno.flatten()}")
        print(f"Vertical (z) forces: {Fnavpicno.flatten()}")
        print(f"Total forces: {Fskupna.flatten()}")

        # Right side conductor forces
        print("\nRight side conductor forces:")
        Fprecno = Fdesno[indF_x, force_column] if len(Fdesno.shape) == 2 else Fdesno[indF_x]
        Fvodoravno = Fdesno[indF_y, force_column] if len(Fdesno.shape) == 2 else Fdesno[indF_y]
        Fnavpicno = Fdesno[indF_z, force_column] if len(Fdesno.shape) == 2 else Fdesno[indF_z]
        Fskupna = np.sqrt(Fprecno ** 2 + Fvodoravno ** 2 + Fnavpicno ** 2)

        print(f"Transverse (x) forces: {Fprecno.flatten()}")
        print(f"Horizontal (y) forces: {Fvodoravno.flatten()}")
        print(f"Vertical (z) forces: {Fnavpicno.flatten()}")
        print(f"Total forces: {Fskupna.flatten()}")
    else:
        print("\nNo conductor force indices found.")

    if len(Fci) > 0:
        # Insulator force indices
        indF_xi = Fci[0::3]
        indF_yi = Fci[1::3]
        indF_zi = Fci[2::3]

        print(f"\nInsulator force index counts:")
        print(f"  indF_xi: {len(indF_xi)} indices")
        print(f"  indF_yi: {len(indF_yi)} indices")
        print(f"  indF_zi: {len(indF_zi)} indices")

        # Left side insulator forces
        print("\nLeft side insulator forces:")
        Fprecno = Flevo[indF_xi, force_column] if len(Flevo.shape) == 2 else Flevo[indF_xi]
        Fvodoravno = Flevo[indF_yi, force_column] if len(Flevo.shape) == 2 else Flevo[indF_yi]
        Fnavpicno = Flevo[indF_zi, force_column] if len(Flevo.shape) == 2 else Flevo[indF_zi]
        Fskupna = np.sqrt(Fprecno ** 2 + Fvodoravno ** 2 + Fnavpicno ** 2)

        print(f"Transverse (x) forces: {Fprecno.flatten()}")
        print(f"Horizontal (y) forces: {Fvodoravno.flatten()}")
        print(f"Vertical (z) forces: {Fnavpicno.flatten()}")
        print(f"Total forces: {Fskupna.flatten()}")

        # Right side insulator forces
        print("\nRight side insulator forces:")
        Fprecno = Fdesno[indF_xi, force_column] if len(Fdesno.shape) == 2 else Fdesno[indF_xi]
        Fvodoravno = Fdesno[indF_yi, force_column] if len(Fdesno.shape) == 2 else Fdesno[indF_yi]
        Fnavpicno = Fdesno[indF_zi, force_column] if len(Fdesno.shape) == 2 else Fdesno[indF_zi]
        Fskupna = np.sqrt(Fprecno ** 2 + Fvodoravno ** 2 + Fnavpicno ** 2)

        print(f"Transverse (x) forces (N): {Fprecno.flatten()}")
        print(f"Horizontal (y) forces (N): {Fvodoravno.flatten()}")
        print(f"Vertical (z) forces (N) {Fnavpicno.flatten()}")
        #print(f"Total forces: {Fskupna.flatten():.4f} ")
        print(f"Total forces (N): {[f'{x:.2f}' for x in Fskupna.flatten()]}")
    else:
        print("\nNo insulator force indices found.")

    #

    print("\n" + "=" * 60)
    print("END OF ANALYSIS")
    print("=" * 60)

    # Return results for further processing if needed
    results = {
        'ec': ec,
        'ci': ci,
        'tocke_x': tocke_x,
        'tocke_y': tocke_y,
        'tocke_z': tocke_z,
        'tocke_xi': tocke_xi,
        'tocke_yi': tocke_yi,
        'tocke_zi': tocke_zi
    }

    return results




# Add main execution block
if __name__ == "__main__":

    # without file extension '.txt' !!
    File_ToRead = f'TestFileCleavesX'  # Change this to your file
    #File_ToRead = 'primer_1a'

    File_ToRead=f'ySR_{File_ToRead}.h5'
    F2_drawStatic(File_ToRead)
