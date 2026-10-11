import numpy as np
from pathlib import Path
from SaveAndLoad.SaveLoad import load_menuForDynamic
from postprocess.mergedProcessFiles import mergeFiles_and_backup
from prestatic.reactionNodes import sileMomenti
from postprocess.systemEquationPostProcess import systemEquationPostProcess
from SaveAndLoad.SaveLoad import save_dynamic_postprocess




def FD32_postprocessing(FileToOpen):
    """
    Post-processing function for calculating reactions with basic order
    """
    filename = Path(FileToOpen).stem  # only file name without extension

    # Combine interval coordinate files (due to memory limitations in dynamics) into one
    Z, t = mergeFiles_and_backup(filename, keep_individual_files=True)


    # Load rest of the data from static results
    FileToRead = f'ySR_{filename}.h5'  # Should be file name+extension, to read and write new results
    cab, _, _, hh, ne_c, tipverige, le, zacpog, Ce_active, Bdi, U, V, Cdi, Ced, M, Qe_g = load_menuForDynamic(FileToRead)
    tipverige = np.array([x.decode('utf-8') for x in tipverige])  # because h5py saving method,  to avoid error with prefix b letter

    # General constants
    g = 9.81  # acceleration m/s^2

    # Result in basic order
    #Qeg_0 = Qe_g  # Just to keep through program same variables
    e_a = zacpog.copy() # for first

    # Calculate ALL coordinates through constraint equations
    nt = Z.shape[0]  # number of time steps


    # Initialize output arrays
    izhod_e_a = np.zeros((nt, e_a.size))
    izhod_e_dot = np.zeros((nt, e_a.size))
    izhod_edotdot = np.zeros((nt, e_a.size))
    Qc_izhod_a = np.zeros((nt, Ce_active.shape[1]))
    #e_a_out=np.zeros_like(e_a)

    for ia, z_row in enumerate(Z): # Pass t[ia] and z_row;
        #e_a=e_a_out
        e_a, e_dot, edotdot, Qc = systemEquationPostProcess(t[ia], z_row, zacpog, e_a, ne_c, V, U, le,
                                 Qe_g, M, tipverige, hh, Cdi, Ce_active,cab)

        izhod_e_a[ia, :] = e_a  # assume 1‑D
        izhod_e_dot[ia, :] = e_dot
        izhod_edotdot[ia, :] = edotdot
        Qc_izhod_a[ia, :] = Qc

        #izhod_e_a[ia, :] = e_a.T
        #izhod_e_dot[ia, :] = e_dot.T
        #izhod_edotdot[ia, :] = edotdot.T
        #Qc_izhod_a[ia, :] = Qc_izhod.T

    ne=int(Qc_izhod_a.shape[1]/4)
    # Calculate forces and moments
    izhod_Flevo = np.zeros((nt, ne))
    izhod_Fdesno = np.zeros((nt, ne))


    for s, row in enumerate(Qc_izhod_a):
        Flevo, Fdesno = sileMomenti(row)
        izhod_Flevo[s, :len(Flevo)] = Flevo.ravel()
        izhod_Fdesno[s, :len(Fdesno)] = Fdesno.ravel()


    # Save results to file , add them to existing file

    save_dynamic_postprocess(FileToRead, t, izhod_e_a, izhod_e_dot, izhod_Flevo, izhod_Fdesno)

    print('END POST processing')



# ============================================================================
# USAGE
# ============================================================================
if __name__ == "__main__":

    # A. Calling external file with input data about power line field
    FileToOpen = 'TestFileCleavesX.txt'
    FileToOpen = ('primer_1.txt')

    FD32_postprocessing(FileToOpen)
    print(f"Completed Postprocessing")
