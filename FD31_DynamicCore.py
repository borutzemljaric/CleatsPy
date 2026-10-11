import numpy as np
from pathlib import Path
from SaveAndLoad.SaveLoad import load_menuForDynamic, update_cab_optional
from dynamic.d_open_MenuDynamic import d_open_menuDynamic
from dynamic.d1_mainLoop import d1_mainLoop
#from create_classes import cable_class


def FD31_DynamicCore(FileToOpen):
    """
    Preprocessing dynamic and processing
    """
    # A. Prepare Simulation data

    # Load data files
    # Extract filename without extension
    filename = Path(FileToOpen).stem
    FileToRead = f'ySR_{filename}.h5'  # Should be file name+extension, to read and write new results


    # Load data from static calculation
    cab, _, _, hh, ne_c, tipverige, le,zacpog, Ce_active, Bdi, U, V, Cdi, Ced, M, Qe_g= load_menuForDynamic(FileToRead)
    tipverige = np.array([x.decode('utf-8') for x in tipverige])  # because h5py saving method,  to avoid error with prefix b letter

    # Copy conductor data
    g = 9.81  # gravity acceleration

    # Define separately time and electric currents
    d_open_menuDynamic(cab)  # add time in current values
    update_cab_optional(FileToRead, cab) #update cab class with additional attributes


    # B. Go to computation
    d1_mainLoop(zacpog, U, V, hh, ne_c, tipverige, Cdi, Qe_g, Ce_active, M,
                cab, le,filename)



    print("Simulation completed!")
    return cab.ste_delitev  # Or return whatever you need


# ============================================================================
# USAGE
# ============================================================================
if __name__ == "__main__":

    # A. Calling external file with input data about power line field
    FileToOpen = 'TestFileCleavesX.txt'
    FileToOpen = 'primer_1.txt'

    result = FD31_DynamicCore(FileToOpen)  # calculate independent coordinates, and store them into more file intervals

    print(f"Completed {result} intervals")