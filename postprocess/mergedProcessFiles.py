import h5py
import numpy as np
import os
from datetime import datetime


# Version that merge and creates a backup of individual files before merging
def mergeFiles_and_backup(ime_datoteke1, keep_individual_files=True):
    """
    Function that merges multiple output files into one file with backup option.

    Parameters:
    -----------
    ime_datoteke1 : str
        Base filename (e.g., 'test')
    keep_individual_files : bool
        If True, keeps individual files after merging
        If False, deletes individual files after successful merge
    """

    # First open the first file
    ime_datoteke_del = f'0e_process_{ime_datoteke1}.h5'

    with h5py.File(ime_datoteke_del, 'r') as f:
        t = f['t'][:]  # Load time data
        Z = f['Z'][:]  # Load Z data
        ste_delitev = f.attrs['ste_delitev']  # Get number of divisions from attributes

    s = len(t)  # Number of entries
    # Delete the last entry (remove last row)
    Z = Z[:-1, :]  # Remove last row of Z
    t = t[:-1]  # Remove last element of t

    td = t.copy()
    Zd = Z.copy()

    # Collect all individual filenames for potential deletion
    individual_files = [f'0e_process_{ime_datoteke1}.h5']

    # Now open individual files and write to new file
    for i in range(1, ste_delitev ):
        ime_datoteke_del = f'{i}e_process_{ime_datoteke1}.h5'
        individual_files.append(ime_datoteke_del)

        with h5py.File(ime_datoteke_del, 'r') as f:
            t = f['t'][:]
            Z = f['Z'][:]

        if i < ste_delitev:
            s = len(t)  # Number of entries
            Z = Z[:-1, :]  # Remove last row
            t = t[:-1]  # Remove last element

            td = np.concatenate((td, t))  # Compose time vector
            Zd = np.concatenate((Zd, Z))  # Compose coordinates vector
        else:
            td = np.concatenate((td, t))  # Compose time vector
            Zd = np.concatenate((Zd, Z))  # Compose coordinates vector

    # Save merged results
    merged_filename = f'y_merged_file_{ime_datoteke1}.h5'

    with h5py.File(merged_filename, 'w') as f:
        # Save main arrays
        f.create_dataset('Z', data=Zd, compression='gzip')
        f.create_dataset('t', data=td, compression='gzip')

        # Add comprehensive metadata
        f.attrs['ste_delitev'] = int(ste_delitev)
        #f.attrs['merged_from'] = f'{ste_delitev} individual files'
        f.attrs['creation_date'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S').encode('ascii', 'ignore')
        #f.attrs['Z_dtype'] = str(Zd.dtype)
        #f.attrs['t_dtype'] = str(td.dtype)
        f.attrs['Z_shape'] = list(Zd.shape)
        f.attrs['t_shape'] = [len(td)]

    print(f"✓ Successfully merged {ste_delitev} files into: {merged_filename}")
    print(f"  Total entries: {len(td)}, Z dimensions: {Zd.shape}")

    # Delete individual files if requested
    if not keep_individual_files:
        import os
        deleted_count = 0
        for file in individual_files:
            if os.path.exists(file):
                os.remove(file)
                deleted_count += 1
        print(f"  Deleted {deleted_count} individual files")

    return Zd, td

