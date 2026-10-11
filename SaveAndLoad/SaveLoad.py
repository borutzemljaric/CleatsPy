import h5py
from datetime import datetime
import numpy as np
import pandas as pd
from dataclasses import dataclass
import scipy.sparse as sp

from prestatic.create_classes import cable_class

def save_open_menu(filename, N, ime_datoteke1, cab):
    with pd.HDFStore(filename, mode='w') as store:
        store.put('N', N, format='fixed')
        store.root._v_attrs['ime_datoteke1'] = ime_datoteke1
        for key, value in vars(cab).items():
            if key.startswith('__'):
                continue
            store.root._v_attrs[f'cab_{key}'] = value
        store.root._v_attrs['save_timestamp'] = pd.Timestamp.now().isoformat()



def update_cab_optional(filename, cab):
    """
    Overwrite only the optional attributes of cab in an existing HDF5 file.
    The file must have been created with save_open_menu (which stores cab_* attributes).
    """
    # List of optional fields (those that can be None)
    optional_fields = ['tfinal', 'tsim', 'trelay','interval', 'ste_delitev', 'note', 'Ik', 'phi','Nseg','Ncon','a','b']

    # Open in append mode (creates if missing, allows updating)
    with pd.HDFStore(filename, mode='a') as store:
        attrs = store.root._v_attrs
        for field in optional_fields:
            value = getattr(cab, field, None)
            if value is not None:
                # Update the attribute (create or overwrite)
                setattr(attrs, f'cab_{field}', value)



def save_segmentation(filename, e_zp, e_zp_s ,hh,ne_c,tipverige,le):
    """
    Save segmentation results from segmentacija_ANCF
    """
    with h5py.File(filename, 'a') as f:
        # Save main arrays
        f.create_dataset('e_zp', data=e_zp, compression='gzip')
        f.create_dataset('e_zp_s', data=e_zp_s, compression='gzip')
        f.create_dataset('ne_c', data=ne_c)
        f.create_dataset('tipverige', data=tipverige)
        f.create_dataset('le', data=le, compression='gzip')

        # Save scalars as attributes
        f.attrs['hh'] = int(hh)
        f.attrs['saved_date'] = datetime.now().isoformat()

        # Add shape information for easy loading
        f.attrs['e_zp_shape'] = str(e_zp.shape)
        f.attrs['e_zp_s_shape'] = str(e_zp_s.shape)

    #print(f"Segmentation saved to {filename}")
    #print(f"  - e_zp shape: {e_zp.shape}")
    #print(f"  - e_zp_s shape: {e_zp_s.shape}")

def save_CeActive(filename,Ce_active,  konstante, neod_fix ):
    """
    Save segmentation results from segmentacija_ANCF
    """
    with h5py.File(filename, 'a') as f:
        # Save main arrays
        f.create_dataset('Ce_active', data=Ce_active, compression='gzip')

        f.create_dataset('konstante', data=konstante, compression='gzip')
        f.create_dataset('neod_fix', data=neod_fix, compression='gzip')

    #print(f"Segmentation saved to {filename}")

def save_Bmatrix(filename,Bdi, U, V, pari, Cdi, Ced, Ce ):
    """
    Save segmentation results from segmentacija_ANCF
    """
    with h5py.File(filename, 'a') as f:
        # Save main arrays
        f.create_dataset('Bdi', data=Bdi, compression='gzip')
        f.create_dataset('U', data=U, compression='gzip')
        f.create_dataset('V', data=V, compression='gzip')
        f.create_dataset('pari', data=pari, compression='gzip')
        f.create_dataset('Cdi', data=Cdi, compression='gzip')
        f.create_dataset('Ced', data=Ced, compression='gzip')
        f.create_dataset('Ce', data=Ce, compression='gzip')

    #print(f"Segmentation saved to {filename}")


def save_staticANCF(filename,zacpog,  Qe_g):
    """
    Save segmentation results from segmentacija_ANCF
    """
    with h5py.File(filename, 'a') as f:
        # Save main arrays
        f.create_dataset('zacpog', data=zacpog, compression='gzip')
        f.create_dataset('Qe_g', data=Qe_g, compression='gzip')

    #print(f"Segmentation saved to {filename}")


# Add this function definition after
def save_H_M(filename,H_ll, H_ee, H_el, H_le, M):
    """
    Save  results from H_M_matrike
    """

    with h5py.File(filename, 'a') as f:
        # Save main arrays
        f.create_dataset('H_ll', data=H_ll, compression='gzip')
        f.create_dataset('H_ee', data=H_ee, compression='gzip')
        f.create_dataset('H_el', data=H_el, compression='gzip')
        f.create_dataset('H_le', data=H_le, compression='gzip')
        f.create_dataset('M',    data=M,    compression='gzip')


        # Add shape information for easy loading
        f.attrs['M_shape'] = str(M.shape)


    #print(f"Segmentation saved to {filename}")

    # Save automatically
    #ime_datoteke1,H_ll, H_ee, H_el, H_le, M

# Add this function definition after
def save_reactions(filename, Flevo, Fdesno):
    """
    Save  results from reactions
    """
    with h5py.File(filename, 'a') as f:
        # Save main arrays
        f.create_dataset('Flevo', data=Flevo, compression='gzip')
        f.create_dataset('Fdesno', data=Fdesno, compression='gzip')




"""
load
"""


def load_saved_variables(filename):
    """
    Load all saved variables from the HDF5 file and return them in original format.

    Returns:
    - zacpog, zacpog_db, zacpog_db_W: numpy arrays
    - Flevo, Fdesno, Folevo, Fodesno: numpy arrays
    - izo: dictionary of numpy arrays
    - hh: integer
    - ne_c: integer
    - tip_verige: array
    - e_zp: array
    - le: array
    """


    def load_object(h5_obj):
        """Recursively load objects from HDF5."""
        if isinstance(h5_obj, h5py.Group):
            # Check if it's a DataFrame (has 'data' dataset and 'columns' attribute)
            if 'data' in h5_obj and 'columns' in h5_obj.attrs:
                # Reconstruct DataFrame
                data = h5_obj['data'][()]
                columns = list(h5_obj.attrs['columns'])

                # Create DataFrame from structured array
                if data.dtype.names:  # It's a structured array
                    df = pd.DataFrame.from_records(data)
                else:  # Regular array
                    df = pd.DataFrame(data, columns=columns)

                # Restore dtypes if saved
                if 'dtypes' in h5_obj.attrs:
                    dtypes = h5_obj.attrs['dtypes']
                    for col, dtype_str in zip(columns, dtypes):
                        if 'float' in dtype_str:
                            df[col] = df[col].astype(float)
                        elif 'int' in dtype_str:
                            df[col] = df[col].astype(int)
                        elif 'object' in dtype_str or 'string' in dtype_str:
                            df[col] = df[col].astype(str)

                return df
            else:
                # Regular group -> dictionary
                result = {}
                for key in h5_obj.keys():
                    result[key] = load_object(h5_obj[key])

                # Also include attributes that aren't datasets
                for attr_name in h5_obj.attrs:
                    # Skip shape attributes (handled separately)
                    if '_shape' not in attr_name:
                        result[attr_name] = h5_obj.attrs[attr_name]

                return result
        elif isinstance(h5_obj, h5py.Dataset):
            data = h5_obj[()]

            # Convert bytes to string if needed
            if isinstance(data, bytes):
                return data.decode('utf-8')
            elif isinstance(data, np.ndarray) and data.dtype.kind == 'S':
                # Convert string arrays back to Python strings
                return data.astype(str)

            return data
        return None

    with h5py.File(filename, 'r') as f:
        # 1. Load izo (dictionary of numpy arrays)
        izo = {}
        if 'izo' in f:
            izo_group = f['izo']
            # Check if izo was saved as a group (dictionary)
            if isinstance(izo_group, h5py.Group):
                izo = load_object(izo_group)
                # Convert any non-array values to numpy arrays
                for key, value in izo.items():
                    if not isinstance(value, np.ndarray):
                        izo[key] = np.array(value)
            else:
                # It might be a dataset (direct array)
                izo = {'izo': izo_group[()]}

        # 2. Load segmentation variables
        hh = None
        if 'hh' in f.attrs:
            hh = int(f.attrs['hh'])

        ne_c = None
        if 'ne_c' in f:
            ne_c_data = f['ne_c'][()]
            # Convert to appropriate type
            if isinstance(ne_c_data, np.ndarray):
                ne_c = ne_c_data
            else:
                ne_c = np.array(ne_c_data)

        tip_verige = None
        if 'tipverige' in f:
            tipverige = f['tipverige'][()]

        e_zp = None
        if 'e_zp' in f:
            e_zp = f['e_zp'][()]

        le = None
        if 'le' in f:
            le = f['le'][()]

        # 3. Load static calculation results
        zacpog = None
        if 'zacpog' in f:
            zacpog = f['zacpog'][()]

        zacpog_db = None
        if 'zacpog_db' in f:
            zacpog_db = f['zacpog_db'][()]

        zacpog_db_W = None
        if 'zacpog_db_W' in f:
            zacpog_db_W = f['zacpog_db_W'][()]

        # 4. Load reaction forces
        Flevo = None
        if 'Flevo' in f:
            Flevo = f['Flevo'][()]

        Fdesno = None
        if 'Fdesno' in f:
            Fdesno = f['Fdesno'][()]

        Folevo = None
        if 'Folevo' in f:
            Folevo = f['Folevo'][()]

        Fodesno = None
        if 'Fodesno' in f:
            Fodesno = f['Fodesno'][()]

        # 5. Load include_wind flag
        include_wind = None
        if 'include_wind' in f.attrs:
            include_wind = bool(f.attrs['include_wind'])

        # Print loaded info
        print("Loaded variables:")
        print(f"  izo keys: {list(izo.keys()) if izo else 'None'}")
        print(f"  hh: {hh}")
        print(f"  ne_c shape: {ne_c.shape if hasattr(ne_c, 'shape') else ne_c}")
        print(f"  tip_verige shape: {tipverige.shape if hasattr(tipverige, 'shape') else tip_verige}")
        print(f"  e_zp shape: {e_zp.shape if hasattr(e_zp, 'shape') else e_zp}")
        print(f"  le shape: {le.shape if hasattr(le, 'shape') else le}")
        print(f"  zacpog shape: {zacpog.shape if hasattr(zacpog, 'shape') else zacpog}")
        print(f"  zacpog_db shape: {zacpog_db.shape if hasattr(zacpog_db, 'shape') else zacpog_db}")
        print(f"  zacpog_db_W shape: {zacpog_db_W.shape if hasattr(zacpog_db_W, 'shape') else zacpog_db_W}")
        print(
            f"  Reaction forces loaded: {Flevo is not None}, {Fdesno is not None}, {Folevo is not None}, {Fodesno is not None}")

        return (zacpog, zacpog_db, zacpog_db_W, Flevo, Fdesno, Folevo, Fodesno,
                izo, hh, ne_c, tipverige, e_zp, le)

def load_StaticResults(filename):
    """
    Load data saved by save_open_menu.
    Returns: (cab, e_zp, e_zp_s, hh, ne_c, tipverige, le)
    """
    with h5py.File(filename, mode='r') as store:

        # 2. Reconstruct cab object from attributes
        cab_dict = {}
        for key, value in store.attrs.items():
            if key.startswith('cab_'):
                attr_name = key[4:]          # remove 'cab_' prefix
                cab_dict[attr_name] = value

        # Recreate the cable_class instance
        # Make sure cable_class is defined in your namespace
        cab = cable_class(**cab_dict)

        # 3. Load all arrays and scalars
        e_zp = store['e_zp'][:]
        e_zp_s = store['e_zp_s'][:]
        hh = store.attrs['hh']               # scalar attribute
        ne_c = store['ne_c'][:]
        tipverige = store['tipverige'][:]
        le = store['le'][:]
        zacpog = store['zacpog'][:]

    return cab, e_zp, e_zp_s, hh, ne_c, tipverige, le,zacpog


def load_menuForDynamic(filename):
    """
    Load data saved in static program part.
    Returns: (...)
    """
    with h5py.File(filename, mode='r') as store:

        # 2. Reconstruct cab object from attributes
        cab_dict = {}
        for key, value in store.attrs.items():
            if key.startswith('cab_'):
                attr_name = key[4:]          # remove 'cab_' prefix
                cab_dict[attr_name] = value

        # Recreate the cable_class instance
        # Make sure cable_class is defined in your namespace
        cab = cable_class(**cab_dict)

        # 3. Load all arrays and scalars
        e_zp = store['e_zp'][:]
        e_zp_s = store['e_zp_s'][:]
        hh = store.attrs['hh']               # scalar attribute
        ne_c = store['ne_c'][:]
        tipverige = store['tipverige'][:]
        le = store['le'][:]
        zacpog = store['zacpog'][:]
        Ce_active = store['Ce_active'][:]
        Bdi = store['Bdi'][:]
        U = store['U'][:]
        V = store['V'][:]
        Cdi = store['Cdi'][:]
        Ced = store['Ced'][:]
        M = store['M'][:]
        Qe_g = store['Qe_g'][:]


    return cab, e_zp, e_zp_s, hh, ne_c, tipverige, le,zacpog, Ce_active, Bdi, U, V, Cdi, Ced, M, Qe_g




"""
save dynamic calculation
"""
def save_dynamic_process(filename,Z, t, ste_delitev):
    """
    Save  results from H_M_matrike
    """

    with h5py.File(filename, 'w') as f:
        # Save main arrays
        f.create_dataset('Z', data=Z, compression='gzip')
        f.create_dataset('t', data=t, compression='gzip')

        # Add shape information for easy loading
        f.attrs['ste_delitev'] = int(ste_delitev)

    #print(f" Dynamical processing saved to {filename}")


def save_dynamic_postprocess(filename, t, izhod_e_a, izhod_e_dot, izhod_Flevo, izhod_Fdesno):
    """
    Save segmentation results from postprocessing
    """
    with h5py.File(filename, 'a') as f:
        # Before creating your new datasets, delete any existing ones with the same names
        # In the case of re-lunching postprocessing
        for name in ['t', 'izhod_e_a', 'izhod_e_dot','izhod_Flevo','izhod_Fdesno']:
            if name in f:
                del f[name]
        f.create_dataset('t', data=t, compression='gzip')
        f.create_dataset('izhod_e_a', data=izhod_e_a, compression='gzip')
        f.create_dataset('izhod_e_dot', data=izhod_e_dot, compression='gzip')
        f.create_dataset('izhod_Flevo', data=izhod_Flevo, compression='gzip')
        f.create_dataset('izhod_Fdesno', data=izhod_Fdesno, compression='gzip')


        # Add shape information for easy loading
        #f.attrs['e_zp_shape'] = str(e_zp.shape)
    #print(f"Postprocess saved to {filename}")




def load_for_plot(filename):
    """
    Load  results from HDF5 file saved by postprocess.

    Parameters:
    -----------
    filename : str
        Path to the HDF5 file

    Returns:
    --------
    tuple : (hh, ne_c, tipverige, zacpog, zacpog_db, t,
             izhod_e_a, izhod_e_dot, izhod_Flevo, izhod_Fdesno)
    """


    with h5py.File(filename, 'r') as f:
        # Load scalar attributes
        hh = int(f.attrs['hh'])
        ime_datoteke1 = f.attrs.get('ime_datoteke1', 'unknown')
        saved_date = f.attrs.get('saved_date', 'unknown')

        # Load datasets
        zacpog = f['zacpog'][()]
        t = f['t'][()]
        ne_c = f['ne_c'][()]
        tipverige = f['tipverige'][()]
        izhod_e_a = f['izhod_e_a'][()]
        izhod_e_dot = f['izhod_e_dot'][()]
        izhod_Flevo = f['izhod_Flevo'][()]
        izhod_Fdesno = f['izhod_Fdesno'][()]

        # Print summary
        # print(f"✓ Loaded postprocessing data from: {filename}")
        # print(f"  File saved: {saved_date}")
        # print(f"  Original file name: {ime_datoteke1}")
        # print(f"  Number of spans (hh): {hh}")
        # print(f"  Time steps: {len(t)}")
        # print(f"  Time range: {t[0]:.2f} to {t[-1]:.2f} s")
        # print(f"  Displacement matrix shape: {izhod_e_a.shape}")
        # print(f"  Force matrices shapes: {izhod_Flevo.shape}, {izhod_Fdesno.shape}")

        return hh, ne_c, tipverige, zacpog, t, izhod_e_a, izhod_e_dot, izhod_Flevo, izhod_Fdesno




