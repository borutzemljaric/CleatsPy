from SaveAndLoad.SaveLoad import load_menuForDynamic
from plots.readCoordinate import readCoordinate
from plots.readForces import readForces
from plots.readForcesMidspan import readForcesMidspan
from plots.readCoordinateALL import readCoordinateALL

def load_curve_data(n_krivulj,file_names):
    """
    Load data for n_curves
    """
    #print(n_krivulj,file_names)

    results = {'t': None}

    # Validate input
    if n_krivulj < 1 or n_krivulj > 5:
        raise ValueError("number of curves must be between 1 and 5")

    # Load displacement data in MIDSPAN
    for i in range(n_krivulj):
        try:
            t_curr, X_i, Y_i, Z_i = readCoordinate(file_names[i])
            cab, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _ = load_menuForDynamic(file_names[i])

            if i == 0:
                results['t'] = t_curr
                results['cab_d'] = cab.d

            # Store all data
            for name, value in [
                ('X', X_i), ('Y', Y_i), ('Z', Z_i)
            ]:
                results[f'{name}_{i + 1}'] = value

        except IndexError:
            print(f"Warning: File {file_names[i]} not available")
            print(f"Warning: File {i} not available")
            break


    # Load force data in CLEATES for function 2- g_cComponentForce_t(gr.k, data, relative=False)
    for i in range(n_krivulj):

        try:
            t_curr, FXcleave_i,FYcleave_i,FZcleave_i= readForces(file_names[i])

            #print('izhod', t_curr, Flevo_i, Fdesno_i)
            # Fleve are absolute values, FXlevo, FYleco and FZlevo are the component forces in given direction

            for name, value in [
                ('FXcleave', FXcleave_i),('FYcleave', FYcleave_i),('FZcleave', FZcleave_i)
            ]:
                results[f'{name}_{i + 1}'] = value

        except IndexError:
            print(f"Warning: Force data for {file_names[i]} not available")
            break

    # Load force data in MIDSPAN for function 3- g_cCompForceMid_t(gr.k, data, relative=False)
    for i in range(n_krivulj):

        try:
            t_curr, FXmid_i, FYmid_i, FZmid_i = readForcesMidspan(file_names[i])

            # print('izhod', t_curr, Flevo_i, Fdesno_i)
            # Fleve are absolute values, FXlevo, FYleco and FZlevo are the component forces in given direction

            for name, value in [
                ('FX_mid', FXmid_i), ('FY_mid', FYmid_i), ('FZ_mid', FZmid_i)
            ]:
                results[f'{name}_{i + 1}'] = value

        except IndexError:
            print(f"Warning: Force data for {file_names[i]} not available")
            break

    # Load displacement data along cable line- for function 4, g_xy_at_given_t(data)
    for i in range(n_krivulj):
        try:
            t_curr, X_i, Y_i, Z_i = readCoordinateALL(file_names[i])

            if i == 0:
                results['t'] = t_curr

            # Store all data
            for name, value in [
                ('allX', X_i), ('allY', Y_i), ('allZ', Z_i)
            ]:
                results[f'{name}_{i + 1}'] = value

        except IndexError:
            print(f"Warning: File {file_names[i]} not available")
            print(f"Warning: File {i} not available")
            break

    return results
