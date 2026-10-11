import numpy as np
import pandas as pd
from prestatic.create_classes import cable_class

def open_menu(FileToOpen):
    """
    Read data about route and determine conductor data
    """

    ## A. Read cable route data from file
    # Format: label, cleave_type (SP-start point, EP-endpoint, IC-intermediate cleave, CL-cleave), a (distance)x,  element division,
    try:
        # Using pandas to read the data
        df = pd.read_csv(FileToOpen, sep=r'\s+', header=None,
                         names=['label','cleave_type','a','element_division'])
        N = df
    except FileNotFoundError:
        print(f"File {FileToOpen} not found. Using empty DataFrame.")
        N = pd.DataFrame()

    #N['cleave_type'] = N['cleave_type'].astype('string')

    # Get filename without extension
    from pathlib import Path
    ime_datoteke1 = Path(FileToOpen).stem






    # predefined material properties
    materials = {
        'Aluminum': {
            'E_y': 69000,  # N/mm² (69 Gpa aluminium) for solid conductor
            'rho': 2.70e-3,  # g/mm³ (2.70 g/cm³)
            'color': 'blue'
        },
        'XLPE': {
            'E_y': 43,  # N/mm² (43 MPa) for annular section
            'rho': 0.92e-3,  # g/mm³ (0.92 g/cm³)
            'color': 'green'
        },
        'Copper': {
            'E_y': 110000,  # N/mm² (110 GPa Copper) for solid conductor
            'rho': 8.96e-3,  # g/mm³ (8.96 g/cm³)
            'color': 'orange'
        }
    }

    # Reading cable data

    #vrv_c = ['CableName',vrv_material, CrossSection nominal,inner diameter, outer diameter]
    vrv_c = ['CableName', 'Aluminum', 1000, 35.7, 84]  #



    ## Calculate EI, EA and mass properties for each conductor size and wall thickness
    # Flexural rigidity
    #diameters = 2 * np.sqrt(vrv_c[3] / np.pi)  # mm no need because it is given
    I_solid = (np.pi * vrv_c[3] ** 4) / 64  # mm⁴
    FR_solid = materials[vrv_c [1]]['E_y']* I_solid  # N·mm²
    # Second moment for annular section, insulation
    I_pipe = (np.pi * (vrv_c [4] ** 4 - vrv_c [3] ** 4)) / 64  # mm⁴
    FR_pipe = materials['XLPE']['E_y'] * I_pipe  # N·mm²
    # Combined flexural rigidity (conductor + pipe)
    EI_total = FR_solid+ FR_pipe

    # Axial rigidity
    #diameters = 2 * np.sqrt(vrv_c[3] / np.pi)  # mm no need because it is given
    A_solid = (np.pi * vrv_c[3] ** 2) / 4  # mm2
    EA_solid = materials[vrv_c [1]]['E_y']* A_solid  # N
    # Second moment for annular section, insulation
    A_pipe = (np.pi * (vrv_c [4] ** 2 - vrv_c [3] ** 2)) / 4  # mm2
    EA_pipe = materials['XLPE']['E_y'] * A_pipe  # N
    # Combined flexural rigidity (conductor + pipe)
    EA_total = EA_solid+ EA_pipe


    # Calculate masses
    core_masses = A_solid  * materials[vrv_c [1]]['rho']   # g/mm=kg/m
    xlpe_vol = (np.pi / 4 * (vrv_c[4] ** 2 - vrv_c[3] ** 2))   # mm³
    xlpe_mass = xlpe_vol * materials['XLPE']['rho']  # g/mm=kg/m
    mass_total = core_masses + xlpe_mass

    # Prepare for contact
    # EE is given in N/mm^2, converted to N/m^2 (Pascals) by multiplying by 1e6
    EE = EA_total/(A_solid+A_pipe) * 1e6
    # Poisson's ratio for aluminum
    poisson = 0.33
    # Diameter (mm) converted to meters, then divided by 4 to get R*.
    # (d_mm * 1e-3 converts mm to m. /2 gives radius. /2 again gives effective radius R* for crossed cylinders)
    Rstar = float(vrv_c[4])  * 1e-3 / 2 / 2
    # Effective elastic modulus ( Contact theory)
    Estar = EE / (2 * (1 - poisson ** 2))


    """
    'pass conductor data'
    """

    cab = cable_class
    cab.name=vrv_c[0]  # str
    cab.material=vrv_c[1]  # str
    cab.A=A_solid + A_pipe  # float mm2
    cab.d=float(vrv_c[4])  # float mm
    cab.m=mass_total  # float kg/m
    cab.EI=EI_total*1e-6  # float conversion to Nm2 (from Nmm2)
    cab.EA=EA_total  # float N
    cab.kn=1e6#4 * Rstar**0.5 * Estar / 3 #contact stiffness (sphere)
    cab.mi=0.01#0.1 # Contact damping


    return N, ime_datoteke1, cab