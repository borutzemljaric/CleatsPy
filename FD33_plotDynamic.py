
import numpy as np
from prestatic.create_classes import plot_class
from pathlib import Path
from plots.load_curve_data import load_curve_data
from plots.g_xyz_t import g_xyz_t
from plots.g_ABSdisplacement_t import g_ABSdisplacement_t
from plots.g_cComponentForce_t import g_cComponentForce_t
from plots.g_cCompForceMid_t import g_cCompForceMid_t
from plots.g_xy_at_given_t import g_xy_at_given_t
from plots.g_cInCleatForces import g_cIncleatForces
from plots.g_cInCleatForces2D import g_cIncleatForces2D



def FD33_plotDynamic(file_names,gr):
    """
    Function for plotting dynamic simulation results.
    """

    new_files = []
    for f in file_names:
        if f is not None:  # skip None
            stem = Path(f).stem  # gets filenames without extension
            new_name = f'ySR_{stem}.h5'  # creates 'ySR_yDR_FILENAME.h5', etc.
            new_files.append(new_name)

    #print(new_files)


    # Usage
    data = load_curve_data(gr.n_files, new_files)




    # Switch case for different plot types
    if gr.plot_type == 0:
        relative=True  # if True show relative distance from initial position, otherwise absolute position
        g_xyz_t(gr.k,data,relative)


    elif gr.plot_type == 1:
        relative = True  # if True show relative distance from initial position, otherwise absolute position
        drawFigure = True#True  # control to plot or not results
        g_ABSdisplacement_t(gr.k,data,relative,drawFigure)


    elif gr.plot_type == 2:
        g_cComponentForce_t(gr.k, data, relative=False)


    elif gr.plot_type == 3:
        g_cCompForceMid_t(gr.k, data, relative=False)


    elif gr.plot_type == 4:
        g_xy_at_given_t(data)


    elif gr.plot_type == 5:
        #g_cIncleatForces(gr.k,data)
        g_cIncleatForces2D(gr.k,data)  # testing 2D results



# ============================================================================
# USAGE
# ============================================================================
if __name__ == "__main__":

    # Defining and Calling external files with  data to be presented
    # it is possible to define up to 5 files

    # File names
    file_names = ['primer_4.txt', 'TestFileCleavesX.txt', 'primer_5.txt',
                  None, None]

    gr=plot_class
    # User parameters
    gr.k = 3  # span or cleave index to plot  (first span is 1, internally to read python arrays become index 0)
    gr.n_files= 1 # number of files to be opened and compared, Temporary only one!
    gr.n_curves = 1  # number of curves to plot, that is used if only one file is opened
    gr.plot_type =5  # type of graph

    # plot type
    # 0- xyz by component displacement, cable in the middle of span, dimensional
    # 1- normalised ABSOLUTE displacement in the middle of span
    # 2 Forces by components x,y,z in Cleats
    # 3 Forces internal in cable by components x,y,z in Midspan
    # 4 Displacements x,z along cable line y by time steps
    # 5 Forces Axial InsideCleat


    # still unused
    # 6
    # 7
    # 8

    FD33_plotDynamic(file_names,gr)
