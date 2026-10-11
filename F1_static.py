import numpy as np
from prestatic.open_menu import open_menu
from prestatic.segmentation import segmentation_trifoil
from prestatic.matrix_Ce_B import matrix_Ce,B_matrix
from prestatic.calculateStaticCoordinates import calculateStaticCoordinates
from prestatic.H_M_matrike import H_M_matrike
from prestatic.reactionNodes import reakctionNodes
from SaveAndLoad.SaveLoad import save_open_menu, save_segmentation,save_CeActive,save_Bmatrix
from SaveAndLoad.SaveLoad import save_staticANCF, save_H_M,save_reactions


def F1_static(FileToOpen):
    """
    Main function for static analysis of cleaves
    """
    # 1. Open data file with input data
    print("Reading data")
    N, ime_datoteke1, cab = open_menu(FileToOpen)
    filename = f'ySR_{ime_datoteke1}.h5' #Save  results into file with partially predefined name
    save_open_menu(filename, N, ime_datoteke1, cab) #save results

    # # 2. Execute segmentation of conductors in spans
    print("Line segmentation")
    e_zp, e_zp_s ,hh,ne_c,tipverige,le=segmentation_trifoil(N, cab)
    save_segmentation(filename,e_zp, e_zp_s ,hh,ne_c,tipverige,le) #save results


    # # 3. Calculate Jacobian determinant Ce
    print("Prepare Jacobian")
    Ce_active, konstante, neod_fix = matrix_Ce(hh, ne_c, tipverige)
    save_CeActive(filename, Ce_active,  konstante, neod_fix) #save results


    # # 4. Calculate dependent/independent coordinates
    print("Determine dependant/independent coordinates")
    Bdi, U, V, pari, Cdi, Ced, Ce= B_matrix(Ce_active,konstante,neod_fix)
    save_Bmatrix(filename, Bdi, U, V, pari, Cdi, Ced, Ce) #save results


    # # 5. Choose whether to calculate with or without wind
    print("Calculating coordinates")
    zacpog, Qe_g=calculateStaticCoordinates( cab, hh, ne_c, le,
                            tipverige, e_zp, Bdi,U,V,konstante,pari)  # static coordinates calculation
    save_staticANCF(filename, zacpog, Qe_g) #save results


    # # 6. Prepare mass matrix and H matrices for dynamic calculation
    print("Prepare mass matrix and H matrices")
    H_ll, H_ee, H_le, H_el, M=H_M_matrike(hh, ne_c, cab, le , Ce_active)
    save_H_M(filename, H_ll, H_ee, H_el, H_le, M) #save results


    # # 7. Calculate  reactions
    print("Calculate  reactions forces")
    Flevo, Fdesno=reakctionNodes(zacpog, e_zp,Qe_g, hh, ne_c, le, tipverige,cab)
    save_reactions(filename, Flevo, Fdesno) #save results

    print("End of static program")






# ============================================================================
# Usage- static calculation
# ============================================================================
if __name__ == "__main__":

    # A. Calling external file with input data about power line field
    FileToOpen = 'TestFileCleavesX.txt'
    FileToOpen = 'primer_1.txt'
    F1_static(FileToOpen)