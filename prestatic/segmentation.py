import numpy as np
from numpy.ma.core import arctan

from prestatic.draw3DtestCoordinates import draw3DtestCoordinates

def segmentation_trifoil(N, cab):

    # A. Reminder conductor data
    # cab.name = vrv_c[0]  # str
    # cab.material = vrv_c[1]  # str
    # cab.A = A_solid + A_pipe  # float mm2
    # cab.d = float(vrv_c[4])  # float mm
    # cab.m = mass_total  # float kg/m
    # cab.EI = EI_total  # float N/mm2
    # cab.EA = EA_total  # float N
    global s
    g = 9.81  # gravitational constant m/s2

    # Convert route data N from vectors to arrays
    # Format of N: label, cleave_type (SP-start point, EP-endpoint, IC-intermediate cleave, CL-cleave), a (distance)x,  element division,

    pom = N['a']     # distance between cleaves
    hh = len(pom)
    #print('hh',hh)  hh is number of cleaves, that is hh=spans+1
    a = np.array(pom)

    # Cleave type definition
    pom = N['cleave_type']  #chain_type
    tipverige = np.array(pom)
    
    # Number of divisions per spans
    pom = N['element_division']
    ne_c = np.array(pom)

    # # Initial horizontal force for conductor
    # pom = N['tension']
    # S1_c = np.array(pom)
    #
    # # # Convert to force in N
    # # H1_c = np.zeros(hh)
    # # for i in range(hh-1):
    # #     H1_c[i] = S1_c[i] * 1e3
    #
    # # # Breaking angles (convert from degrees to radians)
    # # pom = N['route_angle']
    # # atd = np.array(pom)
    # # for i in range(hh-1):
    # #     atd[i] = atd[i] * PI / 180
    
    # Height points , the cleaves are in this height, only for more control
    h=np.zeros(hh)
    for i in range(hh):
        h[i] = 0
    


    # transform values to construct trefoil configuration, rows are conductors from 1 to 3, columns x and z movement,
    # y is always zero
    add=0.001 # for better control distance between cable,
    tm=np.array([[0,0],[cab.d+add,0],[(cab.d+add)/2,np.sqrt(3)*(cab.d+add)/2]])
    tm=tm*1e-3  # conversion from mm to m!



    # CALCULATION OF LOCAL COORDINATES
    # Find maximum number of divisions
    max_ne_c = int(np.max(ne_c[:hh-1])) + 1
    
    # Initialize matrices
    y_js = np.zeros((hh-1, max_ne_c))
    z_js = np.zeros((hh-1, max_ne_c))
    x_js = np.zeros((hh-1, max_ne_c))
    by_js = np.zeros((hh-1, max_ne_c))
    bz_js = np.zeros((hh-1, max_ne_c))
    bx_js = np.zeros((hh-1, max_ne_c))
    le = np.zeros((hh-1, max_ne_c))


    # Determine points along spans
    lo=a # the length of curved beam is almost equal to span, so to simplify
    for k in range(hh-1):
        pp = int(ne_c[k])
        ss = np.linspace(0, lo[k], pp + 1)
        for s in range(pp + 1):
            zz, by, bz = tocke(ss[s],lo[k],cab.EI,cab.m,g)
            y_js[k, s] =ss[s]
            z_js[k, s] = zz
            by_js[k, s] = by
            bz_js[k, s] = bz
            le[k, s] = lo[k] / ne_c[k] # it gives one value per row more than needed, but it is ok later

    # PREPARATION OF INITIAL CONDITIONS FOR ANCF COORDINATES

    e_zp2 = np.array([])  # vector containing initial ( estimated) ANCF coordinates

    for m in range (3):  # three conductors loop
        e_zp = np.array([])  # vector containing initial ( estimated) ANCF coordinates
        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # Convert z coordinate to final value (adjust z for his position in trefoil configuration )
        #zz_js = np.zeros(z_js.size)
        zz_js= -z_js + tm[m, 1]


        yy_js = y_js.copy() # just rewriting into new variable

        xx_js=np.zeros(zz_js.shape)+tm[m,0]

        # TRANSITION FROM LOCAL COORDINATES Y TO GLOBAL COORDINATES
        #  in local the spans are calculated separately in local ks, now assemble into line

        y0 = 0  # starting point of the line

        for k in range(hh-1):
            last_idx = int(ne_c[k])
            for s in range(last_idx+1):
                yy_js[k, s] += y0

            y0 = yy_js[k, last_idx]

        # print(zz_js)
        # print(yy_js)
        # print(xx_js)
        #breakpoint()
        # print('zz_js', zz_js.shape)
        ep=np.zeros(12) # helper

        for k in range(hh-1):
            for i in range(int(ne_c[k])):
                if k == 0 and i == 0:  # A. left fixing point - conductor connection

                    ep[0] = xx_js[k, i]
                    ep[1] = yy_js[k, i]
                    ep[2] = zz_js[k, i]
                    ep[3] = bx_js[k, i]
                    ep[4] = by_js[k, i]
                    ep[5] = bz_js[k, i]
                    ep[6] = xx_js[k, i+1]
                    ep[7] = yy_js[k, i+1]
                    ep[8] = zz_js[k, i+1]
                    ep[9] = bx_js[k, i+1]
                    ep[10] = by_js[k, i+1]
                    ep[11] = bz_js[k, i+1]
                    e_zp = np.concatenate((e_zp, ep))

                elif i >= 1 and i != int(ne_c[k]) - 1:  # B. conductor-conductor
                    ep[0] = xx_js[k, i]
                    ep[1] = yy_js[k, i]
                    ep[2] = zz_js[k, i]
                    ep[3] = bx_js[k, i]
                    ep[4] = by_js[k, i]
                    ep[5] = bz_js[k, i]
                    ep[6] = xx_js[k, i+1]
                    ep[7] = yy_js[k, i+1]
                    ep[8] = zz_js[k, i+1]
                    ep[9] = bx_js[k, i+1]
                    ep[10] = by_js[k, i+1]
                    ep[11] = bz_js[k, i+1]
                    e_zp = np.concatenate((e_zp, ep))

                elif i >= 1 and i == int(ne_c[k]) - 1 and k != hh-2:  # C. conductor-conductor-Cleave connection
                    # chain_type = int(tipverige[k+1])
                    #
                    # if chain_type == 0:
                    #     pass
                    # elif chain_type == 1:  # cleave
                        # Conductor element before insulator
                    ep[0] = xx_js[k, i]
                    ep[1] = yy_js[k, i]
                    ep[2] = zz_js[k, i]
                    ep[3] = bx_js[k, i]
                    ep[4] = by_js[k, i]
                    ep[5] = bz_js[k, i]
                    ep[6] = xx_js[k, i + 1]
                    ep[7] = yy_js[k, i + 1]
                    ep[8] = zz_js[k, i + 1]
                    ep[9] = bx_js[k, i + 1]
                    ep[10] = by_js[k, i + 1]
                    ep[11] = bz_js[k, i + 1]

                    e_zp = np.concatenate((e_zp, ep))


                    # Conductor element after insulator
                    ep[0] = xx_js[k+1, 0]
                    ep[1] = yy_js[k+1, 0]
                    ep[2] = zz_js[k+1, 0]
                    ep[3] = bx_js[k+1, 0]
                    ep[4] = by_js[k+1, 0]
                    ep[5] = bz_js[k+1, 0]
                    ep[6] = xx_js[k+1, 1]
                    ep[7] = yy_js[k+1, 1]
                    ep[8] = zz_js[k+1, 1]
                    ep[9] = bx_js[k+1, 1]
                    ep[10] = by_js[k+1, 1]
                    ep[11] = bz_js[k+1, 1]

                    e_zp = np.concatenate((e_zp, ep))

                elif i >= 1 and i == int(ne_c[k]) - 1 and k == hh-2:  # C. conductor - right cleve connection
                    ep[0] = xx_js[k, i]
                    ep[1] = yy_js[k, i]
                    ep[2] = zz_js[k, i]
                    ep[3] = bx_js[k, i]
                    ep[4] = by_js[k, i]
                    ep[5] = bz_js[k, i]
                    ep[6] = xx_js[k, i + 1]
                    ep[7] = yy_js[k, i + 1]
                    ep[8] = zz_js[k, i + 1]
                    ep[9] = bx_js[k, i + 1]
                    ep[10] = by_js[k, i + 1]
                    ep[11] = bz_js[k, i + 1]

                    e_zp = np.concatenate((e_zp, ep))
                    break

        # Create vertical vector
        e_zp2 = np.concatenate((e_zp2, e_zp))
    e_zp=e_zp2.copy()  # form out matrix
    #print(e_zp)
    #print('e_zp.shape',e_zp.shape)

    # Create initial conditions vector for dynamic simulation
    e_zp_s = np.zeros(2 * e_zp.size)
    locations = np.arange(0, 2 * e_zp.size, 2)
    e_zp_s[locations[:len(e_zp)]] = e_zp
    #print('e_zp_s',e_zp_s)

    #draw3DtestCoordinates(e_zp) help-draw e-zp for visualisation
    
    return e_zp, e_zp_s ,hh,ne_c,tipverige,le




def tocke(x,Lo,EI,m,g):
    """Calculate point coordinates along simply supported beam"""
    #zz = m*g/(24*EI)*(Lo-x)**2*x**2  # sag fixed-fixed beam under own self weight
    zz= m * g * x * (Lo**3 - 2*Lo*x**2 + x**3) / (24 * EI)  # simply supported
    #dzzdx = m*g/(24*EI)*(Lo-x)*(Lo-2*x) # inclination dzz/dx
    dzzdx= m * g * (6*Lo*x**2 - 4*x**3 - Lo**3) / (24 * EI)
    by = np.cos(arctan(dzzdx)) # derivation in y-axis
    bz = np.sin(arctan(dzzdx)) # derivation in z-axis
    # note that zz is calculated as positive curve so later is subtracted
    return zz, by, bz


