#!/usr/bin/env python
"""
Andrew's mapping to nearest image borrowed from DL_Poly

06/04/16
"""
import numpy as np
import sys

def invertcellmatrix(cell):
    """
    invert the matrix corresponding to the 1-d array cell for nearest image mapping
    """
    invcell = np.zeros(9, dtype=float)

    # adjoint
    invcell[0] = cell[4]*cell[8] - cell[5]*cell[7]
    invcell[1] = cell[2]*cell[7] - cell[1]*cell[8]
    invcell[2] = cell[1]*cell[5] - cell[2]*cell[4]
    invcell[3] = cell[5]*cell[6] - cell[3]*cell[8]
    invcell[4] = cell[0]*cell[8] - cell[2]*cell[6]
    invcell[5] = cell[2]*cell[3] - cell[0]*cell[5]
    invcell[6] = cell[3]*cell[7] - cell[4]*cell[6]
    invcell[7] = cell[1]*cell[6] - cell[0]*cell[7]
    invcell[8] = cell[0]*cell[4] - cell[1]*cell[3]

    # determinant
    det = cell[0]*invcell[0] + cell[3]*invcell[1] + cell[6]*invcell[2]

    if det == 0:
        raise ValueError('*** error *** : indeterminant matrix')
    else:
        tmp = 1.0/det
        for ia in xrange(9):
            invcell[ia] = invcell[ia]*tmp
        return invcell
######################

def mapdisplacement(dr,cell,invcell):
    """
    take the 1-d array corresponding to matric of cell vectors and its inverse, 
    to return the displacement dr to it's unique nearest-image value
    """

    newdr = np.zeros(3, dtype=float)

    alpha = invcell[0]*dr[0] + invcell[3]*dr[1] + invcell[6]*dr[2]
    beta = invcell[1]*dr[0] + invcell[4]*dr[1] + invcell[7]*dr[2] 
    gamma = invcell[2]*dr[0] + invcell[5]*dr[1] + invcell[8]*dr[2] 

    alpha_m = alpha - round(alpha)   
    beta_m = beta - round(beta)      
    gamma_m = gamma - round(gamma)          

    newdr[0] = cell[0]*alpha_m + cell[3]*beta_m + cell[6]*gamma_m       
    newdr[1] = cell[1]*alpha_m + cell[4]*beta_m + cell[7]*gamma_m        
    newdr[2] = cell[2]*alpha_m + cell[5]*beta_m + cell[8]*gamma_m  

    return newdr


#dr = mapdisplacement(dr,cell,invcell)


# Addition from Peter
L1 = np.array([1, 0, 0])
L2 = np.array([0, 1, 0])
L3 = np.array([0, 0, 1])
cell_vecs = np.vstack((L1, L2, L3))
try:
    inv_cell = np.linalg.pinv(cell_vecs)
except LinAlgError:
    print "Singular matrix, aborting."
    sys.exit()

dr = np.random.rand(3)
G = np.dot(inv_cell, dr)
G_n = G - np.round(G)
dr_n = np.dot(cell_vecs, G_n)
print dr, "\n", dr_n



