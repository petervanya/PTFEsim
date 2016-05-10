#!/usr/bin/env python3
"""Usage:
    parse_ms.py <infile>

Parse a Materials Studio xcf file and create two files, one with atoms and position, the other with bonds.

pv278@cam.ac.uk, 10/05/16
"""
import numpy as np
from docopt import docopt
import lmp_lib as ll


if __name__ == "__main__":
    args = docopt(__doc__)
    infile = args["<infile>"]
    beads = {"A": 1, "B": 2, "C": 3, "W": 4, "E": 5}

    print("File name:", infile)
    xmlstring = open(infile, "r").readlines()
    print(len(xmlstring))
    N = len([line for line in xmlstring if "Bead ID=" in line])
    xyz_mat = np.zeros((N, 3))
    cnt = 0

    # Electrodes
    E = [line for line in xmlstring if "Bead ID" in line and "ForcefieldType=\"G\"" in line]
    NE = len(E)
    print("Number of G beads:", NE)
    for line in E:
        key = [field for field in line.split() if "XYZ=" in field][0]
        key = np.array(key[5:-1].split(",")).astype(float)
        cnt += 1
    del E

    # Bead A
    A = [line for line in xmlstring if "Bead ID" in line and "ForcefieldType=\"A\"" in line]
    NA = len(A)
    print("Number of A beads:", NA)
    for line in A:
        key = [field for field in line.split() if "XYZ=" in field][0]
        key = np.array(key[5:-1].split(",")).astype(float)
        cnt += 1
    del A

    # Bead B
    B = [line for line in xmlstring if "Bead ID" in line and "ForcefieldType=\"B11\"" in line]
    NB = len(B)
    print("Number of B beads:", NB)
    for line in B:
        key = [field for field in line.split() if "XYZ=" in field][0]
        key = np.array(key[5:-1].split(",")).astype(float)
        cnt += 1

    # Bead C
    C = [line for line in xmlstring if "Bead ID" in line and "ForcefieldType=\"C\"" in line]
    NC = len(C)
    print("Number of C beads:", NC)
    for line in C:
        key = [field for field in line.split() if "XYZ=" in field][0]
        key = np.array(key[5:-1].split(",")).astype(float)
        cnt += 1

    # Bead W
    W = [line for line in xmlstring if "Bead ID" in line and "ForcefieldType=\"W\"" in line]
    NW = len(W)
    print("Number of W beads:", NW)
    for line in W:
        key = [field for field in line.split() if "XYZ=" in field][0]
        key = np.array(key[5:-1].split(",")).astype(float)
        cnt += 1
    del W

    bead_vec = np.matrix([beads["A"]]*NA + [beads["B"]]*NB + [beads["C"]]*NC + \
               [beads["W"]]*NW + [beads["E"]]*NE).T
    xyz_mat = np.hstack((bead_vec, xyz_mat))
    fname = "nafion_ms.xyz"
    ll.save_xyzfile(fname, xyz_mat)
    print("xyz matrix was saved into", fname)

