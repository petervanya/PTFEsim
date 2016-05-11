#!/usr/bin/env python3
"""Usage:
    parse_ms.py <infile>

Parse a Materials Studio xcf file and create two files, one with atoms and position, the other with bonds.

pv278@cam.ac.uk, 10/05/16
"""
import numpy as np
from docopt import docopt
import lmp_lib as ll


def parse_beads(xmlstring, bead_dict):
    """XML string has a form of lines,
    return:
    * xyz matrix (N, 3)
    * bead types vector (N, 1)
    * bead IDs vector (N, 1)
    * molecule IDs, vector (N, 1)
    """
    # Examples of lines:
    # <Bead ID="60669" Mapping="60486" Parent="6" UserID="93" XYZ="0.335008132054088,0.56406164179294,0.35715898635349"
    #    Connections="60670,60672,60746" ForcefieldType="A" BeadType="388551"/>
    # <Molecule ID="6" Mapping="60486" Parent="2" Children="ci(149):149+60636" Name="Nafion_EW1244_MW11040"/>
    N = len([line for line in xmlstring if "Bead ID=" in line])
    xyz_mat = np.zeros((N, 3), dtype=float)
    bead_IDs = np.zeros(N, dtype=int)
    user_IDs = np.zeros(N, dtype=int)
    mol_IDs = np.zeros(N, dtype=int)
    cnt = 0
    Nb = {}

    for b in bead_dict.keys():
        E = [line for line in xmlstring if "<Bead ID" in line and "ForcefieldType=\"%s\"" % b in line]
        Nb[b] = len(E)
        print("Number of %s beads: %i" % (b, Nb[b]))
        for line in E:
            key = [field for field in line.split() if "XYZ=" in field][0]
            bead_ID = [field for field in line.split() if "ID=" in field][0]
            user_ID = [field for field in line.split() if "UserID=" in field][0]
            mol_ID = [field for field in line.split() if "Parent=" in field][0]
            key = np.array(key[5:-1].split(",")).astype(float)
 
            xyz_mat[cnt] = key
            bead_IDs[cnt] = int(bead_ID[4:-1])
            user_IDs[cnt] = int(user_ID[8:-1])
            mol_IDs[cnt] = int(mol_ID[8:-1])
            cnt += 1
        del E

    bead_types = []
    for b in bead_dict.keys():
        bead_types += [bead_dict[b]]*Nb[b]
    bead_types = np.array(bead_types)
#    parsed_mat = np.hstack((bead_types, xyz_mat))
    return xyz_mat, bead_types, bead_IDs, user_IDs, mol_IDs


def parse_bonds(xmlstring):
    """Extract bonds from a XML string in the form of lines,
    return (N, 3) matrix with columns:
    * bond ID
    * bead 1 ID
    * bead 2 ID
    """
    # <BeadConnector ID="60670" Mapping="60486" Parent="6" Connects="60667,60669"/>
    N = len([line for line in xmlstring if "<BeadConnector" in line])
    parsed_mat = np.zeros((N, 3), dtype=int)
    cnt = 0

    good = [line for line in xmlstring if "<BeadConnector" in line]
    for line in good:
        bond_id = [field for field in line.split() if "ID=" in field][0]
        beads = [field for field in line.split() if "Connects=" in field][0]
        beads = np.array(beads[10:-3].split(",")).astype(int)

        parsed_mat[cnt, 0] = int(bond_id[4:-1])
        parsed_mat[cnt, 1:] = beads
        cnt += 1

    return parsed_mat


if __name__ == "__main__":
    args = docopt(__doc__)
    infile = args["<infile>"]
    bead_dict = {"A": 1, "B11": 2, "C": 3, "W": 4, "G": 5}

    xmlstring = open(infile, "r").readlines()
    print("File name: %s, length: %i" % (infile, len(xmlstring)))
    
    print("Parsing data...")
    xyz_mat, bead_types, bead_IDs, user_IDs, mol_IDs = parse_beads(xmlstring, bead_dict)
    bond_mat = parse_bonds(xmlstring)
    N = xyz_mat.shape[0]

    print("Cleaning data...")
    # Molecule IDs to start with 1
#    mol_IDs -= np.min(mol_IDs) - 1


    print("Saving data...")
    fname_beads = "nafion_ms.xyz"
#    ll.save_xyzfile(fname_beads, np.hstack((np.matrix(bead_types).T, np.matrix(mol_IDs).T, xyz_mat)))
    with open(fname_beads, "w") as f:
         f.write("# bead_types, mol_IDs, xyz_mat\n")
         for i in range(N):
             f.write("%i   %i   %i   %i   %.6f   %.6f   %.6f\n" % \
                    (bead_types[i], bead_IDs[i], user_IDs[i], mol_IDs[i],\
                     xyz_mat[i, 0], xyz_mat[i, 1], xyz_mat[i, 2]))
    print("xyz matrix was saved into", fname_beads)

    fname_bonds = "nafion_ms.bonds"
    np.savetxt(fname_bonds, bond_mat, fmt="%i")
    print("Bond matrix saved into", fname_bonds)


