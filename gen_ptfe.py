#!/usr/bin/env python
"""Usage:
    gen_ptfe.py <input> [--save <fname>]

* 5 beads in a monomer
* 15 monomers in a string 
* number of strings to fit the density from papers (dry or wet)
* finally, add water beads
Beads:
1. A: (CF2)6
2. B: O CF2 C(CF3)F O CF2
3. C: CF3 SO3H
4. W: (H2O)6

pv278@cam.ac.uk, 09/11/15
"""
import numpy as np
import os, sys, yaml
from docopt import docopt


def bead_wt(bead, arg="dry"):
    """Return atomic weigths of of beads A, B, C or W"""
    if bead == 1:          # (CF2)6
        return 6*12 + 12*19
    elif bead == 2:        # O CF2 C(CF3)F O CF2
        return 2*16 + 4*12 + 8*19
    elif bead == 3:
        if arg == "dry":   # CF3 SO3H
            return 12 + 3*19 + 32 + 3*16 + 1
        else:              # CF3 SO3H 3H2O
            return 12 + 3*19 + 32 + 3*16 + 1 + 6*1 + 3*16
    elif bead == 4:        # (H2O)6
        return 12*1 + 6*16
    else:
        print "No such bead."
        return 0


def num_poly_strings(rho, V, m=15, arg="dry"):
    """Calculate number of strings in the given volume
    and using density from papers"""
    n = 17 + 1  # DISCUSS neglecting one F (+1 term)
    tot_mass = (rho*1000) * V*(1e-10)**3 / 1.67e-27
    one_string_mass = (3*bead_wt(1) + bead_wt(2) + \
                      bead_wt(3, arg)) * m
    return tot_mass/one_string_mass


def grow_string(m, L, mol_id=1, mu=6.0, sigma=0.6):
    """Return xyz matrix of one Nafion string"""
#    names = ["B", "C", "A", "A", "A"]*m
    types = [2, 3, 1, 1, 1]*m
    types = np.matrix(types).T
    mol_ids = np.matrix([mol_id]*5*m).T
    pos = np.zeros((5*m, 3))
    pos[0] = np.random.rand(3)*L
    for i in range(1, 5*m):
        pos[i] = pos[i-1] + np.random.randn(3)*sigma + mu
    return np.hstack((mol_ids, types, pos))


def get_backbone_atoms(Ns, m, L, mu, sigma):
    """Generate xyz matrix from a given number of strings"""
    Nbs = 5*m
    mol_ids = range(1, Ns+1)
    xyz = np.zeros((Ns*Nbs, 5))
    for i in range(Ns):
        xyz[i*Nbs : (i+1)*Nbs] = grow_string(m, L, mol_ids[i], mu, sigma)
    return xyz


def get_wb_atoms(Nwb, L, count=1):
    """Generate xyz matrix from a given number of water beads.
    count -- where to start molecular id counting"""
    xyz = np.zeros((Nwb, 5))
    xyz[:, 2:5] = np.random.rand(Nwb, 3)*L
    xyz[:, 1] = 1
    xyz[:, 0] = range(count, count+Nwb)
    return xyz


def bond_mat(Ns):
    """
    Ns -- num strings
    Nm -- num monomers in a string
    """
    Nbonds = Ns * (Nm-1 + Nm*4) # between m'mers + 4 in each m'mer
    bonds = np.zeros((Nbonds, 3), dtype=int)
    bonds[:, 2] = 1  # ONE BOND TYPE NOW
    for i in range(Ns):
        pass
        # FILL
    return bonds


def get_electrodes():
    pass

# ===== printing to string
# TO DO: FIND BETTER NAMES
def get_header(N, Nbonds, L):
    """Generate LAMMPS header"""
    s = "#blabla\n"
    s += str(N) + " atoms\n"
    s += str(Nbonds) + " bonds\n"
    s += "4 atom types\n"  # MAKE THIS AUTOMATICALLY FIND OUT
    s += "1 bond types\n"  # SAME
    s += "\n"
    s += "0.0 " + str(L) + " xlo xhi\n"
    s += "0.0 " + str(L) + " ylo yhi\n"
    s += "0.0 " + str(L) + " zlo zhi\n\n"
    return s


def mass2str(masses):
    s = "Masses\n\n"
    for k, v in masses.iteritems():
        s += str(k) + " " + str(v)  + "\n"
    s += "\n"
    return s


def get_pair_coeffs(coeffs, gamma, cutoff):   # part1 part2 force gamma cutoff
    s = "PairIJ Coeffs\n\n"
    for k, v in coeffs.iteritems():
        s += "%s %s %s %s\n" % (str(k), str(v), str(gamma), str(cutoff))
    s += "\n"
    return s


def get_bond_coeffs(id_bondcoeffs):
    """Save bond coefficients into a string"""
    s = "Bond Coeffs\n\n"
    for k, v in id_bondcoeffs.iteritems():
        str_v = ""
        for i in len(v):
            str_v += str(i) + "\t"
        s += str(k) + "\t" + str_v + "\n"
    return s


def atoms2str(atom_mat):
    """Convert atomic matrix to string, atom_type molecular
    xyz_mat[:, 0] are atom ids"""
    M = len(atom_mat)
    s = ""
    for i in range(M):
        s += "%i\t%i\t%i\t%f\t%f\t%f\n" % \
             (i+1, atom_mat[i, 0], atom_mat[i, 1], atom_mat[i, 2], atom_mat[i, 3], atom_mat[i, 4])
    return s


def bonds2str(bond_mat):
    """Convert bond_matrix to string"""
    M, N = bond_mat.shape
    s = ""
    for i in range(M):
        s += str(i+1) + "\t"
        for j in range(N):
            s += str(bond_mat[i, j]) + "\t"
        s += "\n"
    return s


if __name__ == "__main__":
    args = docopt(__doc__)
    try:
        data = yaml.load(open(args["<input>"]))
    except IOError:
        print "File does not exist:", args["<input>"]
    
    # ===== general parameters
    L = data["box-size"]
    rho_dry = 1.95   # g/cm^3
    rho_wet = 1.68
    elem_wts = {"C": 12, "F": 19, "O": 16, "H": 1, "S": 32}
    masses = dict((i, bead_wt(i)) for i in range(1, 5))
    gamma = 1.0
    cutoff = 1.0

    # ===== pair parameters
    coeff2num = dict((coeff, num) for (coeff, num) in zip("ABCW", range(1, 5)))
    a_ij = dict()
    for k,v in data["int-params"].iteritems():
        a_ij[" ".join([str(coeff2num[i]) for i in k.split()])] = v

    # ===== polymer parameters
    m = 15           # num of monomers in one polymer string
    Ns = int(round(num_poly_strings(rho_wet, L**3, arg="wet")))
    Nbs = 5*m        # num beads per string
    Nbm = 5          # num beads per monomer
    Nb = Nbs*Ns      # tot num beads
    mean_bead_dist = L/Nb**(1./3)
    mu = mean_bead_dist
    sigma = mean_bead_dist/10  # 10 is arbitrary

    # ===== atoms
    poly_xyz = get_backbone_atoms(Ns, m, L, mu, sigma)
    lmbda = data["membrane"]["lambda"]
    Nwb = (Ns*m*(lmbda-3))/6   # number of water beads

    wb_xyz = get_wb_atoms(Nwb, L, count=Ns+1)
    xyz = np.vstack((poly_xyz, wb_xyz))
    xyz_str = atoms2str(xyz)

    # ===== bonds
    # FILL

    # ===== putting it together
    final_string = get_header(len(xyz), 0, L) + \
                   mass2str(masses) + \
                   get_pair_coeffs(a_ij, gamma, cutoff) + \
                   "Atoms\n\n" + xyz_str #+ \
#                   "Bonds\n\n" + final_bonds_str

    if args["--save"]:
        fname = args["<fname>"]
        open(fname, "w").write(final_string)
        print "Data file written in", fname
    else:
        print final_string



