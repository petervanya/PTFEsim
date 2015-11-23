#!/usr/bin/env python
"""Usage:
    gen_ptfe.py <input> [--save <fname>]

* 5 beads in a monomer
* 15 monomers in a chain
* number of chains to fit the density from papers (dry or wet)
* finally, add water beads

pv278@cam.ac.uk, 09/11/15
"""
import numpy as np
import os, sys, yaml
from docopt import docopt
from parse_topo import *


def nafion_bead_wt(bead, arg="dry"):
    """Return weigths in atomic units of beads A, B, C or W"""
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


def num_poly_chains(rho, V, Nm=15, arg="dry"):
    """Calculate number of chains in the given volume
    and using density from papers, using SI units
    TODO: customise the numbers of beads, not just for Nafion"""
    tot_mass = rho * V  
    one_chain_mass = (3*nafion_bead_wt(1) + nafion_bead_wt(2) + \
                     nafion_bead_wt(3, arg)) * Nm * 1.67e-27
    return tot_mass/one_chain_mass


def grow_chain(Nm, L, mol_id=1, mu=6.0, sigma=0.6):
    """Return xyz matrix of one Nafion chain"""
#    names = ["A", "A", "A", "B", "C"]*m
    types = [1, 1, 1, 2, 3]*Nm
    types = np.matrix(types).T
    mol_ids = np.matrix([mol_id]*5*Nm).T
    pos = np.zeros((5*Nm, 3))    # MAKE THIS MORE GENERAL
    pos[0] = np.random.rand(3)*L
    for i in range(1, 5*Nm):
        pos[i] = pos[i-1] + np.random.randn(3)*sigma + mu
    return np.hstack((mol_ids, types, pos))


def get_polymer_atoms(Ns, Nm, L, mu, sigma):
    """Generate xyz matrix from a given number of chains"""
    Nbs = 5*Nm
    mol_ids = range(1, Ns+1)
    xyz = np.zeros((Ns*Nbs, 5))
    for i in range(Ns):
        xyz[i*Nbs : (i+1)*Nbs] = grow_chain(Nm, L, mol_ids[i], mu, sigma)
    return xyz


def get_wb_atoms(Nwb, L, count=1):
    """Generate xyz matrix from a given number of water beads.
    count -- where to start molecular id counting"""
    xyz = np.zeros((Nwb, 5))
    xyz[:, 2:5] = np.random.rand(Nwb, 3)*L
    xyz[:, 1] = 1
    xyz[:, 0] = range(count, count+Nwb)
    return xyz


def bonds_mat(Na, Nm, Nbm=5):
    """
    Ns -- num chains
    Nm -- num monomers in a chain
    Nbm -- num beads in monomer
    Types:
    * 1: AA bond
    * 2: AB bond
    * 3: AC bond
    Order:
    count, id, part1, part2
    """
    Nbonds = Ns * (Nm-1 + Nm*4)        # between m'mers + 4 in each m'mer
    bonds = np.zeros((Nbonds, 3), dtype=int)
    cnt = 0
    for i in range(Ns):
        for j in range(Nm):     # within a monomer
            first = Ns*i+j+1    # first bead in a monomer
            bonds[cnt]   = [1, first,   first+1]
            bonds[cnt+1] = [1, first+1, first+2]
            bonds[cnt+2] = [2, first,   first+3]
            bonds[cnt+3] = [3, first+3, first+4]
            cnt += 4
    for i in range(Ns):
        firstm = Ns*i+1         # first monomer in str
        for j in range(Nm-1):   # between monomers
            bonds[cnt] = [1, firstm+j*Nbm+2, firstm+j*Nbm+5]
            cnt += 1
    return bonds


def get_electrodes():
    pass

# =====
# ===== printing to str
# =====
# TO DO: FIND BETTER NAMES
def header2str(N, Nbonds, atomtypes, bondtypes, L):
    """Generate LAMMPS header"""
    s = "#blabla\n"
    s += str(N) + " atoms\n"
    s += str(Nbonds) + " bonds\n"
    s += str(atomtypes) + " atom types\n"
    s += str(bondtypes) + " bond types\n"
    s += "\n"
    s += "0.0 " + str(L) + " xlo xhi\n"
    s += "0.0 " + str(L) + " ylo yhi\n"
    s += "0.0 " + str(L) + " zlo zhi\n\n"
    return s


def mass2str(masses):
    s = "Masses\n\n"
    for k, v in masses.iteritems():
        s += str(k) + " " + str(v)  + "\n"
    return s + "\n"


def paircoeffs2str(coeffs, gamma, cutoff):   # part1 part2 force gamma cutoff
    s = "PairIJ Coeffs\n\n"
    for k, v in coeffs.iteritems():
        s += "%s %s %s %s\n" % (str(k), str(v), str(gamma), str(cutoff))
    return s + "\n"


def bondcoeffs2str(id_bondcoeffs):
    """Save bond coefficients into a str"""
    s = "Bond Coeffs\n\n"
    for k, v in id_bondcoeffs.iteritems():
        s += str(k) + "\t" + str(v[0]) + "\t" + str(v[1]) + "\n"
    return s + "\n"


def atoms2str(atom_mat):
    """Convert atomic matrix to str, atom_type molecular
    xyz_mat[:, 0] are atom ids"""
    M = len(atom_mat)
    s = ""
    for i in range(M):
        s += "%i\t%i\t%i\t%e\t%e\t%e\n" % \
             (i+1, atom_mat[i, 0], atom_mat[i, 1], atom_mat[i, 2], atom_mat[i, 3], atom_mat[i, 4])
    return s + "\n"


def bonds2str(bond_mat):
    """Convert bond_matrix to str"""
    M, N = bond_mat.shape
    s = ""
    for i in range(M):
        s += str(i+1) + "\t"
        for j in range(N):
            s += str(bond_mat[i, j]) + "\t"
        s += "\n"
    return s + "\n"


if __name__ == "__main__":
    args = docopt(__doc__)
    try:
        data = yaml.load(open(args["<input>"]))
    except IOError:
        print "File does not exist:", args["<input>"]
    
    # ===== general parameters
    elem_wts = {"C": 12, "F": 19, "O": 16, "H": 1, "S": 32}
    kB = 1.38e-23
    Maw = 1.67e-27
    d_DPD = 8.14e-10        # DPD distance unit
    T = data["temperature"]
    L = data["box-size"]*d_DPD
    rho_dry = 1950          # kg/m^3
    rho_wet = 1680
    masses = dict( (i, nafion_bead_wt(i)*Maw) for i in range(1, 5) )  # SI units
    gamma = 1.0 * d_DPD**3  # in SI units, CHECK
    rc = d_DPD              # cutoff distance

    # ===== pair and bond parameters
    atomtypes = 4
    bondtypes = 3
    coeff2num = dict((coeff, num) for (coeff, num) in zip("ABCW", range(1, 5)))
    a_ij = {}
    for k, v in data["int-params"].iteritems():
        a_ij[" ".join([str(coeff2num[i]) for i in k.split()])] \
            = v * kB*T/rc
    for i in range(1, 5):
        a_ij["%i %i" % (i, i)] = 25 * kB*T/rc

    k_ij = {}
    r0 = 0.0
    for k, v in data["spring-const"].iteritems():
        k_ij[" ".join([str(coeff2num[i]) for i in k.split()])] \
            = ["%e" % (v * kB*T/rc**2), r0]
#    for i in range(1, bondtypes+1):
#        k_ij[i] = [1.0 * kB*T/rc**2, 0.0]   # spring const, eq. dist

    # ===== polymer parameters
    Nm = 15                 # num of monomers in one polymer chain
    Ns = int(round(num_poly_chains(rho_wet, L**3, Nm, arg="wet")))
    Nbs = 5*Nm              # num beads per chain
    Nbm = 5                 # num beads per monomer
    Nb = Nbs*Ns             # tot num beads
    mean_bead_dist = L/Nb**(1./3)
    mu = mean_bead_dist
    sigma = mean_bead_dist/10  # 10 is arbitrary
    print Ns, "polymer chains in a given volume"

    # ===== atoms
    poly_xyz = get_polymer_atoms(Ns, Nm, L, mu, sigma)
    lmbda = data["membrane"]["lambda"]
    Nwb = (Ns*Nm*(lmbda-3))/6   # number of water beads
    wb_xyz = get_wb_atoms(Nwb, L, count=Ns+1)

    xyz = np.vstack((poly_xyz, wb_xyz))
    xyz_str = atoms2str(xyz)
    print len(xyz), "beads created"

    # ===== bonds
    bonds = bonds_mat(Ns, Nm, Nbm)
    bonds_str = bonds2str(bonds)
    print len(bonds), "bonds created"

    # ===== putting it together
    final_string = header2str(len(xyz), len(bonds), atomtypes, bondtypes, L) + \
                   mass2str(masses) + \
                   paircoeffs2str(a_ij, gamma, rc) + \
                   bondcoeffs2str(k_ij) + \
                   "Atoms\n\n" + xyz_str + \
                   "Bonds\n\n" + bonds_str

    if args["--save"]:
        fname = args["<fname>"]
        open(fname, "w").write(final_string)
        print "Data file written in", fname
    else:
        print final_string



