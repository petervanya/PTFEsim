#!/usr/bin/env python
"""Usage:
    gen_ptfe.py <input> [--save <fname>] [--parsetopo]

Arguments:
    <input>                  Input yaml file

Options:
    --save <fname>           Save the data file
    --parsetopo              Use if parameters should be obtained from topo string

* 5 beads in a monomer
* 15 monomers in a chain
* number of chains to fit the density from the papers (dry or wet)
* finally, add water beads

pv278@cam.ac.uk, 09/11/15
"""
import numpy as np
import os, sys, yaml
from docopt import docopt
import parse_topo as pt

kB = 1.38e-23
Maw = 1.67e-27
d_DPD = 8.14e-10             # DPD distance unit
elem_wts = {"C": 12, "F": 19, "O": 16, "H": 1, "S": 32}
rho_dry = 1950               # kg/m^3
rho_wet = 1680
rho_DPD = 3

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


def get_polymer_atoms(Nc, Nm, L, mu, sigma):
    """Generate xyz matrix from a given number of chains"""
    Nbc = 5*Nm
    mol_ids = range(1, Nc+1)
    xyz = np.zeros((Nc*Nbc, 5))
    for i in range(Nc):
        xyz[i*Nbc : (i+1)*Nbc] = grow_chain(Nm, L, mol_ids[i], mu, sigma)
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
    Nc -- num chains
    Nm -- num monomers in a chain
    Nbm -- num beads in monomer
    Types:
    * 1: AA bond
    * 2: AB bond
    * 3: AC bond
    Order:
    count, id, part1, part2
    """
    Nbonds = Nc * (Nm-1 + Nm*4)        # between m'mers + 4 in each m'mer
    bonds = np.zeros((Nbonds, 3), dtype=int)
    cnt = 0
    for i in range(Nc):
        for j in range(Nm):            # within a monomer
            first = Nc*i+j+1           # first bead in a monomer
            bonds[cnt]   = [1, first,   first+1]
            bonds[cnt+1] = [1, first+1, first+2]
            bonds[cnt+2] = [2, first,   first+3]
            bonds[cnt+3] = [3, first+3, first+4]
            cnt += 4
    for i in range(Nc):
        firstm = Nc*i + 1              # first monomer in str
        for j in range(Nm-1):          # between monomers
            bonds[cnt] = [1, firstm+j*Nbm+2, firstm+j*Nbm+5]
            cnt += 1
    return bonds


def bonds_mat2(raw_topo, Nc):
    """Create bond matrix from topology string and number of chains"""
    bead_list, Nm = pt.parse_beads(raw_topo)
    Nbm = len(bead_list)
    bonds = pt.construct_bonds(bead_list, Nm, 0)
    for i in range(1, Nc):
        bonds = np.vstack( (bonds, pt.construct_bonds(bead_list, Nm, i*Nm*Nbm)) ) 
    return bonds


def gen_pair_coeffs(atom_types, atoms_yaml, gamma, rc):
    """
    Generate atomic params a_ij for all possible combinations 
    given number of atom types. Read custom bonds from input.yaml file
    Nafion atom_types = "ABCW"
    """
    a_ij = {}
    Nat = len(atom_types)
    num2coeff = dict((num, coeff) for (num, coeff) in zip(range(1, Nat+1), atom_types))
    for i in range(1, Nat+1):
        for j in range(1, i+1):
            key = "%i %i" % (j, i)
            lkey = "%s %s" % (num2coeff[j], num2coeff[i])
            if lkey in atoms_yaml.keys() or lkey[::-1] in atoms_yaml.keys():
                try:
                    a_ij[key] = [atoms_yaml[lkey] * kB*T/rc, gamma, rc]
#                    a_ij[key] = [atoms_yaml[lkey], 1.0, 1.0]
                except KeyError: # "B A" -> "A B"
                    a_ij[key] = [atoms_yaml[lkey[::-1]] * kB*T/rc, gamma, rc]
#                    a_ij[key] = [atoms_yaml[lkey[::-1]], 1.0, 1.0]
            else:
                a_ij[key] = [25 * kB*T/rc, gamma, rc]
#                a_ij[key] = [25.0, 1.0, 1.0]
    return a_ij


def gen_bond_coeffs(atom_types, bonds_yaml, r0):
    """
    Generate bond coeffs k_ij and r0 for all possible combinations
    given number of atom types. Read custom bonds from input.yaml file
    Nafion atom types = "ABCW"
    """
    k_ij = {}
    bmap = pt.bond_map(atom_types)   # "ABCW"
    Nat = len(atom_types)
    for i in range(Nat):
        for j in range(i+1): 
            key = atom_types[i] + " " + atom_types[j]
            if key in bonds_yaml.keys():
                k_ij[bmap[key]] = [bonds_yaml[key] * kB*T/(rc**2), r0]
#                k_ij[bmap[key]] = [bonds_yaml[key], r0]
            else:
                k_ij[bmap[key]] = [4 * kB*T/(rc**2), r0]
#                k_ij[bmap[key]] = [4, r0]
    return k_ij


def get_electrodes():
    pass


# =====
# ===== functions to convert stuff to str
# =====
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


def paircoeffs2str(coeffs):   # part1 part2 force gamma cutoff
    s = "PairIJ Coeffs\n\n"
    for k, v in coeffs.iteritems():
        s += "%s %s %s %s\n" % (str(k), str(v[0]), str(v[1]), str(v[2]))
    return s + "\n"


def bondcoeffs2str(k_ij):
    """Save bond coefficients into a str"""
    s = "Bond Coeffs\n\n"
    for k, v in k_ij.iteritems():
        s += "%s %s %s\n" % (str(k), "%e" % v[0], "%e" % v[1])
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
    """Convert bond matrix to str"""
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
    T = data["temperature"]
    L = data["box-size"]*d_DPD
    atom_types = "ABCW"      # MAKE THIS GENERAL
    Nat = len(atom_types)    # num atom types
    masses = dict( (i, nafion_bead_wt(i)*Maw) for i in range(1, 5) )  # SI units
    gamma = 1.0 * d_DPD**3   # in SI units, CHECK
    rc = d_DPD               # cutoff distance
    r0 = data["equilibrium-dist"] * d_DPD

    # ===== polymer parameters, MAKE THIS PARSE TOPOLOGY
    if args["--parsetopo"]:
        raw_topo = data["topology"]
        bead_list, Nm = pt.parse_beads(raw_topo)
        Nbm = len(bead_list)
        Nbc = Nbm*Nm
    else:
        Nm = 15                 # num monomers in one polymer chain
        Nbm = 5                 # num beads per monomer
        Nbc = Nbm*Nm            # num beads per chain
    print "Nm: %i, Nbm: %i" % (Nm, Nbm)
    Nc = int(round(num_poly_chains(rho_wet, L**3, Nm, arg="wet")))
    Nb = Nbc*Nc                 # tot num beads
    mean_bead_dist = L/Nb**(1./3)
    mu = mean_bead_dist
    sigma = mean_bead_dist/10   # 10 is arbitrary
    print Nc, "polymer chains in a given volume"

    # ===== pair and bond parameters
    a_ij = gen_pair_coeffs("ABCW", data["pair-coeffs"], gamma, rc)
    k_ij = gen_bond_coeffs("ABCW", data["bond-coeffs"], r0)
    Nbt = len(k_ij)                  # num bond types

#    coeff2num = dict((coeff, num) for (coeff, num) in zip(atom_types, range(1, Nat+1)))
#    a_ij = {}
#    for k, v in data["pair-coeffs"].iteritems():
#        a_ij[" ".join([str(coeff2num[i]) for i in k.split()])] \
#            = [v * kB*T/rc, gamma, rc]
#    for i in range(1, Nat+1):   # same bonds
#        a_ij["%i %i" % (i, i)] = [25 * kB*T/rc, gamma, rc]
#
#    k_ij = {}
#    r0 = 0.0
#    for k, v in data["bond-coeffs"].iteritems():
#        k_ij[" ".join([str(coeff2num[i]) for i in k.split()])] \
#            = ["%e" % (v * kB*T/rc**2), r0]
#    for i in range(1, Nat+1):   # same bonds
#        k_ij["%i %i" % (i, i)] = ["%e" % (4 * kB*T/rc**2), r0]

    # ===== atoms
    poly_xyz = get_polymer_atoms(Nc, Nm, L, mu, sigma)
    lmbda = data["water-uptake"]
    Nwb = (Nc*Nm*(lmbda-3))/6   # number of water beads
    wb_xyz = get_wb_atoms(Nwb, L, count=Nc+1)

    xyz = np.vstack((poly_xyz, wb_xyz))
    xyz_str = atoms2str(xyz)
    print len(xyz), "beads created"

    # ===== bonds
    bonds = bonds_mat2(data["topology"], Nc)  # extracting Nm, Nbm
    bonds_str = bonds2str(bonds)
    print len(bonds), "bonds created"

    # ===== putting it together
    final_string = header2str(len(xyz), len(bonds), Nat, Nbt, L) + \
                   mass2str(masses) + \
                   paircoeffs2str(a_ij) + \
                   bondcoeffs2str(k_ij) + \
                   "Atoms\n\n" + xyz_str + \
                   "Bonds\n\n" + bonds_str

    if args["--save"]:
        fname = args["--save"]
        open(fname, "w").write(final_string)
        print "Data file written in", fname
    else:
        print final_string


