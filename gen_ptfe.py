#!/usr/bin/env python
"""Usage:
    gen_ptfe.py <input> [--save <fname>] [--parsetopo] [--expt-rho]

Arguments:
    <input>                  Input yaml file

Options:
    --save <fname>           Save the data file
    --parsetopo              Use if parameters should be obtained from topo string
    --expt-rho               Use experimental density value instead of DPD (rho=3)

* Nbm beads in a monomer
* Nmc monomers in a chain
* number of chains (Nc) to fit the density from the papers (dry or wet) or DPD rho=3
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


def num_poly_chains(rho, V, Nmc=15, arg="dry"):
    """Calculate number of chains in the given volume
    and using density from papers, using SI units
    TODO: customise the numbers of beads, not just for Nafion"""
    tot_mass = rho * V  
    one_chain_mass = (3*nafion_bead_wt(1) + nafion_bead_wt(2) + \
                     nafion_bead_wt(3, arg)) * Nmc * 1.67e-27
    return tot_mass/one_chain_mass


def calc_nc_nw(N, Nmc, Nbm, lmbda):
    """Return number of water beads and chains based on 
    universal DPD density, number of beads and water uptake"""
    const = (lmbda-3)/6.
    Nc = N/(Nmc*(const + Nbm))
    Nw = const*Nmc*Nc
    return int(Nc), int(Nw)


def grow_one_chain(beads, Nmc, L, mol_id=1, mu=6.0, sigma=0.6):
    """Return xyz matrix of one Nafion chain"""
    Nbm = len(beads)
    types = beads*Nmc
    types = np.matrix(types).T
    mol_ids = np.matrix([mol_id]*Nbm*Nmc).T
    pos = np.zeros((Nbm*Nmc, 3))
    pos[0] = np.random.rand(3)*L
    for i in range(1, len(types)):
        pos[i] = pos[i-1] + np.random.randn(3)*sigma + mu
    return np.hstack((mol_ids, types, pos))


def get_polymer_atoms(beads, Nc, Nmc, L, mu, sigma):
    """Generate xyz matrix from a given number of chains"""
    Nbm = len(beads)
    Nbc = Nbm*Nmc
    mol_ids = range(1, Nc+1)
    xyz = np.zeros((Nc*Nbc, Nbm))
    for i in range(Nc):
        xyz[i*Nbc : (i+1)*Nbc] = grow_one_chain(beads, Nmc, L, mol_ids[i], mu, sigma)
    return xyz


def gen_wb_atoms(Nw, L, count=1):
    """Generate xyz matrix from a given number of water beads.
    count -- where to start molecular id counting"""
    xyz = np.zeros((Nw, 5))
    xyz[:, 2:5] = np.random.rand(Nw, 3)*L
    xyz[:, 1] = 1
    xyz[:, 0] = range(count, count+Nw)
    return xyz


def bonds_mat(Na, Nmc, Nbm=5):
    """
    NOT USED NOW
    Nc -- num chains
    Nmc -- num monomers in a chain
    Nbm -- num beads in monomer
    Types:
    * 1: AA bond
    * 2: AB bond
    * 3: AC bond
    Order:
    count, id, part1, part2
    """
    Nbonds = Nc * (Nmc-1 + Nmc*4)        # between m'mers + 4 in each m'mer
    bonds = np.zeros((Nbonds, 3), dtype=int)
    cnt = 0
    for i in range(Nc):
        for j in range(Nmc):            # within a monomer
            first = Nc*i+j+1           # first bead in a monomer
            bonds[cnt]   = [1, first,   first+1]
            bonds[cnt+1] = [1, first+1, first+2]
            bonds[cnt+2] = [2, first,   first+3]
            bonds[cnt+3] = [3, first+3, first+4]
            cnt += 4
    for i in range(Nc):
        firstm = Nc*i + 1              # first monomer in str
        for j in range(Nmc-1):          # between monomers
            bonds[cnt] = [1, firstm+j*Nbm+2, firstm+j*Nbm+5]
            cnt += 1
    return bonds


def bonds_mat2(raw_topo, Nc):
    """Create bond matrix from topology string and number of chains"""
    bead_list, Nmc = pt.parse_beads(raw_topo)
    Nbm = len(bead_list)
    bonds = pt.construct_bonds(bead_list, Nmc, 0)
    for i in range(1, Nc):
        bonds = np.vstack( (bonds, pt.construct_bonds(bead_list, Nmc, i*Nmc*Nbm)) ) 
    return bonds


def gen_pair_coeffs(bead_types, atoms_yaml, gamma, rc):
    """
    Generate atomic params a_ij for all possible combinations 
    given number of atom types. Read custom bonds from input.yaml file
    Nafion bead_types = "ABCW"
    """
    a_ij = {}
    Nbt = len(bead_types)
    num2coeff = dict((num, coeff) for (num, coeff) in zip(range(1, Nbt+1), bead_types))
    for i in range(1, Nbt+1):
        for j in range(1, i+1):
            key = "%i %i" % (j, i)
            lkey = "%s %s" % (num2coeff[j], num2coeff[i])
            if lkey in atoms_yaml.keys() or lkey[::-1] in atoms_yaml.keys():
                try:
                    a_ij[key] = [atoms_yaml[lkey] * kB*T/rc, gamma, rc]
                except KeyError: # "B A" -> "A B"
                    a_ij[key] = [atoms_yaml[lkey[::-1]] * kB*T/rc, gamma, rc]
            else:
                a_ij[key] = [25 * kB*T/rc, gamma, rc]
    return a_ij


def gen_bond_coeffs(bead_types, bonds_yaml, r0):
    """
    Generate bond coeffs k_ij and r0 for all possible combinations
    given number of atom types. Read custom bonds from input.yaml file
    Nafion atom types = "ABCW"
    """
    k_ij = {}
    bmap = pt.bond_map(bead_types)   # "ABCW"
    Nbt = len(bead_types)
    for i in range(Nbt):
        for j in range(i+1): 
            key = bead_types[i] + " " + bead_types[j]
            if key in bonds_yaml.keys():
                k_ij[bmap[key]] = [bonds_yaml[key] * kB*T/(rc**2), r0]
            else:
                k_ij[bmap[key]] = [4 * kB*T/(rc**2), r0]
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
    
    # ===== set general parameters
    T = data["temperature"]
    L = data["box-size"]*d_DPD
    rc = d_DPD                         # cutoff distance
    r0 = data["equilibrium-dist"] * d_DPD
    gamma = data["gamma"] * d_DPD**3   # SI units
    lmbda = data["water-uptake"]

    # ===== set polymer parameters
    if args["--parsetopo"]:
        raw_topo = data["topology"]
        print "Topology:", raw_topo
        bead_list, Nmc = pt.parse_beads(raw_topo)   # generate list of information about connectivity
        Nbm = len(bead_list)
        bead_dict = pt.gen_bead_dict(raw_topo)      # generate dict of beads in one monomer
        bead_types = sorted("".join(bead_dict.keys()) + "W")
        Nbt = len(bead_types)        # num atom types
        coeff2num = dict((coeff, num) for num, coeff in zip(range(1, Nbt+1), bead_types))
        beads = []
        [[beads.append(coeff2num[k]) for i in range(v)] for k, v in sorted(bead_dict.items())]
    else:
        Nmc = 15                 # num monomers in one polymer chain
        Nbm = 5                  # num beads per monomer
        bead_types = "ABCW"
        beads = [1, 1, 1, 2, 3]
        Nbt = len(bead_types)    # num atom types
    print "Nmc: %i, Nbm: %i" % (Nmc, Nbm)

    if args["--expt-rho"]:
        Nc = int(round(num_poly_chains(rho_wet, L**3, Nmc, arg="wet")))
        Nw = (Nc*Nmc*(lmbda-3))/6   # number of water beads
        N = Nc*Nmc*Nbm + Nw
    else:
        N = rho_DPD*int(L/d_DPD)**3
        Nc, Nw = calc_nc_nw(N, Nmc, Nbm, lmbda)
    print Nc, "polymer chains in a given volume"

    Nbc = Nbm*Nmc               # num beads per chain
    Nb = Nbc*Nc                 # tot num beads
    mean_bead_dist = L/Nb**(1./3)
    mu = mean_bead_dist
    sigma = mean_bead_dist/10   # 10 is arbitrary
    masses = dict( (i, nafion_bead_wt(i)*Maw) for i in range(1, Nbt) )  # SI units

    # ===== pair and bond parameters
    a_ij = gen_pair_coeffs(bead_types, data["pair-coeffs"], gamma, rc)
    k_ij = gen_bond_coeffs(bead_types, data["bond-coeffs"], r0)
    Nbondt = len(k_ij)          # num bond types

    # ===== atoms
    poly_xyz = get_polymer_atoms(beads, Nc, Nmc, L, mu, sigma)
    wb_xyz = gen_wb_atoms(Nw, L, count=Nc+1)
    xyz = np.vstack((poly_xyz, wb_xyz))
    xyz_str = atoms2str(xyz)
    print len(xyz), "beads created"

    # ===== bonds
    bonds = bonds_mat2(data["topology"], Nc)  # extracting Nmc, Nbm
    bonds_str = bonds2str(bonds)
    print len(bonds), "bonds created"

    # ===== putting it together
    final_string = header2str(len(xyz), len(bonds), Nbt, Nbondt, L) + \
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


