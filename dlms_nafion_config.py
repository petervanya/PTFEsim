#!/usr/bin/env python
"""Usage:
    dlms_nafion_config.py <input> [--el <el> --xyz <xyz>]

Generate DL_MESO input CONFIG and FIELDS files.
Electrode support: carbon, quartz or (amorphous) silica.

Options:
    --el <el>      Electrode support
    --xyz <xyz>    Save coordinates into xyz file <xyz>

pv278@cam.ac.uk, 03/06/16
"""
import numpy as np
from numpy import pi, cos, sin
import sys
import yaml
from docopt import docopt
import lmp_lib as ll

NA = 6.022e23
AMU = 1.66e-27
rho_DPD = 3.0
a_DPD = 25.0
rc = 8.14e-10               # DPD distance unit
k0 = 4.0
rho_el = {"carbon": 2000.0, "quartz": 2648.0, "silica": 2196.0}
rho_Pt = 21450.0
elem_wts = {"C": 12.01, "O": 15.99, "Si": 28.08}


def set_electrodes(data, L):
    """From yaml file extract:
    * Lcl: electrode layer width
    * Lpt: Pt layer width on electrodes
    * Nsup: number of support beads
    * Npt: number of platinum beads
    * Scheme:
        *-----------*
        |C |     |C | L-Lpt |
        |--|PTFE |--|        } L
        |Pt|     |Pt| Lpt   |
        *-----------*
         Lcl
    """
    Lcl = data["electrodes"]["width"]  # electrode layer width
    if Lcl > L/2:
        print("Electrode layer thicker than membrane, aborting.")
        sys.exit()

    if Lcl > 1e-10:
        Pt_ratio = data["electrodes"]["Pt-ratio"]
        elmat = data["electrodes"]["material"]
        Lpt = L * Pt_ratio
        Vcl = 2 * Lcl * L**2
        print("Electrodes on | Width: %.2f | Material: %s" % (Lcl, elmat))

        Nelb = int(rho_DPD * int(Vcl))     # not (4*pi/3*rc**3) ), must be cubes
        Nsup = int((1-Pt_ratio) * Nelb)
        Npt = int(Pt_ratio * Nelb)
        if Nsup%2 != 0: Nsup += 1
        if Npt%2 != 0: Npt += 1
        print("Beads: Electrode: %i | Support: %i | Pt: %i" % (Nelb, Nsup, Npt))
    else:
        Lcl, Lpt = 0.0, 0.0
        Nsup, Npt = 0, 0
        print("Electrodes off.")
    return Lcl, Lpt, Nsup, Npt


def calc_nc_nw(N, Nmc, Nbm, lmbda):
    """Return number of water beads and chains based on 
    universal DPD density, number of beads and water uptake"""
    const = (lmbda-3)/6.0
    Nc = N/(Nmc*(const + Nbm))
    Nw = const*Nmc*Nc
    return int(Nc), int(Nw)


def grow_polymer(Nbm, Nc, Nmc, L, Lcl, mu):
    """Generate xyz matrix from a given number of chains,
    by stacking one bead next to another in distance mu"""
    Nbc = Nbm*Nmc       # num beads per chain
    xyz = np.zeros((Nc*Nbc, 3))
    for i in range(Nc):
        xyz[i*Nbc : (i+1)*Nbc] = grow_one_chain(Nbm, Nmc, L, Lcl, mu)
    return xyz


def grow_one_chain(Nbm, Nmc, L, Lcl, mu=1.0):
    """Return xyz matrix of one polymer chain.
    Return (Nbc*len(beads_in_m), 3) xyz matrix
    """
    N = Nbm*Nmc   # beads per chain
    xyz = np.zeros((N, 3))
    xyz[0] = np.random.rand(3)*L
    for i in range(1, N):
        theta = np.random.rand()*pi
        phi = np.random.rand()*2*pi
        r = mu
        new_bead_pos = [r*cos(theta), r*sin(theta)*cos(phi), r*sin(theta)*sin(phi)]
        xyz[i] = xyz[i-1] + new_bead_pos
    xyz[:, 0] = xyz[:, 0]*(L-2*Lcl)/L + Lcl        # fit into proper volume
    return xyz


def gen_water_beads(Nw, L, Lcl):
    """Generate xyz matrix from a given number of water beads.
    Return (Nw, 3) xyz matrix
    Nw: number of water beads"""
    if Nw == 0:
        return np.empty((0, 3))
    xyz = np.zeros((Nw, 3))
    xyz = np.random.rand(Nw, 3)
    xyz[:, 0] *= L - 2*Lcl
    xyz[:, 0] += Lcl
    xyz[:, 1:3] *= L
    return xyz


def gen_el_support(Nsup, L, Lcl, Lpt):
    """Generate electrode support on both sides of the membrane,
    squeeze into Lsup in z-dir and shift by Lpt in z-dir.
    Return (Nsup, 3) xyz matrix
    Nsup: number of electrode beads"""
    Lsup = L - Lpt
    if Nsup == 0:
        return np.empty((0, 3), dtype=float)
    xyz1 = np.random.rand(Nsup/2, 3)
    xyz1[:, 0] *= Lcl
    xyz1[:, 1] *= Lsup
    xyz1[:, 1] += Lpt
    xyz1[:, 2] *= L
    xyz2 = np.random.rand(Nsup/2, 3)
    xyz2[:, 0] *= Lcl
    xyz2[:, 0] += L - Lcl
    xyz2[:, 1] *= Lsup
    xyz2[:, 1] += Lpt
    xyz2[:, 2] *= L
    return np.vstack((xyz1, xyz2))


def gen_platinum(Npt, L, Lcl, Lpt):
    """Randomly generate Pt beads in catalyst layer
    Return (Npt, 3) xyz matrix
    Npt: number of platinum beads"""
    if Npt == 0:
        return np.empty((0, 3), dtype=float)
    xyz1 = np.zeros((Npt/2, 3))
    xyz1 = np.random.rand(Npt/2, 3)
    xyz1[:, 0] *= Lcl
    xyz1[:, 1] *= Lpt
    xyz1[:, 2] *= L
    xyz2 = np.zeros((Npt - Npt/2, 3))
    xyz2 = np.random.rand(Npt - Npt/2, 3)
    xyz2[:, 0] *= Lcl
    xyz2[:, 0] += L - Lcl
    xyz2[:, 1] *= Lpt
    xyz2[:, 2] *= L
    xyz = np.vstack((xyz1, xyz2))
    return xyz


def gen_bonds(Nmc, bead_types):
    """
    Generate bonds for one chain.
    * bead_types: e.g. "AAABC", length Nbm (Number of beads per monomer)
    * Nmc: num monomers in a chain
    Return (Nmc * Nbm, 2) matrix: [bead1, bead2]
    """
    Nbm = len(bead_types)  # number of beads per monomer
    mono_bond_block = np.array([[1, -2],\
                                [2, 1],\
                                [3, 2],\
                                [4, 3],\
                                [5, 4]], dtype=int)
    bond_mat = np.zeros((Nbm*Nmc, 2))
    for i in range(Nmc):
        bond_mat[Nbm*i : Nbm*(i+1)] = mono_bond_block + Nbm
    return bond_mat[1:]    # reject 1st dummy bond


def gen_inter_coeffs(atoms_yaml, bead_types, gamma=4.5, rc=1.0):
    """
    Generate atomic params a_ij for all possible combinations 
    given number of atom types. Read custom bonds from input.yaml file
    * bead_types, e.g. "ABCW"
    * 3.27 = coeff from Groot-Warren, JCP, 1997 for rho=3
    """
    Nbt = len(bead_types)
    a_ij = {}
    for i in range(Nbt):
        for j in range(i+1):
            pair = "%s %s" % (bead_types[j], bead_types[i])
            if pair in atoms_yaml.keys():
                a_ij[pair] = [a_DPD + 3.27*atoms_yaml[pair], rc, gamma]
            elif pair[::-1] in atoms_yaml.keys():
                a_ij[pair] = [a_DPD + 3.27*atoms_yaml[pair[::-1]], rc, gamma]
            else:
                a_ij[pair] = [a_DPD, rc, gamma]
    return a_ij


def species2str(bead_types, bead_pop):
    """
    * bead_types: e.g. "ABCWEP"
    * bead_pop: dict of bead population
    """
    s = "SPECIES %i\n" % len(bead_types)
    for b in bead_types:
        s += b + "    1.0 0.0 " + str(bead_pop[b])
        if b == "E" or b == "P":
            s += " 1\n"
        else:
            s += " 0\n"
    return s + "\n"


def inter2str(a_ij, method="dpd"):
    s = "INTERACTIONS %i\n" % len(a_ij)
    for k, v in a_ij.items():
        s += k + "    %s   %.3f  %.1f  %.2f\n" % (method, v[0], v[1], v[2])
    return s + "\n"


def mol2str(molname, Nmols, bead_list, bond_mat, bond_type="harm", k0=4.0, r0=0.1):
    """Input:
    * molecule name
    * Nmols: number of molecules of this type
    * bead_list: list of bead types in one molecule, each one char
    * bond_mat: (Nbonds, 2) matrix of connected beads
    * bond_type
    * k0: spring constant
    * r0: equilibrium distance
    """
    s = molname + "\n"
    s += "nummols %s \n" % str(Nmols)
    s += "beads %i\n" % len(bead_list)
    for n in bead_list:
        s += n + "\n"
    s += "bonds %i\n" % len(bond_mat)
    for i in range(len(bond_mat)):
        s += "%s  %.2i %.2i %.3f %.3f\n" % \
             (bond_type, bond_mat[i, 0], bond_mat[i, 1], k0, r0)
    s += "finish\n"
    return s


def save_config(fname, names, xyz, imcon=0):
    """Save positions into file
    * imcom: include box coords (0 or 1)"""
    N = len(xyz)
    conf_str = "pokus\n" + "0\t%i\n" % imcon
    if imcon == 1:
        box = L*np.eye(3)
        for i in range(len(box)):
            conf_str += "%f\t%f\t%f\n" % (box[i, 0], box[i, 1], box[i, 2])

    for i in range(N):
        conf_str += "%s        %i\n" % (names[i], i+1)
        # careful about the spaces
        conf_str += "    %.10f    %.10f    %.10f\n" % (xyz[i, 0], xyz[i, 1], xyz[i, 2])

    open(fname, "w").write(conf_str)
    print("Initial configuration saved in %s" % fname)


if __name__ == "__main__":
    args = docopt(__doc__)
    try:
        data = yaml.load(open(args["<input>"]))
    except IOError:
        print("Input file not found:", args["<input>"])
        sys.exit()
    np.random.seed(1234)
    gamma = data["gamma"] 
    r0 = data["equilibrium-dist"]
    lmbda = data["water-uptake"]
    L = data["box-size"]
    if lmbda < 3:
        print("Water uptake should be more than 3, aborting.")
        sys.exit()

    print("=== Creating DL_MESO input file for Nafion ===")
    print("Box size: %.1f | Water uptake: %.1f" % (L, lmbda))

    # ===== setting numbers
    Lcl, Lpt, Nsup, Npt = set_electrodes(data, L)
    Nbm, Nmc = 5, data["mono-per-chain"]
    N = int(rho_DPD * L**2*(L-2*Lcl))   # number of polymer beads, must be cubes
    Nc, Nw = calc_nc_nw(N, Nmc, Nbm, lmbda)
    Nbc = Nbm*Nmc              
    Nb = Nbc*Nc                
    mu, sigma = 1.0, 0.01   # arbitrary sigma for generating chains
    print("Monomers per chain: %i, Beads per monomer: %i" % (Nmc, Nbm))
    print(Nc, "polymer chains created")
    bead_types = "ABCWEP"
    bead_pop = {"A": 0, "B": 0, "C": 0, "W": Nw, "E": Nsup, "P": Npt}

    # ===== beads
    mono_beads = "AAABC"
    poly_xyz = grow_polymer(Nbm, Nc, Nmc, L, Lcl, mu)
    wb_xyz = gen_water_beads(Nw, L, Lcl)
    el_xyz = gen_el_support(Nsup, L, Lcl, Lpt)
    pt_xyz = gen_platinum(Npt, L, Lcl, Lpt)
    xyz = np.vstack((wb_xyz, el_xyz, pt_xyz, poly_xyz))
    config_names = ["W"]*Nw + ["E"]*Nsup + ["P"]*Npt + list(mono_beads)*Nmc*Nc
    print(len(xyz), "beads created, density:", float(len(xyz)) / L**3)
    
    fname = "CONFIG"
    save_config(fname, config_names, xyz)

    if args["--xyz"]:
        fname = args["--xyz"]
        bt2num = {}
        for i, bt in enumerate(bead_types):
            bt2num[bt] = i+1
        names = [bt2num[bt] for bt in config_names]
        xyz = np.hstack((np.matrix(names).T, xyz))
        ll.save_xyzfile(fname, xyz)
        print("xyz file saved in", fname)


    # ===== bonds
    bond_mat = gen_bonds(Nmc, mono_beads)
    print("FIELD: %i bonds in a chain" % len(bond_mat))
    bead_list = list("AAABC"*Nmc)
    nafion_mol_str = mol2str("nafion", Nc, bead_list, bond_mat)

    # ===== interaction and bond parameters
    a_ij = gen_inter_coeffs(data["chi-params"], bead_types, gamma)

    # ===== final FIELDS string
    field_string = "pokus\n\n" +\
                   species2str(bead_types, bead_pop) +\
                   inter2str(a_ij, method="dpd") + \
                   "MOLECULES 1\n" + \
                   nafion_mol_str + "\n" + \
                   "close\n"

    fname = "FIELD"
    open(fname, "w").write(field_string)
    print("FIELD file saved in %s" % fname)


