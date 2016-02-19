#!/usr/bin/env python
"""Usage:
    rdf_water.py (--bead <b> | water) <fnames> [--bins <nbins>] [--binalg <b>]

Read LAMMPS data files using regex and compute radial distribution function
for any DPD beads or separetely water beads (that could be found in beads C and W):
* C: 3 molecules (SO3 + 3H2O) (bead 3)
* W: 6 molecules (bead 4)

Arguments:
    --beadtype <b>   Number from 1 to num. bead types, typically 4
    <fnames>         Regex for all the required xyz files

Options:
    --bins <nbins>   Number of bins
    --binalg <b>    'numpy' or 'fortran [Default: fortran]

pv278@cam.ac.uk, 09/01/16
"""
import numpy as np
from math import *
import glob, sys
from docopt import docopt
import lmp_lib as ll
from f_rdf import f_rdf      # Fortran module


def set_num_bins(N, method="sqrt"):
    """Set number of histogram bins, various recipes.
    Available methods: rice, sturges, sqrt (default)"""
    if method == "rice":
        return int(2*N**(1./3)) + 1   # Rice rule
    elif method == "sturges":
        return int(log(N, 2)+1) + 1   # Sturges' formula
    else:
        return int(sqrt(N)) + 1       # most primitive


def read_outfile(outfile):
    """Read one xyz outfile into a numpy matrix"""
    A = open(outfile, "r").readlines()[2:]
    A = [line.split() for line in A]
    A = np.array(A, order="F").astype(float)
    return A


def save_data(outfile, *args):
    """Save two vectors into file"""
    m, n = len(args), len(args[0])   # cols, rows
    args = zip(*args)
    with open(outfile, "w") as f:
        for i in range(n):
            line = ""
            for j in range(m):
                line += str(args[i][j]) + "\t"
            line += "\n"
            f.write(line)


def compute_rdf_water_np(outfile, nbins=30):
    """Compute radial dist'n fcn from the xyz frame using 
    Fortran routine for pair distances and numpy binning
    """
    A = read_outfile(outfile)
    xyz_C = A[A[:, 0] == 3][:, 1:]
    xyz_W = A[A[:, 0] == 4][:, 1:]

    d_C = f_rdf.pair_dist_arr(xyz_C)         # rdf for beads C
#    print "  Distance matrix for C beads done."
    rdf_raw_C, r = np.histogram(d_C, nbins)
#    print "  Binning for C beads done."
    del d_C
    d_W = f_rdf.pair_dist_arr(xyz_W)         # rdf for beads W
#    print "  Distance matrix for W beads done."
    rdf_raw_W, r = np.histogram(d_W, nbins)
#    print "  Binning for W beads done."
    del d_W
    d_CW = f_rdf.pair_dist_arr2(xyz_C, xyz_W) # rdf for combined beads C and W
#    print "  Distance matrix for CW beads done."
    rdf_raw_CW, r = np.histogram(d_CW, nbins)
#    print "  Binning for CW beads done."
    del d_CW

    rdf_raw = rdf_raw_C * 3**2 + rdf_raw_W * 6**2 + rdf_raw_CW * 3*6
    r = r[:-1] + np.diff(r)/2.0
    dr = r[1] - r[0]
    rdf = rdf_raw/(4*pi*r**2 * dr)
    return r, rdf


def compute_rdf_water(outfile, nbins=30):
    """Compute radial dist'n fcn from the xyz frame
    using Fortran routine, both distance matrix and binning"""
    A = read_outfile(outfile)
    xyz_C = A[A[:, 0] == 3][:, 1:]
    xyz_W = A[A[:, 0] == 4][:, 1:]

    rdf_raw_C, r = f_rdf.pair_dist_hist(xyz_C, nbins)
#    print "  Bead C pair dist and binning beads done"
    rdf_raw_W, r = f_rdf.pair_dist_hist(xyz_W, nbins)
#    print "  Bead W pair dist and binning beads done"
    rdf_raw_CW, r = f_rdf.pair_dist_hist2(xyz_C, xyz_W, nbins)
#    print "  Beads C and W pair dist and binning beads done"
    
    rdf_raw = rdf_raw_C * 3**2 + rdf_raw_W * 6**2 + rdf_raw_CW * 3*6
    r = r[:-1] + np.diff(r)/2.0
    dr = r[1] - r[0]
    rdf = rdf_raw/(4*pi*r**2 * dr)
    return r, rdf


def master_rdf_water(outfiles, nbins=30, method="fortran"):
    """Construct an rdf for water beads
    from all the available xyz files"""
    rdf_mat = []
    if method == "fortran":
        for outfile in outfiles:
            r, rdf_i = compute_rdf_water(outfile, nbins)
            rdf_mat.append(rdf_i)
            print outfile, "done."
    if method == "numpy":
        for outfile in outfiles:
            r, rdf_i = compute_rdf_water_np(outfile, nbins)
            rdf_mat.append(rdf_i)
            print outfile, "done."
    rdf_mat = np.array(rdf_mat).T
    np.savetxt("rdf_mat.out", rdf_mat)
    print "rdf matrix saved in rdf_mat.out"
    rdf = np.array(np.sum(rdf_mat, 1) / len(outfiles))
    return r, rdf


def compute_rdf(outfile, beadtype, nbins=30):
    """Compute RDF from the xyz frame uusing Fortran routine, distance matrix and binning"""
    A = read_outfile(outfile)
    xyz = A[A[:, 0] == beadtype][:, 1:]
    rdf_raw, r = f_rdf.pair_dist_hist(xyz, nbins)  # key routine
    r = r[:-1] + np.diff(r)/2.0
    dr = r[1] - r[0]
    rdf = rdf_raw/(4*pi*r**2 * dr)
    return r, rdf


def master_rdf(outfiles, beadtype, nbins=30):
    """Construct an rdf for given bead type 
    from all the available xyz files"""
    rdf_mat = []
    for outfile in outfiles:
        r, rdf_i = compute_rdf(outfile, beadtype, nbins)
        rdf_mat.append(rdf_i)
        print outfile, "done."
    rdf_mat = np.array(rdf_mat).T
    np.savetxt("rdf_mat.out", rdf_mat)
    print "rdf matrix saved in rdf_mat.out"
    rdf = np.array(np.sum(rdf_mat, 1) / len(outfiles))
    return r, rdf


if __name__ == "__main__":
    args = docopt(__doc__)
#    print args
    outfiles = glob.glob(args["<fnames>"])
    beadtype = "water" if args["water"] else args["<b>"]
    if len(outfiles) == 0:
        raise ValueError("No xyz files captured, aborting.")
    print outfiles
    Nfiles = len(outfiles)
    N = int(open(outfiles[0], "r").readline())
    
    if args["--bins"]:
        Nbins = int(args["--bins"])
    else:
        Nbins = set_num_bins(N, method="sturges")
    method = args["--binalg"]

    A = ll.read_xyzfile(outfiles[0])
    Nbt = len(set(A[:, 0]))
    Nb = len(A[A[:, 0] == int(beadtype)])

    print "Total beads:", N, "| Num. beadtypes:", Nbt, "| Bins:", Nbins
    print "Bead type:", beadtype, "| Beads of this type:", Nb
    
    if beadtype == "water":
        r, vals = master_rdf_water(outfiles, Nbins, method)
        fname = "rdf_water.out"
    else:
        r, vals = master_rdf(outfiles, int(beadtype), Nbins)
        fname = "rdf.out"
    save_data(fname, r, vals)
    print "rdf saved in", fname


