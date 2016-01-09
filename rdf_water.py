#!/usr/bin/env python
"""Usage:
    rdf.py <fnames> [--binsize <bs>]

Read LAMMPS data files using regex and compute radial distribution function
for water beads:
* C: 3 molecules (SO3 + 3H2O) (bead 3)
* W: 6 molecules (bead 4)

Arguments:
    <fnames>         Regex for all the required xyz files

Options:
    --binsize <bs>   Size of the bins for the histogram, CAREFUL WITH UNITS [default: 0.05]

pv278@cam.ac.uk, 09/01/16
"""
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
from math import pi
import glob, sys
from docopt import docopt
import f_rdf      # Fortran module

bins = np.arange(0)

def set_num_bins(N):
    """Set number of histogram bins, various recipes"""
    nbins = int(2*N**(1./3)) + 1   # Rice rule
    nbins = int(sqrt(N)) + 1       # most primitive
    nbins = int(log(N, 2)+1) + 1   # Sturges' formula, from wiki
    return nbins


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


def compute_rdf(outfile, nbins=nbins):
    """Compute radial dist'n fcn from the xyz frame using Fortran routine"""
    global bins
    A = read_outfile(outfile)
    xyz_C = A[A[:, 0] == 3][:, 1:]
    xyz_W = A[A[:, 0] == 4][:, 1:]

    d_C = f_rdf.f_rdf.pair_dist_mat(xyz_C)         # rdf for beads C
    rdf_raw_C, r = np.histogram(d, bins=nbins)
    del d_C
    d_W = f_rdf.f_rdf.pair_dist_mat(xyz_W)         # rdf for beads W
    rdf_raw_W, r = np.histogram(d_W, bins=nbins)
    del d_W
    d_CW = f_rdf.f_rdf.pair_dist_mat(xyz_C, xyz_W) # rdf for combined beads C and W
    rdf_raw_CW, r = np.histogram(d_CW, bins=nbins)
    del d_CW
    
    rdf_raw = rdf_raw_C * 3**2 + rdf_raw_W * 6**2 + rdf_raw_CW * 3*6
    r = r[:-1] + np.diff(r)/2.0
    dr = r[1] - r[0]
    rdf = rdf_raw/(4*pi*r**2 * dr)
    return r, rdf


def master_rdf(outfiles, nbins=nbins):
    """Construct an rdf from all the available xyz files"""
    rdf_mat = []
    for outfile in outfiles:
        r, rdf_i = compute_rdf(outfile, nbins)
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
    dr = float(args["--binsize"])
    outfiles = glob.glob(args["<fnames>"])
    if len(outfiles) == 0:
        print "ERROR: No xyz files found. Aborting."
        sys.exit()
    print outfiles
    Nfiles = len(outfiles)
    N = int(open(outfiles[0], "r").readline())
    print "Total particles:", N
    Nbins = set_num_bins(N)
    print "Number of bins:", Nbins

    r, vals = master_rdf(outfiles, Nbins)
    fname = "rdf_water.out"
    save_data(fname, r, vals)
    print "rdf saved in", fname

