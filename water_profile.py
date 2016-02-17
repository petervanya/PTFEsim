#!/usr/bin/env python
"""Usage:
    water_profile.py <files> (1d <axis> | 2d <plane>) [--bins <bins>]
                     [--subst <s>] [--boxsize <L>]

Create 1d or 2d profile of water molecules in Nafion w or wo electrodes.

Arguments:
    <files>         Dump files from LAMMPS
    1d              1d profile of water
    2d              2d profile of water
    --axis <axis>   Select axis to profile along ("x", "y", "z")
    --plane <plane> Which plane to profile on ("xy", "yz", "xz")

Options:
    --boxsize <L>   Box size [default: 40]
    --subst <s>     Requested substance: "water", "sulfonic", "backbone" [default: water]

pv278@cam.ac.uk, 03/01/15
"""
import numpy as np
from math import log, sqrt
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import glob, sys
from docopt import docopt
import lmp_lib as ll

rc = 8.14e-10

def read_outfile(outfile):
    """Read one xyz outfile into a numpy matrix"""
    A = open(outfile, "r").readlines()[2:]
    A = np.array([line.split() for line in A], order="F").astype(float)
    return A


def set_num_bins(N, method="sqrt"):
    """Set number of histogram bins, various recipes.
    Available methods: rice, sturges, sqrt (default)"""
    if method == "rice":
        return int(2*N**(1./3)) + 1   # Rice rule
    elif method == "sturges":
        return int(log(N, 2)+1) + 1   # Sturges' formula
    else:
        return int(sqrt(N)) + 1       # most primitive


def create_1d_profile(dumpfiles, axis, nbins):
    """DEPRECATED"""
    Nfiles = len(dumpfiles)
    res = np.zeros(nbins)
    for dumpfile in dumpfiles:
        A = read_outfile(dumpfile)
        A = A[A[:, 0] == 4][:, 1:]   # water beads, SHOULD CONSIDER C AS WELL?
        profile, bins = np.histogram(A[:, axis], bins=nbins)
        res += profile/float(Nfiles)
    bins = bins[:-1] + np.diff(bins)/2.0
    return res, bins


def create_1d_profile2(dumpfiles, axis, subst, bins):
    """NEW function with pre-set bins and fixed water plot"""
    Nfiles = len(dumpfiles)
    res = np.zeros(len(bins)-1)
    if subst == "water":
        for dumpfile in dumpfiles:
            A = read_outfile(dumpfile)
            beads3 = A[A[:, 0] == 3][:, 1:]   # bead C, 3 H2O molecules
            beads4 = A[A[:, 0] == 4][:, 1:]   # bead W, 6 H2O molecules
            profile3, nic = np.histogram(beads3[:, axis], bins=bins)
            profile4, nic = np.histogram(beads4[:, axis], bins=bins)
            res += 3*profile3/float(Nfiles)
            res += 6*profile4/float(Nfiles)
    elif subst == "sulfonic":
        for dumpfile in dumpfiles:
            A = read_outfile(dumpfile)
            beads = A[A[:, 0] == 3][:, 1:]
            profile, nic = np.histogram(beads[:, axis], bins=bins)
            res += profile/float(Nfiles)
    elif subst == "backbone":
        for dumpfile in dumpfiles:
            A = read_outfile(dumpfile)
            beads1 = A[A[:, 0] == 1][:, 1:]   # bead A, 6 CF2 groups
            beads2 = A[A[:, 0] == 2][:, 1:]   # bead B, 5 CF2 groups
            beads3 = A[A[:, 0] == 3][:, 1:]   # bead B, 1 CF3 group
            profile1, nic = np.histogram(beads1[:, axis], bins=bins)
            profile2, nic = np.histogram(beads2[:, axis], bins=bins)
            profile3, nic = np.histogram(beads3[:, axis], bins=bins)
            res += 6*profile1/float(Nfiles)
            res += 5*profile2/float(Nfiles)
            res += profile3/float(Nfiles)
    return res, bins


def create_2d_profile(dumpfiles, plane, nbins):
    """TO DO: rethink the binning process"""
    Nfiles = len(dumpfiles)
    res = np.zeros((nbins, nbins))
    for dumpfile in dumpfiles:
        A = read_outfile(dumpfile)
        A = A[A[:, 0] == 4][:, 1:]    # water beads, SHOULD CONSIDER C AS WELL?
        profile, x, y = np.histogram2d(A[:, plane[0]], A[:, plane[1]], bins=nbins)
        res += profile/float(Nfiles)
    return res


if __name__ == "__main__":
    args = docopt(__doc__)
#    print args
    L = float(args["--boxsize"])*rc
    axes = {"x": 0, "y": 1, "z": 2}
    planes = {"xy": (0, 1), "yz": (1, 2), "xz": (0, 2)}
    dumpfiles = glob.glob(args["<files>"])
    if not dumpfiles:
        print "No files captured, aborting."
        sys.exit()
    A = ll.read_xyzfile(dumpfiles[0])
    if args["--bins"]:
        nbins = int(args["<bins>"])
    else:
        nbins = set_num_bins(len(A))
    subst = args["--subst"]
    if subst not in ["water", "sulfonic", "backbone"]:
        print "Only 'water', 'sulfonic', 'backbone' available, aborting."
        sys.exit()
    print "Substance:", subst, "| Number of bins:", nbins, "| Box size:", L/rc
    print dumpfiles
    
    if args["1d"]:
        try:
            axis = axes[args["<axis>"]]
        except KeyError:
            print "Wrong axis, choose from 'x', 'y', 'z'. Aborting."
            sys.exit()
        bins = np.linspace(0, L, nbins)
#        profile, bins = create_1d_profile(dumpfiles, axis, nbins)
        profile, bins = create_1d_profile2(dumpfiles, axis, subst, bins)
        bins = bins[:-1] + np.diff(bins)/2.0
        outname = "profile_1d.out"
        np.savetxt(outname, zip(bins, profile))
        print "Array saved in", outname

        plt.plot(bins, profile)
        plt.xlim([0, bins[-1]])
        plotname = "profile_1d.png"
        plt.savefig(plotname)
        print "Plot saved in", plotname

    elif args["2d"]:
        try:
            plane = planes[args["<plane>"]]
        except KeyError:
            print "Wrong plane, choose from 'xy', 'yz', 'xz'."
            sys.exit()
        profile = create_2d_profile(dumpfiles, plane, nbins)
        outname = "profile_2d.out"
        np.savetxt(outname, profile)
        print "Matrix saved in", outname

        plt.imshow(profile, cmap = cm.Greys_r)
        plt.axis("off")
        plotname = "profile_2d.png"
        plt.savefig(plotname)
        print "Plot saved in", plotname



