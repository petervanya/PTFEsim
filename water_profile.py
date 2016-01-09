#!/usr/bin/env python
"""Usage:
    water_profile.py <files> (1d --axis <axis> | 2d --plane <plane>) [--plot]

Arguments:
    <files>         Dump files from LAMMPS
    1d              1d profile of water
    2d              2d profile of water
    --axis <axis>   Select axis to profile along ("x", "y", "z")
    --plane <plane> Which plane to profile on ("xy", "yz", "xz")

pv278@cam.ac.uk, 03/01/15
"""
import numpy as np
from math import log, sqrt
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import glob, sys
from docopt import docopt


def read_outfile(outfile):
    """Read one xyz outfile into a numpy matrix"""
    A = open(outfile, "r").readlines()[2:]
    A = np.array([line.split() for line in A], order="F").astype(float)
    return A


def set_num_bins(A):
    """Set number of histogram bins, various recipes"""
    N = len(A)
    nbins = int(log(N, 2)+1) + 1   # Sturges' formula, from wiki
    nbins = int(2*N**(1./3)) + 1   # Rice rule
    nbins = int(sqrt(N)) + 1
    return nbins


def create_1d_profile(dumpfiles, axis):
    Nfiles = len(dumpfiles)
    nbins = 0
    for dumpfile in dumpfiles:
        A = read_outfile(dumpfile)
        A = A[A[:, 0] == 4][:, 1:]   # water beads, SHOULD CONSIDER C AS WELL?
        if nbins == 0:
            nbins = set_num_bins(A)
            print "Number of bins =", nbins
            res = np.zeros(nbins)
        profile, bins = np.histogram(A[:, axis], bins=nbins)
        res += profile/float(Nfiles)
    return res


def create_2d_profile(dumpfiles, plane):
    Nfiles = len(dumpfiles)
    nbins = 0
    for dumpfile in dumpfiles:
        A = read_outfile(dumpfile)
        A = A[A[:, 0] == 4][:, 1:]    # water beads, SHOULD CONSIDER C AS WELL?
        if nbins == 0:
            nbins = set_num_bins(A) 
            print "Number of bins =", nbins
            res = np.zeros((nbins, nbins))
        profile, x, y = np.histogram2d(A[:, plane[0]], A[:, plane[1]], bins=nbins)
        res += profile/float(Nfiles)
    return res


if __name__ == "__main__":
    args = docopt(__doc__)
#    print args
    axes = {"x": 0, "y": 1, "z": 2}
    planes = {"xy": (0, 1), "yz": (1, 2), "xz": (0, 2)}
    dumpfiles = glob.glob(args["<files>"])
    if not dumpfiles:
        print "No files captured, aborting."
        sys.exit()
    else:
        print dumpfiles

    if args["1d"]:
        try:
            axis = axes[args["<axis>"]]
        except KeyError:
            print "Wrong axis, choose from 'x', 'y', 'z'."
            sys.exit()
        profile = create_1d_profile(dumpfiles, axis)
        outname = "profile_1d.out"
        np.savetxt(outname, profile)
        print "Array saved in", outname
    elif args["2d"]:
        try:
            plane = planes[args["<plane>"]]
        except KeyError:
            print "Wrong plane, choose from 'xy', 'yz', 'xz'."
            sys.exit()
        profile = create_2d_profile(dumpfiles, plane)
        outname = "profile_2d.out"
        np.savetxt(outname, profile)
        print "Matrix saved in", outname
        if args["--plot"]:
            plt.imshow(profile, cmap = cm.Greys_r)
            plotname = "plot_2d.png"
            plt.savefig(plotname)
            print "Plot saved in", plotname
#            plt.show()



