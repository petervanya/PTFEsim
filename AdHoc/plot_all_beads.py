#!/usr/bin/env python
"""Usage:
    plot_all_beads.py <fileW> <fileB> [--boxsize <L>]

[AD HOC] script to format plotting of normalised water profiles
for water and PTFE backbone.
Previously also plotted sulfonic acid groups.

Options:
    --boxsize <L>       Box size [default: 20.0]

25/02/16
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import glob, sys
from docopt import docopt

rc = 8.14e-10
matplotlib.rcParams.update({'font.size': 16})

args = docopt(__doc__)
L = float(args["--boxsize"])
fileW = args["<fileW>"]
fileB = args["<fileB>"]
lw = 3.0    # linewidth
print("Box size: %.2f" % L)

A = np.loadtxt(fileB)
A[:, 0] *= rc*1e9   # convert to nm
A[:, 1] /= sum(A[:, 1])
plt.plot(A[:, 0], A[:, 1], "red", label="backbone", linewidth=lw)

A = np.loadtxt(fileW)
A[:, 0] *= rc*1e9
A[:, 1] /= sum(A[:, 1])
plt.plot(A[:, 0], A[:, 1], "blue", label="water", linewidth=lw)

plt.xlabel("$x$ (nm)", fontsize=20)
plt.xlim([0.5, np.max(A[:, 0]) - 0.5])
plt.legend(loc="best")

imgname = "plot.png"
plt.savefig(imgname)
print("Plot saved in", imgname)


