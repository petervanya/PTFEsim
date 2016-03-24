#!/usr/bin/env python
"""Usage:
    plot_all_beads.py <fileW> <fileS> <fileB> <savefile> [--title <ttl>]

[AD HOC] script to format plotting of normalised water profiles
for water, sulfonic acid groups and PTFE backbone

25/02/16
"""
import numpy as np
import matplotlib.pyplot as plt
import glob, sys
from docopt import docopt


args = docopt(__doc__)
#print args
fileW = args["<fileW>"]
fileS = args["<fileS>"]
fileB = args["<fileB>"]
lw = 2.0    # linewidth

A = np.loadtxt(fileB)
A[:, 0] *= 1e9
A[:, 1] /= sum(A[:, 1])
plt.plot(A[:, 0], A[:, 1], "red", label="backbone", linewidth=lw)

A = np.loadtxt(fileS)
A[:, 0] *= 1e9
A[:, 1] /= sum(A[:, 1])
plt.plot(A[:, 0], A[:, 1], "green", label="sulfonic", linewidth=lw)

A = np.loadtxt(fileW)
A[:, 0] *= 1e9
A[:, 1] /= sum(A[:, 1])
plt.plot(A[:, 0], A[:, 1], "blue", label="water", linewidth=lw)

plt.xlabel("$x$ (nm)")
plt.xlim([0.0, 32.0])
plt.yticks([])

if args["--title"]:
   plt.title(args["<ttl>"])
plt.savefig(args["<savefile>"])
print "Plot saved in", args["<savefile>"]

