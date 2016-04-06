#!/usr/bin/env python
"""
[AD HOC] Test how length of the distance vector changes
with limiting cutoff rc.

Usage:
   test_rdf_scaling.py <N> [--L <L> --bins <nb> --plot]

Arguments:
    <N>            Number of particles [default: 1000]

Options:
    --L <L>        Box size [default: 1.0]
    --bins <nb>    Number of bins [default: 100]

05/04/16
"""
import numpy as np
import matplotlib.pyplot as plt
from docopt import docopt


def dist_vec_HH(xyz):
    """No PBCs"""
    dist_mat = np.sqrt(np.sum((xyz[None, :]-xyz[:, None])**2, 2))
    dist_vec = np.tril(dist_mat).reshape(-1)
    return dist_vec[dist_vec != 0.0]


def dist_vec_HH2(xyz):
    """Consider PBCs by replicating the system 6x
    GIVES WRONG RESULTS, UNFINISHED"""
    N = len(xyz)
    dist_mat = xyz[None, :] - xyz[:, None]  # FINISH
    dist_vec = np.tril(dist_mat).reshape(-1)
    return dist_vec[dist_vec != 0.0]


def dist_vec_naive(xyz, L):
    """With PBCs, WRONG"""
    N = xyz.shape[0]
    dist_vec = np.zeros(N*(N-1)/2)
    d_imag = np.zeros(6)
    imag_mat = L * np.vstack((np.zeros((1, 3)), np.eye(3), -np.eye(3)))
    cnt = 0
    for i in range(N):
        for j in range(i):
            for k in range(6):
                d_imag[k] = np.sqrt( np.sum((xyz[i] - xyz[j] + imag_mat[k])**2) )
            dist_vec[cnt] = np.min(d_imag)
            cnt += 1
    return dist_vec


def dist_vec_naive2(xyz, cell):
    """With PBCs"""
    N = xyz.shape[0]
    inv_cell = np.linalg.inv(cell)
    dist_vec = np.zeros(N*(N-1)/2)
    d_imag = np.zeros(6)
    imag_mat = L * np.vstack((np.zeros((1, 3)), np.eye(3), -np.eye(3)))
    cnt = 0
    for i in range(N):
        for j in range(i):
            dr = xyz[i] - xyz[j]
            G = np.dot(inv_cell, dr)   # coords of dr in basis of the cell, range [-1,1]
            G_n = G- np.round(G)       # clever bit
            dr_n = np.dot(cell, G_n)
            dist_vec[cnt] = np.linalg.norm(dr_n)
            cnt += 1
    return dist_vec


if __name__ == "__main__":
    args = docopt(__doc__)
    N = int(args["<N>"])
    L = float(args["--L"])
    Nbins = int(args["--bins"])

    xyz = np.random.rand(N, 3) * L
    cell = L * np.eye(3)
    dist_vec = dist_vec_naive2(xyz, cell)
#    dist_vec = dist_vec_HH(xyz)
#    dist_vec = dist_vec_HH2(xyz)
    
    start = 0.1
    end = 1.5
    pts = np.round((end-start)/start) + 1
    print pts
    rc = np.linspace(start, end, pts)
    data = np.zeros((len(rc), 2))
    for i in range(len(rc)):
        l = len(dist_vec[dist_vec < rc[i]])
        print rc[i], l
        data[i] = [rc[i], l]
    
    if args["--plot"]:
        plt.plot(data[:, 0], data[:, 1])
        plt.show()
    
    bins = np.linspace(0.0, 1.0, Nbins)
    h, r = np.histogram(dist_vec, bins)
    r = r[:-1] + np.diff(r)/2.0
    dr = r[1] - r[0]
    print len(h), len(r)
     
    # Correct normalisation, from researchgate.net/bla
    h = h * L**3/len(dist_vec) / (4 * np.pi * r**2 * dr)
    plt.plot(r, h)
    plt.ylim([0.0, 2.0])
    plt.savefig("rdf_test.png")


