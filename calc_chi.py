#!/usr/bin/env python
"""Usage:
    calc_chi.py <data> [--temp <T>]

Calculate Flory-Huggins chi params for all beads based on
solubility and molar volume

Options: 
    --temp <T>          Temperature [default: 300]

pv278@cam.ac.uk, 17/02/16
"""
import numpy as np
import pandas as pd
from docopt import docopt

R = 8.314

def chi(b1, b2, T):
    """Return chi based from lines of DataFrame"""
    return (b1.Vm + b2.Vm)/(2*R*T) * (b1.deltaT - b2.deltaT)**2


def gen_chi_mat(tbl, T):
    """Return matrix of chi values"""
    Nidx = len(tbl)
    beads = tbl.index.values
    chi_mat = np.zeros((Nidx, Nidx))
    for i in range(Nidx):
        for j in range(i):
            chi_mat[i, j] = chi(tbl.ix[beads[i]], tbl.ix[beads[j]], T)
    return chi_mat + chi_mat.T


if __name__ == "__main__":
    args = docopt(__doc__)
    T = float(args["--temp"])
    tbl = pd.read_csv(args["<data>"], index_col=0)
#    print tbl
    
    cols = tbl.index.values
#    print tbl.index.values
#    print tbl.ix["A"].Vm
#    print type(tbl.ix["A"])
#    print len(tbl)

    mat = gen_chi_mat(tbl, T)
    print pd.DataFrame(mat, columns=cols, index=cols)

