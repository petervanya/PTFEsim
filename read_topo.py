#!/usr/bin/env python
"""Usage:
    read_topo.py [<topostr>] [<Nm>]

Playing with reading topologies to create a polymer chain from 
the given topology.
Two types of beads: A1, and [A1]
* Nm -- number of monomers in chain (polymerisation)

Options:
    <topostr>     The topology string
    <Nm>          Polymerisation
"""
import numpy as np
import re
from docopt import docopt


def bead_parser(mono):
    """From a string parse the beads and numbers in one monomer"""
    topo_entries = [s.strip(" ") for s in mono.split(",")]
    print "Topology:", topo_entries
    Nbm = sum([int(i) for i in re.findall(r"\d+", mono)])
    bead_list = []
#    for i in range(Nbm):
    cnt = 1
    for topo in topo_entries:
         if topo[0] == "[" and topo[-1] == "]":   # sides
             typ, num = topo.strip("[]").split()[0], int(topo.strip("[]").split()[1])
             for i in range(num):
                 bead_list.append([cnt, typ, "side", cnt-1])
                 cnt += 1
         else:                                    # backbone
             typ, num = topo.split()[0], int(topo.split()[1])
             for i in range(num):
                 if i == 0:
                     bond = cnt - 3
                 else:
                     bond = cnt - 1
                 bead_list.append([cnt, typ, "backbone", bond])
                 cnt += 1
    return bead_list


def construct_chain(bead_list, Nm):
    """From a list of beads in one monomer construct whole chain 
    with numbers and bonds"""
    pass


args = docopt(__doc__)
if args["<topostr>"]:
    s = args["<topostr>"]
else:
    s = "(A 3, [B 1], [C 1]) 15"

mono, Nm = s.split(")")[0], int(s.split(")")[-1])
mono = mono.lstrip("(")
Nbm = sum([int(i) for i in re.findall(r"\d+", mono)]) # or re.findall(r"[0-9]", mono)
print "Number of beads in monomer:", Nbm

mat = bead_parser(mono)
print np.array(mat)



