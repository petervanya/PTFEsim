#!/usr/bin/env python
"""Usage:
    default_ptfe_input.py [--L <L> --lmbda <l> --gamma <g>]
                          [--el <el> --width <w> --rpt <rpt>]

Generate input.yaml file that will serve to generate
a LAMMPS data file by being read by gen_ptfe.py script.
This script supersedes the default_ptfe_input*.sh scripts.

Options:
    --L <L>             Box size [default: 40]
    --lmbda <l>         Water uptake [default: 9]
    --gamma <g>         DPD friction gamma [default: 4.5]
    --el <el>           Electrode type: carbon, silica, quartz [default: carbon]
    --width <w>         Electrode width [default: 5.0]
    --rpt <rpt>         Ratio of width of Pt segment on electrode [default: 0.0]

pv278@cam.ac.uk, 16/06/16
"""
import yaml
import sys
from docopt import docopt


chi = {}
chi["carbon"] = \
"""
A B: 1.23
A C: 7.44
B C: 2.70
A W: 3.36
B W: 1.53
C W: 1.48
A E: 0.29
B E: 0.06
C E: 0.06
W E: 2.96
"""
chi["quartz"] = \
"""
A B: 1.23
A C: 7.44
B C: 2.70
A W: 3.36
B W: 1.53
C W: 1.48
A E: 1.82
B E: 0.25
C E: 0.25
W E: 0.94
"""
chi["silica"] = \
"""
A B: 1.23
A C: 7.44
B C: 2.70
A W: 3.36
B W: 1.53
C W: 1.48
A E: 1.51
B E: 0.10
C E: 0.10
W E: 1.56
"""


args = docopt(__doc__)
L = float(args["--L"])
T = 1.0
Nmc = 15
lmbda = float(args["--lmbda"])
gamma = float(args["--gamma"])
k0 = 4.0
r0 = 0.1
el = args["--el"]
w = float(args["--width"])
rPt = float(args["--rpt"])

if el not in ["carbon", "silica", "quartz"]:
    print("ERROR. Choose from electrodes: carbon, silica, quartz.")
    sys.exit()

if w < 0.0 or w > L/2:
    print("ERROR: Electrode width is at most one half of box size.")
    sys.exit()

if rPt < 0.0 or rPt > 1.0:
    print("ERROR: Platinum segment must be between 0 and 1.")
    sys.exit()

s = \
"""# ===== Bead types:
# * A, B, C
# * W: water bead 6 H2O
# * E: electrodes from carbon/quartz/silica
# * P: platinum
# =====\n"""

s += "box-size:          %.0f        # DPD units, 1 = 8.14 AA\n" % L
s += "temperature:       %.1f       # Units of kB T\n" % T
s += "mono-per-chain:    %i\n" % Nmc
s += "water-uptake:      %.0f         # number of H2O/SO3H\n" % lmbda
s += "gamma:             %.1f       # DPD drag coefficient\n" % gamma
s += "topology:          (A 3, [B 1], [C 1])15\n\n"

s += "chi-params:\n"
for k, v in yaml.load(chi[el]).items():
    s += "    %s: %.2f\n" % (k, v)

s += "\n"
s += "bond-coeffs:\n    A A: %.1f\n\n" % k0
s += "equilibrium-dist:  %.1f\n\n" % r0

s += "electrodes:\n"
s += "    width:         %.1f       # electrode width\n" % w 
s += "    Pt-amount:     %.1f       # ratio of Pt/C segment\n" % rPt

fname = "input.yaml"
open(fname, "w").write(s)
print("Parameter file saved in", fname)



