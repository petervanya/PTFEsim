#!/bin/bash
# Default input.yaml serving for creation of LAMMPS data file to simulate Nafion
# 09/11/15

cat >input.yaml <<EOF
# ===== notation:
# W: water bead 6 H2O
# E: electrodes from carbon
# P: platinum on electrodes
# =====
box-size:       5        # DPD units, 1 = 8.14 AA
temperature:    300       # K
PTFE:           Nafion    # defines beads
topology:       (A 3, [B 1], [C 1])15
water-uptake:   9         # water uptake, num H2O/SO3H
gamma:          1.0       # DPD drag coefficient

pair-coeffs:           # DPD pair coeffs, Wu etal., EES (2008), DPD units of kBT
    A B: 29.0
    A C: 49.3
    A W: 36.0
    B C: 33.8
    B W: 30.0
    C W: 29.9

bond-coeffs:           # spring const param k_ij
    A A: 4.0
    A B: 4.0
    A C: 4.0

equilibrium-dist: 1.0  # r0 in k(r-r0)^2

electrodes:
    width: 1           # width of one electrode, same on the other side
    Pt-amount: 0.01    # number of Pt beads per number of carbon black beads
EOF
echo "Parameter file saved in input.yaml"
