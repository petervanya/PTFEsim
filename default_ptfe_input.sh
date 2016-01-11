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

ksi-params:           # DPD Flory-Huggins params, Wu etal., EES (2008), DPD units of kBT
    A B: 1.23
    A C: 7.44
    A W: 3.36
    B C: 2.70
    B W: 1.53
    C W: 1.48
    A E: 0.12
    B E: 0.71
    C E: 0.74
    W E: 3.03
    A P: 63.75
    B P: 60.23
    C P: 62.80
    W P: 33.10

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
