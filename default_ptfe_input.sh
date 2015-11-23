#!/bin/bash
# Default input.yaml serving for creation of LAMMPS data file to simulate Nafion
# 09/11/15

cat >input.yaml <<EOF
box-size:    20        # DPD units, 1 = 8.14 AA
temperature: 300       # K
PTFE:        Nafion    # defines beads
topology:    (A 3, [B 1], [C 1])15

membrane:              # Nafion: polymerisation m, num CF2 n
    m: 15
    n: 17
    lambda: 9          # water uptake, num H2O/SO3H

int-params:            # from Wu etal., EES (2008), DPD units of kBT
    A B: 29.0
    A C: 49.3
    A W: 36.0
    B C: 33.8
    B W: 30.0
    C W: 29.9

spring-const:
    A A: 4.0
    A B: 4.0
    A C: 4.0

electrodes:
    width: FILL
    Pt-amount: FILL
    Pt-size: FILL
EOF
