#!/bin/bash
# Default input.yaml serving for creation of LAMMPS data file to simulate Nafion
# 09/11/15

cat >input.yaml <<EOF
box-size: 200
temperature: 300

membrane:   # Nafion values for m, n
    m: 15
    n: 17
    lambda: 9

int-params:  # from We et. al., EES (2008)
    A B: 29.0
    A C: 49.3
    A W: 36.0
    B C: 33.8
    B W: 30.0
    C W: 29.9

electrodes:
    width: FILL
    Pt-amount: FILL
    Pt-size: FILL
EOF
