#!/usr/bin/env bash

mkdir -p seismograms
mkdir -p snapshots

FDFK2D ./input inpar.dat ./seismograms seisx.su seisz.su
