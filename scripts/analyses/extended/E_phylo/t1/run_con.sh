#!/bin/bash
h=$1
iqtree2 -s red.faa -m LG+G4 -g con_$h.nwk -t start_$h.nwk -T 1 --prefix cs_$h -redo -quiet > cs_$h.stdout 2>&1
