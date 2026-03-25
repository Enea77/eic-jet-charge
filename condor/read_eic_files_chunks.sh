#!/bin/bash
echo "Starting job on $(hostname)"

cd /afs/cern.ch/user/e/eprifti/EIC_DIS/CHJetsReCluster

CHUNK=$1
OUTNAME=$(basename $CHUNK .list).root
echo "Processing chunk: $CHUNK -> $OUTNAME"

root -l -b -q "process_eic_trees_list.C(\"${CHUNK}\", \"/eos/home-e/eprifti/EIC_files/${OUTNAME}\", 9999)"

echo "Job finished!"