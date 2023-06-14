#!/bin/bash

# Homopolish database install
mkdir data
cd data
wget http://140.123.104.107/bioinfo/mash_sketches/bacteria.msh.gz
gunzip bacteria.msh.gz
