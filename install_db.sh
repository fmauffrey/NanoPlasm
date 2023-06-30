#!/bin/bash

# Homopolish database install
mkdir data
cd data
wget http://140.123.104.107/bioinfo/mash_sketches/bacteria.msh.gz
gunzip bacteria.msh.gz

# KARGA
cd ..
git clone https://github.com/DataIntellSystLab/KARGA.git
cd data
wget https://www.meglab.org/downloads/megares_v3.00/megares_database_v3.00.fasta