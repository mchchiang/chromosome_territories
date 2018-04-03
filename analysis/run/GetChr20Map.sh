#!/bin/bash

# A script to extra chr20 map data from original contact map
growing_mat="growing_50000.matrix"
growing_chr20="contact_hic_growing_chr_20.dat"
senescence_mat="senescence_50000.matrix"
senescence_chr20="contact_hic_senescence_chr_20.dat"

awk '{ if (($1>=54383&&$1<=55643)&&($2>=54383&&$2<=55643)) {print $1-54383 FS $2-54383 FS $3} }' $growing_mat > $growing_chr20

awk '{ if (($1>=54383&&$1<=55643)&&($2>=54383&&$2<=55643)) {print $1-54383 FS $2-54383 FS $3} }' $senescence_mat > $senescence_chr20
