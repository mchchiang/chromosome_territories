# ContactFraction.py
#
# A program that reads the lammpstrj file and determine
# the fraction of beads in contact with the wall
#

import sys
import math

args = sys.argv
args.pop(0) # Ignore self

# Input arguments
cut_off = float(args.pop(0))
pos_file = args.pop(0)
frac_file = args.pop(0)

# Position file data columns
xcol = 1
ycol = 2
zcol = 3
ixcol = 4
iycol = 5
izcol = 6
typecol = 7

# Simulation parameters
lz = 40.0
nbeads = 6303
wall_pos = 20.0
time_inc = 1000

writer = open(frac_file, "w")

time = 0
index = 1
wall_bead_count = 0

with open(pos_file, "r") as f:
    for line in f:
        data = line.strip().split()
        if (len(data) == 8):
            z = float(data[zcol])
            iz = float(data[izcol])
            z += iz * lz

            if (z > (wall_pos-cut_off)):
                wall_bead_count += 1

            index += 1

        if (index > nbeads):
            # Compute fraction of beads touching the wall
            frac = wall_bead_count / float(nbeads)
            output = "%d %.5f\n" % (time, frac)
            writer.write(output)

            wall_bead_count = 0
            index = 1
            time += time_inc

writer.close()
