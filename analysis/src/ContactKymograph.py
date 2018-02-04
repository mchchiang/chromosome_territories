# ContactKemograph.py
#
# A program that reads the lammpstrj file and determine
# when beads are in contact with the wall
#

import sys
import math

args = sys.argv
args.pop(0) # Ignore self

# Input arguments
cut_off = float(args.pop(0))
pos_file = args.pop(0)
kymo_file = args.pop(0)

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

writer = open(kymo_file, "w")

time = 0
index = 1

with open(pos_file, "r") as f:
    for line in f:
        data = line.strip().split()
        if (len(data) == 8):
            z = float(data[zcol])
            iz = float(data[izcol])
            z += iz * lz
            colour = int(data[typecol])
            if (z <= (wall_pos-cut_off)):
                colour = 0

            output = "%d %d %d\n" % (time, index, colour)
            writer.write(output)

            index += 1

        if (index > nbeads):
            writer.write("\n\n")
            
            index = 1
            time += time_inc

writer.close()

