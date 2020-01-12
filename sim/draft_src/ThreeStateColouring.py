# ThreeStateColouring.py
# A program that assign beads to one of the 3 possible states 
# (EU, H3K9me3, H3K27me3) based on ChIP-seq and LAD data

import sys

args = sys.argv

args.pop(0) # Ignore self
frac_file = args.pop(0) # The fractional state file
out_file = args.pop(0) # Output file

reader = open(frac_file, "r")
writer = open(out_file, "w")

def get_type(LAD, K9, K27):
    if (K27 > K9):
        return 3
    elif (LAD > 0 or K9 > 0):
        return 2
    else:
        return 1

for line in reader:
    # line format:
    # index fracLAD fracH3K9me3 fracH3K27me3
    data = line.split()
    index = int(data[0])
    LAD = float(data[1])
    H3K9 = float(data[2])
    H3K27 = float(data[3])
    bead_type = get_type(LAD, H3K9, H3K27)
    writer.write("{:d} {:d}\n".format(index, bead_type))

reader.close()
writer.close()

