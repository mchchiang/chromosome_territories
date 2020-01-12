# rebin_chromatin_state.py
# A program which rebins the chromatin state to suit lower resolutions
# Only works for a two-state (1 or 2) system!

import sys
import math

args = sys.argv
args.pop(0) # Ignore self

if (len(args) != 4):
    print "Usage: bin_chromatin_state.py [nbeads] [binsize] [infile] [outfile]"
    sys.exit(1)

nbeads = int(args.pop(0))
binsize = int(args.pop(0))
infile = args.pop(0)
outfile = args.pop(0)

nbins = int(math.ceil(float(nbeads) / float(binsize)))
state = [0 for i in xrange(nbins)]

with open(infile, "r") as f:
    for line in f:
        if (len(line) == 0 or line.startswith("#")):
            continue
        data = line.strip().split()
        if (len(data) == 0):
            continue
        index = int(data[0])
        value = int(data[1])
        bin_index = index / binsize
        if (value == 2):
            state[bin_index] += 1

# Normalise
for i in xrange(nbins):
    state[i] /= float(binsize)
    if (state[i] >= 0.5):
        state[i] = 2
    else:
        state[i] = 1

# Output
with open(outfile, "w") as writer:
    for i in xrange(nbins):
        writer.write("{:d} {:d}\n".format(i, state[i]))
