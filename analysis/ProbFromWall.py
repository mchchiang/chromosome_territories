# ProbFromWall.py
# Compute the probability of each bead from z wall

import sys
import math

args = sys.argv

args.pop(0) # ignore self

z_wall_pos = float(args.pop(0))
min_val = float(args.pop(0))
max_val = float(args.pop(0))
binsize = float(args.pop(0))
output_file = args.pop(0)

value_col = 2 # Measuring the z position
max_col = 2

files = [open(i, "r") for i in args]

n = int(math.ceil((max_val - min_val) / binsize))
distrb = [0 for i in xrange(n)]
count = 0

for f in files:
    for line in f:
        if (not line.startswith("#")):
            data = line.strip().split()

            # Skip blank lines or if the line does not contain enough columns
            if (data == [] or len(data) < (max_col+1)):
                continue
            
            value = abs(z_wall_pos - float(data[value_col]))
            col = int(math.floor((value - min_val) / binsize))
            if (col < 0):
                print "The value %.5f is less than min = %.5f" \
                    % (value, min_val)
            elif (col >= n):
                print "The value %.5f is greater than or equal to max = %.5f" \
                    % (value, max_val)
            else:
                distrb[col] += 1
                count += 1
    f.close()

for i in xrange(n):
    distrb[i] /= float(count)

writer = open(output_file, "w")

for i in xrange(n):
    left = i * binsize + min_val
    right = (i+1) * binsize + min_val
    centre = (right - left) / 2.0 + left
    output = "%.5f %.5f %.5f %.5f\n" % (left, right, centre, distrb[i])
    writer.write(output)

writer.close()
