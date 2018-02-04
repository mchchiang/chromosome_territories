# VanillaNorm.py
#
# A program that reads in a contact matrix in the format
# (bin_i bin_j count_ij) (full or upper triangle) and
# performs the square root vanilla normalisation on the matrix.
# 
# The square root vanilla normalisation procedure is as follows:
# M_ij* = M_ij/sqrt(sum_k M_kj * sum_l M_il)
# For references, see Rao et al. Cell (1509) 1665-1680, 2014
#

import sys
import math

args = sys.argv
args.pop(0) # Ignore self

# Check there is enough number of arguments
if (len(args) != 5):
    print "Usage: VanillaNorm.py [nx] [ny] " \
          "[mode=full,upper] [input] [output]"
    sys.exit(1)

nx = int(args.pop(0)) # Number of columns in matrix
ny = int(args.pop(0)) # Number of rows in matrix
mode = args.pop(0) # Either full or upper
matrix_file = args.pop(0) # Input contact matrix
norm_file = args.pop(0) # Output normalised matrix

# Read the contact matrix
contact = [[0.0 for j in xrange(ny)] for i in xrange(nx)]

with open(matrix_file, 'r') as f:
    for line in f:
        # Igore any empty lines or lines beginning with '#'
        if (line.startswith('#')):
            continue
        data = line.strip().split()
        if (len(data) != 3):
            continue
        x = int(data[0])
        y = int(data[1])
        count = float(data[2])
        
        if (x < 0 or x >= nx):
            print "Index x out of range: %d" % x
        elif (y < 0 or y >= ny):
            print "Index y out of range: %d" % y
        elif (mode == "full"):
            contact[x][y] = count
        else:
            contact[x][y] = count
            contact[y][x] = count

# Compute row and column sums
x_sum = [0.0 for i in xrange(nx)]
y_sum = [0.0 for i in xrange(ny)]

for i in xrange(nx):
    for j in xrange(ny):
        x_sum[i] += contact[i][j]
        y_sum[j] += contact[i][j]

# Perform vanilla normalisation
for i in xrange(nx):
    for j in xrange(ny):
        if (x_sum[i] == 0 or y_sum[j] == 0):
            contact[i][j] = 0
        else:
            contact[i][j] /= math.sqrt(x_sum[i] * y_sum[j])

# Output result
with open(norm_file, 'w') as f:
    for i in xrange(nx):
        for j in xrange(ny):
            f.write("%d %d %.10f\n" % (i, j, contact[i][j]))
        f.write("\n")
