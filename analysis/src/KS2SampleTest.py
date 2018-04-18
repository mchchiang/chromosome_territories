# KS2SampleTest.py
#
# Perform the Kolmogorov-Smmirnov two sample test
#

import sys
import scipy
import scipy.stats
import math

def read_array(filename, value_col):
    values= []
    with open(filename, 'r') as f:
        for line in f:
            if (line.startswith("#")): continue
            data = line.strip().split()
            if (data == []): continue
            values.append(float(data[value_col]))
    return scipy.array(values)


args = sys.argv
args.pop(0) # Ignore self

if (len(args) < 3):
    print "Usage: [value_col] [array1_file] [array2_file] [out_file]"
    sys.exit(1)

value_col = int(args.pop(0))
array1_file = args.pop(0)
array2_file = args.pop(0)
out_file = args.pop(0)

# Read in arrays 
array1 = read_array(array1_file, value_col)
array2 = read_array(array2_file, value_col)

array1 = scipy.nan_to_num(array1)
array2 = scipy.nan_to_num(array2)

# Compute KS 2-Sample Test
result = scipy.stats.ks_2samp(array1, array2)

with open(out_file, 'w') as writer:
    writer.write("{:.5g} {:.5g}\n".format(*result))


