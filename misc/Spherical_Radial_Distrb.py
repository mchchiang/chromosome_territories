# Get_Radial_Distrb.py
# Obtain the radial distribution of each chromosome ball

import sys
import math

args = sys.argv
args.pop(0) # Ignore self

outFile = args.pop(0)

xCol = 0
yCol = 1
zCol = 2

# For storing distances for each ball
numOfBins = 50
maxBinNum = numOfBins-1
minBinNum = 0
min = 0.0
max = 10.0
binWidth = (max - min) / numOfBins
distrb = [0 for i in xrange(numOfBins)] 

count = 0

for file in args:
    print "Analysing " + file
    
    # Skip the first 29 lines (header stuff)
    with open(file, "r") as f:
        for line in f:
            data = line.strip().split()
            x = float(data[xCol])
            y = float(data[yCol])
            z = float(data[zCol])
            r = math.sqrt(x*x+y*y+z*z)
                
            # Bin the data - find bin number
            binNum = int(r / binWidth)
            if (binNum > maxBinNum):
                print ("binNum = %d > maxBinNum = %d") % (binNum, maxBinNum)
                binNum = maxBinNum
                
            elif (binNum < minBinNum):
                print ("binNum = %d < minBinNum = %d") % (binNum, minBinNum)
                binNum = minBinNum
                
            distrb[binNum] += 1
                
            count += 1

# Normalise distribution
for i in xrange(numOfBins):
    distrb[i] /= float(count)

# Output distribution
writer = open(outFile, "w")

for i in xrange(numOfBins):
    r = binWidth * (i + 0.5)
    output = "%.3f %.5f" % (r, distrb[i])
    output += "\n"
    writer.write(output)

writer.close()
