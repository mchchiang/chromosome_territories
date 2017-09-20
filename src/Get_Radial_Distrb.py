# Get_Radial_Distrb.py
# Obtain the radial distribution of each chromosome ball

import sys
import math

args = sys.argv
args.pop(0) # Ignore self

outFile = args.pop(0)

numOfChromos = 46
idCol = 0
xCol = 3
yCol = 4
zCol = 5
nxCol = 6
nyCol = 7
nzCol = 8
boxSize = 200.0

# For storing distances for each ball
numOfBins = 10
maxBinNum = numOfBins-1
minBinNum = 0
min = 0.0
max = 70.0
binWidth = (max - min) / numOfBins
distrb = [[0 for j in xrange(numOfBins)] for i in xrange(numOfChromos)]

fileCount = 0

for file in args:
    print "Analysing " + file
    
    lineCount = 0
    chromoCount = 0
    
    # Skip the first 29 lines (header stuff)
    with open(file, "r") as f:
        for line in f:
            if (lineCount < 29):
                lineCount += 1
                continue

            data = line.strip().split()
            id = int(data[idCol])

            if (id <= numOfChromos): # Is a chromosome
                x = float(data[xCol]) + float(data[nxCol])*boxSize
                y = float(data[yCol]) + float(data[nyCol])*boxSize
                z = float(data[zCol]) + float(data[nzCol])*boxSize
                r = math.sqrt(x*x+y*y+z*z)
                chromoCount += 1
                
                # Bin the data - find bin number
                binNum = int(r / binWidth)
                if (binNum > maxBinNum):
                    print ("binNum = %d > maxBinNum = %d") % (binNum, maxBinNum)
                    binNum = maxBinNum

                elif (binNum < minBinNum):
                    print ("binNum = %d < minBinNum = %d") % (binNum, minBinNum)
                    binNum = minBinNum
                    
                distrb[id-1][binNum] += 1
                
            if (chromoCount == numOfChromos):
                break

            lineCount += 1
        
    fileCount += 1

# Normalise distribution

print fileCount

for i in xrange(numOfChromos):
    for j in xrange(numOfBins):
        distrb[i][j] / float(fileCount)

# Output distribution
writer = open(outFile, "w")

for j in xrange(numOfBins):
    r = binWidth * (j + 0.5)
    output = "%.3f" % r
    for i in xrange(numOfChromos):
        output += " %.5f" % (distrb[i][j])
    output += "\n"
    writer.write(output)

writer.close()
