# distance.py 

import sys
import math

args = sys.argv
args.pop(0) # Ignore self

nbeads = int(args.pop(0))
lx = float(args.pop(0))
ly = float(args.pop(0))
lz = float(args.pop(0))
outfile = args.pop(0)

# compute distance
def dist2(a, b):
    dx = a[0] - b[0]
    dy = a[1] - b[1]
    dz = a[2] - b[2]
    return dx*dx+dy*dy+dz*dz

pos = [[0.0 for j in xrange(3)] for i in xrange(nbeads)]
r2_avg = [0.0 for i in xrange(nbeads)]
#rg2_avg = [0.0 for i in xrange(nbeads)]
r2_count = [0 for i in xrange(nbeads)]

# read x,y,z
for infile in args:
    with open(infile, "r") as f:
        for line in f:
            data = line.strip().split();
            if (len(data) == 9):
                index = int(data[0])
                x = float(data[3])
                y = float(data[4])
                z = float(data[5])
                ix = int(data[6])
                iy = int(data[7])
                iz = int(data[8])
                x += (lx*ix)
                y += (ly*iy)
                z += (lz*iz)
                pos[index-1][0] = x
                pos[index-1][1] = y
                pos[index-1][2] = z

    print "Computing distances for ", infile

    for i in xrange(nbeads):
        cm = [0.0, 0.0, 0.0]
        for j in xrange(i,nbeads):
            sep = abs(j-i)
            r2 = dist2(pos[i],pos[j])
            r2_avg[sep] += r2
            r2_count[sep] += 1

            # Compute cm
"""            cm = [0.0, 0.0, 0.0]
            count = j-i+1
            
            for k in xrange(i,j+1):
                cm[0] += pos[i][0]
                cm[1] += pos[i][1]
                cm[2] += pos[i][2]
            cm[0] /= float(count)
            cm[1] /= float(count)
            cm[2] /= float(count)
            rg2 = 0.0
            for k in xrange(i,j+1):
                rg2 += dist2(pos[k],cm)
            rg2 /= float(count)
            rg2_avg[sep] += rg2
"""                    

# Normalise
for i in xrange(nbeads):
    r2_avg[i] /= float(r2_count[i])
#    rg2_avg[i] /= float(r2_count[i])

with open(outfile, "w") as w:
    for i in xrange(nbeads):
        w.write("{:d} {:.5f}\n".format(i, r2_avg[i]))
#        w.write("{:d} {:.5f} {:.5f}\n".format(i, r2_avg[i], rg2_avg[i]))
