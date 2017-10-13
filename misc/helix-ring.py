import sys
from math import *

phi_min = 0.0
phi_max = 4.0*pi
phi_inc = 0.1
phi = phi_min

a = 0.38
k = 6.0

r_chr = 12.0
p = 1.0

writer = open("helix.dat", "w")

while (phi < phi_max):
    x = r_chr*(a+(1-a)*(cos(k*phi))**2*cos(phi))
    y = r_chr*(a+(1-a)*(cos(k*phi))**2*sin(phi))
    z = p*phi/(2*pi)
    writer.write("%.5f %.5f %.5f\n" % (x, y, z))

    phi += phi_inc

writer.close()

writer = open("ring.dat", "w")


r_t = 42.0
T = 378.0
phi = phi_min
phi_max = 2.0*pi*T

while (phi < phi_max):
    x = r_t/(r_chr+r_t)*(r_chr*(a+(1-a)*cos(k*phi)**2)*cos(phi)+r_t*cos(phi/T))
    y = r_chr*(a+(1-a)*cos(k*phi)**2*sin(phi))
    z = r_t*sin(phi/T)
    writer.write("%.5f %.5f %.5f\n" % (x, y, z))

    phi += phi_inc

writer.close()
