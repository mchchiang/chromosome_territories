DIPLOID NUCLEUS SIMULATIONS 

STEP 1.
We generate an input file for LAMMPS with 46 spheres - one for chromosomes.
The file is called LammpsInput.46Chr.1BallEach;
the program is Gen_46Chr_1BallEach.c++

Each sphere has a "type" (from 1 to 9) associated with the total LAD
fraction (read from the file "LAD.Pk.genome.full.dat"). In the genome browser 
LADs are stretches of chromatin often found associated with the lamina, and 
extracted from Tig3 cells. Tig3 cells are human diploid embryonic lung 
fibroblasts.

STEP 2.
We generate equilibrated positions of the 46 spheres modelling chromosomes.
To do this, the LAMMPS input script
Run.AllChr.1Ball.onlyLAM.lam
is used. This takes the output of STEP 1 as a start file; then it:
(i) performs an initial run to separate out spheres (with soft potential);
(ii) creates a shell and fills it randomly with beads of type 10 (these
are models for sticky regions for LADs); the interactions between spheres
(chromosomes) and shell (lamina) are set in such a way that higher types
are attracted to beads in the shell more strongly. 
The LAMMPS simulations involves random motion of the chromosome spheres with
repeated warm-up/cool-down to create uncorrelated snapshots.

As an output we get 10 uncorrelated snapshots of the 46 spheres, 
Eq.46Chr1Ball.LAMonly.X; with X=1,...,10
to check these with gnuplot, e.g.
sp'Eq.46Chr1Ball.LAMonly.2' u 4:5:6:3 palette

3)Gen_46Chr_fromBalls.c++
This file reads the position of the spheres in "Eq.46Chr1Ball.LAMonly.X"
(see code) and builds some cylinder-like configurations (approximations 
for the mitotic structure). It also occasionally checks to see whether 
it is doing overlap with other mitotic-like cylinders. The cylinder beads 
have the epigenetic which they read from genome ENCODE files. 
NB: "H3K4me3.Pk..dat" = this is the histone mark for promoters
"PCG.Pk...dat" = this is the histone mark for polycomb (H3K27me3)
these marks are for GM12878 cells - initially we shall assume they are
OK for Tig3 cells (ideally IMR90 would be a better choice though)
"LAD.Pk...dat" files give the positions of the LADs as used in STEP 1
The output of STEP 3 is a file called "LammpsInput.46Chr" to be plotted with
gnuplot as follows,
sp "LammpsInput.22Chr" u 4:5:6:3 palette -- to see the epigenetics
sp "LammpsInput.22Chr" u 4:5:6:2 palette -- to see the chromosomes

4) Now the following task (relaxation of mitotic cylinders) are split
in parts as follows:

4a) Run.AllChrFromBalls.LAM.lam1: do a warm up run of 1000 steps
and generates a "WarmedUP" file in the directory "WarmUP". 

4b) Now have to remove (by hand) the lines corresponding to pair&bond&angle
coefficients in the WarmedUP.. data file.

4c)Run.AllChrFromBalls.LAM.lam2: Takes the "warmedup" file and uses an indenter
to decrease confinement radius to 50\sigma. This increases volume fraction to
desired value of Nbeads*0.5^3/R^3 \simeq 0.3. It then creates the outer membrane again + some short relaxation. One can plot & visualise territories by using
sp "WarmUPwLAM/WarmedUPwLAM.22ChrFromBalls.1" u 4:5:6:($3<10?$2:1/0) palette

4d)Run.AllChrFromBalls.LAM.lam3: Takes the file
"WarmUPwLAM/WarmedUPwLAM.22ChrFromBalls.1" and does the long relaxation dynamics
+ Polymerase (through external file WritePolymerase.c++). The external file
finds looped promoters and recolours a region around them.
To see epigenetics: p "newdataRestart.22Chr.Dt6000E1.0" i 0 u 1:3:3 palette

The first thing is to see whether the distribution of chromosomes from
STEP 1 is reasonable (wrt the Cremers paper).
