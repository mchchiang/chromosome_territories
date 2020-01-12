/* distance.cpp
 * A code to calculate the end-to-end distance using the 
 * final lammps output trajectory file
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <memory>
#include <cmath>
#include "Bead.hpp"
#include "Polymer.hpp"
#include "LAMMPS.hpp"

using std::cout;
using std::endl;
using std::string;
using std::ofstream;
using std::shared_ptr;
using std::make_shared;

double dist2(double x, double y, double z);

int main(int argc, char* argv[]) {
  
  if (argc < 6) {
    cout << "Usage: distance [nbeads] [polymerid] " 
	 << "[outfile] [mapfile] [lammpsfile]" << endl;
    return 1;
  }

  int argi {};
  int nbeads {stoi(string(argv[++argi]), nullptr, 10)};
  int polymerid {stoi(string(argv[++argi]), nullptr, 10)};
  string outfile {argv[++argi]};
  string mapfile {argv[++argi]};
  string lammpsfile;
  
  shared_ptr<LAMMPS> lammps = make_shared<LAMMPS>();

  shared_ptr<Bead> bead1, bead2;
  shared_ptr<Polymer> polymer;
  vector<double> r2avg (nbeads, 0.0);
  vector<int> count (nbeads, 0);
  argi++;

  for (; argi < argc; argi++) {
    lammpsfile = string(argv[argi]);
    lammps->importData(lammpsfile, mapfile);
    polymer = lammps->getPolymer(polymerid);

    double lx {lammps->getLx()};
    double ly {lammps->getLy()};
    double lz {lammps->getLz()};

    // Compute distances
    int sep;
    double r2;
    double x, y, z, dx, dy, dz;
    for (int j {}; j < nbeads; j++) {
      bead1 = polymer->getBead(j);
      x = bead1->getPosition(0)+(bead1->getBoundaryCount(0)*lx);
      y = bead1->getPosition(1)+(bead1->getBoundaryCount(1)*ly);
      z = bead1->getPosition(2)+(bead1->getBoundaryCount(2)*lz);
      for (int k {}; k <= j; k++) {
	bead2 = polymer->getBead(k);
	sep = abs(j-k);
	dx = x-(bead2->getPosition(0)+(bead2->getBoundaryCount(0)*lx));
	dy = y-(bead2->getPosition(1)+(bead2->getBoundaryCount(1)*ly));
	dz = z-(bead2->getPosition(2)+(bead2->getBoundaryCount(2)*lz));
	r2 = dist2(dx,dy,dz);
	r2avg[sep] += r2;
	count[sep]++;
      }
    }
  }

  // Normalise results
  for (int i {}; i < nbeads; i++) {
    r2avg[i] /= static_cast<double>(count[i]);
  }
  cout << "Outputting results ..." << endl;
  // Output results
  ofstream writer;
  writer.open(outfile);
  if (!writer) {
    cout << "Problem with opening output file" << endl;
    return 1;
  }
  writer << std::fixed << std::setprecision(5);
  for (int i {}; i < nbeads; i++) {
    writer << i << " " << r2avg[i] << endl;
  }
  writer.close();
}

double dist2(double x, double y, double z) {
  return x*x+y*y+z*z;
}
