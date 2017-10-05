// Polymer.h

#ifndef POLYMER_H
#define POLYMER_H

#include <vector>
#include <memory>
#include "Bead.h"

using std::vector;
using std::shared_ptr;

class Polymer {

private:
  vector< shared_ptr<Bead> > beads {};

public:
  // Constructor
  Polymer(); 
  Polymer(int numOfBeads);

  // Accessor methods
  shared_ptr<Bead> getBead(int id);
  vector< shared_ptr<Bead> >& getBeads();

  int getNumOfBeads();

  // For changing the polymer
  void addBead(int id, shared_ptr<Bead> bead);
  void removeBead(int id);
  
  // Static factory methods
  static shared_ptr<Polymer> createRandomWalkPolymer(int numOfBeads,
						     double lx, double ly, 
						     double lz);

};

#endif
