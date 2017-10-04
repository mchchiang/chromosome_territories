// Polymer.h

#ifndef POLYMER_H
#define POLYMER_H

#include <vector>
#include <memory>

using std::vector;
using std::shared_ptr;

class Polymer {

private:
  Polymer(); // Private constructor - use static method to create a polymer
  vector< shared_ptr<Bead> > beads;

public:

  // Accessor methods
  shared_ptr<Bead> getBead(int id);
  vector< shared_ptr<Bead> >& getBeads();

  int getNumOfBeads();

  // For changing the polymer
  int addBead(shared_ptr<Bead> bead);
  void removeBead(int id);
  
  int getNumOfBeads();
  
  // Static factory methods
  static shared_ptr<Polymer> createRandomWalkPolymer(int numOfBeads);

};

#endif
