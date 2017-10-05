// Bond.h

#ifndef BOND_H
#define BOND_H

#include <memory>
#include "Bead.h"

using std::shared_ptr;

class Bead;

class Bond : std::enable_shared_from_this<Bond> {

private:
  static const int numOfBeads;
  int type;
  shared_ptr<Bead> beads[2];

public:
  
  // Constructor
  Bond(int type, shared_ptr<Bead> bead1, shared_ptr<Bead> bead2);

  // Accessor methods
  shared_ptr<Bead> getBead(int id);
  void setType(int type);
  int getType();
  
  void unbond();
};

#endif
