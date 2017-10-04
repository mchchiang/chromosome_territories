// Bond.h

#ifndef BOND_H
#define BOND_H

#include <memory>
#include "Bead.h"

using std::shared_ptr;

class Bond {

private:
  shared_ptr<Bead>[2];
  Bond(); // Private constructor - must use static method to create a bond

public:
  shared_ptr<Bead> getBead(int id);
  void unbond();
  static void createBond(int type, 
			 shared_ptr<Bead> bead1, shared_ptr<Bead> bead2);
  static void removeBond(shared_ptr<Bead> bead1, shared_ptr<Bead> bead2);
};

#endif
