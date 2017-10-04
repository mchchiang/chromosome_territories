// Angle.h

#ifndef ANGLE_H
#define ANGLE_H

#include <memory>
#include "Bead.h"

using std::shared_ptr;

class Angle {

private:
  shared_ptr<Bead>[3];
  Angle(); // Private constructor - must use static method to create a bond

public:
  shared_ptr<Bead> getBead(int id);
  void unbond();
  static void createAngle(int type, shared_ptr<Bead> bead1, 
			 shared_ptr<Bead> bead2, shared_ptr<Bead> bead3);
  static void removeAngle(shared_ptr<Bead> bead1, shared_ptr<Bead> bead2,
			  shared_ptr<Bead> bead3);
};

#endif
