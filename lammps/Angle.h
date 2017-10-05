// Angle.h

#ifndef ANGLE_H
#define ANGLE_H

#include <memory>
#include "Bead.h"

using std::shared_ptr;

class Angle : std::enable_shared_from_this<Angle> {

private:
  static const int numOfBeads;
  int type;
  shared_ptr<Bead> beads[3];

public:

  // Constructor
  Angle(int type, shared_ptr<Bead> bead1, 
		shared_ptr<Bead> bead2, shared_ptr<Bead> bead3);

  // Accessor methods
  shared_ptr<Bead> getBead(int id);
  void setType(int type);
  int getType();

  void unbond();
};

#endif
