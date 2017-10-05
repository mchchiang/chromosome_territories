
// Bond.cpp

#include "Bond.h"

const int Bond::numOfBeads {2};

// Constructor
Bond::Bond(int t, shared_ptr<Bead> b1, shared_ptr<Bead> b2) :
  type {t}, beads {b1, b2} {
	b1->addBond(shared_from_this());
	b2->addBond(shared_from_this());
}

// Accessor methods
shared_ptr<Bead> Bond::getBead(int id){
  return beads[id];
}

void Bond::setType(int t){
  type = t;
}

int Bond::getType(){
  return type;
}

void Bond::unbond(){
  for (int i {}; i < numOfBeads; i++){
	beads[i]->removeBond(shared_from_this());
  }
}
