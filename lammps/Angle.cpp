// Angle.cpp

#include "Angle.h"

const int Angle::numOfBeads {3};

// Constructor
Angle::Angle(int t, shared_ptr<Bead> b1, 
			 shared_ptr<Bead> b2, shared_ptr<Bead> b3) :
  type {t}, beads {b1, b2, b3} {
	b1->addAngle(shared_from_this());
	b2->addAngle(shared_from_this());
	b3->addAngle(shared_from_this());
}

// Accessor methods
shared_ptr<Bead> Angle::getBead(int id){
  return beads[id];
}

void Angle::setType(int t){
  type = t;
}

int Angle::getType(){
  return type;
}

void Angle::unbond(){
  for (int i {}; i < numOfBeads; i++){
	beads[i]->removeAngle(shared_from_this());
  }
}

