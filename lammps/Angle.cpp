// Angle.cpp

#include <memory>
#include "Bead.hpp"
#include "Angle.hpp"

using std::shared_ptr;

// Constructor
Angle::Angle(int t, shared_ptr<Bead> b1, 
	     shared_ptr<Bead> b2, shared_ptr<Bead> b3) :
  type {t}, beads {b1, b2, b3} {}

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
