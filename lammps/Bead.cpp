// Bead.cpp

#include <map>
#include <memory>
#include "Bead.h"

using std::map;
using std::pair;
using std::shared_ptr;
using Bond = Bead::Bond;
using Angle = Bead::Angle;
using BondMap = map< shared_ptr<Bond>, shared_ptr<Bond> >;
using AngleMap = map< shared_ptr<Angle>, shared_ptr<Angle> >;

// Constructors
Bead::Bead(double x, double y, double z,
		   double vx, double vy, double vz,
		   int nx, int ny, int nz, int t, int l) :
  position {x, y, z}, velocity {vx, vy, vz},
  boundaryCount {nx, ny, nz}, type {t}, label {l} {
}

Bead::Bead(double x, double y, double z) :
  Bead {x, y, z, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0} {}

Bead::Bead() :
  Bead {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0} {}


// Accessor methods
void Bead::setPosition(int dim, double value){
  position[dim] = value;
}

double Bead::getPosition(int dim){
  return position[dim];
}

void Bead::setVelocity(int dim, double value){
  velocity[dim] = value;
}

double Bead::getVelocity(int dim){
  return velocity[dim];
}

void Bead::setBoundaryCount(int dim, int count){
  boundaryCount[dim] = count;
}

int Bead::getBoundaryCount(int dim){
  return boundaryCount[dim];
}

void Bead::setType(int t){
  type = t;
}

int Bead::getType(){
  return type;
}

void Bead::setLabel(int l){
  label = l;
}

int Bead::getLabel(){
  return label;
}

void Bead::addBond(shared_ptr<Bond> bond){
  bondList[bond] = bond;
}

void Bead::addBondWith(int t, shared_ptr<Bead> bead){
  Bond(t, shared_from_this(), bead);
}

void Bead::removeBond(shared_ptr<Bond> bond){
  bondList.erase(bond);
}

void Bead::removeBondWith(shared_ptr<Bead> bead){
  for (BondMap::iterator it = bondList.begin(); it != bondList.end(); it++){
	shared_ptr<Bond> bond = it->second;
	if (bond->getBead(0) == bead || bond->getBead(1) == bead){
	  bead->removeBond(bond);
	  bondList.erase(it);
	  break;
	}
  }
}

void Bead::addAngle(shared_ptr<Angle> angle){
  angleList[angle] = angle;
}

void Bead::addAngleWith(int t, 
						shared_ptr<Bead> bead1, shared_ptr<Bead> bead2){
  Angle(t, shared_from_this(), bead1, bead2);
}

void Bead::removeAngle(shared_ptr<Angle> angle){
  angleList.erase(angle);
}


void Bead::removeAngleWith(shared_ptr<Bead> bead1, shared_ptr<Bead> bead2){
  for (AngleMap::iterator it = angleList.begin(); it != angleList.end(); it++){
	shared_ptr<Angle> angle = it->second;
	shared_ptr<Bead> b1 = angle->getBead(0);
	shared_ptr<Bead> b2 = angle->getBead(1);
	shared_ptr<Bead> b3 = angle->getBead(2);
	if ((bead1 == b1 && (bead2 == b2 || bead2 == b3)) ||
		(bead1 == b2 && (bead2 == b1 || bead2 == b3)) ||
		(bead1 == b3 && (bead2 == b1 || bead2 == b2))){
	  bead1->removeAngle(angle);
	  bead2->removeAngle(angle);
	  angleList.erase(it);
	  break;
	}
  }
}

void Bead::removeAllBonds(){
  for (BondMap::iterator it = bondList.begin(); it != bondList.end(); it++){
	shared_ptr<Bond> bond = it->second;
	if (bond->getBead(0) != shared_from_this())
	  bond->getBead(0)->removeBond(bond);
	else
	  bond->getBead(1)->removeBond(bond);
  }
  bondList.clear();
}

void Bead::removeAllAngles(){
  for (AngleMap::iterator it = angleList.begin(); it != angleList.end(); it++){
	shared_ptr<Angle> angle = it->second;
	if (angle->getBead(0) != shared_from_this()){
	  angle->getBead(1)->removeAngle(angle);
	  angle->getBead(2)->removeAngle(angle);
	} else if (angle->getBead(1) != shared_from_this()){
	  angle->getBead(0)->removeAngle(angle);
	  angle->getBead(2)->removeAngle(angle);
	} else {
	  angle->getBead(0)->removeAngle(angle);
	  angle->getBead(1)->removeAngle(angle);
	}
  }
  angleList.clear();
}

//=================================================================

// Methods for Bond class

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

//=================================================================

// Methods for Angle class

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
