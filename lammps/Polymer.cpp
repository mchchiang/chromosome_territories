// Polymer.cpp

#include <iostream>
#include <vector>
#include <memory>
#include <cstdlib>
#include <cmath>
#include "Bead.h"
#include "Polymer.h"

using std::cout;
using std::endl;
using std::vector;
using std::shared_ptr;
using std::make_shared;

// Constructors
Polymer::Polymer(int nBeads){

  beads.reserve(nBeads);

  int i {};
  int type {};
  
  beads.push_back( make_shared<Bead>());
  i++;

  if (i < nBeads){
    beads.push_back(make_shared<Bead>());
    beads[i-1]->addBondWith(type, beads[i]);
	i++;
  }

  for (; i < nBeads; i++){
    beads.push_back(make_shared<Bead>());
    beads[i-1]->addBondWith(type, beads[i]);
    beads[i-2]->addAngleWith(type, beads[i-1], beads[i]);
  }
}

Polymer::Polymer() : Polymer {0} {}

// Accessor methods
shared_ptr<Bead> Polymer::getBead(int id){
  return beads[id];
}

vector< shared_ptr<Bead> >& Polymer::getBeads(){
  return beads;
}

int Polymer::getNumOfBeads(){
  return beads.size();
}

// Adding or removing beads
void Polymer::addBead(int id, shared_ptr<Bead> bead){
  beads.insert(beads.begin()+id, bead);
}

void Polymer::removeBead(int id){
  // Find the bond and angle type with neighbouring bead
  int bondType {};
  int angleType {};
  shared_ptr<Bead::Bond> bond {}; 
  shared_ptr<Bead::Angle> angle {};
  int bead1Index {id+1};
  int bead2Index {id+2};

  bond = beads[id]->getBondWith(beads[bead1Index]);
  angle = beads[id]->getAngleWith(beads[bead1Index], beads[bead2Index]);

  if (bond != nullptr)
	bondType = bond->getType();
  if (angle != nullptr)
	angleType = angle->getType();

  // Remove all bonds and angles
  beads[id]->removeAllBonds();
  beads[id]->removeAllAngles();
  
  // Make sure the polymer remains connected
  beads[id-1]->addBondWith(bondType, beads[id+1]);
  beads[id-2]->addAngleWith(angleType, beads[id-1], beads[id+1]);
  beads[id-1]->addAngleWith(angleType, beads[id+1], beads[id+2]);

  // Erase the current bead
  beads.erase(beads.begin()+id);
}

void Polymer::removeAllBeads(){
  for (auto const& bead : beads){
	bead->removeAllBonds();
	bead->removeAllAngles();
  }
  beads.clear();
}

// Static factory methods for creating polymers
shared_ptr<Polymer> Polymer::createRandomWalkPolymer(int nBeads,
						     double lx, double ly,
						     double lz){
  // Initialise random number generator
  srand(time(NULL));
  double pi {M_PI};
  shared_ptr<Polymer> polymer = make_shared<Polymer>(nBeads);
  double x, y, z, r, costheta, sintheta, phi;
  shared_ptr<Bead> previous {};
  shared_ptr<Bead> current {};

  // Set the first bead to be centred at the origin
  previous = polymer->getBead(0);
  for (int i {}; i < 3; i++){
    previous->setPosition(i, 0.0);
  }

  for (int i {1}; i < nBeads; i++){
    current = polymer->getBead(i);
    do {
      r = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
      costheta = 1.0-2.0*r;
      sintheta = sqrt(1-costheta*costheta);
      r = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
      phi = 2.0*pi*r;
      x = previous->getPosition(0) + sintheta * cos(phi);
      y = previous->getPosition(1) + sintheta * sin(phi);
      z = previous->getPosition(2) + costheta;
    } while (fabs(x) > lx/2.0 || fabs(y) > ly/2.0 || fabs(z) > lz/2.0);
    current->setPosition(0, x);
    current->setPosition(1, y);
    current->setPosition(2, z);
	previous = current;
  }
  return polymer;
}
