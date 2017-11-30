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
Polymer::Polymer(int nBeads, int beadType, 
		 int bondType, int angleType, bool createBead){

  if (nBeads > 0){
    beads.reserve(nBeads);

    if (createBead){

      int i {};
      shared_ptr<Bead> bead {};

      bead = make_shared<Bead>();
      bead->setType(beadType);
      beads.push_back(bead);
      i++;
	
      if (i < nBeads){
	bead = make_shared<Bead>();
	bead->setType(beadType);
	beads.push_back(bead);
	beads[i-1]->addBondWith(bondType, beads[i]);
	i++;
      }
	
      for (; i < nBeads; i++){
	bead = make_shared<Bead>();
	bead->setType(beadType);
	beads.push_back(bead);
	beads[i-1]->addBondWith(bondType, beads[i]);
	beads[i-2]->addAngleWith(angleType, beads[i-1], beads[i]);
      }
    }
  }
}

Polymer::Polymer() : Polymer {0, false} {} 
Polymer::Polymer(int nBeads) : Polymer {nBeads, true} {}
Polymer::Polymer(int nBeads, bool createBead) : 
  Polymer {nBeads, 1, 1, 1, createBead} {}
Polymer::Polymer(int nBeads, int beadType) : 
  Polymer {nBeads, beadType, 1, 1} {}


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
void Polymer::addBead(shared_ptr<Bead> bead){
  beads.push_back(bead);
}

void Polymer::addBead(int id, shared_ptr<Bead> bead){
  beads.insert(beads.begin()+id, bead);
}

void Polymer::removeBead(int id){
  // Find the bond and angle type with neighbouring bead
  int bondType {1};
  int angleType {1};
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

// Statistics of polymer
vector<double> Polymer::getCentreOfMass(double lx, double ly, double lz){
  double x {}, y {}, z {};
  for (auto const& b : beads){
    x += (b->getPosition(0) + lx*b->getBoundaryCount(0));
    y += (b->getPosition(1) + ly*b->getBoundaryCount(1));
    z += (b->getPosition(2) + lz*b->getBoundaryCount(2));
  }
  double numOfBeads = beads.size();
  x /= numOfBeads;
  y /= numOfBeads;
  z /= numOfBeads;
  return {x, y, z};
}

double Polymer::getGyrationRadius(double lx, double ly, double lz){
  vector<double> cm = getCentreOfMass(lx, ly, lz);
  double dx {}, dy {}, dz {}, sum {};
  for (auto const& b : beads){
    dx = b->getPosition(0) + lx*b->getBoundaryCount(0) - cm[0];
    dy = b->getPosition(1) + ly*b->getBoundaryCount(1) - cm[1];
    dz = b->getPosition(2) + lz*b->getBoundaryCount(2) - cm[2];
    sum += dx*dx+dy*dy+dz*dz;
  }
  sum /= beads.size();
  sum = sqrt(sum);
  return sum;
}

// Static factory methods for creating polymers
shared_ptr<Polymer> Polymer::createRandomWalkPolymer(int nBeads, int beadType,
						     double x0, double y0, 
						     double z0, double lx, 
						     double ly, double lz){
  // Initialise random number generator
  srand(time(NULL));
  double pi {M_PI};
  shared_ptr<Polymer> polymer = make_shared<Polymer>(nBeads, beadType);
  double x, y, z, r, costheta, sintheta, phi;
  shared_ptr<Bead> previous {};
  shared_ptr<Bead> current {};

  // Set the first bead to be centred at the origin
  previous = polymer->getBead(0);
  previous->setPosition(0, x0);
  previous->setPosition(1, y0);
  previous->setPosition(2, z0);

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
    } while (fabs(x-x0) > lx/2.0 || fabs(y-y0) > ly/2.0 || fabs(z-z0) > lz/2.0);
    current->setPosition(0, x);
    current->setPosition(1, y);
    current->setPosition(2, z);
	previous = current;
  }
  return polymer;
}
