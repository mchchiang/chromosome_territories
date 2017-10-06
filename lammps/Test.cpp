// Test.cpp

#include <iostream>
#include <memory>
#include "LAMMPS.h"
#include "Polymer.h"
#include "Bead.h"

using std::cout;
using std::endl;
using std::shared_ptr;
using std::make_shared;

int main(){
  double lx {100.0};
  double ly {100.0};
  double lz {100.0};
  shared_ptr<LAMMPS> lammps = make_shared<LAMMPS>(lx, ly, lz);
  shared_ptr<Polymer> p1 = 
	Polymer::createRandomWalkPolymer(1000, lx, ly, lz);
  shared_ptr<Polymer> p2 = 
	Polymer::createRandomWalkPolymer(1000, lx, ly, lz);
  shared_ptr<Bead> b1 = make_shared<Bead>(20,30,40);
  lammps->setTypesOfBeads(1);
  lammps->setTypesOfBonds(1);
  lammps->setTypesOfAngles(1);
  lammps->addPolymer(0,p1);
  lammps->addPolymer(1,p2);
  lammps->removePolymer(0);
  lammps->addBead(0,b1);
  shared_ptr<Bead> b2 = lammps->getPolymer(0)->getBead(324);
  shared_ptr<Bead> b3 = lammps->getPolymer(0)->getBead(483);
  b1->addBondWith(0, b2);
  b1->addAngleWith(0, b2, b3);
  b2->removeBondWith(b1);
  lammps->getPolymer(0)->removeBead(765);
  lammps->exportData("mydata.dat", "map.dat");
}
