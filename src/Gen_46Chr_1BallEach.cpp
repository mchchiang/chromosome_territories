/* Gen_46Chr_1BallEach.cpp
 * This code generates a bead to model each chromosome of 
 * a human cell. It classifies the beads into nine types (1-9)
 * depending on their degree of contact with the nucleus laminar
 */

#include <iostream>
#include <vector>
#include <string>
#include <cmath>

using std::cout;
using std::cin;
using std::endl;
using std::string;
using std::vector;
using std::cos;
using std::sin;

int main(int argc, char * argv[])
{
  const int numOfFibre {46};
  const long chromoLength [numOfFibre]
     {249250621, 243199373, 198022430, 191154276, 180915260,
	  171115067, 159138663, 146364022, 141213431, 135534747,
	  135006516, 133851895, 115169878, 107349540, 102531392,
	  90354753,  81195210,  78077248,  59128983,  63025520,
	  48129895,  51304566,  155270560,
	  249250621, 243199373, 198022430, 191154276, 180915260,
	  171115067, 159138663, 146364022, 141213431, 135534747,
	  135006516, 133851895, 115169878, 107349540, 102531392,
	  90354753,  81195210,  78077248,  59128983,  63025520,
	  48129895,  51304566,  155270560
	  };
  const int numOfBeads {1};

  const int totalBeads {numOfFibre * numOfBeads};

  // Box size of simulation
  double lx {120};
  double ly {120};
  double lz {120};
  
  // Initialise the laminar contact and colour of each bead
  vector<vector<double> > laminar 
	(numOfFibre, vector<double>(numOfBeads, 0.0));
  vector<vector<double> > colour 
	(numOfFibre, vector<double>(numOfBeads, 0.0));
  
  // Initialise positions of beads, bonds, and angles
  vector<vector<vector<double> > > position
	(numOfFibre, vector<vector<double> > 
	 (numOfBeads, vector<double> (3,0.0)));
  vector<vector<int> > bond (totalBeads, vector<int>(2,0));
  vector<vector<int> > angle (totalBeads, vector<int>(3,0));
  
  // Initialise the number of beads in each fibre
  vector<int> numOfBeadInChromo (numOfFibres, numOfBeads);

  
}


