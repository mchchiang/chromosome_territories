/* Gen_46Chr_1BallEach.cpp
 * This code generates a bead to model each chromosome of 
 * a human cell. It classifies the beads into nine types (1-9)
 * depending on their degree of contact with the nucleus laminar
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>

using namespace std;

int main(int argc, char * argv[])
{
  const int haploidNum {23};
  const int numOfFibres {haploidNum*2};
  const long chromoLength [numOfFibres]
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
  const int totalBeads {numOfFibres * numOfBeads};

  // Box size of simulation
  double lx {120};
  double ly {120};
  double lz {120};
  
  // Initialise the laminar contact and colour of each bead
  vector<vector<double> > laminar 
	(numOfFibres, vector<double>(numOfBeads, 0.0));
  vector<vector<double> > colour 
	(numOfFibres, vector<double>(numOfBeads, 0.0));
  
  // Initialise positions of beads, bonds, and angles
  vector<vector<vector<double> > > position
	(numOfFibres, vector<vector<double> > 
	 (numOfBeads, vector<double> (3,0.0)));
  vector<vector<int> > bond (totalBeads-1, vector<int>(2,0));
  vector<vector<int> > angle (totalBeads-2, vector<int>(3,0));
  
  // Initialise the number of beads in each fibre
  vector<int> numOfBeadsInChromo (numOfFibres, numOfBeads);
  
  // Reading laminar file
  cout << "Reading laminar file " << endl;

  ifstream readLAM;
  string readLAMFile {"LAD.Pk.genome.full.dat"};
  readLAM.open(readLAMFile);

  if (!readLAM)
    {
      cout << "Unable to read the LAM file \"" 
	   << readLAMFile << "\"\n"
	   << "Aborting reading process" << endl;
      return 1;
    }

  // Remove header line
  string header;
  getline(readLAM, header); 

  const string chrPrefix {"chr"};
  double lamStart {};
  double lamEnd {};
  double lamContent {};
  string ch {};
  int bin {};
  int chromo {};

  while (!readLAM.eof())
    {
      readLAM >> bin >> ch >> lamStart >> lamEnd;
      ch.erase(0, chrPrefix.length());

      if (ch != "" && ch != "Y") // Ignore chrY data
	{
	  if (ch == "X")
	    chromo = haploidNum;
	  else 
	    chromo = stoi(ch, nullptr, 10);

	  lamContent = (lamEnd - lamStart) / chromoLength[chromo-1];
	  laminar[chromo-1][0] += lamContent;
	  laminar[chromo+haploidNum-1][0] += lamContent;
	}
    }

  readLAM.close();

  cout << "Finish reading laminar file" << endl;

  // Assigning colour depending on laminar content
  // Binning contact strength into 10 categories (1-10)
  cout << "Assigning contact strength type to each chromosome" << endl;
  int numOfBins {10};
  double min {0.1};
  double max {0.6};
  double binWidth = (max - min) / numOfBins;
  int binNum {};

  for (int i {}; i < numOfFibres; i++)
    for (int j {}; j < numOfBeadsInChromo[i]; j++)
      {
	binNum = static_cast<int>((laminar[i][j] - min) / binWidth + 1);
	if (binNum < 0) binNum = 0;
	colour[i][0] = binNum;
      }

  // Output total laminar content for each chromosome
  cout << "The laminar content for each chromosome is as follows:" << endl;
  cout << setw(15) << "Chromosome #"
       << setw(15) << "Length [Mbp]"
       << setw(15) << "LAM content"
       << setw(15) << "Bin/Type"
       << endl;

  for (int i {}; i < numOfFibres; i++)
    {
      cout << setw(15) << (i+1)
	   << setw(15) << (chromoLength[i] / 1000000.0)
	   << setw(15) << laminar[i][0]
	   << setw(15) << colour[i][0]
	   << endl;
    }
  
  // Creating the polymer
  // Create all the bonds and angles
  int startBead {};
  int localNumOfBeads {};
  int angleIndex {};
  int bondIndex {};

  for (int i {}; i < numOfFibres; i++)
    {
      localNumOfBeads = numOfBeadsInChromo[i];
      // Store bonds
      for (int j {}; j < localNumOfBeads-1; j++)
	{
	  bond[bondIndex][0] = j+startBead;
	  bond[bondIndex][1] = j+1+startBead;
	  bondIndex++;
	}

      // Store angles
      for (int j {}; j < localNumOfBeads-2; j++)
	{
	  angle[angleIndex][0] = j+startBead;
	  angle[angleIndex][1] = j+1+startBead;
	  angle[angleIndex][2] = j+2+startBead;
	  angleIndex++;
	}

      startBead += numOfBeadsInChromo[i];
    }

  // Generate position of each bead
  srand(time(NULL)); // Initialise the random generator
  
  double theta[2];
  double non {static_cast<double>(rand())/static_cast<double>(RAND_MAX)};
  double tilt {0};
  double shift[3] {0.0, 0.0, 0.0};
  double sh[2] {cos(non*92*M_PI/46.0), sin(non*92*M_PI/46.0)};
  double R {5.0};
  
  
}


