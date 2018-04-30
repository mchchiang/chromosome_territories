/* ContactPlot.cpp
 *
 * A program that creates the gnuplot script for making 
 * an arc diagram of the contacts for a specific genome region
 */

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include "ContactMapLib.hpp"

using std::cout;
using std::endl;
using std::string;
using std::ofstream;

class RGB {
public:
  unsigned char R;
  unsigned char G;
  unsigned char B;

  RGB(unsigned char r, unsigned char g, unsigned char b){
    R = r;
    G = g;
    B = b;
  }

  bool Equals(RGB rgb){
    return (R == rgb.R) && (G == rgb.G) && (B == rgb.B);
  }
};

class HSV{
public:
  double H;
  double S;
  double V;

  HSV(double h, double s, double v){
    H = h;
    S = s;
    V = v;
  }

  bool Equals(HSV hsv){
    return (H == hsv.H) && (S == hsv.S) && (V == hsv.V);
  }
};

static RGB HSVToRGB(HSV hsv){
  double r = 0, g = 0, b = 0;

  if (hsv.S == 0){
    r = hsv.V;
    g = hsv.V;
    b = hsv.V;
  } else {
    int i;
    double f, p, q, t;
    
    if (hsv.H == 360)
      hsv.H = 0;
    else
      hsv.H = hsv.H / 60;
    
    i = (int)trunc(hsv.H);
    f = hsv.H - i;
    
    p = hsv.V * (1.0 - hsv.S);
    q = hsv.V * (1.0 - (hsv.S * f));
    t = hsv.V * (1.0 - (hsv.S * (1.0 - f)));
    
    switch (i){
    case 0:
      r = hsv.V;
      g = t;
      b = p;
      break;
      
    case 1:
      r = q;
      g = hsv.V;
      b = p;
      break;
      
    case 2:
      r = p;
      g = hsv.V;
      b = t;
      break;
      
    case 3:
      r = p;
      g = q;
      b = hsv.V;
      break;
      
    case 4:
      r = t;
      g = p;
      b = hsv.V;
      break;
      
    default:
      r = hsv.V;
      g = p;
      b = q;
      break;
    }
    
  }
  
  return RGB((unsigned char)(r * 255), (unsigned char)(g * 255), (unsigned char)(b * 255));
}

struct contact {
  int x, y;
  double score;
  contact(int i, int j, double s) : x {i}, y {j}, score {s} {}
};

double getParaConst(double sep, double height, int numOfBeads);
double getColour(double score, double max);

int main(int argc, char* argv[]){

  if (argc < 9){
    cout << "Usage: [numOfBeads] [startBead] [endBead] [localDist] [height]"
	 << "[mode=full/upper] [matrixFile] [outFile]" << endl;
    return 1;
  }

  int argi {};
  int numOfBeads {stoi(string(argv[++argi]), nullptr, 10)};
  int startBead {stoi(string(argv[++argi]), nullptr, 10)};
  int endBead {stoi(string(argv[++argi]), nullptr, 10)};
  int localDist {stoi(string(argv[++argi]), nullptr, 10)};
  double height {stod(string(argv[++argi]), nullptr)};
  string mode (argv[++argi]);
  string matrixFile (argv[++argi]);
  string outFile (argv[++argi]);

  bool full {true};
  if (mode != "full") full = false;

  CMap map = ContactMap::createFromMatrixFile(numOfBeads, full, matrixFile);

  ofstream writer;
  writer.open(outFile);
  if (!writer){
    cout << "Problem with opening the plot file!" << endl;
  }

  // Find the local and distal contacts
  vector<contact> local {};
  vector<contact> distal {};
  double score;
  double maxScore {};
  for (int i {startBead}; i <= endBead; i++){
    for (int j {}; j < i; j++){
      score = map->get(i,j);
      if (score > maxScore) maxScore = score;
      if (score < 1e-10) continue;
      if (abs(i-j) <= localDist){ // Local contact
	local.push_back(contact(i,j,score));
      } else { // Distal contact
	distal.push_back(contact(i,j,score));
      }
    }
    for (int j {endBead+1}; j < numOfBeads; j++){
      score = map->get(i,j);
      if (score > maxScore) maxScore = score;
      if (score > 1e-10) continue;
      if (abs(i-j) <= localDist){ // Local contact
	local.push_back(contact(i,j,score));
      } else { // Distal contact
	distal.push_back(contact(i,j,score));
      }
    }
  }

  // Create the gnuplot script
  writer << "# gnuplot plotting file\n" << endl;
  writer << "f(x,a,b,c) = c*(x-a)*(x-b)\n" << endl;

  writer << "unset key" << endl;
  writer << "set xrange [0:" << numOfBeads-1 << "]" << endl;
  writer << "set yrange [0:" << height*1.1 << "]\n" << endl;

  writer << "p ";
  // Write distal contacts first
  for (size_t i {}; i < distal.size(); i++){
    writer << "f(x," << distal[i].x << "," << distal[i].y << ","
	   << getParaConst(abs(distal[i].x-distal[i].y), height, numOfBeads) 
	   << ") w l lw 1 lc rgb '#" 
	   << std::hex 
	   << getColour(distal[i].score, maxScore) 
	   << std::dec << "', \\" << endl; 
  }

  // Write local contacts
  for (size_t i {}; i < local.size(); i++){
    writer << "f(x," << local[i].x << "," << local[i].y << ","
	   << getParaConst(abs(local[i].x-local[i].y), height, numOfBeads) 
	   << ") w l lw 1 lc rgb '#"
	   << std::hex 
	   << getColour(local[i].score, maxScore) 
	   << std::dec << "', \\" << endl; 
  }
  writer.close();
}

double getParaConst(double sep, double max, int numOfBeads){
  double c {max/numOfBeads*sep};
  return -4.0*c/(sep*sep);
}

double getColour(double score, double max){
  // Use 10 bins
  const int numOfBins {10};
  double inc {max/numOfBins};
  int binIndex {static_cast<int>(score/inc)+1};
  if (binIndex > numOfBins) binIndex = numOfBins;
  //  HSV hsv (180, static_cast<double>(binIndex)/numOfBins, 1.0);
  //RGB rgb = HSVToRGB(hsv);
  //  return rgb.R+256*256+rgb.G+256+rgb.B;
  return static_cast<int>(static_cast<double>(binIndex)/numOfBins*256);
}
