//Here I generate one ball per chromosome. The quantity of LAMIN in each chromosome
//is the indicator of its affinity to the nuclear envelope (different types depending on that). 
//For the final position of the ball I create the chromosomes.

using namespace std;
#include<iostream>
#include<math.h>
#include<stdlib.h>
#include<fstream>
#include<sstream>
#include<vector> 
#include <unistd.h>
#include<ctime>

int InitialConfg=0;
//0 == RANDOM WALK; 1 == Asymmetric RW (stretched); 2 == CRUMPLED inside SPHERE; 3==Fractal globule
// 5 == Straight line
int bond1=2;
int ang1=3;
int c=1;
int nb=1;
int na=1;

const int nfiber=46;
const long int Nbeads=1;
const long int BasePairsChr=6072610000; //chr 1- ... -22 X times 2 // diploid
string Chromosome;
string CHR = "46chr";
const int cg=1000000;

//LENGTHs
long int Length[46]={249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560,249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560};
int nbeadChr[46];

//STATES
//double PCG[nfiber][Nbeads];
//double LAM[nfiber][Nbeads];
//double PROM[nfiber][Nbeads];
vector<double>beadsn(Nbeads,0);
vector<vector<double> >PCG(nfiber,beadsn);
vector<vector<double> >LAM(nfiber,beadsn);
vector<vector<double> >PROM(nfiber,beadsn);

    
//double AbundancePROM[nfiber][Nbeads];
//double AbundancePCG[nfiber][Nbeads];
//double AbundanceLAM[nfiber][Nbeads];
vector<vector<double> >AbundancePCG(nfiber,beadsn);
vector<vector<double> >AbundanceLAM(nfiber,beadsn);
vector<vector<double> >AbundancePROM(nfiber,beadsn);

//int color[nfiber][Nbeads];
vector<vector<double> >color(nfiber,beadsn);

int reps;
int nbond=0;
int nangle=0;

double center[3]={0.0,0.0,0.0};
double distance(double v1[3], double v2[3]);
void InitialiseColors();
void Paint();

//MAIN BODY OF FILE
int main(int argc, char* argv[]){

double Lx=120;
double Ly=120;
double Lz=120;
double LxL=1000;
double LyL=1000;
double LzL=1000;
long int nbead;
int npoly,ntype,type;
int nbeads,ntypepol;

ntypepol=10;

int coarsegrain;
coarsegrain=cg;

//beads per polymer -- just spheres!
nbeads = 1;
cout << "Nbeads " << nbeads<<endl;


string Type="RW";
InitialConfg=1;//RW
///////////////////////////////////////////////////////////////////////////////////////
InitialiseColors();
///////////////////////////////////////////////////////////////////////////////////////
    
///////////////////////////////////////////////////////////////////////////////////////
int nmaxx=nfiber*Nbeads;
vector<int>bond1(2,0);
vector<int>ang1(3,0);
vector<vector<int> >bond(nmaxx,bond1);
vector<vector<int> >angle(nmaxx,ang1);
    
vector<double>pos(3,0);
vector<vector<double> >position1(Nbeads,pos);
vector<vector<vector<double> > >position(nfiber,position1);

for(int r=0; r< 46; r++)for(int cc=0; cc< Nbeads; cc++)color[r][cc]=1;
nbead =0;
   
//LENGHTS
for(int p=0;p<46;p++)nbeadChr[p]=1;
///////////////////////////////////////////////////////////////////////////////////////

    
///////////////////////////////////////////////////////////////////////////////////////
ofstream writeW;
stringstream writeFileW;
writeFileW << "LammpsInput.46Chr.1BallEach";
writeW.open(writeFileW.str().c_str());
///////////////////////////////////////////////////////////////////////////////////////

    
///////////////////////////////////////////////////////////////////////////////////////
//READ EXTERNAL FILES
///////////////////////////////////////////////////////////////////////////////////////
int Rstart=0;
int chrnum=0;
int count=0;

for(reps=Rstart;reps<46;reps++){
sleep(1);
string ChrID;
stringstream convert;

chrnum = reps%23 +1;
convert << "chr"<<chrnum;
ChrID = convert.str();
cout << "ChrID " <<ChrID <<endl;
 if(chrnum==23) ChrID="chrX";
   
srand(time(NULL));
cout << reps <<endl;

ifstream readProms;
stringstream readPromsFile;
readPromsFile << "H3K4me3.Pk.genome.full.dat";
readProms.open(readPromsFile.str().c_str());

ifstream readLAM;
stringstream readLAMFile;
readLAMFile << "LAD.Pk.genome.full.dat";
readLAM.open(readLAMFile.str().c_str());

ifstream readPCG;
stringstream readPCGFile;
readPCGFile << "PCG.genome.full.dat";
readPCG.open(readPCGFile.str().c_str());

//////////////////////////////////////////////////////
string dummy;
//////////////////////////////////////////////////////
//LAMINB1
cout << "READING " <<readLAMFile.str().c_str() <<endl; //cin.get();
//cout << dummy;
string aux;
double score=0;  
double th=1;
int dir=1;
double ladstart=0;double ladend=0;
double dum;
string ch;
string cdum;double startold=0;
int bin;
for(int i=0;i<1;i++){getline(readLAM,cdum);cout<<cdum<<endl;}
while(!readLAM.eof()){
	readLAM >>bin >> ch>>ladstart >> ladend;
	int nnb=0;
	 if(1==1 && ch==ChrID){
		
        	LAM[reps][nnb]+=(ladend-ladstart)*1.0/Length[reps];

	//cout << "LAM " << reps << " " << nnb << " " << endl;
	//cin.get();
        }
}
readLAM.close();
//////////////////////////////////////////////////////
cout <<"Chromosome " <<  reps+1 << " with length [Mbp] " << Length[reps]*1.0/1000000 << " & LAM content " << LAM[reps][0] << endl;//cin.get();
cout << "press a key " << endl; cin.get();
//////////////////////////////////////////////////////////////
////////		PICK COLOR	//////////////////////
//////////////////////////////////////////////////////////////
cout << "I'M PAINTING..." <<endl;
Paint();
cout << "I'VE FINISHED PAINTING..." <<endl;
//////////////////////////////////////////////////////////////
////////		MAKE POLYMER     /////////////////////
//////////////////////////////////////////////////////////////
cout << "BONDS" <<endl;
int cumulative=0;
for(int rr=Rstart;rr<reps;rr++)cumulative+=nbeadChr[rr];
int nbeadT=0;
for(int i=0;i<nbeadChr[reps]-1;i++){
nbeadT++;
    //cout << " if " <<endl;
    bond[nbond][0] = i+cumulative;
    bond[nbond][1] = i+1+cumulative;
    nbond++;
    //cout << nbond << " " << i << endl;
}
    
for(int i=0;i<nbeadChr[reps]-2;i++){
    angle[nangle][0] = i+cumulative;
    angle[nangle][1] = i+1+cumulative;
    angle[nangle][2] = i+2+cumulative;
    nangle++;
}

cout << "MAKING THE POLYMER ..." <<endl;
///////////////////////////////////////////////////////////
//make polymer
double theta[2];
double phi[4];
double r = 1.1;
double non1=((double)(rand())/((double)(RAND_MAX)));
double non=((double)(rand())/((double)(RAND_MAX)));
    
///////////////////////////////////////////////////////////////////////////
// GENERAL WAY OF MAKING CHROMOSOMES
// HERE only 1 bead per poly -> Spheres
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// TILTED AND NOISY HELIX
///////////////////////////////////////////////////////////////////////////
if(InitialConfg==1){
    double Tilt=0;
    double shift[3]={0.0,0.0,0.0};
    //cout << " sh " << non << endl; cin.get();
    //here i'm doing a shift in the quadrants of a plane
    double sh[2]={cos(non*92*M_PI/46.),sin(non*92*M_PI/46.)};
    //
    double noise1=0.7;
    double R=5.;
    int Nturns=20; //this is because there are ~200 beads per turn
    double p=4./20.; //in units of sigma -- jump 1sigma per turn
    int s=0;
    int m=0;
int checknear=1;
while(checknear==1){
cout << "doing " <<  reps <<endl;
    for(int m=0;m<nbeadChr[reps];m++){
        if(m==0){
        shift[0]=((double)(rand())/((double)(RAND_MAX)))*35+sh[0]*25;
        shift[1]=((double)(rand())/((double)(RAND_MAX)))*35+sh[1]*25;
        shift[2]=((double)(rand())/((double)(RAND_MAX))-0.5)*20+((reps%4)-2)*10;//*20+(reps%4)*5
        if(reps%2==0)Tilt=(((double)(rand())/((double)(RAND_MAX)))-0.5)*M_PI/4.;
        if(reps%2==1)Tilt=-(((double)(rand())/((double)(RAND_MAX)))-0.5)*M_PI/4.;
            //cout << Tilt << " T " << endl; cin.get();
        }
        double noise=((double)(rand())/(double)(RAND_MAX)-0.5)*0.25;
        theta[1]=m*1.0/Nturns*M_PI;
        double Rv=5+5*cos(8*theta[1]);
        position[reps][m][0]=+R*cos(theta[1])+noise;
        position[reps][m][1]=+R*sin(theta[1])+noise;
        position[reps][m][2]=+m*p;
        //
        position[reps][m][0]=position[reps][m][0]*cos(Tilt)+position[reps][m][2]*sin(Tilt)-shift[0]+20;
        position[reps][m][1]=position[reps][m][1]-shift[1]+20;
        position[reps][m][2]=-position[reps][m][0]*sin(Tilt)+position[reps][m][2]*cos(Tilt)-shift[2]-10;
       //if(m>0)cout << "distance " << m <<" " << sqrt((position[reps][m][0]-position[reps][m-1][0])*(position[reps][m][0]-position[reps][m-1][0])+(position[reps][m][1]-position[reps][m-1][1])*(position[reps][m][1]-position[reps][m-1][1]+(position[reps][m][2]-position[reps][m-1][2])*(position[reps][m][2]-position[reps][m-1][2])))<<endl;
        //cin.get();
        //if(m==nbeadChr[reps]-1){cout << c << " " << position[reps][m][0]<< " " << position[reps][m][1] << " " << position[reps][m][2]<<endl;cin.get();}
       // if(m==0){cout << c << " " << position[reps][m][0]<< " " << position[reps][m][1] << " " << position[reps][m][2] << endl;cin.get();}
        //writeW << c << " " << reps+1 << " " << color[reps][m] << " " << position[reps][m][0]-shift[0]<<" " << position[reps][m][1]-shift[1] << " " << position[reps][m][2]-shift[2] << " " << 0 << " " << 0 << " " << 0 << endl;
checknear=0;
if(m%500==0)for(int nn=Rstart;nn<reps;nn++)for(int mm=0;mm<nbeadChr[nn];mm++)if(pow(position[nn][mm][0]-position[reps][m][0],2.0)+pow(position[nn][mm][1]-position[reps][m][1],2.0)+pow(position[nn][mm][2]-position[reps][m][2],2.0)<2.0){cout << "AIA " ; checknear=1; break;}
if(checknear==1){cout << reps << endl; break;}//cin.get(); break;}
        c++;
        }
    }
}
        cout << "DONE RW "<<reps << endl;
} //close loop over reps
////////////////////////////////////////////////////
////////		WRITE FILE     /////////////////////
////////////////////////////////////////////////////
int totb=0;
for(int rr=Rstart;rr<nfiber;rr++)totb+=nbeadChr[rr];
writeW<< "LAMMPS data file from restart file: timestep = 0,\t procs = 1"<<endl;
writeW << totb << " atoms "<<endl;
writeW << nbond << " bonds "<<endl;
writeW << nangle << " angles "<<endl;
writeW << "\n";
writeW << ntypepol << " atom types "<<endl;
writeW << 1 << " bond types "<<endl;
writeW << 1 << " angle types "<<endl;
writeW << "\n";
writeW << -Lx/2.0 << " " << (Lx-Lx/2.0) << " xlo xhi"<<endl;
writeW << -Ly/2.0 << " " << (Ly-Ly/2.0) << " ylo yhi"<<endl;
writeW << -Lz/2.0 << " " << (Lz-Lz/2.0) << " zlo zhi"<<endl; ///TO BE CHANGED!!!
//
writeW << "\nMasses \n"<<endl;
for(int j=1; j<=ntypepol;j++) writeW << j << " " << 1 << endl; 
 //
int cc=1;
writeW << "\nAtoms \n"<<endl;
for(int nn=Rstart;nn<46; nn++){
    cout << "nl " << nbeadChr[nn] <<endl;
    for(int m=0;m<nbeadChr[nn]; m++){
        writeW << cc << " " << nn+1 << " " << color[nn][m] << " " << position[nn][m][0]<<" " << position[nn][m][1] << " " << position[nn][m][2] << " " << 0 << " " << 0 << " " << 0 << endl;
        //cout << "cc " << cc << " " << nn+1 << " "  << color[nn][m]<<endl;
        //cin.get();
        cc++;
    }
}
///////////////////////////////////////////////////////////
//FINISHED WRITING ATOMS
///////////////////////////////////////////////////////////
writeW << endl;
writeW << endl;
writeW << "\n Velocities \n" <<endl;
for(int j=0;j<totb; j++) writeW<<j+1 << " "<<0 << " "<<0<<" "<<0 <<endl;

writeW << endl;

/*
writeW << "\n Bonds \n"<<endl;
for(int i=0;i<nbond;i++) {
writeW << i+1 <<" "<< 1 <<" "<< bond[i][0]+1<<" "<<bond[i][1]+1 << endl;
}

writeW << "\n Angles \n"<<endl; 
for(int i=0;i<nangle;i++) writeW << i+1 <<" "<< 1 <<" "<< angle[i][0]+1<<" "<<angle[i][1]+1 <<" "<< angle[i][2]+1<< endl;
*/

writeW.close();
stringstream finalcomm; 

return 0;
}
//END BODY OF FILE

//FUNCTIONS
double distance(double v1[3], double v2[3]){
	double d=0;
	for(int n=0;n<3;n++) d+=(v1[n]-v2[n])*(v1[n]-v2[n]);
	return sqrt(d);
}

void InitialiseColors(){
    for(int n=0;n<nfiber;n++){
	for(int m=0;m<Nbeads;m++){
		color[n][m]=0;
	}
}
}

void Paint(){
int bins=10;
double db=0.05;
	for(int n=0;n<Nbeads;n++){
		for(int i=0;i<bins;i++)	if(LAM[reps][n]>i*db+0.1 && LAM[reps][n]<=(i+1)*db+0.1){color[reps][n]=i+1;cout << color[reps][n]<<endl;}
	}
}

