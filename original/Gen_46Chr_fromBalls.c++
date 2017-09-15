//This program "COLORS" Human Chromosomes 19-22  in such a way that
//there are proms-LAM-PCG

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
int bond1=2;
int ang1=3;
int c=1;
int nb=1;
int na=1;

const int nfiber=46;
const long int Nbeads=50000;
const long int BasePairsChr=6072610000; //3079843747; //chr 1- ... -22 X times 2 // diploid
string Chromosome;
string CHR = "46chr";
const int cg=10000; //coarsegraining level 10kb

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

//BALLS POSITIONS
vector<double>CMBall1(3,0);
vector<vector<double> >CMBall(nfiber,CMBall1);
    
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

//MAIN BODY
int main(int argc, char* argv[]){
//
double Lx=200;
double Ly=200;
double Lz=200;
double LxL=1000;
double LyL=1000;
double LzL=1000;
long int nbead;
int npoly,ntype,type;
int nbeads,ntypepol;
//
ntypepol=10;
//
int coarsegrain;
coarsegrain=cg;
//
nbeads = ceil(BasePairsChr*1.0/cg);
cout << "Nbeads in total " << nbeads<<endl;
//
string Type="RW";
InitialConfg=1;//RW

////////////////////////////////////////////////////////////////////////////////////
InitialiseColors();
//////////////////////////////////////////////////////////////////////////////////
    
int nmaxx=nfiber*Nbeads;
vector<int>bond1(2,0);
vector<int>ang1(3,0);
vector<vector<int> >bond(nmaxx,bond1);
vector<vector<int> >angle(nmaxx,ang1);
    
//double position[nfiber][Nbeads][3];
vector<double>pos(3,0);
vector<vector<double> >position1(Nbeads,pos);
vector<vector<vector<double> > >position(nfiber,position1);
//
for(int r=0; r< 46; r++)for(int cc=0; cc< Nbeads; cc++)color[r][cc]=1;
nbead =0;
//
//LEngths
for(int p=0;p<46;p++)nbeadChr[p]=round(Length[p]/cg);
//This time the lengths are the real lengths after coarse graining
//
ofstream writeW;
stringstream writeFileW;
writeFileW << "LammpsInput.46Chr";
writeW.open(writeFileW.str().c_str());
///////////////////////////////////////////////////////////////////
int Rstart=0;
int chrnum=0;

// LOOP OVER CHROMOSOMES
for(reps=Rstart;reps<46;reps++){
sleep(1);
string ChrID;
stringstream convert;

chrnum = reps%23 +1;
convert << "chr"<<chrnum;
ChrID = convert.str();
cout << "ChrID " <<ChrID <<endl;
if(chrnum==23) ChrID="chrX";

//
srand(time(NULL));
cout << reps <<endl;
//
//READ FILES
ifstream readPos;
stringstream readPosFile;
readPosFile << "Eq.46Chr1Ball.LAMonly.1";
readPos.open(readPosFile.str().c_str());
//
ifstream readProms;
stringstream readPromsFile;
readPromsFile << "H3K4me3.Pk.genome.full.dat";
readProms.open(readPromsFile.str().c_str());
//
ifstream readLAM;
stringstream readLAMFile;
readLAMFile << "LAD.Pk.genome.full.dat";
readLAM.open(readLAMFile.str().c_str());
//
ifstream readPCG;
stringstream readPCGFile;
readPCGFile << "PCG.genome.full.dat";
readPCG.open(readPCGFile.str().c_str());
//
////////////////////////////////////////////////
////////	READ POSITIONS OF BALLS		/////////
////////////////////////////////////////////////
string dummm; double xball,yball,zball, idball, molball,typeball;int ixx,iyy,izz;
if(!readPos){cout << "where's the Eq.46Chr... file?"<<endl;cin.get();}
for(int n=0;n<36;n++)getline(readPos,dummm);
for(int i=0;i<5046;i++){
	readPos >> idball>> molball>> typeball>>xball>>yball>>zball>>ixx>>iyy>>izz;
	if(molball>0){
	CMBall[idball-1][0]=xball;
	CMBall[idball-1][1]=yball;
	CMBall[idball-1][2]=zball;
	}
//if(i==5045){cout << idball << " " << xball<<endl;cin.get();}
}
readPos.close();

///////////////////////////////////////////
////////	SET	CHRMM STATE         //////////
///////////////////////////////////////////
string dummy;

//getline(readState,dummy);
//cout << dummy <<endl;
long int start,end;
double Awild;
int nb=0;
double sscore;
string ch;
string dumm;
for(int i=0;i<1;i++)getline(readProms,ch);
cout << "READING " <<readPromsFile.str().c_str() <<endl; //cin.get();
if(!readProms){cout << "where's the promoters file?"<<endl;cin.get();}
int countbeads=0;
int bin;
double bo;
while(!readProms.eof()){
	readProms >>bin >> ch>> start >> end >> dumm >>sscore>>dumm>>bo >> bo >>bo;
	//cout << ch << " " <<  start << " " << end << " " << sscore << endl;
		int nstart = floor(start*1.0/coarsegrain);
		int nend = ceil(end*1.0/coarsegrain);
		for(int nnb=nstart; nnb<nend; nnb++){
			if(ch==ChrID && sscore>0){
			// HI STICKY - PROM
			PROM[reps][nnb]=1;
	      		if(nnb==nstart)AbundancePROM[reps][nnb]+=min(1-(start*1.0/cg-nstart),(end-start)*1.0/cg);
	       		if(nnb==nend-1 && nnb!=nstart)AbundancePROM[reps][nnb]+=(end*1.0/cg-(nend-1));
			if(nnb!=nstart && nnb!=nend-1) AbundancePROM[reps][nnb]=1;

          //     cout << "PROM " << reps+1 << " " << nnb << " " << AbundancePROM[reps][nnb]<<endl; cin.get();
	//cout << 1-(start*1.0/cg-nstart) << " " << (end-start)*1.0/cg << " " << (end*1.0/cg-(nend-1)) << " " << AbundancePROM[reps][nnb]<<endl; cin.get();
			}
        }
    countbeads=nend;
}
readProms.close();
    
//LAMINB1
cout << "READING " <<readLAMFile.str().c_str() <<endl; //cin.get();
if(!readLAM){cout << "where's the LAD file?"<<endl;cin.get();}
//cout << dummy;
string aux;
double score=0;  
double th=1;
int dir=1;
double ladstart=0;double ladend=0;
double dum;
string cdum;double startold=0;
for(int i=0;i<1;i++){getline(readLAM,cdum);}//cout<<cdum<<endl;}
while(!readLAM.eof()){
	//cout <<  Chromosome << start << end << aux << score << dir << dum << dum;
	readLAM >>bin >> ch>>ladstart >> ladend;
    //cout <<  bin << " " << ch << " " << ladstart << " "<<ladend<<endl; cin.get();
    int nstart = floor(ladstart*1.0/coarsegrain);
    int nend = ceil(ladend*1.0/coarsegrain);
    if(1==1 && ch==ChrID){
		for(int nnb=nstart; nnb<nend; nnb++){
        LAM[reps][nnb]=1;
	if(nnb==nstart)AbundanceLAM[reps][nnb]+=min(1-(ladstart*1.0/cg-nstart),(ladend-ladstart)*1.0/cg);
	if(nnb==nend-1 && nnb!=nstart)AbundanceLAM[reps][nnb]+=(ladend*1.0/cg-(nend-1));
	if(nnb!=nstart && nnb!=nend-1) AbundanceLAM[reps][nnb]=1;
	//cout << "LAM " << reps << " " << nnb << " " << endl;
	//cin.get();
        }
	}
    countbeads=max(countbeads,nend);
}
readLAM.close();

    
//PCG
//cout << dummy;
cout << "READING " <<readPCGFile.str().c_str() <<endl; //cin.get();
if(!readPCG){cout << "where's the PCG file?"<<endl;cin.get();}
string pcgaux;
double pcgscore=0;
double pcgth=40;
for(int i=0;i<1;i++){getline(readPCG,pcgaux);}//cout<<pcgaux<<endl;}
while(!readPCG.eof()){
    readPCG >>ch>> start >> end >> pcgscore;
    //cout << ch << " " <<  start << " " << end <<endl;
    int nstart = floor(start*1.0/coarsegrain);
    int nend = ceil(end*1.0/coarsegrain);
    for(int nnb=nstart; nnb<nend; nnb++){
        if(ch==ChrID && pcgscore>15){
        // PCG
        PCG[reps][nnb]=1;
 	if(nnb==nstart)AbundancePCG[reps][nnb]+=min(1-(start*1.0/cg-nstart),(end-start)*1.0/cg);
	if(nnb==nend-1 && nnb!=nstart)AbundancePCG[reps][nnb]+=(end*1.0/cg-(nend-1));
	if(nnb!=nstart && nnb!=nend-1) AbundancePCG[reps][nnb]=1;
          // cout << "PCG " << reps << " " << nnb << " " << AbundancePCG[reps][nnb]<< " " << AbundancePROM[reps][nnb]<< " " << AbundanceLAM[reps][nnb] << endl; //cin.get();
        }
    }
    countbeads=max(countbeads,nend);
}
readPCG.close();
//nbeadChr[reps]=countbeads;
cout << "chromosome " << reps+1 << " with length [Mbp] " << nbeadChr[reps]*1.0/1000000 <<endl;//cin.get();
//19706 -- 20989 -- 16040 -- 17075
//////////////////////////////////////////////////////////////////////
////////						//////////////
//////////////////////////////////////////////////////////////////////
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
// USE BALLS CM + RANDOM TILT
///////////////////////////////////////////////////////////////////////////
if(InitialConfg==1){
    double shift[3]={CMBall[reps][0],CMBall[reps][1],CMBall[reps][2]};
    //
    double R=10.;
    int Nturns=250; //this is because there are ~200 beads per turn
    double p=1./500.; //in units of sigma -- jump 1sigma per turn
    int s=0;
    int m=0;
    double Tilt;
int checknear=1;
while(checknear==1){
cout << "doing chromosome " <<  reps+1 <<endl;
int whichtilt=int(rand())%2;
    for(int m=0;m<nbeadChr[reps];m++){
        if(m==0){
        Tilt=(((double)(rand())/((double)(RAND_MAX)))-0.5)*M_PI;
        //Tilt2=(((double)(rand())/((double)(RAND_MAX)))-0.5)*1.0*M_PI;
        }
        double noise=((double)(rand())/(double)(RAND_MAX)-0.5)*0.25;
        theta[1]=m*1.0/Nturns*M_PI;
        double Rv=6+4*cos(10*theta[1]);
        position[reps][m][0]=+Rv*cos(theta[1])+noise;
        position[reps][m][1]=+Rv*sin(theta[1])+noise;
        position[reps][m][2]=+m*p-(nbeadChr[reps]*p)/2.;
        // MAKE ROTATION
	double oldpos[3]={position[reps][m][0],position[reps][m][1],position[reps][m][2]};
       if(whichtilt==0){
       position[reps][m][0]=oldpos[0]*cos(Tilt)-oldpos[2]*sin(Tilt)+shift[0];
       position[reps][m][1]=oldpos[1]+shift[1];
       position[reps][m][2]=oldpos[2]*cos(Tilt)+oldpos[0]*sin(Tilt)+shift[2];
       }
	else{
       position[reps][m][0]=oldpos[0]*cos(Tilt)-oldpos[1]*sin(Tilt)+shift[0];
       position[reps][m][1]=oldpos[1]*cos(Tilt)+oldpos[0]*sin(Tilt)+shift[1];
       position[reps][m][2]=oldpos[2]+shift[2];
	}
        //if(m==nbeadChr[reps]-1){cout << c << " " << position[reps][m][0]<< " " << position[reps][m][1] << " " << position[reps][m][2]<<endl;cin.get();}
       // if(m==0){cout << c << " " << position[reps][m][0]<< " " << position[reps][m][1] << " " << position[reps][m][2] << endl;cin.get();}
        //writeW << c << " " << reps+1 << " " << color[reps][m] << " " << position[reps][m][0]-shift[0]<<" " << position[reps][m][1]-shift[1] << " " << position[reps][m][2]-shift[2] << " " << 0 << " " << 0 << " " << 0 << endl;
checknear=0;
if(m%1000==0)for(int nn=Rstart;nn<reps;nn++)for(int mm=0;mm<nbeadChr[nn];mm++)if(pow(position[nn][mm][0]-position[reps][m][0],2.0)+pow(position[nn][mm][1]-position[reps][m][1],2.0)+pow(position[nn][mm][2]-position[reps][m][2],2.0)<3){cout << "AIA OVERLAP! " ; checknear=1; break;}
if(checknear==1){cout << reps << endl; break;}//cin.get(); break;}
        c++;
        }
    }
}
        cout << "DONE Chromosome "<<reps+1 << endl;
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
    cout << "#beads in Chromosome " << nn+1<<" = " << nbeadChr[nn] << " @ " << cg <<"bp resolution" <<endl;
    for(int m=0;m<nbeadChr[nn]; m++){
        writeW << cc << " " << nn+1 << " " << color[nn][m] << " " << position[nn][m][0]<<" " << position[nn][m][1] << " " << position[nn][m][2] << " " << 0 << " " << 0 << " " << 0 << endl;
        //cout << "cc " << cc << " " << nn+1 << " "  << color[nn][m]<<endl;
        //cin.get();
        
        cc++;
    }
}
///////////////////////////////////////////////////////////
//FINISHED POLYMER
///////////////////////////////////////////////////////////

writeW << endl;
writeW << endl;
writeW << "\n Velocities \n" <<endl;
for(int j=0;j<totb; j++) writeW<<j+1 << " "<<0 << " "<<0<<" "<<0 <<endl;

writeW << endl;
writeW << "\n Bonds \n"<<endl; 
for(int i=0;i<nbond;i++) {
writeW << i+1 <<" "<< 1 <<" "<< bond[i][0]+1<<" "<<bond[i][1]+1 << endl;
}

writeW << "\n Angles \n"<<endl; 
for(int i=0;i<nangle;i++) writeW << i+1 <<" "<< 1 <<" "<< angle[i][0]+1<<" "<<angle[i][1]+1 <<" "<< angle[i][2]+1<< endl;

writeW.close();
stringstream finalcomm; 

return 0;
}

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
	for(int n=0;n<Nbeads;n++){
    //SOLVE CONFLICT
    //cout << reps << " " << n << " " << AbundancePCG[reps][n]<< " " << AbundancePROM[reps][n]<< " " << AbundanceLAM[reps][n] << endl; cin.get();

    if(PROM[reps][n]==1 && (LAM[reps][n]==1 or PCG[reps][n]==1)){
        if(AbundancePROM[reps][n]>=3*AbundanceLAM[reps][n])LAM[reps][n]=0;else PROM[reps][n]=0;
        if(AbundancePROM[reps][n]>=3*AbundancePCG[reps][n])PCG[reps][n]=0;else PROM[reps][n]=0;
//        if(AbundancePROM[reps][n]<AbundancePCG[reps][n])PROM[reps][n]=0;
//        if(AbundancePROM[reps][n]<AbundanceLAM[reps][n])PROM[reps][n]=0;
    }
	//SOLVE CONFLICT
    if(PCG[reps][n]==1 && LAM[reps][n]==1){
 	 if(AbundancePCG[reps][n]>=2*AbundanceLAM[reps][n])LAM[reps][n]=0;else PCG[reps][n]=0;
    }

    //IF PROM LEVEL NOT IHGH ENOUGH PUT NEUTR
    if(AbundancePROM[reps][n]<=0.5)PROM[reps][n]=0;
    if(AbundanceLAM[reps][n]>0.8){PROM[reps][n]=0;LAM[reps][n]=1;PCG[reps][n]=0;}

    color[reps][n]=1;
        if(PROM[reps][n]==1){color[reps][n]=2; cout << n << " is promoter "<<endl;} //prom
        if(LAM[reps][n]==1 && PCG[reps][n]==0){color[reps][n]=3;  cout << n << " is lam" <<endl;} //LAM
        if(PCG[reps][n]==1 && LAM[reps][n]==0){color[reps][n]=4; cout << n << " is pcg" <<endl;} //PCG
        if(LAM[reps][n]==1 && PCG[reps][n]==1){ color[reps][n]=5;  cout << "aiaiai" <<endl; cin.get();}//mixed
	//if(color[reps][n]>1)cin.get();
        
	}
}

