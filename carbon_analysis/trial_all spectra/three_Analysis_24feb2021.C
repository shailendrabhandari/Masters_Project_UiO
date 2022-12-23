//Example how to run Analysis
//./runAnalysis 0.5604 0.5573 0.5518 -13 4 15   beta values from Shailendra{GoldTarget}
//./runAnalysis 0.5338 0.5305 0.5242 -13 4 15   beta values from Shailendra{Carbon Target}
// with ./runAnalsysis beta1 beta2 beta3 timeCutLow timeCutHigh maxaddbackdistance betaDiffLow betaDiffHigh
#include "Event.h"
#include "TChain.h"

#include <vector>
#include "TTree.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"

#include "TROOT.h"
#include "TRint.h"
#include "TVector3.h"

#include "TMath.h"

#include <fstream>
#include <iostream>
#include <strstream>
#include <stdlib.h>

#include <TFile.h>
#include "TCutG.h"
#include <search.h>

using namespace std;

const int NUMBEROFDALICRYSTALS = 186;
const int NUMBEROFDALIADDBACKCRYSTALS = 30;
//float ADDBACKDISTANCE = 10;

//For addback
int fAddbackTable[NUMBEROFDALICRYSTALS][NUMBEROFDALIADDBACKCRYSTALS];
int fNumberOfAddbackPartners[NUMBEROFDALICRYSTALS];
bool crystalUsedForAddback[NUMBEROFDALICRYSTALS] = {false};

float fTheta[NUMBEROFDALICRYSTALS];
float fPosX[NUMBEROFDALICRYSTALS];
float fPosY[NUMBEROFDALICRYSTALS];
float fPosZ[NUMBEROFDALICRYSTALS];

int fDaliFold   = 0;//Fold
int fDaliFoldTa = 0;//Fold at Target time
int fDaliMult   = 0;//multiplicity
int fDaliMultTa = 0;//multiplicity at target time

//This is only for Add-back based analysis
struct dali{
    int id;
    int layer;
    float theta;   //Angle for Doppler Correction
    float x;
    float y;
    float z;
    float e;       //Energy
    float t;       //Time
    bool  ttrue;   //Bool if time is true
    float dopp1, dopp2, dopp3; //Doppler energy. Three doppler corrections for the three betas
    float doppwa1, doppwa2, doppwa3; //Doppler energy with true multiplicity and addback.
    float idwa;
    //float twa;
    //float ttruewa;
}fDali[NUMBEROFDALICRYSTALS];

///Doppler correctionformula
//__________________________________________________________
Double_t DopplerCorrect(Double_t beta, Double_t theta, Double_t energy) {
    
    if(energy <= 0.) return -1;
    else return energy * (1-beta*TMath::Cos(theta))/TMath::Sqrt((1.0-beta*beta));
}


//__________________________________________________________
bool IncludeAddbackTable(int detid[2],Double_t maxDistance, TVector3 fPos[]) {
    
    float distance = TMath::Sqrt(TMath::Power(fPos[detid[0]].X()-fPos[detid[1]].X(),2) +
                                 TMath::Power(fPos[detid[0]].Y()-fPos[detid[1]].Y(),2) +
                                 TMath::Power(fPos[detid[0]].Z()-fPos[detid[1]].Z(),2));
    
    //cout<<"Distance: "<<distance<<endl;
    if( distance > maxDistance ) return false;
    else return true;
}

//__________________________________________________________
void CreateAddBackTable(Double_t maxDistance) {
    
    cout<<"Creating the add-back table"<<endl;
    
    TVector3 detPos[NUMBEROFDALICRYSTALS];
    float x,y,z;
    int number;
    float angle;
    
    FILE *fdetPos  = fopen("./input/AverageInteractionPoint.out","r");
    char dummy[100];
    
    for(int i=0;i<NUMBEROFDALICRYSTALS;i++) {
        fgets(dummy,100,fdetPos);
        sscanf(dummy,"%f %f %f %i %i %f",&x,&y,&z,&number,&number,&angle);
        detPos[i].SetX(x);
        detPos[i].SetY(y);
        detPos[i].SetZ(z);
        fTheta[i] = angle*3.14159/180.;
        fPosX[i] = x;
        fPosY[i] = y;
        fPosZ[i] = z;
        //cout<<detPos[i].X()<<" "<<detPos[i].Y()<<" "<<detPos[i].Z()<<endl;
    }
    FILE *fAddbackTableOut = fopen("./input/AddbackTable.out","w");
    
    int detid[2];
    bool inTable;
    int counter = 0;
    for(int i=0;i<NUMBEROFDALICRYSTALS;i++) {
        fNumberOfAddbackPartners[i] = 0;
        fprintf(fAddbackTableOut," %i",i);
        
        for(int j=i+1;j<NUMBEROFDALICRYSTALS;j++) {
            detid[0] = i;
            detid[1] = j;
            
            inTable = IncludeAddbackTable(detid,maxDistance,detPos);
            
            if(inTable && counter< NUMBEROFDALIADDBACKCRYSTALS) {
                fprintf(fAddbackTableOut," %i",j);
                fAddbackTable[i][counter] = j;
                fNumberOfAddbackPartners[i]++;
                counter++;
            }
            if(counter == NUMBEROFDALIADDBACKCRYSTALS)  { //Too many detectors
                cout<<"You have to increase the variable NUMBEROFDALI2ADDBACKCRYSTALS!!!"<<endl;
                //STD::abort();
            }
        }
        counter = 0;
        fprintf(fAddbackTableOut," \n");
    }
    fclose(fAddbackTableOut);
}

//__________________________________________________________

int CompareByEnergy(const void *element1, const void *element2) {
    return((dali*)element1)->e < ((dali*)element2)->e? 1: -1;
    /////CompareByEnergy()>0: el1 goes after el2 --- CompareByEnergy()<0: el1 goes before el2 ---
}

void SortData(int count, dali fDali[]) {
    qsort(fDali,count,sizeof(dali),CompareByEnergy);
}

//__________________________________________________________
int main(int argc, char **argv) {
    
    if(argc!=7) {
        std::cout << "Invalid arguments..." << endl;
        return 1;
    }
    
    double beta1 = atof(argv[1]);
    double beta2 = atof(argv[2]);
    double beta3 = atof(argv[3]);
    
    //cout<<"The beta values are: "<<beta1<<" "<<endl;
    //"<<beta2<<" "<<beta3<<"
    
    double timeCutLow = atof(argv[4]);
    double timeCutHigh = atof(argv[5]);
    
    //Create the add back table
    CreateAddBackTable(atof(argv[6]));
    
    //double betaDiffLow = atof(argv[7]);
    //double betaDiffHigh = atof(argv[8]);
    
    // Create interactive interface
    TRint *theApp = new TRint("ROOT example", &argc, argv, NULL, 0);  //close the line "theApp->Run();" at the end of the code
    //__________________________________________________________
   
    
    //PID cuts:
    TFile *BRcut[3];
    TFile *ZDcut[3];
    
    //Cuts
    
    char name[100];
    
    for(int i=0;i<3;i++) {
        sprintf(name,"./cut/brcut%i.root",i);
        BRcut[i] = new TFile(name,"READ");
    }
    
    for(int i=0;i<3;i++) {
        sprintf(name,"./cut/zdcut%i.root",i);
        ZDcut[i] = new TFile(name,"READ");
    }
    
    TCutG *brcut[3];
    TCutG *zdcut[3];
    
    for(int i=0;i<3;i++) {
        sprintf(name,"brcut%i",i);
        BRcut[i]->GetObject("CUTG",brcut[i]);
    }
    
    for(int i=0;i<3;i++) {
        sprintf(name,"zdcut%i",i);
        ZDcut[i]->GetObject("CUTG",zdcut[i]);
    }
    
    Event *Tree = new Event();
    
    Long64_t nentries = Tree->fChain->GetEntriesFast();
    
    TFile *rootfile = new TFile("three_Analysis_carbonTaget.root","RECREATE");
    TTree *tree = new TTree("tree","tree");
    rootfile->cd();
    
    //________________________________________________________________________
    

    //Specific variables;
    int daliMult;
    int daliTimeTrueMult;
    int daliFold;
    int daliTimeTrueFold;
    
    Double_t DopplerAngle;
    
    int countingFilling = 0;
    //________________________________________________________________________
    
    //Define spectra:
    int minBin = 0;
    int maxBin = 6000;
    int binning = 25;
    int numBin = (maxBin-minBin)/binning;
    double maxE = 6000.;
    double minE = 0.;
    double binE = 10.;
    double numbinE = (maxE - minE)/binE;
    //________________________________________________________________________
    
    //Spectra:

//////TH1F *hEnergy_mult1_wa_77Zn = new TH1F("hEnergy_mult1_wa_77Zn","hEnergy_mult1_wa_77Zn",n1umBin,1minBin,1maxBin); //////change binning on 09 March 2021.
    
TH1F *hEnergy_mult1_wa_77Zn = new TH1F("hEnergy_mult1_wa_77Zn","hEnergy_mult1_wa_77Zn",240,0,6000);
TH1F *hEnergy_mult12_wa_77Zn = new TH1F("hEnergy_mult12_wa_77Zn","hEnergy_mult12_wa_77Zn",240,0,6000);
TH1F *hEnergy_mult123_wa_77Zn = new TH1F("hEnergy_mult123_wa_77Zn","hEnergy_mult123_wa_77Zn",240,0,6000);
TH1F *hEnergy_mult1_wa_78Zn = new TH1F("hEnergy_mult1_wa_78Zn","hEnergy_mult1_wa_78Zn",240,0,6000);
TH1F *hEnergy_mult2_wa_78Zn = new TH1F("hEnergy_mult2_wa_78Zn","hEnergy_mult2_wa_78Zn",240,0,6000);
TH1F *hEnergy_mult3_wa_78Zn = new TH1F("hEnergy_mult3_wa_78Zn","hEnergy_mult3_wa_78Zn",240,0,6000);
TH1F *hEnergy_mult4_wa_78Zn = new TH1F("hEnergy_mult4_wa_78Zn","hEnergy_mult4_wa_78Zn",240,0,6000);
TH1F *hEnergy_mult12_wa_78Zn = new TH1F("hEnergy_mult12_wa_78Zn","hEnergy_mult12_wa_78Zn",240,0,6000);
TH1F *hEnergy_mult123_wa_78Zn = new TH1F("hEnergy_mult123_wa_78Zn","hEnergy_mult123_wa_78Zn",240,0,6000);
TH1F *hEnergy_mult1234_wa_78Zn = new TH1F("hEnergy_mult1234_wa_78Zn","hEnergy_mult1234_wa_78Zn",240,0,6000);
TH1F *hEnergy_mult12345_wa_78Zn = new TH1F("hEnergy_mult12345_wa_78Zn","hEnergy_mult12345_wa_78Zn",240,0,6000);
TH1F *hEnergy_mult123456_wa_78Zn = new TH1F("hEnergy_mult123456_wa_78Zn","hEnergy_mult123456_wa_78Zn",240,0,6000);
TH1F *hEnergy_mult1_wa_79Zn = new TH1F("hEnergy_mult1_wa_79Zn","hEnergy_mult1_wa_79Zn",240,0,6000);
TH1F *hEnergy_mult12_wa_79Zn = new TH1F("hEnergy_mult12_wa_79Zn","hEnergy_mult12_wa_79Zn",240,0,6000);
TH1F *hEnergy_mult123_wa_79Zn = new TH1F("hEnergy_mult123_wa_79Zn","hEnergy_mult123_wa_79Zn",240,0,6000);

    TH1F *hEnergy_mult1_77Zn = new TH1F("hEnergy_mult1_77Zn","hEnergy_mult1_77Zn",240,0,6000);
    TH1F *hEnergy_mult12_77Zn = new TH1F("hEnergy_mult12_77Zn","hEnergy_mult12_77Zn",240,0,6000);
    TH1F *hEnergy_mult123_77Zn = new TH1F("hEnergy_mult123_77Zn","hEnergy_mult123_77Zn",240,0,6000);
    TH1F *hEnergy_mult1_78Zn = new TH1F("hEnergy_mult1_78Zn","hEnergy_mult1_78Zn",240,0,6000);
    TH1F *hEnergy_mult2_78Zn = new TH1F("hEnergy_mult2_78Zn","hEnergy_mult2_78Zn",240,0,6000);
    TH1F *hEnergy_mult3_78Zn = new TH1F("hEnergy_mult3_78Zn","hEnergy_mult3_78Zn",240,0,6000);
    TH1F *hEnergy_mult4_78Zn = new TH1F("hEnergy_mult4_78Zn","hEnergy_mult4_78Zn",240,0,6000);
    TH1F *hEnergy_mult12_78Zn = new TH1F("hEnergy_mult12_78Zn","hEnergy_mult12_78Zn",240,0,6000);
    TH1F *hEnergy_mult123_78Zn = new TH1F("hEnergy_mult123_78Zn","hEnergy_mult123_78Zn",240,0,6000);
    TH1F *hEnergy_mult1234_78Zn = new TH1F("hEnergy_mult1234_78Zn","hEnergy_mult1234_78Zn",240,0,6000);
    TH1F *hEnergy_mult12345_78Zn = new TH1F("hEnergy_mult12345_78Zn","hEnergy_mult12345_78Zn",240,0,6000);
    TH1F *hEnergy_mult123456_78Zn = new TH1F("hEnergy_mult123456_78Zn","hEnergy_mult123456_78Zn",240,0,6000);
    TH1F *hEnergy_mult1_79Zn = new TH1F("hEnergy_mult1_79Zn","hEnergy_mult1_79Zn",240,0,6000);
    TH1F *hEnergy_mult12_79Zn = new TH1F("hEnergy_mult12_79Zn","hEnergy_mult12_79Zn",240,0,6000);
    TH1F *hEnergy_mult123_79Zn = new TH1F("hEnergy_mult123_79Zn","hEnergy_mult123_79Zn",240,0,6000);
    
    
    
    TH2F *h2D_Energy_mult1_wa_79Zn = new TH2F("h2D_Energy_mult1_wa_79Zn","h2D_Energy_mult1_wa_79Zn",186,0,186,240,0,6000);
    TH2F *h2D_Energy_mult12_wa_79Zn = new TH2F("h2D_Energy_mult12_wa_79Zn","h2D_Energy_mult12_wa_79Zn",186,0,186,240,0,6000);
    TH2F *h2D_Energy_mult123_wa_79Zn = new TH2F("h2D_Energy_mult123_wa_79Zn","h2D_Energy_mult123_wa_79Zn",186,0,186,240,0,6000);
    TH2F *h2D_Energy_mult1_wa_77Zn = new TH2F("h2D_Energy_mult1_wa_77Zn","h2D_Energy_mult1_wa_77Zn",186,0,186,240,0,6000);
    TH2F *h2D_Energy_mult12_wa_77Zn = new TH2F("h2D_Energy_mult12_wa_77Zn","h2D_Energy_mult12_wa_77Zn",186,0,186,240,0,6000);
    TH2F *h2D_Energy_mult123_wa_77Zn = new TH2F("h2D_Energy_mult123_wa_77Zn","h2D_Energy_mult123_wa_77Zn",186,0,186,240,0,6000);
    TH2F *h2D_Energy_mult1_wa_78Zn = new TH2F("h2D_Energy_mult1_wa_78Zn","h2D_Energy_mult1_wa_78Zn",186,0,186,240,0,6000);
    TH2F *h2D_Energy_mult2_wa_78Zn = new TH2F("h2D_Energy_mult2_wa_78Zn","h2D_Energy_mult2_wa_78Zn",186,0,186,240,0,6000);
    TH2F *h2D_Energy_mult3_wa_78Zn = new TH2F("h2D_Energy_mult3_wa_78Zn","h2D_Energy_mult3_wa_78Zn",186,0,186,240,0,6000);
    TH2F *h2D_Energy_mult4_wa_78Zn = new TH2F("h2D_Energy_mult4_wa_78Zn","h2D_Energy_mult4_wa_78Zn",186,0,186,240,0,6000);
    TH2F *h2D_Energy_mult12_wa_78Zn = new TH2F("h2D_Energy_mult12_wa_78Zn","h2D_Energy_mult12_wa_78Zn",186,0,186,240,0,6000);
    TH2F *h2D_Energy_mult123_wa_78Zn = new TH2F("h2D_Energy_mult123_wa_78Zn","h2D_Energy_mult123_wa_78Zn",186,0,186,240,0,6000);
    TH2F *h2D_Energy_mult1234_wa_78Zn = new TH2F("h2D_Energy_mult1234_wa_78Zn","h2D_Energy_mult1234_wa_78Zn",186,0,186,240,0,6000);
    TH2F *h2D_Energy_mult12345_wa_78Zn = new TH2F("h2D_Energy_mult12345_wa_78Zn","h2D_Energy_mult12345_wa_78Zn",186,0,186,240,0,6000);
    TH2F *h2D_Energy_mult123456_wa_78Zn = new TH2F("h2D_Energy_mult123456_wa_78Zn","h2D_Energy_mult123456_wa_78Zn",186,0,186,240,0,6000);
    
    TH2F *h2D_Energy_mult1_79Zn = new TH2F("h2D_Energy_mult1_79Zn","h2D_Energy_mult1_79Zn",186,0,186,240,0,6000);
    TH2F *h2D_Energy_mult12_79Zn = new TH2F("h2D_Energy_mult12_79Zn","h2D_Energy_mult12_79Zn",186,0,186,240,0,6000);
    TH2F *h2D_Energy_mult123_79Zn = new TH2F("h2D_Energy_mult123_79Zn","h2D_Energy_mult123_79Zn",186,0,186,240,0,6000);
    TH2F *h2D_Energy_mult1_77Zn = new TH2F("h2D_Energy_mult1_77Zn","h2D_Energy_mult1_77Zn",186,0,186,240,0,6000);
    TH2F *h2D_Energy_mult12_77Zn = new TH2F("h2D_Energy_mult12_77Zn","h2D_Energy_mult12_77Zn",186,0,186,240,0,6000);
    TH2F *h2D_Energy_mult123_77Zn = new TH2F("h2D_Energy_mult123_77Zn","h2D_Energy_mult123_77Zn",186,0,186,240,0,6000);
    TH2F *h2D_Energy_mult1_78Zn = new TH2F("h2D_Energy_mult1_78Zn","h2D_Energy_mult1_78Zn",186,0,186,240,0,6000);
    TH2F *h2D_Energy_mult2_78Zn = new TH2F("h2D_Energy_mult2_78Zn","h2D_Energy_mult2_78Zn",186,0,186,240,0,6000);
    TH2F *h2D_Energy_mult3_78Zn = new TH2F("h2D_Energy_mult3_78Zn","h2D_Energy_mult3_78Zn",186,0,186,240,0,6000);
    TH2F *h2D_Energy_mult4_78Zn = new TH2F("h2D_Energy_mult4_78Zn","h2D_Energy_mult4_78Zn",186,0,186,240,0,6000);
    TH2F *h2D_Energy_mult12_78Zn = new TH2F("h2D_Energy_mult12_78Zn","h2D_Energy_mult12_78Zn",186,0,186,240,0,6000);
    TH2F *h2D_Energy_mult123_78Zn = new TH2F("h2D_Energy_mult123_78Zn","h2D_Energy_mult123_78Zn",186,0,186,240,0,6000);
    TH2F *h2D_Energy_mult1234_78Zn = new TH2F("h2D_Energy_mult1234_78Zn","h2D_Energy_mult1234_78Zn",186,0,186,240,0,6000);
    TH2F *h2D_Energy_mult12345_78Zn = new TH2F("h2D_Energy_mult12345_78Zn","h2D_Energy_mult12345_78Zn",186,0,186,240,0,6000);
    TH2F *h2D_Energy_mult123456_78Zn = new TH2F("h2D_Energy_mult123456_78Zn","h2D_Energy_mult123456_78Zn",186,0,186,240,0,6000);
    
    
     TH1F *hEnergy_mult_All_77Zn = new TH1F("hEnergy_mult_All_77Zn","hEnergy_mult_All_77Zn",240,0,6000);
     TH1F *hEnergy_fold12_77Zn = new TH1F("hEnergy_fold12_77Zn","hEnergy_fold12_77Zn",240,0,6000);
     TH1F *hEnergy_fold1_77Zn = new TH1F("hEnergy_fold1_77Zn","hEnergy_fold1_77Zn",240,0,6000);
     TH1F *hEnergy_mult_All_78Zn = new TH1F("hEnergy_mult_All_78Zn","hEnergy_mult_All_78Zn",240,0,6000);
     TH1F *hEnergy_fold12_78Zn = new TH1F("hEnergy_fold12_78Zn","hEnergy_fold12_78Zn",240,0,6000);
     TH1F *hEnergy_fold1_78Zn = new TH1F("hEnergy_fold1_78Zn","hEnergy_fold1_78Zn",240,0,6000);
     TH1F *hEnergy_mult_All_79Zn = new TH1F("hEnergy_mult_All_79Zn","hEnergy_mult_All_79Zn",240,0,6000);
     TH1F *hEnergy_fold12_79Zn = new TH1F("hEnergy_fold12_79Zn","hEnergy_fold12_79Zn",240,0,6000);
     TH1F *hEnergy_fold1_79Zn = new TH1F("hEnergy_fold1_79Znnew","hEnergy_fold1_79Znnew",240,0,6000);
    
    TH2F *hBR_aoq_z_full = new TH2F("hBR_aoq_z_full","hBR_aoq_z_full",2000,2.45,2.75,2000,27,34);
    TH2F *hZD_aoq_z_full = new TH2F("hZD_aoq_z_full","hZD_aoq_z_full",2000,2.35,2.75,2000,26,34);
    TH2F *hBR_aoq_z_DS1_77Zn = new TH2F("hBR_aoq_z_DS1_77Zn","hBR_aoq_z_DS1_77Zn",2000,2.45,2.75,2000,27,34);
    TH2F *hBR_aoq_z_DS1_78Zn = new TH2F("hBR_aoq_z_DS1_78Zn","hBR_aoq_z_DS1_78Zn",2000,2.45,2.75,2000,27,34);
    TH2F *hBR_aoq_z_DS1_79Zn = new TH2F("hBR_aoq_z_DS1_79Zn","hBR_aoq_z_DS1_79Zn",2000,2.45,2.75,2000,27,34);
    
    TH2F *hZD_aoq_z_F11_77Zn = new TH2F("hZD_aoq_z_F11_77Zn","hZD_aoq_z_F11_77Zn",2000,2.35,2.75,2000,26,34);
    TH2F *hZD_aoq_z_F11_78Zn = new TH2F("hZD_aoq_z_F11_78Zn","hZD_aoq_z_F11_78Zn",2000,2.35,2.75,2000,26,34);
    TH2F *hZD_aoq_z_F11_79Zn = new TH2F("hZD_aoq_z_F11_79Zn","hZD_aoq_z_F11_79Zn",2000,2.35,2.75,2000,26,34);
  
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //Start looping through data;
    
    Long64_t i=0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        
        if(jentry%10000 ==0) cout << jentry <<"/"<<nentries<<" Events DONE!"<<endl;
        
        Tree->fChain->GetEvent(jentry);
    
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        //Sorting the DALI2 data
        fDaliFold   = 0;//Fold
        fDaliFoldTa = 0;//Fold
        fDaliMultTa = 0;//multiplicity
        for(int j=0;j<NUMBEROFDALICRYSTALS;j++) {
            crystalUsedForAddback[j] = false;
        }
        
        //cout<<"DALINaI: "<<Tree->DALINaI_<<endl;
        
        for(int j=0;j<Tree->DALINaI_;j++){
            fDali[j].id         = Tree->DALINaI_id[j]-1;
            fDali[j].layer      = Tree->DALINaI_layer[j];
            fDali[j].theta      = fTheta[fDali[j].id];
            fDali[j].x          = fPosX[fDali[j].id];
            fDali[j].y          = fPosY[fDali[j].id];
            fDali[j].z          = fPosZ[fDali[j].id];
            fDali[j].e          = Tree->DALINaI_fEnergy[j];
            fDali[j].t          = Tree->DALINaI_fTimeOffseted[j];
            if(fDali[j].e>0){
                
                
                
                //excluded detectors
                
                //if (fDali[j].id < 52) fDali[j].e = -999.; this cut is for the physical runs to supress atomic bg
                //if (fDali[j].id == 97) fDali[j].e = -999.;Victor excluded 97 and 98  too but I could not find a problem with them
                if (fDali[j].id < 60) fDali[j].e = -999.;
                if (fDali[j].id == 10) fDali[j].e = -999.;
                if (fDali[j].id == 23) fDali[j].e = -999.;
                if (fDali[j].id == 31) fDali[j].e = -999.;
                if (fDali[j].id == 48) fDali[j].e = -999.;
                if (fDali[j].id == 63) fDali[j].e = -999.;
                if (fDali[j].id == 127) fDali[j].e = -999.;
                if (fDali[j].id == 128) fDali[j].e = -999.;
                if (fDali[j].id == 141) fDali[j].e = -999.;
                if (fDali[j].id == 167) fDali[j].e = -999.;
                
                
                fDali[j].dopp1 = DopplerCorrect(beta1,fDali[j].theta,fDali[j].e);
                fDali[j].dopp2 = DopplerCorrect(beta2,fDali[j].theta,fDali[j].e);
                fDali[j].dopp3 = DopplerCorrect(beta3,fDali[j].theta,fDali[j].e);
                
                if(fDali[j].t>timeCutLow-500&&fDali[j].t<timeCutHigh+500)fDaliFold++;
                
                if(fDali[j].t>timeCutLow && fDali[j].t<timeCutHigh){
                    fDali[j].ttrue = true;
                    fDaliFoldTa++;
                }
                else fDali[j].ttrue = false;
            }
            else {
                fDali[j].dopp1 = -999.;
                fDali[j].dopp2 = -999.;
                fDali[j].dopp3 = -999.;
                //fDali[j].dopp[3] = -999.;
                //fDali[j].ttrue   = false;
            }
        }
        
        for(int j=Tree->DALINaI_;j<NUMBEROFDALICRYSTALS;j++){
            fDali[j].id         = -1;
            fDali[j].layer      = -1;
            fDali[j].theta      = -1;
            fDali[j].x          = -999;
            fDali[j].y          = -999;
            fDali[j].z          = -999;
            fDali[j].e          = -999;
            fDali[j].t          = -999;
            fDali[j].ttrue      = false;
            fDali[j].dopp1    = -999;
            fDali[j].dopp2    = -999;
            fDali[j].dopp3    = -999;
            //if (fDali[j].id <= 61) fDali[j].e = -999.;
            
        }
           
        //if(fDali[0].e>0) /////Pieter asked me to close this line March 2019
        SortData(fDaliFold,fDali); /// A. Junglaus burayi da kapattirmisti alttaki satiri actirmisti.
        //SortData(Tree->DALINaI_,fDali);
        
        //Going to the add-back:
        float dummyEnergy[NUMBEROFDALICRYSTALS][6] = {{0.}};
        //Making add-back and true multiplicity:
        //The Energy must be sorted already according to the highest detected one.
        //cout<<"Starting addback"<<endl;
        for(int i=0;i<fDaliFold;i++){
            
            if(crystalUsedForAddback[fDali[i].id] == true || fDali[i].ttrue == false) continue;
            
            dummyEnergy[fDaliMultTa][0] = DopplerCorrect(beta1,fDali[i].theta,fDali[i].e);
            dummyEnergy[fDaliMultTa][1] = DopplerCorrect(beta2,fDali[i].theta,fDali[i].e);
            dummyEnergy[fDaliMultTa][2] = DopplerCorrect(beta3,fDali[i].theta,fDali[i].e);
            
            crystalUsedForAddback[fDali[i].id]=true;
            fDali[fDaliMultTa].idwa = fDali[i].id;
            
            for(int j = i+1;j<fDaliFold;j++)  {
                if(crystalUsedForAddback[fDali[j].id]==false && fDali[j].ttrue==true)  {
                    for(int k = 0;k<fNumberOfAddbackPartners[fDali[i].id] ;k++) {
                        if(fDali[j].id == fAddbackTable[fDali[i].id][k+1])  {
                            
                            crystalUsedForAddback[fDali[j].id]=true;
                            
                            dummyEnergy[fDaliMultTa][0] += DopplerCorrect(beta1,fDali[i].theta,fDali[j].e);
                            dummyEnergy[fDaliMultTa][1] += DopplerCorrect(beta2,fDali[i].theta,fDali[j].e);
                            dummyEnergy[fDaliMultTa][2] += DopplerCorrect(beta3,fDali[i].theta,fDali[j].e);
                            
                        }
                    }
                }
            }
            fDaliMultTa++;
        }
        for(int i = 0;i<fDaliMultTa;i++) {
            fDali[i].doppwa1 = dummyEnergy[i][0];
            fDali[i].doppwa2 = dummyEnergy[i][1];
            fDali[i].doppwa3 = dummyEnergy[i][2];
        }
        for(int i = fDaliMultTa;i<NUMBEROFDALICRYSTALS;i++) {
            fDali[i].doppwa1 = -999;
             fDali[i].doppwa2 = -999;
             fDali[i].doppwa3 = -999;
            fDali[i].idwa      = -999;
        }
        
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        //Starting with the analysis conditions
        bool brzn77 = false;
        bool brzn78 = false;
        bool brzn79 = false;
        
        bool zdzn77 = false;
        bool zdzn78 = false;
        bool zdzn79 = false;
        
        //BRS
        if(brcut[0]->IsInside(Tree->aoqc[2],Tree->zetc[2])) brzn77 = true;
        if(brcut[1]->IsInside(Tree->aoqc[2],Tree->zetc[2])) brzn78 = true;
        if(brcut[2]->IsInside(Tree->aoqc[2],Tree->zetc[2])) brzn79 = true;
        
        //ZDS
        if(zdcut[0]->IsInside(Tree->aoqc[5],Tree->zetc[5])) zdzn77 = true;
        if(zdcut[1]->IsInside(Tree->aoqc[5],Tree->zetc[5])) zdzn78 = true;
        if(zdcut[2]->IsInside(Tree->aoqc[5],Tree->zetc[5])) zdzn79 = true;
        
        /*//Beta diff
         float beta_diff = Tree->BigRIPSBeam_beta[0]-Tree->BigRIPSBeam_beta[3];
         if(brzn77&&zdzn77) h_beta_diff[0]->Fill(beta_diff);
         if(brzn78&&zdzn78) h_beta_diff[1]->Fill(beta_diff);
         if(brzn79&&zdzn79) h_beta_diff[2]->Fill(beta_diff);*/
        
        //Trigger register information
        bool DSB = false;
        bool F11 = false;
        bool DaliTrigger = false;
        
     /*   //Need to check the settings
        if(Tree->fbit==1||Tree->fbit==3||Tree->fbit==7) DSB = true;
        if(Tree->fbit==7) DaliTrigger = true;
        if(Tree->fbit==3) F11 = true;  */

        //Need to check the settings
        if(Tree->EventInfo_fBit[0]==1||Tree->EventInfo_fBit[0]==3||Tree->EventInfo_fBit[0]==7) DSB = true;
        if(Tree->EventInfo_fBit[0]==7) DaliTrigger = true;
        if(Tree->EventInfo_fBit[0]==3) F11 = true;



        
       
        //Getting the statistics for the cross-section:
        //Beam
        if(DSB && brzn77)
            hBR_aoq_z_DS1_77Zn->Fill(Tree->aoqc[2],Tree->zetc[2]);
        if(DSB && brzn78)
            hBR_aoq_z_DS1_78Zn->Fill(Tree->aoqc[2],Tree->zetc[2]);
        if(DSB && brzn79)
            hBR_aoq_z_DS1_79Zn->Fill(Tree->aoqc[2],Tree->zetc[2]);


////////////////////////Gamma-ray sepctra /////////////////////
       
        
        /////With Addback//////
        
        //////77Zn//////////////
        for(int j=0;j<fDaliMultTa;j++){
            if((fDali[j].t)>-13 && (fDali[j].t)<4) {
                if((fDaliMultTa == 1 )&& brzn77 && zdzn77 && fDali[j].idwa >= 60 )
                    hEnergy_mult1_wa_77Zn->Fill(fDali[j].doppwa1);
                
                if((fDaliMultTa<=2) && brzn77 && zdzn77 && fDali[j].idwa >= 60)
                    hEnergy_mult12_wa_77Zn->Fill(fDali[j].doppwa1);
                
                if((fDaliMultTa<=3) && brzn77 && zdzn77 && fDali[j].idwa >= 60)
                    hEnergy_mult123_wa_77Zn->Fill(fDali[j].doppwa1);
            }
      //////78Zn//////////////
            if((fDali[j].t)>-13 && (fDali[j].t)<4) {
        if((fDaliMultTa == 1 )&& brzn78 && zdzn78 && fDali[j].idwa >= 60 )
            hEnergy_mult1_wa_78Zn->Fill(fDali[j].doppwa2);
        
        if((fDaliMultTa<=2)&& brzn78 && zdzn78 && fDali[j].idwa >= 60)
            hEnergy_mult12_wa_78Zn->Fill(fDali[j].doppwa2);
        
        if((fDaliMultTa<=3) && brzn78 && zdzn78 && fDali[j].idwa >= 60)
            hEnergy_mult123_wa_78Zn->Fill(fDali[j].doppwa2);

        if((fDaliMultTa<=4) && brzn78 && zdzn78 && fDali[j].idwa >= 60)
            hEnergy_mult1234_wa_78Zn->Fill(fDali[j].doppwa2);

        if((fDaliMultTa<=5) && brzn78 && zdzn78 && fDali[j].idwa >= 60)
            hEnergy_mult12345_wa_78Zn->Fill(fDali[j].doppwa2);

        if((fDaliMultTa<=6) && brzn78 && zdzn78 && fDali[j].idwa >= 60)
            hEnergy_mult123456_wa_78Zn->Fill(fDali[j].doppwa2);




        if((fDaliMultTa == 2 )&& brzn78 && zdzn78 && fDali[j].idwa >= 60 )
            hEnergy_mult2_wa_78Zn->Fill(fDali[j].doppwa2);

        if((fDaliMultTa == 3 )&& brzn78 && zdzn78 && fDali[j].idwa >= 60 )
            hEnergy_mult3_wa_78Zn->Fill(fDali[j].doppwa2);

        if((fDaliMultTa == 4 )&& brzn78 && zdzn78 && fDali[j].idwa >= 60 )
            hEnergy_mult4_wa_78Zn->Fill(fDali[j].doppwa2);




            }
        //////79Zn//////////////
            if((fDali[j].t)>-13 && (fDali[j].t)<4) {
                if((fDaliMultTa == 1 )&& brzn79 && zdzn79 && fDali[j].idwa >= 60 )
                    hEnergy_mult1_wa_79Zn->Fill(fDali[j].doppwa3);
                
                if((fDaliMultTa<=2) && brzn79 && zdzn79 && fDali[j].idwa >= 60)
                    hEnergy_mult12_wa_79Zn->Fill(fDali[j].doppwa3);
                
                if((fDaliMultTa<=3) && brzn79 && zdzn79 && fDali[j].idwa >= 60)
                    hEnergy_mult123_wa_79Zn->Fill(fDali[j].doppwa3);
            }
        }
    
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       ////Without Addback
            //////77Zn//////////////
        for(int j=0;j<fDaliFold;j++){
            if(fDali[j].ttrue) {
      
            if((fDali[j].t)>-13 && (fDali[j].t)<4) {
                        if((fDaliMultTa == 1) && brzn77 && zdzn77 && fDali[j].id >= 60)
                            hEnergy_mult1_77Zn->Fill(fDali[j].dopp1);
                        if((fDaliMultTa<=2)&& brzn77 && zdzn77 && fDali[j].id >= 60)
                            hEnergy_mult12_77Zn->Fill(fDali[j].dopp1);
                        if((fDaliMultTa<=3) && brzn77 && zdzn77 && fDali[j].id >= 60)
                            hEnergy_mult123_77Zn->Fill(fDali[j].dopp1);
                    }
        //////78Zn//////////////
        if((fDali[j].t)>-13 && (fDali[j].t)<4) {
                    if((fDaliMultTa == 1) && brzn78 && zdzn78 && fDali[j].id >= 60)
                        hEnergy_mult1_78Zn->Fill(fDali[j].dopp2);
            
                    if((fDaliMultTa<=2) && brzn78 && zdzn78 && fDali[j].id >= 60)
                        hEnergy_mult12_78Zn->Fill(fDali[j].dopp2);
                    
                    if((fDaliMultTa<=3) && brzn78 && zdzn78 && fDali[j].id >= 60)
                        hEnergy_mult123_78Zn->Fill(fDali[j].dopp2);
                    
                    if((fDaliMultTa<=4) && brzn78 && zdzn78 && fDali[j].id >= 60)
                        hEnergy_mult1234_78Zn->Fill(fDali[j].dopp2);
                    
                    if((fDaliMultTa<=5) && brzn78 && zdzn78 && fDali[j].id >= 60)
                        hEnergy_mult12345_78Zn->Fill(fDali[j].dopp2);
                    
                    if((fDaliMultTa<=6) && brzn78 && zdzn78 && fDali[j].id >= 60)
                        hEnergy_mult123456_78Zn->Fill(fDali[j].dopp2);


                    if((fDaliMultTa == 2) && brzn78 && zdzn78 && fDali[j].id >= 60)
                        hEnergy_mult2_78Zn->Fill(fDali[j].dopp2);

                    if((fDaliMultTa == 3) && brzn78 && zdzn78 && fDali[j].id >= 60)
                        hEnergy_mult3_78Zn->Fill(fDali[j].dopp2);

                    if((fDaliMultTa == 4) && brzn78 && zdzn78 && fDali[j].id >= 60)
                        hEnergy_mult4_78Zn->Fill(fDali[j].dopp2);





                }
        //////79Zn//////////////
    if((fDali[j].t)>-13 && (fDali[j].t)<4) {
                    if((fDaliMultTa == 1) && brzn79 && zdzn79 && fDali[j].id >= 60)
                        hEnergy_mult1_79Zn->Fill(fDali[j].dopp3);
        
                    if((fDaliMultTa<=2) && brzn79 && zdzn79 && fDali[j].id >= 60)
                        hEnergy_mult12_79Zn->Fill(fDali[j].dopp3);
                    
                    if((fDaliMultTa<=3) && brzn79 && zdzn79 && fDali[j].id >= 60)
                        hEnergy_mult123_79Zn->Fill(fDali[j].dopp3);
                }
            }
        }
        
        
        
        /////2D histograms/////
        
        /////With Addback//////
        
        for(int j=0;j<fDaliMultTa;j++){
            //////77Zn//////////////
            if((fDali[j].t)>=-13 && (fDali[j].t)<=4) {
                if((fDaliMultTa == 1 )&& brzn77 && zdzn77)
                    h2D_Energy_mult1_wa_77Zn->Fill(fDali[j].idwa,fDali[j].doppwa2);
                
                if((fDaliMultTa<=2)&& brzn77 && zdzn77)
                    h2D_Energy_mult12_wa_77Zn->Fill(fDali[j].idwa,fDali[j].doppwa2);
                
                if((fDaliMultTa<=3) && brzn77 && zdzn77)
                    h2D_Energy_mult123_wa_77Zn->Fill(fDali[j].idwa,fDali[j].doppwa2);
            }
            //////78Zn//////////////
            if((fDali[j].t)>=-13 && (fDali[j].t)<=4) {
                if((fDaliMultTa == 1 )&& brzn78 && zdzn78)
                    h2D_Energy_mult1_wa_78Zn->Fill(fDali[j].idwa,fDali[j].doppwa3);
                
                if((fDaliMultTa<=2) && brzn78 && zdzn78)
                    h2D_Energy_mult12_wa_78Zn->Fill(fDali[j].idwa,fDali[j].doppwa3);
                
                if((fDaliMultTa<=3) && brzn78 && zdzn78)
                    h2D_Energy_mult123_wa_78Zn->Fill(fDali[j].idwa,fDali[j].doppwa3);

                if((fDaliMultTa<=4) && brzn78 && zdzn78)
                    h2D_Energy_mult1234_wa_78Zn->Fill(fDali[j].idwa,fDali[j].doppwa3);

                if((fDaliMultTa<=5) && brzn78 && zdzn78)
                    h2D_Energy_mult12345_wa_78Zn->Fill(fDali[j].idwa,fDali[j].doppwa3);

                if((fDaliMultTa<=6) && brzn78 && zdzn78)
                    h2D_Energy_mult123456_wa_78Zn->Fill(fDali[j].idwa,fDali[j].doppwa3);



                if((fDaliMultTa == 2 )&& brzn78 && zdzn78)
                    h2D_Energy_mult2_wa_78Zn->Fill(fDali[j].idwa,fDali[j].doppwa3);

                if((fDaliMultTa == 3 )&& brzn78 && zdzn78)
                    h2D_Energy_mult3_wa_78Zn->Fill(fDali[j].idwa,fDali[j].doppwa3);

                if((fDaliMultTa == 4 )&& brzn78 && zdzn78)
                    h2D_Energy_mult4_wa_78Zn->Fill(fDali[j].idwa,fDali[j].doppwa3);




                
            }
            //////79Zn//////////////
            if((fDali[j].t)>=-13 && (fDali[j].t)<=4) {
                if((fDaliMultTa == 1 )&& brzn79 && zdzn79 )
                    h2D_Energy_mult1_wa_79Zn->Fill(fDali[j].idwa,fDali[j].doppwa1);
                
                if((fDaliMultTa<=2) && brzn79 && zdzn79)
                    h2D_Energy_mult12_wa_79Zn->Fill(fDali[j].idwa,fDali[j].doppwa1);
                
                if((fDaliMultTa<=3) && brzn79 && zdzn79)
                    h2D_Energy_mult123_wa_79Zn->Fill(fDali[j].idwa,fDali[j].doppwa1);
            }
        }
        
        
        
        ///without Addback
        for(int j=0;j<fDaliFold;j++){
            if(fDali[j].ttrue) {
                //////77Zn//////////////
                if((fDali[j].t)>=-10 && (fDali[j].t)<=4) {
                    if((fDaliMultTa == 1) && brzn77 && zdzn77)
                        h2D_Energy_mult1_77Zn->Fill(fDali[j].id,fDali[j].dopp2);
                    
                    if((fDaliMultTa<=2) && brzn77 && zdzn77)
                        h2D_Energy_mult12_77Zn->Fill(fDali[j].id,fDali[j].dopp2);
                    
                    if((fDaliMultTa<=3) && brzn77 && zdzn77)
                        h2D_Energy_mult123_77Zn->Fill(fDali[j].id,fDali[j].dopp2);

                }
                //////78Zn//////////////
                if((fDali[j].t)>=-13 && (fDali[j].t)<=4) {
                    if((fDaliMultTa == 1) && brzn78 && zdzn78)
                        h2D_Energy_mult1_78Zn->Fill(fDali[j].id,fDali[j].dopp3);
                    
                    if((fDaliMultTa<=2) && brzn78 && zdzn78)
                        h2D_Energy_mult12_78Zn->Fill(fDali[j].id,fDali[j].dopp3);
                    
                    if((fDaliMultTa<=3) && brzn78 && zdzn78)
                        h2D_Energy_mult123_78Zn->Fill(fDali[j].id,fDali[j].dopp3);
                    
                    if((fDaliMultTa<=4) && brzn78 && zdzn78)
                        h2D_Energy_mult1234_78Zn->Fill(fDali[j].id,fDali[j].dopp3);
                    
                    if((fDaliMultTa<=5) && brzn78 && zdzn78)
                        h2D_Energy_mult12345_78Zn->Fill(fDali[j].id,fDali[j].dopp3);
                    
                    if((fDaliMultTa<=6) && brzn78 && zdzn78)
                        h2D_Energy_mult123456_78Zn->Fill(fDali[j].id,fDali[j].dopp3);

                    if((fDaliMultTa == 2) && brzn78 && zdzn78)
                        h2D_Energy_mult2_78Zn->Fill(fDali[j].id,fDali[j].dopp3);

                    if((fDaliMultTa == 3) && brzn78 && zdzn78)
                        h2D_Energy_mult3_78Zn->Fill(fDali[j].id,fDali[j].dopp3);

                    if((fDaliMultTa == 4) && brzn78 && zdzn78)
                        h2D_Energy_mult4_78Zn->Fill(fDali[j].id,fDali[j].dopp3);












                }
                //////79Zn//////////////
                if((fDali[j].t)>=-9 && (fDali[j].t)<=5) {
                    if((fDaliMultTa == 1) && brzn79 && zdzn79)
                        h2D_Energy_mult1_79Zn->Fill(fDali[j].id,fDali[j].dopp1);
                    
                    if((fDaliMultTa<=2)&& brzn79 && zdzn79)
                        h2D_Energy_mult12_79Zn->Fill(fDali[j].id,fDali[j].dopp1);
                    
                    if((fDaliMultTa<=3) && brzn79 && zdzn79)
                        h2D_Energy_mult123_79Zn->Fill(fDali[j].id,fDali[j].dopp1);
                }
            }
        }
        
        
        
     ///////   Fold spectra
        for(int j=0;j<fDaliFold;j++){
            if(fDali[j].ttrue) {
               //////77Zn//////////////
                if((fDali[j].t)>-13 && (fDali[j].t)<4 && fDali[j].id >= 60) {
         hEnergy_mult_All_77Zn->Fill(fDali[j].dopp1);
         if(fDaliFold<3)  hEnergy_fold12_77Zn->Fill(fDali[j].dopp1);
         if(fDaliFold==1) hEnergy_fold1_77Zn->Fill(fDali[j].dopp1);
                }
       //////78Zn//////////////
            if((fDali[j].t)>-13 && (fDali[j].t)<4 && fDali[j].id >= 60) {
                hEnergy_mult_All_78Zn->Fill(fDali[j].dopp2);
                if(fDaliFold<3)  hEnergy_fold12_78Zn->Fill(fDali[j].dopp2);
                if(fDaliFold==1) hEnergy_fold1_78Zn->Fill(fDali[j].dopp2);
                }
        //////79Zn//////////////
        if((fDali[j].t)>-13 && (fDali[j].t)<4 && fDali[j].id >= 60) {
            hEnergy_mult_All_79Zn->Fill(fDali[j].dopp3);
            if(fDaliFold<3)  hEnergy_fold12_79Zn->Fill(fDali[j].dopp3);
            if(fDaliFold==1) hEnergy_fold1_79Zn->Fill(fDali[j].dopp3);
                }
            }
        }
        //ZDS
        // hZD_aoq_z_full->Fill(Tree->aoqc[5],Tree->zetc[5]);
        if(F11 & zdzn77)
            hZD_aoq_z_F11_77Zn->Fill(Tree->aoqc[5],Tree->zetc[5]);
        if(F11 & zdzn78)
            hZD_aoq_z_F11_78Zn->Fill(Tree->aoqc[5],Tree->zetc[5]);
        if(F11 & zdzn79)
            hZD_aoq_z_F11_79Zn->Fill(Tree->aoqc[5],Tree->zetc[5]);
        
    }
    
hEnergy_mult1_wa_77Zn->Write();
hEnergy_mult12_wa_77Zn->Write();
hEnergy_mult123_wa_77Zn->Write();
hEnergy_mult1_77Zn->Write();
hEnergy_mult12_77Zn->Write();
hEnergy_mult123_77Zn->Write();
    
hEnergy_mult1_wa_78Zn->Write();
hEnergy_mult2_wa_78Zn->Write();
hEnergy_mult3_wa_78Zn->Write();
hEnergy_mult4_wa_78Zn->Write();
hEnergy_mult12_wa_78Zn->Write();
hEnergy_mult123_wa_78Zn->Write();
hEnergy_mult1234_wa_78Zn->Write();
hEnergy_mult12345_wa_78Zn->Write();
hEnergy_mult123456_wa_78Zn->Write();
hEnergy_mult1_78Zn->Write();
hEnergy_mult2_78Zn->Write();
hEnergy_mult3_78Zn->Write();
hEnergy_mult4_78Zn->Write();
hEnergy_mult12_78Zn->Write();
hEnergy_mult123_78Zn->Write();
hEnergy_mult1234_78Zn->Write();
hEnergy_mult12345_78Zn->Write();
hEnergy_mult123456_78Zn->Write();
    
hEnergy_mult1_wa_79Zn->Write();
hEnergy_mult12_wa_79Zn->Write();
hEnergy_mult123_wa_79Zn->Write();
hEnergy_mult1_79Zn->Write();
hEnergy_mult12_79Zn->Write();
hEnergy_mult123_79Zn->Write();
    
    
    h2D_Energy_mult1_77Zn->Write();
    h2D_Energy_mult12_77Zn->Write();
    h2D_Energy_mult123_77Zn->Write();
    h2D_Energy_mult1_78Zn->Write();
    h2D_Energy_mult2_78Zn->Write();
    h2D_Energy_mult3_78Zn->Write();
    h2D_Energy_mult4_78Zn->Write();
    h2D_Energy_mult12_78Zn->Write();
    h2D_Energy_mult123_78Zn->Write();
    h2D_Energy_mult1234_78Zn->Write();
    h2D_Energy_mult12345_78Zn->Write();
    h2D_Energy_mult123456_78Zn->Write();
    h2D_Energy_mult1_79Zn->Write();
    h2D_Energy_mult12_79Zn->Write();
    h2D_Energy_mult123_79Zn->Write();
    
    h2D_Energy_mult1_wa_77Zn->Write();
    h2D_Energy_mult12_wa_77Zn->Write();
    h2D_Energy_mult123_wa_77Zn->Write();
    h2D_Energy_mult1_wa_78Zn->Write();
    h2D_Energy_mult2_wa_78Zn->Write();
    h2D_Energy_mult3_wa_78Zn->Write();
    h2D_Energy_mult4_wa_78Zn->Write();
    h2D_Energy_mult12_wa_78Zn->Write();
    h2D_Energy_mult123_wa_78Zn->Write();
    h2D_Energy_mult1234_wa_78Zn->Write();
    h2D_Energy_mult12345_wa_78Zn->Write();
    h2D_Energy_mult123456_wa_78Zn->Write();
    h2D_Energy_mult1_wa_79Zn->Write();
    h2D_Energy_mult12_wa_79Zn->Write();
    h2D_Energy_mult123_wa_79Zn->Write();
    
    
hEnergy_mult_All_77Zn->Write();
hEnergy_fold12_77Zn->Write();
hEnergy_fold1_77Zn->Write();
hEnergy_mult_All_78Zn->Write();
hEnergy_fold12_78Zn->Write();
hEnergy_fold1_78Zn->Write();
hEnergy_mult_All_79Zn->Write();
hEnergy_fold12_79Zn->Write();
hEnergy_fold1_79Zn->Write();
            
    rootfile->Write();
    rootfile->Close();
     theApp->Run();
    return 0;
}

