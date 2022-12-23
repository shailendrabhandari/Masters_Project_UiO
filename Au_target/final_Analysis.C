//Example how to run Analysis
//./runAnalysis 0.5604 0.5573 0.5518 -13 4 15   beta values from Shailendra
// with ./runAnalsysis beta1 beta2 beta3 timeCutLow timeCutHigh maxaddbackdistance betaDiffLow betaDiffHigh


//./runAnalysis 0.5604 0.5573 0.5518 -13 4 15   beta values from Shailendra{Gold_Target}

//./runAnalysis 0.5338 0.5305 0.5242 -13 4 15   beta values from Shailendra{Carbon_Target}

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
#include "TXMLParser.h"


using namespace std;


int neve = 0;
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

double brho_35, brho_57, beta_35, beta_57, ratio_beta_35_brho_35, ratio_beta_57_brho_57;
double brho_89, brho_911, beta_89, beta_911, ratio_beta_89_brho_89, ratio_beta_911_brho_911;

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
    
    /*TFile *F3cut = new TFile("../cut_Zn/F3_plastic.root","READ");
    TCutG *F3bgcut;
    F3cut->GetObject("CUTG",F3bgcut);
    
    TFile *F7cut = new TFile("../cut_Zn/F7_plastic.root","READ");
    TCutG *F7bgcut;
    F7cut->GetObject("CUTG",F7bgcut);
    
    TFile *F8cut = new TFile("../cut_Zn/F8_plastic.root","READ");
    TCutG *F8bgcut;
    F8cut->GetObject("CUTG",F8bgcut);
    
    TFile *F11cut = new TFile("../cut_Zn/F11_plastic.root","READ");
    TCutG *F11bgcut;
    F11cut->GetObject("CUTG",F11bgcut);
    
    TFile *FPF3cut = new TFile("../cut_Zn/F3_fp.root","READ");
    TCutG *FPF3bgcut;
    FPF3cut->GetObject("CUTG",FPF3bgcut);
    
    TFile *FPF7cut = new TFile("../cut_Zn/F7_fp.root","READ");
    TCutG *FPF7bgcut;
    FPF7cut->GetObject("CUTG",FPF7bgcut);
    
    TFile *FPF8cut = new TFile("../cut_Zn/F8_fp.root","READ");
    TCutG *FPF8bgcut;
    FPF8cut->GetObject("CUTG",FPF8bgcut);
    
    TFile *FPF11cut = new TFile("../cut_Zn/F11_fp.root","READ");
    TCutG *FPF11bgcut;
    FPF11cut->GetObject("CUTG",FPF11bgcut);
    
    
    TFile *CSBRcut = new TFile("../cut_Zn/ChargeState_BR.root","READ");
    TCutG *CSBRbgcut;
    CSBRcut->GetObject("CUTG",CSBRbgcut);
    
    TFile *CSZDcut = new TFile("../cut_Zn/ChargeState_ZD.root","READ");
    TCutG *CSZDbgcut;
    CSZDcut->GetObject("CUTG",CSZDbgcut);
    
    
    TFile *betaBR_betaZD_cut = new TFile("../cut_Zn/betaBRZD_cut.root","READ");
    TCutG *betaBR_betaZD_bgcut;
    betaBR_betaZD_cut->GetObject("CUTG",betaBR_betaZD_bgcut);
    */
    
    
    Event *Tree = new Event();
    
    Long64_t nentries = Tree->fChain->GetEntriesFast();
    cout<<"nentries = "<<nentries <<endl;
    
    TFile *rootfile = new TFile("three_Analysis_Au_20May21_100dets_off.root","RECREATE");
    TTree *tree = new TTree("tree","tree");
    rootfile->cd();
    
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
    
    //Specific variables;
    int daliMult;
    int daliTimeTrueMult;
    int daliFold;
    int daliTimeTrueFold;
    
    Double_t DopplerAngle;
    
    int countingFilling = 0;
    //________________________________________________________________________
    
   
    //Spectra:

    ///------------------------------------------------------------------------------------------------------------------
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







    TH2F *gg_mult2_78Zn = new TH2F("gg_mult2_78Zn","",6000,0,5999,6000,0,5999);
    TH2F *gg_mult3_78Zn = new TH2F("gg_mult3_78Zn","",6000,0,5999,6000,0,5999);
    TH2F *gg_mult23_78Zn = new TH2F("gg_mult23_78Zn","",6000,0,5999,6000,0,5999);
    TH2F *gg_mult234_78Zn = new TH2F("gg_mult234_78Zn","",6000,0,5999,6000,0,5999);
    TH2F *gg_multAll_78Zn = new TH2F("gg_multAll_78Zn","",6000,0,5999,6000,0,5999);
    
    TH2F *gg_mult2_79Zn = new TH2F("gg_mult2_79Zn","",6000,0,5999,6000,0,5999);
    TH2F *gg_mult3_79Zn = new TH2F("gg_mult3_79Zn","",6000,0,5999,6000,0,5999);
    TH2F *gg_mult23_79Zn = new TH2F("gg_mult23_79Zn","",6000,0,5999,6000,0,5999);
    TH2F *gg_mult234_79Zn = new TH2F("gg_mult234_79Zn","",6000,0,5999,6000,0,5999);
    TH2F *gg_multAll_79Zn = new TH2F("gg_multAll_79Zn","",6000,0,5999,6000,0,5999);
    
    
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
        
        
        
       /* // Plastic cuts
        
        bool F3cutBool = false;
        bool F7cutBool = false;
        bool F8cutBool = false;
        bool F11cutBool = false;
        
        
        if(F3bgcut->IsInside(Tree->F3PLA_TR-Tree->F3PLA_TL,log(Tree->F3PLA_QL_raw/Tree->F3PLA_QR_raw))) F3cutBool = true;
        if(F7bgcut->IsInside(Tree->F7PLA_TR-Tree->F7PLA_TL,log(Tree->F7PLA_QL_raw/Tree->F7PLA_QR_raw))) F7cutBool = true;
        if(F8bgcut->IsInside(Tree->F8PLA_TR-Tree->F8PLA_TL,log(Tree->F8PLA_QL_raw/Tree->F8PLA_QR_raw))) F8cutBool = true;
        if(F11bgcut->IsInside(Tree->F11PLA_TR-Tree->F11PLA_TL,log(Tree->F11PLA_QL_raw/Tree->F11PLA_QR_raw))) F11cutBool = true;
        
        
        // Focal plane cuts
        
        bool FPF3cutBool = false;
        bool FPF7cutBool = false;
        bool FPF8cutBool = false;
        bool FPF11cutBool = false;
        
        if(FPF3bgcut->IsInside(Tree->F3PLA_TR-Tree->F3PLA_TL,Tree->F3X)) FPF3cutBool = true;
        if(FPF7bgcut->IsInside(Tree->F7PLA_TR-Tree->F7PLA_TL,Tree->F7X)) FPF7cutBool = true;
        if(FPF8bgcut->IsInside(Tree->F8PLA_TR-Tree->F8PLA_TL,Tree->F8X)) FPF8cutBool = true;
        if(FPF11bgcut->IsInside(Tree->F11PLA_TR-Tree->F11PLA_TL,Tree->F11X)) FPF11cutBool = true;
        
        
        // Charge state cuts
        bool CSBRcutBool = false;
        bool CSZDcutBool = false;
        
        
        if(CSBRbgcut->IsInside(ratio_beta_57_brho_57/ratio_beta_35_brho_35,brho_35)) CSBRcutBool = true;
        if(CSZDbgcut->IsInside(ratio_beta_911_brho_911/ratio_beta_89_brho_89,brho_89)) CSZDcutBool = true;
        
        
        //Beta cuts
        
        bool betaBR_betaZD_cutBool = false;
        if(betaBR_betaZD_bgcut->IsInside(beta_57,beta_911)) betaBR_betaZD_cutBool = true;
        */
        
        //////////////////////////////////////////////////////////////////
       // if( F3cutBool == false || F7cutBool == false || F8cutBool == false || F11cutBool == false) continue;
        
        
        //if( FPF3cutBool == false || FPF7cutBool == false || FPF8cutBool == false || FPF11cutBool == false) continue;
        
        // if(CSBRcutBool == true && CSZDcutBool == true) continue;
        
       // if(CSBRcutBool == false || CSZDcutBool == false) continue;
        
        //if(betaBR_betaZD_cutBool== false) continue;
        
        //if(F3cutBool == true && F7cutBool == true && F8cutBool == true && F11cutBool == true &&FPF3cutBool == true && FPF7cutBool == true && FPF8cutBool == true && FPF11cutBool == true &&CSBRcutBool == true  && CSZDcutBool == true && betaBR_betaZD_cutBool == true)
        //{
        ////////////////////////////////////////////////////////////////////
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
                
                 if (fDali[j].id < 100) fDali[j].e = -999.;
                if (fDali[j].id == 10) fDali[j].e = -999.;
                if (fDali[j].id == 23) fDali[j].e = -999.;
                if (fDali[j].id == 31) fDali[j].e = -999.;
                if (fDali[j].id == 53) fDali[j].e = -999.;
                if (fDali[j].id == 63) fDali[j].e = -999.;
                if (fDali[j].id == 65) fDali[j].e = -999.;
                if (fDali[j].id == 97) fDali[j].e = -999.; //kill simu
                if (fDali[j].id == 98) fDali[j].e = -999.;//kill simu
                if (fDali[j].id == 127) fDali[j].e = -999.;//kill simu
                if (fDali[j].id == 128) fDali[j].e = -999.;//kill simu
                if (fDali[j].id == 141) fDali[j].e = -999.;//kill simu
                if (fDali[j].id == 167) fDali[j].e = -999.; //kill simu
                
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
        
       
    //////////////////////////////////////////////////////////////////
        /*//Beta diff
         float beta_diff = Tree->BigRIPSBeam_beta[0]-Tree->BigRIPSBeam_beta[3];
         if(brzn77&&zdzn77) h_beta_diff[0]->Fill(beta_diff);
         if(brzn78&&zdzn78) h_beta_diff[1]->Fill(beta_diff);
         if(brzn79&&zdzn79) h_beta_diff[2]->Fill(beta_diff);*/
        
        //Trigger register information
        bool DSB = false;
        bool F11 = false;
        bool DaliTrigger = false;
        
        //Need to check the settings
        if(Tree->fbit==1||Tree->fbit==3||Tree->fbit==7) DSB = true;
        if(Tree->fbit==6||Tree->fbit==7) DaliTrigger = true;
        if(Tree->fbit==3) F11 = true;
        
        ///if( F3cutBool==true  && F7cutBool == true && F8cutBool ==true && F11cutBool==true ){
       
        //Getting the statistics for the cross-section:
        //Beam
        if(DSB && brzn77)
            hBR_aoq_z_DS1_77Zn->Fill(Tree->aoqc[2],Tree->zetc[2]);
        if(DSB && brzn78)
            hBR_aoq_z_DS1_78Zn->Fill(Tree->aoqc[2],Tree->zetc[2]);
        if(DSB && brzn79)
            hBR_aoq_z_DS1_79Zn->Fill(Tree->aoqc[2],Tree->zetc[2]);
      ///  }






////////////////////////Gamma-ray sepctra /////////////////////
      
        /////With Addback//////
        
        //////77Zn//////////////
        for(int j=0;j<fDaliMultTa;j++){
            if((fDali[j].t)>-13 && (fDali[j].t)<4 &&DaliTrigger) {///&&DaliTrigger
                if((fDaliMultTa == 1 )&& brzn77 && zdzn77 && fDali[j].idwa >= 100 )
                    hEnergy_mult1_wa_77Zn->Fill(fDali[j].doppwa1);
                
                if((fDaliMultTa<=2) && brzn77 && zdzn77 && fDali[j].idwa >= 100)
                    hEnergy_mult12_wa_77Zn->Fill(fDali[j].doppwa1);
                
                if((fDaliMultTa<=3) && brzn77 && zdzn77 && fDali[j].idwa >= 100)
                    hEnergy_mult123_wa_77Zn->Fill(fDali[j].doppwa1);
            }
      //////78Zn//////////////
           if((fDali[j].t)>-13 && (fDali[j].t)<4 && DaliTrigger) {
        if((fDaliMultTa == 1 )&& brzn78 && zdzn78 && fDali[j].idwa >= 100 )
            hEnergy_mult1_wa_78Zn->Fill(fDali[j].doppwa2);
        
        if((fDaliMultTa<=2)&& brzn78 && zdzn78 && fDali[j].idwa >= 100)
            hEnergy_mult12_wa_78Zn->Fill(fDali[j].doppwa2);
        
        if((fDaliMultTa<=3) && brzn78 && zdzn78 && fDali[j].idwa >= 100)
            hEnergy_mult123_wa_78Zn->Fill(fDali[j].doppwa2);

        if((fDaliMultTa<=4) && brzn78 && zdzn78 && fDali[j].idwa >= 100)
            hEnergy_mult1234_wa_78Zn->Fill(fDali[j].doppwa2);

        if((fDaliMultTa<=5) && brzn78 && zdzn78 && fDali[j].idwa >= 100)
            hEnergy_mult12345_wa_78Zn->Fill(fDali[j].doppwa2);

        if((fDaliMultTa<=6) && brzn78 && zdzn78 && fDali[j].idwa >= 100)
            hEnergy_mult123456_wa_78Zn->Fill(fDali[j].doppwa2);




        if((fDaliMultTa == 2 )&& brzn78 && zdzn78 && fDali[j].idwa >= 100 )
            hEnergy_mult2_wa_78Zn->Fill(fDali[j].doppwa2);

        if((fDaliMultTa == 3 )&& brzn78 && zdzn78 && fDali[j].idwa >= 100 )
            hEnergy_mult3_wa_78Zn->Fill(fDali[j].doppwa2);

        if((fDaliMultTa == 4 )&& brzn78 && zdzn78 && fDali[j].idwa >= 100 )
            hEnergy_mult4_wa_78Zn->Fill(fDali[j].doppwa2);
            }
        //////79Zn//////////////
            if((fDali[j].t)>-13 && (fDali[j].t)<4 &&DaliTrigger) {
                if((fDaliMultTa == 1 )&& brzn79 && zdzn79 && fDali[j].idwa >= 100 )
                    hEnergy_mult1_wa_79Zn->Fill(fDali[j].doppwa3);
                
                if((fDaliMultTa<=2) && brzn79 && zdzn79 && fDali[j].idwa >= 100)
                    hEnergy_mult12_wa_79Zn->Fill(fDali[j].doppwa3);
                
                if((fDaliMultTa<=3) && brzn79 && zdzn79 && fDali[j].idwa >= 100)
                    hEnergy_mult123_wa_79Zn->Fill(fDali[j].doppwa3);
            }
        }
    
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       ////Without Addback
            //////77Zn//////////////
        for(int j=0;j<fDaliFold;j++){
            if(fDali[j].ttrue) {
      
            if((fDali[j].t)>-13 && (fDali[j].t)<4 &&DaliTrigger) {
                        if((fDaliMultTa == 1) && brzn77 && zdzn77 && fDali[j].id >= 100)
                            hEnergy_mult1_77Zn->Fill(fDali[j].dopp1);
                        if((fDaliMultTa<=2)&& brzn77 && zdzn77 && fDali[j].id >= 100)
                            hEnergy_mult12_77Zn->Fill(fDali[j].dopp1);
                        if((fDaliMultTa<=3) && brzn77 && zdzn77 && fDali[j].id >= 100)
                            hEnergy_mult123_77Zn->Fill(fDali[j].dopp1);
                    }
        //////78Zn//////////////
       if((fDali[j].t)>-13 && (fDali[j].t)<4 && DaliTrigger) {
                    if((fDaliMultTa == 1) && brzn78 && zdzn78 && fDali[j].id >= 100)
                        hEnergy_mult1_78Zn->Fill(fDali[j].dopp2);
            
                    if((fDaliMultTa<=2) && brzn78 && zdzn78 && fDali[j].id >= 100)
                        hEnergy_mult12_78Zn->Fill(fDali[j].dopp2);
                    
                    if((fDaliMultTa<=3) && brzn78 && zdzn78 && fDali[j].id >= 100)
                        hEnergy_mult123_78Zn->Fill(fDali[j].dopp2);
                    
                    if((fDaliMultTa<=4) && brzn78 && zdzn78 && fDali[j].id >= 100)
                        hEnergy_mult1234_78Zn->Fill(fDali[j].dopp2);
                    
                    if((fDaliMultTa<=5) && brzn78 && zdzn78 && fDali[j].id >= 100)
                        hEnergy_mult12345_78Zn->Fill(fDali[j].dopp2);
                    
                    if((fDaliMultTa<=6) && brzn78 && zdzn78 && fDali[j].id >= 100)
                        hEnergy_mult123456_78Zn->Fill(fDali[j].dopp2);


                    if((fDaliMultTa == 2) && brzn78 && zdzn78 && fDali[j].id >= 100)
                        hEnergy_mult2_78Zn->Fill(fDali[j].dopp2);

                    if((fDaliMultTa == 3) && brzn78 && zdzn78 && fDali[j].id >= 100)
                        hEnergy_mult3_78Zn->Fill(fDali[j].dopp2);

                    if((fDaliMultTa == 4) && brzn78 && zdzn78 && fDali[j].id >= 100)
                        hEnergy_mult4_78Zn->Fill(fDali[j].dopp2);
                }
        //////79Zn//////////////
    if((fDali[j].t)>-13 && (fDali[j].t)<4 &&DaliTrigger) {
                    if((fDaliMultTa == 1) && brzn79 && zdzn79 && fDali[j].id >= 100)
                        hEnergy_mult1_79Zn->Fill(fDali[j].dopp3);
        
                    if((fDaliMultTa<=2) && brzn79 && zdzn79 && fDali[j].id >= 100)
                        hEnergy_mult12_79Zn->Fill(fDali[j].dopp3);
                    
                    if((fDaliMultTa<=3) && brzn79 && zdzn79 && fDali[j].id >= 100)
                        hEnergy_mult123_79Zn->Fill(fDali[j].dopp3);
                }
            }
        }
        
        
        
        /////2D histograms/////
        
        /////With Addback//////
        
        for(int j=0;j<fDaliMultTa;j++){
            //////77Zn//////////////
            if((fDali[j].t)>=-13 && (fDali[j].t)<=4 &&DaliTrigger) {
                if((fDaliMultTa == 1 )&& brzn77 && zdzn77)
                    h2D_Energy_mult1_wa_77Zn->Fill(fDali[j].idwa,fDali[j].doppwa2);
                
                if((fDaliMultTa<=2)&& brzn77 && zdzn77)
                    h2D_Energy_mult12_wa_77Zn->Fill(fDali[j].idwa,fDali[j].doppwa2);
                
                if((fDaliMultTa<=3) && brzn77 && zdzn77)
                    h2D_Energy_mult123_wa_77Zn->Fill(fDali[j].idwa,fDali[j].doppwa2);
            }
            //////78Zn//////////////
            if((fDali[j].t)>=-13 && (fDali[j].t)<=4 && DaliTrigger) {
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
            if((fDali[j].t)>=-13 && (fDali[j].t)<=4 &&DaliTrigger) {
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
                if((fDali[j].t)>=-10 && (fDali[j].t)<=4 &&DaliTrigger) {
                    if((fDaliMultTa == 1) && brzn77 && zdzn77)
                        h2D_Energy_mult1_77Zn->Fill(fDali[j].id,fDali[j].dopp2);
                    
                    if((fDaliMultTa<=2) && brzn77 && zdzn77)
                        h2D_Energy_mult12_77Zn->Fill(fDali[j].id,fDali[j].dopp2);
                    
                    if((fDaliMultTa<=3) && brzn77 && zdzn77)
                        h2D_Energy_mult123_77Zn->Fill(fDali[j].id,fDali[j].dopp2);
                }
                //////78Zn//////////////
                 if((fDali[j].t)>=-13 && (fDali[j].t)<=4 && DaliTrigger) {
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
                if((fDali[j].t)>=-9 && (fDali[j].t)<=5 &&DaliTrigger) {
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
                if((fDali[j].t)>-13 && (fDali[j].t)<4 && fDali[j].id >= 100 &&DaliTrigger) {
         hEnergy_mult_All_77Zn->Fill(fDali[j].dopp1);
         if(fDaliFold<3)  hEnergy_fold12_77Zn->Fill(fDali[j].dopp1);
         if(fDaliFold==1) hEnergy_fold1_77Zn->Fill(fDali[j].dopp1);
                }
       //////78Zn//////////////
            if((fDali[j].t)>-13 && (fDali[j].t)<4 && fDali[j].id >= 100 &&DaliTrigger) {
                hEnergy_mult_All_78Zn->Fill(fDali[j].dopp2);
                if(fDaliFold<3)  hEnergy_fold12_78Zn->Fill(fDali[j].dopp2);
                if(fDaliFold==1) hEnergy_fold1_78Zn->Fill(fDali[j].dopp2);
                }
        //////79Zn//////////////
        if((fDali[j].t)>-13 && (fDali[j].t)<4 && fDali[j].id >= 100 &&DaliTrigger) {
            hEnergy_mult_All_79Zn->Fill(fDali[j].dopp3);
            if(fDaliFold<3)  hEnergy_fold12_79Zn->Fill(fDali[j].dopp3);
            if(fDaliFold==1) hEnergy_fold1_79Zn->Fill(fDali[j].dopp3);
                }
            }
        }
        //ZDS
        hBR_aoq_z_full->Fill(Tree->aoqc[2],Tree->zetc[2]);
         hZD_aoq_z_full->Fill(Tree->aoqc[5],Tree->zetc[5]);
        if(F11 & zdzn77)
            hZD_aoq_z_F11_77Zn->Fill(Tree->aoqc[5],Tree->zetc[5]);
        if(F11 & zdzn78)
            hZD_aoq_z_F11_78Zn->Fill(Tree->aoqc[5],Tree->zetc[5]);
        if(F11 & zdzn79)
            hZD_aoq_z_F11_79Zn->Fill(Tree->aoqc[5],Tree->zetc[5]);
     
        
        
        ///////Creating Coincidence data & gg-matrices/////////////////
        //////78Zn/////
        for(int j=0;j<fDaliMultTa;j++){
            if((fDali[j].t)>-13 && (fDali[j].t)<4 && fDali[j].id >= 100 && DaliTrigger && brzn78 && zdzn78){
                for(int v = 0; v < 186; v++)
                {
                    if(fDali[v].doppwa2 > 0.)
                    {
                        for(int w = v+1; w < 186; w++)
                        {
                            if(fDali[w].doppwa2 > 0.)
                                
                                if(fDaliMultTa ==2 ){   //pp' p2p and p3p reactions,
                                    
                                    gg_mult2_78Zn->Fill(fDali[v].doppwa2,fDali[w].doppwa2);
                                    gg_mult2_78Zn->Fill(fDali[w].doppwa2,fDali[v].doppwa2);
                                    
                                    
                                }
                            if(fDaliMultTa==3 ){
                                
                                gg_mult3_78Zn->Fill(fDali[v].doppwa2,fDali[w].doppwa2);
                                gg_mult3_78Zn->Fill(fDali[w].doppwa2,fDali[v].doppwa2);
                                
                                
                            }
                            if((fDaliMultTa ==2) || (fDaliMultTa==3 )){
                                
                                gg_mult23_78Zn->Fill(fDali[v].doppwa2,fDali[w].doppwa2);
                                gg_mult23_78Zn->Fill(fDali[w].doppwa2,fDali[v].doppwa2);
                                
                                
                            }
                            if((fDaliMultTa ==2) || (fDaliMultTa==3) || (fDaliMultTa==4) ){
                                
                                gg_mult234_78Zn->Fill(fDali[v].doppwa2,fDali[w].doppwa2);
                                gg_mult234_78Zn->Fill(fDali[w].doppwa2,fDali[v].doppwa2);
                                
                                
                            }
                            if(fDaliMultTa >=2 ){
                                
                                
                                gg_multAll_78Zn->Fill(fDali[v].doppwa2,fDali[w].doppwa2);
                                gg_multAll_78Zn->Fill(fDali[w].doppwa2,fDali[v].doppwa2);
                                
                                
                                
                            }
                        }
                    }
                }
                
            }
        }
        
        
        
      
        //////79Zn/////
        for(int j=0;j<fDaliMultTa;j++){
            if((fDali[j].t)>-13 && (fDali[j].t)<4 && fDali[j].id >= 100 && DaliTrigger && brzn79 && zdzn79){
                for(int v = 0; v < 186; v++)
                {
                    if(fDali[v].doppwa2 > 0.)
                    {
                        for(int w = v+1; w < 186; w++)
                        {
                            if(fDali[w].doppwa2 > 0.)
                                
                                if(fDaliMultTa ==2 ){   //pp' p2p and p3p reactions,
                                    
                                    gg_mult2_79Zn->Fill(fDali[v].doppwa3,fDali[w].doppwa3);
                                    gg_mult2_79Zn->Fill(fDali[w].doppwa3,fDali[v].doppwa3);
                                    
                                    
                                }
                            if(fDaliMultTa==3 ){
                                
                                gg_mult3_79Zn->Fill(fDali[v].doppwa3,fDali[w].doppwa3);
                                gg_mult3_79Zn->Fill(fDali[w].doppwa3,fDali[v].doppwa3);
                                
                                
                            }
                            if((fDaliMultTa ==2) || (fDaliMultTa==3 )){
                                
                                gg_mult23_79Zn->Fill(fDali[v].doppwa3,fDali[w].doppwa3);
                                gg_mult23_79Zn->Fill(fDali[w].doppwa3,fDali[v].doppwa3);
                                
                                
                            }
                            if((fDaliMultTa ==2) || (fDaliMultTa==3) || (fDaliMultTa==4) ){
                                
                                gg_mult234_79Zn->Fill(fDali[v].doppwa3,fDali[w].doppwa3);
                                gg_mult234_79Zn->Fill(fDali[w].doppwa3,fDali[v].doppwa3);
                                
                                
                            }
                            if(fDaliMultTa >=2 ){
                                
                                
                                gg_multAll_79Zn->Fill(fDali[v].doppwa3,fDali[w].doppwa3);
                                gg_multAll_79Zn->Fill(fDali[w].doppwa3,fDali[v].doppwa3);
                                
                                
                                
                            }
                        }
                    }
                }
                
            }
        }
        
        
        
    }
        neve++;
        tree->Fill();
    //}

/*
hEnergy_mult1_wa_77Zn->Write();
hEnergy_mult12_wa_77Zn->Write();
hEnergy_mult123_wa_77Zn->Write();
hEnergy_mult1_77Zn->Write();
hEnergy_mult12_77Zn->Write();
hEnergy_mult123_77Zn->Write();
    
hEnergy_mult1_wa_78Zn->Write();
hEnergy_mult12_wa_78Zn->Write();
hEnergy_mult123_wa_78Zn->Write();
hEnergy_mult1_78Zn->Write();
hEnergy_mult12_78Zn->Write();
hEnergy_mult123_78Zn->Write();
    
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
    h2D_Energy_mult12_78Zn->Write();
    h2D_Energy_mult123_78Zn->Write();
    h2D_Energy_mult1_79Zn->Write();
    h2D_Energy_mult12_79Zn->Write();
    h2D_Energy_mult123_79Zn->Write();
    
    h2D_Energy_mult1_wa_77Zn->Write();
    h2D_Energy_mult12_wa_77Zn->Write();
    h2D_Energy_mult123_wa_77Zn->Write();
    h2D_Energy_mult1_wa_78Zn->Write();
    h2D_Energy_mult12_wa_78Zn->Write();
    h2D_Energy_mult123_wa_78Zn->Write();
    h2D_Energy_mult1_wa_79Zn->Write();
    h2D_Energy_mult12_wa_79Zn->Write();
    h2D_Energy_mult123_wa_79Zn->Write();
    */


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
    
    
    gg_mult2_78Zn->Write();
    gg_mult3_78Zn->Write();
    gg_mult23_78Zn->Write();
    gg_mult234_78Zn->Write();
    gg_multAll_78Zn->Write();
    
    gg_mult2_79Zn->Write();
    gg_mult3_79Zn->Write();
    gg_mult23_79Zn->Write();
    gg_mult234_79Zn->Write();
    gg_multAll_79Zn->Write();
    rootfile->Write();
    rootfile->Close();
    // theApp->Run();
    return 0;
}

