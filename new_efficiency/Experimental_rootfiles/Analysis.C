//Example how to run Analysis for Calibration runs
//./runAnalysis 0 0 1100 18 // TOfsetted is different for calibration runs and range between 900 and 1100.

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
    float dopp; //Doppler energy. Three doppler corrections for the three betas
    float doppwa; //Doppler energy with true multiplicity and addback.
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
    
    FILE *fdetPos  = fopen("AverageInteractionPoint.out","r");
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
    FILE *fAddbackTableOut = fopen("AddbackTable.out","w");
    
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
    
    if(argc!=5) {
        std::cout << "Invalid arguments..." << endl;
        return 1;
    }
    
    double beta1 = atof(argv[1]);
    
    
    
    double timeCutLow = atof(argv[2]);
    double timeCutHigh = atof(argv[3]);
    
    //Create the add back table
    CreateAddBackTable(atof(argv[4]));
    
   
    
    // Create interactive interface
    TRint *theApp = new TRint("ROOT example", &argc, argv, NULL, 0);
    //__________________________________________________________
   
     char name[100];
   
    
    Event *Tree = new Event();
    
    Long64_t nentries = Tree->fChain->GetEntriesFast();
    
    TFile *rootfile = new TFile("co60_calibrated_histograms.root","RECREATE");
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
    
    /*//Define spectra:
    int minBin = 0;
    int maxBin = 6000;
    int binning = 25;
    int numBin = (maxBin-minBin)/binning;*/
    //________________________________________________________________________
    
    //Spectra:
    
    TH1F *h_beta_diff;
    
    TH2F *h_doppler[10];
    

     h_beta_diff = new TH1F(name,name,1000,0.6,0.12);
    
    for(int i=0;i<10;i++) {
        sprintf(name,"h_doppler[%i]",i);
        //h_doppler[i] = new TH2F(name,name,186,0,186,numBin,minBin,maxBin);
         h_doppler[i] = new TH2F(name,name,186,0,186,1000,0,4000);
    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //Start looping through data;
    
    Long64_t i=0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        
        if(jentry%10000 ==0) cout << jentry <<"/"<<nentries<<" Events DONE!"<<endl;
        Tree->fChain->GetEvent(jentry);
        
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
                
                fDali[j].dopp = DopplerCorrect(beta1,fDali[j].theta,fDali[j].e);
                
                if(fDali[j].t>timeCutLow-500&&fDali[j].t<timeCutHigh+500)fDaliFold++;
                
                if(fDali[j].t>timeCutLow && fDali[j].t<timeCutHigh){
                    fDali[j].ttrue = true;
                    fDaliFoldTa++;
                }
                else fDali[j].ttrue = false;
            }
            else {
                fDali[j].dopp = -999.;
                
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
            fDali[j].dopp    = -999;
            
        }
        
        
        //if(fDali[0].e>0) /////Pieter asked me to closedthis line March 2019
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
            
            crystalUsedForAddback[fDali[i].id]=true;
            fDali[fDaliMultTa].idwa = fDali[i].id;
            
            for(int j = i+1;j<fDaliFold;j++)  {
                if(crystalUsedForAddback[fDali[j].id]==false && fDali[j].ttrue==true)  {
                    for(int k = 0;k<fNumberOfAddbackPartners[fDali[i].id] ;k++) {
                        if(fDali[j].id == fAddbackTable[fDali[i].id][k+1])  {
                            
                            crystalUsedForAddback[fDali[j].id]=true;
                            
                            dummyEnergy[fDaliMultTa][0] += DopplerCorrect(beta1,fDali[i].theta,fDali[j].e);
                            
                            
                        }
                    }
                }
            }
            fDaliMultTa++;
        }
        for(int i = 0;i<fDaliMultTa;i++) {
            fDali[i].doppwa = dummyEnergy[i][0];
           
        }
        for(int i = fDaliMultTa;i<NUMBEROFDALICRYSTALS;i++) {
            fDali[i].doppwa = -999;
            fDali[i].idwa      = -999;
        }
        
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
       

        
        for(int j=0;j<fDaliFold;j++){
            if(fDali[j].ttrue) {
                h_doppler[0]->Fill(fDali[j].id,fDali[j].dopp);
                if(fDaliFold<=3)  h_doppler[1]->Fill(fDali[j].id,fDali[j].dopp);
                if(fDaliFold==1) h_doppler[2]->Fill(fDali[j].id,fDali[j].dopp);
            }
        }
        for(int j=0;j<fDaliMultTa;j++){
            h_doppler[3]->Fill(fDali[j].idwa,fDali[j].doppwa);
            if(fDaliMultTa==1) h_doppler[4]->Fill(fDali[0].idwa,fDali[0].doppwa);
            if(fDaliMultTa==2) h_doppler[5]->Fill(fDali[j].idwa,fDali[j].doppwa);
            if(fDaliMultTa==3) h_doppler[6]->Fill(fDali[j].idwa,fDali[j].doppwa);
            if(fDaliMultTa>3)  h_doppler[7]->Fill(fDali[j].idwa,fDali[j].doppwa);
            if(fDaliMultTa<=2) {h_doppler[8]->Fill(fDali[j].idwa,fDali[j].doppwa);
            }
            if(fDaliMultTa<=3) {h_doppler[9]->Fill(fDali[j].idwa,fDali[j].doppwa);
            }
        }
        
        
    }
    
    for(int i=0;i<10;i++)
        h_doppler[i]->Write();
    
    
     h_beta_diff->Write();
    
    rootfile->Write();
    rootfile->Close();
     theApp->Run();
    return 0;
}

