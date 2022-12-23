#ifndef Event_h
#define Event_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TMatrixT.h>
#include <TVectorT.h>

// Header file for the classes stored in the TTree if any.
#include <TClonesArray.h>
#include <TObject.h>
#include <TNamed.h>

// Fixed size dimensions of array or collections stored in the TTree if any.
const Int_t kMaxEventInfo = 1;
const Int_t kMaxBigRIPSPPAC = 30;
const Int_t kMaxBigRIPSPlastic = 5;
const Int_t kMaxBigRIPSIC = 4;
const Int_t kMaxBigRIPSFocalPlane = 8;
const Int_t kMaxDALINaI = 147;
const Int_t kMaxBigRIPSRIPS = 4;
const Int_t kMaxBigRIPSTOF = 4;
const Int_t kMaxBigRIPSBeam = 4;

class Event {
 public :
  TChain          *fChain;   //!pointer to the analyzed TTree or TChain

  Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           fbit;
   Int_t           EventInfo_;
   UInt_t          EventInfo_fUniqueID[kMaxEventInfo];   //[EventInfo_]
   UInt_t          EventInfo_fBits[kMaxEventInfo];   //[EventInfo_]
   TString         EventInfo_fName[kMaxEventInfo];
   TString         EventInfo_fTitle[kMaxEventInfo];
   Int_t           EventInfo_runnumber[kMaxEventInfo];   //[EventInfo_]
   Int_t           EventInfo_eventnumber[kMaxEventInfo];   //[EventInfo_]
   TString         EventInfo_subsysname[kMaxEventInfo];
   ULong64_t       EventInfo_timestamp[kMaxEventInfo];   //[EventInfo_]
   Int_t           EventInfo_comp_val[kMaxEventInfo];   //[EventInfo_]
   UInt_t          EventInfo_fBit[kMaxEventInfo];   //[EventInfo_]
   Int_t           BigRIPSPPAC_;
   UInt_t          BigRIPSPPAC_fUniqueID[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   UInt_t          BigRIPSPPAC_fBits[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Int_t           BigRIPSPPAC_id[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Int_t           BigRIPSPPAC_fpl[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   TString         BigRIPSPPAC_name[kMaxBigRIPSPPAC];
   Int_t           BigRIPSPPAC_fDataState[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Double_t        BigRIPSPPAC_xzpos[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Double_t        BigRIPSPPAC_yzpos[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Int_t           BigRIPSPPAC_fTX1Raw[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Int_t           BigRIPSPPAC_fTX2Raw[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Int_t           BigRIPSPPAC_fTY1Raw[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Int_t           BigRIPSPPAC_fTY2Raw[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Int_t           BigRIPSPPAC_fTARaw[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Int_t           BigRIPSPPAC_fQX1Raw[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Int_t           BigRIPSPPAC_fQX2Raw[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Int_t           BigRIPSPPAC_fQY1Raw[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Int_t           BigRIPSPPAC_fQY2Raw[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Int_t           BigRIPSPPAC_fQARaw[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Double_t        BigRIPSPPAC_fTSumX[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Double_t        BigRIPSPPAC_fTSumY[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Double_t        BigRIPSPPAC_fTDiffX[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Double_t        BigRIPSPPAC_fTDiffY[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Double_t        BigRIPSPPAC_fX[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Double_t        BigRIPSPPAC_fY[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Bool_t          BigRIPSPPAC_fFiredX[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Bool_t          BigRIPSPPAC_fFiredY[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
   Int_t           BigRIPSPlastic_;
   UInt_t          BigRIPSPlastic_fUniqueID[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
   UInt_t          BigRIPSPlastic_fBits[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
   Int_t           BigRIPSPlastic_id[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
   Int_t           BigRIPSPlastic_fpl[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
   TString         BigRIPSPlastic_name[kMaxBigRIPSPlastic];
   Int_t           BigRIPSPlastic_fDataState[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
   Double_t        BigRIPSPlastic_zposition[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
   Double_t        BigRIPSPlastic_zoffset[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
   Int_t           BigRIPSPlastic_fTLRaw[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
   Int_t           BigRIPSPlastic_fTRRaw[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
   Int_t           BigRIPSPlastic_fQLRaw[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
   Int_t           BigRIPSPlastic_fQRRaw[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
   Double_t        BigRIPSPlastic_fTime[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
   Double_t        BigRIPSPlastic_fTimeL[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
   Double_t        BigRIPSPlastic_fTimeR[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
   Double_t        BigRIPSPlastic_fTimeLSlew[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
   Double_t        BigRIPSPlastic_fTimeRSlew[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
   Double_t        BigRIPSPlastic_fTimeSlew[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
   Int_t           BigRIPSIC_;
   UInt_t          BigRIPSIC_fUniqueID[kMaxBigRIPSIC];   //[BigRIPSIC_]
   UInt_t          BigRIPSIC_fBits[kMaxBigRIPSIC];   //[BigRIPSIC_]
   Int_t           BigRIPSIC_id[kMaxBigRIPSIC];   //[BigRIPSIC_]
   Int_t           BigRIPSIC_fpl[kMaxBigRIPSIC];   //[BigRIPSIC_]
   TString         BigRIPSIC_name[kMaxBigRIPSIC];
   Int_t           BigRIPSIC_fDataState[kMaxBigRIPSIC];   //[BigRIPSIC_]
   Double_t        BigRIPSIC_zcoef[kMaxBigRIPSIC][2];   //[BigRIPSIC_]
   Double_t        BigRIPSIC_ionpair[kMaxBigRIPSIC];   //[BigRIPSIC_]
   Int_t           BigRIPSIC_nhitchannel[kMaxBigRIPSIC];   //[BigRIPSIC_]
   Int_t           BigRIPSIC_fADC[kMaxBigRIPSIC][6];   //[BigRIPSIC_]
   Double_t        BigRIPSIC_fEnergy[kMaxBigRIPSIC][6];   //[BigRIPSIC_]
   Double_t        BigRIPSIC_fRawADCSqSum[kMaxBigRIPSIC];   //[BigRIPSIC_]
   Double_t        BigRIPSIC_fRawADCAvSum[kMaxBigRIPSIC];   //[BigRIPSIC_]
   Double_t        BigRIPSIC_fCalMeVSqSum[kMaxBigRIPSIC];   //[BigRIPSIC_]
   Double_t        BigRIPSIC_fCalMeVAvSum[kMaxBigRIPSIC];   //[BigRIPSIC_]
   Int_t           BigRIPSFocalPlane_;
   UInt_t          BigRIPSFocalPlane_fUniqueID[kMaxBigRIPSFocalPlane];   //[BigRIPSFocalPlane_]
   UInt_t          BigRIPSFocalPlane_fBits[kMaxBigRIPSFocalPlane];   //[BigRIPSFocalPlane_]
   Int_t           BigRIPSFocalPlane_id[kMaxBigRIPSFocalPlane];   //[BigRIPSFocalPlane_]
   Int_t           BigRIPSFocalPlane_fpl[kMaxBigRIPSFocalPlane];   //[BigRIPSFocalPlane_]
   TString         BigRIPSFocalPlane_name[kMaxBigRIPSFocalPlane];
   Int_t           BigRIPSFocalPlane_fDataState[kMaxBigRIPSFocalPlane];   //[BigRIPSFocalPlane_]
   TVectorT<double> BigRIPSFocalPlane_opt_vector[kMaxBigRIPSFocalPlane];
   Int_t           BigRIPSFocalPlane_nfired_ppacx[kMaxBigRIPSFocalPlane];   //[BigRIPSFocalPlane_]
   Int_t           BigRIPSFocalPlane_nfired_ppacy[kMaxBigRIPSFocalPlane];   //[BigRIPSFocalPlane_]
   Double_t        BigRIPSFocalPlane_zpos[kMaxBigRIPSFocalPlane];   //[BigRIPSFocalPlane_]
   Int_t           DALINaI_;
   UInt_t          DALINaI_fUniqueID[kMaxDALINaI];   //[DALINaI_]
   UInt_t          DALINaI_fBits[kMaxDALINaI];   //[DALINaI_]
   Int_t           DALINaI_id[kMaxDALINaI];   //[DALINaI_]
   Int_t           DALINaI_fpl[kMaxDALINaI];   //[DALINaI_]
   TString         DALINaI_name[kMaxDALINaI];
   Int_t           DALINaI_fDataState[kMaxDALINaI];   //[DALINaI_]
   Int_t           DALINaI_fADC[kMaxDALINaI];   //[DALINaI_]
   Int_t           DALINaI_fTDC[kMaxDALINaI];   //[DALINaI_]
   Int_t           DALINaI_layer[kMaxDALINaI];   //[DALINaI_]
   Double_t        DALINaI_theta[kMaxDALINaI];   //[DALINaI_]
   Double_t        DALINaI_fXPos[kMaxDALINaI];   //[DALINaI_]
   Double_t        DALINaI_fYPos[kMaxDALINaI];   //[DALINaI_]
   Double_t        DALINaI_fZPos[kMaxDALINaI];   //[DALINaI_]
   Double_t        DALINaI_costheta[kMaxDALINaI];   //[DALINaI_]
   Double_t        DALINaI_fEnergy[kMaxDALINaI];   //[DALINaI_]
   Double_t        DALINaI_fDoppCorEnergy[kMaxDALINaI];   //[DALINaI_]
   Double_t        DALINaI_fEnergyWithoutT[kMaxDALINaI];   //[DALINaI_]
   Double_t        DALINaI_fTime[kMaxDALINaI];   //[DALINaI_]
   Double_t        DALINaI_fTimeOffseted[kMaxDALINaI];   //[DALINaI_]
   Double_t        DALINaI_fTimeTrueEnergy[kMaxDALINaI];   //[DALINaI_]
   Double_t        DALINaI_fTimeTrueDoppCorEnergy[kMaxDALINaI];   //[DALINaI_]
   Double_t        DALINaI_fTimeTrueTime[kMaxDALINaI];   //[DALINaI_]
   Double_t        DALINaI_fTimeTrueTimeOffseted[kMaxDALINaI];   //[DALINaI_]
   Int_t           BigRIPSRIPS_;
   UInt_t          BigRIPSRIPS_fUniqueID[kMaxBigRIPSRIPS];   //[BigRIPSRIPS_]
   UInt_t          BigRIPSRIPS_fBits[kMaxBigRIPSRIPS];   //[BigRIPSRIPS_]
   Int_t           BigRIPSRIPS_id[kMaxBigRIPSRIPS];   //[BigRIPSRIPS_]
   Int_t           BigRIPSRIPS_fpl[kMaxBigRIPSRIPS];   //[BigRIPSRIPS_]
   TString         BigRIPSRIPS_name[kMaxBigRIPSRIPS];
   Int_t           BigRIPSRIPS_fDataState[kMaxBigRIPSRIPS];   //[BigRIPSRIPS_]
   Int_t           BigRIPSRIPS_upstream_fpl[kMaxBigRIPSRIPS];   //[BigRIPSRIPS_]
   Int_t           BigRIPSRIPS_downstream_fpl[kMaxBigRIPSRIPS];   //[BigRIPSRIPS_]
   Double_t        BigRIPSRIPS_center_brho[kMaxBigRIPSRIPS];   //[BigRIPSRIPS_]
   Double_t        BigRIPSRIPS_brho[kMaxBigRIPSRIPS];   //[BigRIPSRIPS_]
   Double_t        BigRIPSRIPS_length[kMaxBigRIPSRIPS];   //[BigRIPSRIPS_]
   TMatrixT<double> BigRIPSRIPS_matrix[kMaxBigRIPSRIPS];
   Double_t        BigRIPSRIPS_delta[kMaxBigRIPSRIPS];   //[BigRIPSRIPS_]
   Double_t        BigRIPSRIPS_angle[kMaxBigRIPSRIPS];   //[BigRIPSRIPS_]
   Int_t           BigRIPSTOF_;
   UInt_t          BigRIPSTOF_fUniqueID[kMaxBigRIPSTOF];   //[BigRIPSTOF_]
   UInt_t          BigRIPSTOF_fBits[kMaxBigRIPSTOF];   //[BigRIPSTOF_]
   Int_t           BigRIPSTOF_id[kMaxBigRIPSTOF];   //[BigRIPSTOF_]
   Int_t           BigRIPSTOF_fpl[kMaxBigRIPSTOF];   //[BigRIPSTOF_]
   TString         BigRIPSTOF_name[kMaxBigRIPSTOF];
   Int_t           BigRIPSTOF_fDataState[kMaxBigRIPSTOF];   //[BigRIPSTOF_]
   Double_t        BigRIPSTOF_tof[kMaxBigRIPSTOF];   //[BigRIPSTOF_]
   Double_t        BigRIPSTOF_clight[kMaxBigRIPSTOF];   //[BigRIPSTOF_]
   Double_t        BigRIPSTOF_length[kMaxBigRIPSTOF];   //[BigRIPSTOF_]
   TString         BigRIPSTOF_upstream_plname[kMaxBigRIPSTOF];
   TString         BigRIPSTOF_downstream_plname[kMaxBigRIPSTOF];
   Int_t           BigRIPSTOF_upstream_plfpl[kMaxBigRIPSTOF];   //[BigRIPSTOF_]
   Int_t           BigRIPSTOF_downstream_plfpl[kMaxBigRIPSTOF];   //[BigRIPSTOF_]
   Double_t        BigRIPSTOF_time_offset[kMaxBigRIPSTOF];   //[BigRIPSTOF_]
   Int_t           BigRIPSBeam_;
   UInt_t          BigRIPSBeam_fUniqueID[kMaxBigRIPSBeam];   //[BigRIPSBeam_]
   UInt_t          BigRIPSBeam_fBits[kMaxBigRIPSBeam];   //[BigRIPSBeam_]
   Int_t           BigRIPSBeam_id[kMaxBigRIPSBeam];   //[BigRIPSBeam_]
   Int_t           BigRIPSBeam_fpl[kMaxBigRIPSBeam];   //[BigRIPSBeam_]
   TString         BigRIPSBeam_name[kMaxBigRIPSBeam];
   Int_t           BigRIPSBeam_fDataState[kMaxBigRIPSBeam];   //[BigRIPSBeam_]
   Double_t        BigRIPSBeam_brho[kMaxBigRIPSBeam];   //[BigRIPSBeam_]
   Double_t        BigRIPSBeam_aoq[kMaxBigRIPSBeam];   //[BigRIPSBeam_]
   Double_t        BigRIPSBeam_zet[kMaxBigRIPSBeam];   //[BigRIPSBeam_]
   Double_t        BigRIPSBeam_beta[kMaxBigRIPSBeam];   //[BigRIPSBeam_]
   Double_t        BigRIPSBeam_clight[kMaxBigRIPSBeam];   //[BigRIPSBeam_]
   Double_t        BigRIPSBeam_mnucleon[kMaxBigRIPSBeam];   //[BigRIPSBeam_]
   Int_t           BigRIPSBeam_nrips[kMaxBigRIPSBeam];   //[BigRIPSBeam_]
   TString         BigRIPSBeam_ripsname[2][kMaxBigRIPSBeam];
   TString         BigRIPSBeam_tofname[kMaxBigRIPSBeam];
   TString         BigRIPSBeam_icname[kMaxBigRIPSBeam];
   Double_t        xtar;
   Double_t        ytar;
   Double_t        scattangle;
   Double_t        scattphi;
   Double_t        scatta;
   Double_t        scattb;
   Double_t        f8ppacx[6];
   Double_t        f8ppacy[6];
   Double_t        f8ppacz[6];
   Double_t        f8posx[3];
   Double_t        f8posy[3];
   Double_t        f8posz[3];
   Int_t           fgoodppacfocus[12];
   Int_t           fgoodppacfocusor[12];
   Double_t tof[3];
   Double_t beta[3];
 
   Double_t zet[6];
   Double_t aoq[6];
   Double_t zetc[6];
   Double_t aoqc[6];  //Corrected PID

   Double_t delta[4];  

   Int_t           dalimultwotime;
   Int_t           dalimult;
   Int_t           dalitimetruemult;
   Int_t           dalimultthres;
   Int_t           dalitimetruemultthres;

   Double_t        F3X;
   Double_t        F3A;
   Double_t        F3Y;
   Double_t        F3B;
   Double_t        F5X;
   Double_t        F5A;
   Double_t        F5Y;
   Double_t        F5B;
   Double_t        F7X;
   Double_t        F7A;
   Double_t        F7Y;
   Double_t        F7B;
   Double_t        F8X;
   Double_t        F8A;
   Double_t        F8Y;
   Double_t        F8B;
   Double_t        F9X;
   Double_t        F9A;
   Double_t        F9Y;
   Double_t        F9B;
   Double_t        F11X;
   Double_t        F11A;
   Double_t        F11Y;
   Double_t        F11B;
   Double_t        F3PLA_TL_raw;
   Double_t        F3PLA_TR_raw;
   Double_t        F3PLA_TL;
   Double_t        F3PLA_TR;
   Double_t        F3PLA_T;
   Double_t        F7PLA_TL_raw;
   Double_t        F7PLA_TR_raw;
   Double_t        F7PLA_TL;
   Double_t        F7PLA_TR;
   Double_t        F7PLA_T;
   Double_t        F8PLA_TL_raw;
   Double_t        F8PLA_TR_raw;
   Double_t        F8PLA_TL;
   Double_t        F8PLA_TR;
   Double_t        F8PLA_T;
   Double_t        F11PLA_TL_raw;
   Double_t        F11PLA_TR_raw;
   Double_t        F11PLA_TL;
   Double_t        F11PLA_TR;
   Double_t        F11PLA_T;
  Double_t           F3PLA_QL_raw;
   Double_t            F3PLA_QR_raw;
   Double_t           F3PLA_Q_ave;
   Double_t            F7PLA_QL_raw;
  Double_t            F7PLA_QR_raw;
   Double_t           F7PLA_Q_ave;
   Double_t            F8PLA_QL_raw;
   Double_t           F8PLA_QR_raw;
  Double_t            F8PLA_Q_ave;
   Double_t            F11PLA_QL_raw;
   Double_t           F11PLA_QR_raw;
   Double_t            F11PLA_Q_ave;
  Double_t          F7ICNumHit;
   Double_t        F7ICEnergySqSum;
   Double_t        F7ICEnergyAvSum;
   Double_t          F11ICNumHit;
   Double_t        F11ICEnergySqSum;
   Double_t        F11ICEnergyAvSum;
  
    
Double_t 	ratio_beta_57_brho_57;
Double_t 	ratio_beta_35_brho_35;
Double_t 	brho_35;
Double_t 	ratio_beta_911_brho_911;	
Double_t 	ratio_beta_89_brho_89;
Double_t 	brho_89;	
Double_t 	beta_57;	
Double_t	beta_911;

    
   // List of branches
   TBranch        *b_fbit;   //!
   TBranch        *b_EventInfo_;   //!
   TBranch        *b_EventInfo_fUniqueID;   //!
   TBranch        *b_EventInfo_fBits;   //!
   TBranch        *b_EventInfo_fName;   //!
   TBranch        *b_EventInfo_fTitle;   //!
   TBranch        *b_EventInfo_runnumber;   //!
   TBranch        *b_EventInfo_eventnumber;   //!
   TBranch        *b_EventInfo_subsysname;   //!
   TBranch        *b_EventInfo_timestamp;   //!
   TBranch        *b_EventInfo_comp_val;   //!
   TBranch        *b_EventInfo_fBit;   //!
   TBranch        *b_BigRIPSPPAC_;   //!
   TBranch        *b_BigRIPSPPAC_fUniqueID;   //!
   TBranch        *b_BigRIPSPPAC_fBits;   //!
   TBranch        *b_BigRIPSPPAC_id;   //!
   TBranch        *b_BigRIPSPPAC_fpl;   //!
   TBranch        *b_BigRIPSPPAC_name;   //!
   TBranch        *b_BigRIPSPPAC_fDataState;   //!
   TBranch        *b_BigRIPSPPAC_xzpos;   //!
   TBranch        *b_BigRIPSPPAC_yzpos;   //!
   TBranch        *b_BigRIPSPPAC_fTX1Raw;   //!
   TBranch        *b_BigRIPSPPAC_fTX2Raw;   //!
   TBranch        *b_BigRIPSPPAC_fTY1Raw;   //!
   TBranch        *b_BigRIPSPPAC_fTY2Raw;   //!
   TBranch        *b_BigRIPSPPAC_fTARaw;   //!
   TBranch        *b_BigRIPSPPAC_fQX1Raw;   //!
   TBranch        *b_BigRIPSPPAC_fQX2Raw;   //!
   TBranch        *b_BigRIPSPPAC_fQY1Raw;   //!
   TBranch        *b_BigRIPSPPAC_fQY2Raw;   //!
   TBranch        *b_BigRIPSPPAC_fQARaw;   //!
   TBranch        *b_BigRIPSPPAC_fTSumX;   //!
   TBranch        *b_BigRIPSPPAC_fTSumY;   //!
   TBranch        *b_BigRIPSPPAC_fTDiffX;   //!
   TBranch        *b_BigRIPSPPAC_fTDiffY;   //!
   TBranch        *b_BigRIPSPPAC_fX;   //!
   TBranch        *b_BigRIPSPPAC_fY;   //!
   TBranch        *b_BigRIPSPPAC_fFiredX;   //!
   TBranch        *b_BigRIPSPPAC_fFiredY;   //!
   TBranch        *b_BigRIPSPlastic_;   //!
   TBranch        *b_BigRIPSPlastic_fUniqueID;   //!
   TBranch        *b_BigRIPSPlastic_fBits;   //!
   TBranch        *b_BigRIPSPlastic_id;   //!
   TBranch        *b_BigRIPSPlastic_fpl;   //!
   TBranch        *b_BigRIPSPlastic_name;   //!
   TBranch        *b_BigRIPSPlastic_fDataState;   //!
   TBranch        *b_BigRIPSPlastic_zposition;   //!
   TBranch        *b_BigRIPSPlastic_zoffset;   //!
   TBranch        *b_BigRIPSPlastic_fTLRaw;   //!
   TBranch        *b_BigRIPSPlastic_fTRRaw;   //!
   TBranch        *b_BigRIPSPlastic_fQLRaw;   //!
   TBranch        *b_BigRIPSPlastic_fQRRaw;   //!
   TBranch        *b_BigRIPSPlastic_fTime;   //!
   TBranch        *b_BigRIPSPlastic_fTimeL;   //!
   TBranch        *b_BigRIPSPlastic_fTimeR;   //!
   TBranch        *b_BigRIPSPlastic_fTimeLSlew;   //!
   TBranch        *b_BigRIPSPlastic_fTimeRSlew;   //!
   TBranch        *b_BigRIPSPlastic_fTimeSlew;   //!
   TBranch        *b_BigRIPSIC_;   //!
   TBranch        *b_BigRIPSIC_fUniqueID;   //!
   TBranch        *b_BigRIPSIC_fBits;   //!
   TBranch        *b_BigRIPSIC_id;   //!
   TBranch        *b_BigRIPSIC_fpl;   //!
   TBranch        *b_BigRIPSIC_name;   //!
   TBranch        *b_BigRIPSIC_fDataState;   //!
   TBranch        *b_BigRIPSIC_zcoef;   //!
   TBranch        *b_BigRIPSIC_ionpair;   //!
   TBranch        *b_BigRIPSIC_nhitchannel;   //!
   TBranch        *b_BigRIPSIC_fADC;   //!
   TBranch        *b_BigRIPSIC_fEnergy;   //!
   TBranch        *b_BigRIPSIC_fRawADCSqSum;   //!
   TBranch        *b_BigRIPSIC_fRawADCAvSum;   //!
   TBranch        *b_BigRIPSIC_fCalMeVSqSum;   //!
   TBranch        *b_BigRIPSIC_fCalMeVAvSum;   //!
   TBranch        *b_BigRIPSFocalPlane_;   //!
   TBranch        *b_BigRIPSFocalPlane_fUniqueID;   //!
   TBranch        *b_BigRIPSFocalPlane_fBits;   //!
   TBranch        *b_BigRIPSFocalPlane_id;   //!
   TBranch        *b_BigRIPSFocalPlane_fpl;   //!
   TBranch        *b_BigRIPSFocalPlane_name;   //!
   TBranch        *b_BigRIPSFocalPlane_fDataState;   //!
   TBranch        *b_BigRIPSFocalPlane_opt_vector;   //!
   TBranch        *b_BigRIPSFocalPlane_nfired_ppacx;   //!
   TBranch        *b_BigRIPSFocalPlane_nfired_ppacy;   //!
   TBranch        *b_BigRIPSFocalPlane_zpos;   //!
   TBranch        *b_DALINaI_;   //!
   TBranch        *b_DALINaI_fUniqueID;   //!
   TBranch        *b_DALINaI_fBits;   //!
   TBranch        *b_DALINaI_id;   //!
   TBranch        *b_DALINaI_fpl;   //!
   TBranch        *b_DALINaI_name;   //!
   TBranch        *b_DALINaI_fDataState;   //!
   TBranch        *b_DALINaI_fADC;   //!
   TBranch        *b_DALINaI_fTDC;   //!
   TBranch        *b_DALINaI_layer;   //!
   TBranch        *b_DALINaI_fXPos;   //!
   TBranch        *b_DALINaI_fYPos;   //!
   TBranch        *b_DALINaI_fZPos;   //!
   TBranch        *b_DALINaI_theta;   //!
   TBranch        *b_DALINaI_costheta;   //!
   TBranch        *b_DALINaI_fEnergy;   //!
   TBranch        *b_DALINaI_fDoppCorEnergy;   //!
   TBranch        *b_DALINaI_fEnergyWithoutT;   //!
   TBranch        *b_DALINaI_fTime;   //!
   TBranch        *b_DALINaI_fTimeOffseted;   //!
   TBranch        *b_DALINaI_fTimeTrueEnergy;   //!
   TBranch        *b_DALINaI_fTimeTrueDoppCorEnergy;   //!
   TBranch        *b_DALINaI_fTimeTrueTime;   //!
   TBranch        *b_DALINaI_fTimeTrueTimeOffseted;   //!
   TBranch        *b_BigRIPSRIPS_;   //!
   TBranch        *b_BigRIPSRIPS_fUniqueID;   //!
   TBranch        *b_BigRIPSRIPS_fBits;   //!
   TBranch        *b_BigRIPSRIPS_id;   //!
   TBranch        *b_BigRIPSRIPS_fpl;   //!
   TBranch        *b_BigRIPSRIPS_name;   //!
   TBranch        *b_BigRIPSRIPS_fDataState;   //!
   TBranch        *b_BigRIPSRIPS_upstream_fpl;   //!
   TBranch        *b_BigRIPSRIPS_downstream_fpl;   //!
   TBranch        *b_BigRIPSRIPS_center_brho;   //!
   TBranch        *b_BigRIPSRIPS_brho;   //!
   TBranch        *b_BigRIPSRIPS_length;   //!
   TBranch        *b_BigRIPSRIPS_matrix;   //!
   TBranch        *b_BigRIPSRIPS_delta;   //!
   TBranch        *b_BigRIPSRIPS_angle;   //!
   TBranch        *b_BigRIPSTOF_;   //!
   TBranch        *b_BigRIPSTOF_fUniqueID;   //!
   TBranch        *b_BigRIPSTOF_fBits;   //!
   TBranch        *b_BigRIPSTOF_id;   //!
   TBranch        *b_BigRIPSTOF_fpl;   //!
   TBranch        *b_BigRIPSTOF_name;   //!
   TBranch        *b_BigRIPSTOF_fDataState;   //!
   TBranch        *b_BigRIPSTOF_tof;   //!
   TBranch        *b_BigRIPSTOF_clight;   //!
   TBranch        *b_BigRIPSTOF_length;   //!
   TBranch        *b_BigRIPSTOF_upstream_plname;   //!
   TBranch        *b_BigRIPSTOF_downstream_plname;   //!
   TBranch        *b_BigRIPSTOF_upstream_plfpl;   //!
   TBranch        *b_BigRIPSTOF_downstream_plfpl;   //!
   TBranch        *b_BigRIPSTOF_time_offset;   //!
   TBranch        *b_BigRIPSBeam_;   //!
   TBranch        *b_BigRIPSBeam_fUniqueID;   //!
   TBranch        *b_BigRIPSBeam_fBits;   //!
   TBranch        *b_BigRIPSBeam_id;   //!
   TBranch        *b_BigRIPSBeam_fpl;   //!
   TBranch        *b_BigRIPSBeam_name;   //!
   TBranch        *b_BigRIPSBeam_fDataState;   //!
   TBranch        *b_BigRIPSBeam_brho;   //!
   TBranch        *b_BigRIPSBeam_aoq;   //!
   TBranch        *b_BigRIPSBeam_zet;   //!
   TBranch        *b_BigRIPSBeam_beta;   //!
   TBranch        *b_BigRIPSBeam_clight;   //!
   TBranch        *b_BigRIPSBeam_mnucleon;   //!
   TBranch        *b_BigRIPSBeam_nrips;   //!
   TBranch        *b_BigRIPSBeam_ripsname;   //!
   TBranch        *b_BigRIPSBeam_tofname;   //!
   TBranch        *b_BigRIPSBeam_icname;   //!
   TBranch        *b_xtar;   //!
   TBranch        *b_ytar;   //!
   TBranch        *b_scattangle;   //!
   TBranch        *b_scattphi;   //!
   TBranch        *b_scatta;   //!
   TBranch        *b_scattb;   //!
   TBranch        *b_f8ppacx;   //!
   TBranch        *b_f8ppacy;   //!
   TBranch        *b_f8ppacz;   //!
   TBranch        *b_f8posx;   //!
   TBranch        *b_f8posy;   //!
   TBranch        *b_f8posz;   //!
   TBranch        *b_fgoodppacfocus;   //!
   TBranch        *b_fgoodppacfocusor;   //!
   TBranch        *b_tof;   //!
   TBranch        *b_beta;	 //!
   TBranch        *b_zet;	 //!
   TBranch        *b_aoq;	 //!
   TBranch        *b_zetc;	 //!
   TBranch        *b_aoqc;  //!
   TBranch        *b_delta;//!
   TBranch        *b_dalimultwotime;   //!
   TBranch        *b_dalimult;   //!
   TBranch        *b_dalitimetruemult;   //!
   TBranch        *b_dalimultthres;   //!
   TBranch        *b_dalitimetruemultthres;   //!
   

   TBranch        *b_F3X;   //!
   TBranch        *b_F3A;   //!
   TBranch        *b_F3Y;   //!
   TBranch        *b_F3B;   //!
   TBranch        *b_F5X;   //!
   TBranch        *b_F5A;   //!
   TBranch        *b_F5Y;   //!
   TBranch        *b_F5B;   //!
   TBranch        *b_F7X;   //!
   TBranch        *b_F7A;   //!
   TBranch        *b_F7Y;   //!
   TBranch        *b_F7B;   //!
   TBranch        *b_F8X;   //!
   TBranch        *b_F8A;   //!
   TBranch        *b_F8Y;   //!
   TBranch        *b_F8B;   //!
   TBranch        *b_F9X;   //!
   TBranch        *b_F9A;   //!
   TBranch        *b_F9Y;   //!
   TBranch        *b_F9B;   //!
   TBranch        *b_F11X;   //!
   TBranch        *b_F11A;   //!
   TBranch        *b_F11Y;   //!
   TBranch        *b_F11B;   //!
   TBranch        *b_F3PLA_TL_raw;   //!
   TBranch        *b_F3PLA_TR_raw;   //!
   TBranch        *b_F3PLA_TL;   //!
   TBranch        *b_F3PLA_TR;   //!
   TBranch        *b_F3PLA_T;   //!
   TBranch        *b_F7PLA_TL_raw;   //!
   TBranch        *b_F7PLA_TR_raw;   //!
   TBranch        *b_F7PLA_TL;   //!
   TBranch        *b_F7PLA_TR;   //!
   TBranch        *b_F7PLA_T;   //!
   TBranch        *b_F8PLA_TL_raw;   //!
   TBranch        *b_F8PLA_TR_raw;   //!
   TBranch        *b_F8PLA_TL;   //!
   TBranch        *b_F8PLA_TR;   //!
   TBranch        *b_F8PLA_T;   //!
   TBranch        *b_F11PLA_TL_raw;   //!
   TBranch        *b_F11PLA_TR_raw;   //!
   TBranch        *b_F11PLA_TL;   //!
   TBranch        *b_F11PLA_TR;   //!
   TBranch        *b_F11PLA_T;   //!
   TBranch        *b_F3PLA_QL_raw;   //!
   TBranch        *b_F3PLA_QR_raw;   //!
   TBranch        *b_F3PLA_Q_ave;   //!
   TBranch        *b_F7PLA_QL_raw;   //!
   TBranch        *b_F7PLA_QR_raw;   //!
   TBranch        *b_F7PLA_Q_ave;   //!
   TBranch        *b_F8PLA_QL_raw;   //!
   TBranch        *b_F8PLA_QR_raw;   //!
   TBranch        *b_F8PLA_Q_ave;   //!
   TBranch        *b_F11PLA_QL_raw;   //!
   TBranch        *b_F11PLA_QR_raw;   //!
   TBranch        *b_F11PLA_Q_ave;   //!
   TBranch        *b_F7ICNumHit;   //!
   TBranch        *b_F7ICEnergySqSum;   //!
   TBranch        *b_F7ICEnergyAvSum;   //!
   TBranch        *b_F11ICNumHit;   //!
   TBranch        *b_F11ICEnergySqSum;   //!
   TBranch        *b_F11ICEnergyAvSum;   //!
  
TBranch        *b_ratio_beta_57_brho_57;//!
TBranch        *b_ratio_beta_35_brho_35;//!
TBranch        *b_brho_35;//!
TBranch        *b_ratio_beta_911_brho_911;	//!
TBranch        *b_ratio_beta_89_brho_89;//!
TBranch        *b_brho_89;	//!
TBranch        *b_beta_57;	//!
TBranch        *b_beta_911;//!
    

   Event();
   virtual ~Event();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual void     Init();
   virtual Bool_t   Notify();
};


#endif 

