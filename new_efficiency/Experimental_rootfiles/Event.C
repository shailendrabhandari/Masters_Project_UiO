#include "Event.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

using namespace std;

Event::Event() {
  fChain = new TChain();

//fChain->AddFile("calib_0011.root",0,"tree");//60Co
//fChain->AddFile("calib_0009.root",0,"tree");//88Y
fChain->AddFile("calib_0008.root",0,"tree");//137Cs

/* */
    
  Init();
}

Event::~Event(){
}

Int_t Event::GetEntry(Long64_t entry){
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

void Event::Init() {
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   //if (!tree) return;
   //fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);
   fChain->SetBranchAddress("fbit", &fbit, &b_fbit);
   fChain->SetBranchAddress("EventInfo", &EventInfo_, &b_EventInfo_);
   fChain->SetBranchAddress("EventInfo.fUniqueID", EventInfo_fUniqueID, &b_EventInfo_fUniqueID);
   fChain->SetBranchAddress("EventInfo.fBits", EventInfo_fBits, &b_EventInfo_fBits);
   fChain->SetBranchAddress("EventInfo.fName", EventInfo_fName, &b_EventInfo_fName);
   fChain->SetBranchAddress("EventInfo.fTitle", EventInfo_fTitle, &b_EventInfo_fTitle);
   fChain->SetBranchAddress("EventInfo.runnumber", EventInfo_runnumber, &b_EventInfo_runnumber);
   fChain->SetBranchAddress("EventInfo.eventnumber", EventInfo_eventnumber, &b_EventInfo_eventnumber);
   fChain->SetBranchAddress("EventInfo.subsysname", EventInfo_subsysname, &b_EventInfo_subsysname);
   fChain->SetBranchAddress("EventInfo.timestamp", EventInfo_timestamp, &b_EventInfo_timestamp);
   fChain->SetBranchAddress("EventInfo.comp_val", EventInfo_comp_val, &b_EventInfo_comp_val);
   fChain->SetBranchAddress("EventInfo.fBit", EventInfo_fBit, &b_EventInfo_fBit);
   fChain->SetBranchAddress("BigRIPSPPAC", &BigRIPSPPAC_, &b_BigRIPSPPAC_);
   fChain->SetBranchAddress("BigRIPSPPAC.fUniqueID", BigRIPSPPAC_fUniqueID, &b_BigRIPSPPAC_fUniqueID);
   fChain->SetBranchAddress("BigRIPSPPAC.fBits", BigRIPSPPAC_fBits, &b_BigRIPSPPAC_fBits);
   fChain->SetBranchAddress("BigRIPSPPAC.id", BigRIPSPPAC_id, &b_BigRIPSPPAC_id); 
   fChain->SetBranchAddress("BigRIPSPPAC.fpl", BigRIPSPPAC_fpl, &b_BigRIPSPPAC_fpl);
   fChain->SetBranchAddress("BigRIPSPPAC.name", BigRIPSPPAC_name, &b_BigRIPSPPAC_name);
   fChain->SetBranchAddress("BigRIPSPPAC.fDataState", BigRIPSPPAC_fDataState, &b_BigRIPSPPAC_fDataState);
   fChain->SetBranchAddress("BigRIPSPPAC.xzpos", BigRIPSPPAC_xzpos, &b_BigRIPSPPAC_xzpos);
   fChain->SetBranchAddress("BigRIPSPPAC.yzpos", BigRIPSPPAC_yzpos, &b_BigRIPSPPAC_yzpos);
   fChain->SetBranchAddress("BigRIPSPPAC.fTX1Raw", BigRIPSPPAC_fTX1Raw, &b_BigRIPSPPAC_fTX1Raw);
   fChain->SetBranchAddress("BigRIPSPPAC.fTX2Raw", BigRIPSPPAC_fTX2Raw, &b_BigRIPSPPAC_fTX2Raw);
   fChain->SetBranchAddress("BigRIPSPPAC.fTY1Raw", BigRIPSPPAC_fTY1Raw, &b_BigRIPSPPAC_fTY1Raw);
   fChain->SetBranchAddress("BigRIPSPPAC.fTY2Raw", BigRIPSPPAC_fTY2Raw, &b_BigRIPSPPAC_fTY2Raw);
   fChain->SetBranchAddress("BigRIPSPPAC.fTARaw", BigRIPSPPAC_fTARaw, &b_BigRIPSPPAC_fTARaw);
   fChain->SetBranchAddress("BigRIPSPPAC.fQX1Raw", BigRIPSPPAC_fQX1Raw, &b_BigRIPSPPAC_fQX1Raw);
   fChain->SetBranchAddress("BigRIPSPPAC.fQX2Raw", BigRIPSPPAC_fQX2Raw, &b_BigRIPSPPAC_fQX2Raw);
   fChain->SetBranchAddress("BigRIPSPPAC.fQY1Raw", BigRIPSPPAC_fQY1Raw, &b_BigRIPSPPAC_fQY1Raw);
   fChain->SetBranchAddress("BigRIPSPPAC.fQY2Raw", BigRIPSPPAC_fQY2Raw, &b_BigRIPSPPAC_fQY2Raw);
   fChain->SetBranchAddress("BigRIPSPPAC.fQARaw", BigRIPSPPAC_fQARaw, &b_BigRIPSPPAC_fQARaw);
   fChain->SetBranchAddress("BigRIPSPPAC.fTSumX", BigRIPSPPAC_fTSumX, &b_BigRIPSPPAC_fTSumX);
   fChain->SetBranchAddress("BigRIPSPPAC.fTSumY", BigRIPSPPAC_fTSumY, &b_BigRIPSPPAC_fTSumY);
   fChain->SetBranchAddress("BigRIPSPPAC.fTDiffX", BigRIPSPPAC_fTDiffX, &b_BigRIPSPPAC_fTDiffX);
   fChain->SetBranchAddress("BigRIPSPPAC.fTDiffY", BigRIPSPPAC_fTDiffY, &b_BigRIPSPPAC_fTDiffY);
   fChain->SetBranchAddress("BigRIPSPPAC.fX", BigRIPSPPAC_fX, &b_BigRIPSPPAC_fX);
   fChain->SetBranchAddress("BigRIPSPPAC.fY", BigRIPSPPAC_fY, &b_BigRIPSPPAC_fY);
   fChain->SetBranchAddress("BigRIPSPPAC.fFiredX", BigRIPSPPAC_fFiredX, &b_BigRIPSPPAC_fFiredX);
   fChain->SetBranchAddress("BigRIPSPPAC.fFiredY", BigRIPSPPAC_fFiredY, &b_BigRIPSPPAC_fFiredY);
   fChain->SetBranchAddress("BigRIPSPlastic", &BigRIPSPlastic_, &b_BigRIPSPlastic_);
   fChain->SetBranchAddress("BigRIPSPlastic.fUniqueID", BigRIPSPlastic_fUniqueID, &b_BigRIPSPlastic_fUniqueID);
   fChain->SetBranchAddress("BigRIPSPlastic.fBits", BigRIPSPlastic_fBits, &b_BigRIPSPlastic_fBits);
   fChain->SetBranchAddress("BigRIPSPlastic.id", BigRIPSPlastic_id, &b_BigRIPSPlastic_id);
   fChain->SetBranchAddress("BigRIPSPlastic.fpl", BigRIPSPlastic_fpl, &b_BigRIPSPlastic_fpl);
   fChain->SetBranchAddress("BigRIPSPlastic.name", BigRIPSPlastic_name, &b_BigRIPSPlastic_name);
   fChain->SetBranchAddress("BigRIPSPlastic.fDataState", BigRIPSPlastic_fDataState, &b_BigRIPSPlastic_fDataState);
   fChain->SetBranchAddress("BigRIPSPlastic.zposition", BigRIPSPlastic_zposition, &b_BigRIPSPlastic_zposition);
   fChain->SetBranchAddress("BigRIPSPlastic.zoffset", BigRIPSPlastic_zoffset, &b_BigRIPSPlastic_zoffset);
   fChain->SetBranchAddress("BigRIPSPlastic.fTLRaw", BigRIPSPlastic_fTLRaw, &b_BigRIPSPlastic_fTLRaw);
   fChain->SetBranchAddress("BigRIPSPlastic.fTRRaw", BigRIPSPlastic_fTRRaw, &b_BigRIPSPlastic_fTRRaw);
   fChain->SetBranchAddress("BigRIPSPlastic.fQLRaw", BigRIPSPlastic_fQLRaw, &b_BigRIPSPlastic_fQLRaw);
   fChain->SetBranchAddress("BigRIPSPlastic.fQRRaw", BigRIPSPlastic_fQRRaw, &b_BigRIPSPlastic_fQRRaw);
   fChain->SetBranchAddress("BigRIPSPlastic.fTime", BigRIPSPlastic_fTime, &b_BigRIPSPlastic_fTime);
   fChain->SetBranchAddress("BigRIPSPlastic.fTimeL", BigRIPSPlastic_fTimeL, &b_BigRIPSPlastic_fTimeL);
   fChain->SetBranchAddress("BigRIPSPlastic.fTimeR", BigRIPSPlastic_fTimeR, &b_BigRIPSPlastic_fTimeR);
   fChain->SetBranchAddress("BigRIPSPlastic.fTimeLSlew", BigRIPSPlastic_fTimeLSlew, &b_BigRIPSPlastic_fTimeLSlew);
   fChain->SetBranchAddress("BigRIPSPlastic.fTimeRSlew", BigRIPSPlastic_fTimeRSlew, &b_BigRIPSPlastic_fTimeRSlew);
   fChain->SetBranchAddress("BigRIPSPlastic.fTimeSlew", BigRIPSPlastic_fTimeSlew, &b_BigRIPSPlastic_fTimeSlew);
   fChain->SetBranchAddress("BigRIPSIC", &BigRIPSIC_, &b_BigRIPSIC_);
   fChain->SetBranchAddress("BigRIPSIC.fUniqueID", BigRIPSIC_fUniqueID, &b_BigRIPSIC_fUniqueID);
   fChain->SetBranchAddress("BigRIPSIC.fBits", BigRIPSIC_fBits, &b_BigRIPSIC_fBits);
   fChain->SetBranchAddress("BigRIPSIC.id", BigRIPSIC_id, &b_BigRIPSIC_id);
   fChain->SetBranchAddress("BigRIPSIC.fpl", BigRIPSIC_fpl, &b_BigRIPSIC_fpl);
   fChain->SetBranchAddress("BigRIPSIC.name", BigRIPSIC_name, &b_BigRIPSIC_name);
   fChain->SetBranchAddress("BigRIPSIC.fDataState", BigRIPSIC_fDataState, &b_BigRIPSIC_fDataState);
   fChain->SetBranchAddress("BigRIPSIC.zcoef[2]", BigRIPSIC_zcoef, &b_BigRIPSIC_zcoef);
   fChain->SetBranchAddress("BigRIPSIC.ionpair", BigRIPSIC_ionpair, &b_BigRIPSIC_ionpair);
   fChain->SetBranchAddress("BigRIPSIC.nhitchannel", BigRIPSIC_nhitchannel, &b_BigRIPSIC_nhitchannel);
   //fChain->SetBranchAddress("BigRIPSIC.fADC[6]", BigRIPSIC_fADC, &b_BigRIPSIC_fADC);
   //fChain->SetBranchAddress("BigRIPSIC.fEnergy[6]", BigRIPSIC_fEnergy, &b_BigRIPSIC_fEnergy);
   fChain->SetBranchAddress("BigRIPSIC.fRawADCSqSum", BigRIPSIC_fRawADCSqSum, &b_BigRIPSIC_fRawADCSqSum);
   fChain->SetBranchAddress("BigRIPSIC.fRawADCAvSum", BigRIPSIC_fRawADCAvSum, &b_BigRIPSIC_fRawADCAvSum);
   fChain->SetBranchAddress("BigRIPSIC.fCalMeVSqSum", BigRIPSIC_fCalMeVSqSum, &b_BigRIPSIC_fCalMeVSqSum);
   fChain->SetBranchAddress("BigRIPSIC.fCalMeVAvSum", BigRIPSIC_fCalMeVAvSum, &b_BigRIPSIC_fCalMeVAvSum);
   fChain->SetBranchAddress("BigRIPSFocalPlane", &BigRIPSFocalPlane_, &b_BigRIPSFocalPlane_);
   fChain->SetBranchAddress("BigRIPSFocalPlane.fUniqueID", BigRIPSFocalPlane_fUniqueID, &b_BigRIPSFocalPlane_fUniqueID);
   fChain->SetBranchAddress("BigRIPSFocalPlane.fBits", BigRIPSFocalPlane_fBits, &b_BigRIPSFocalPlane_fBits);
   fChain->SetBranchAddress("BigRIPSFocalPlane.id", BigRIPSFocalPlane_id, &b_BigRIPSFocalPlane_id);
   fChain->SetBranchAddress("BigRIPSFocalPlane.fpl", BigRIPSFocalPlane_fpl, &b_BigRIPSFocalPlane_fpl);
   fChain->SetBranchAddress("BigRIPSFocalPlane.name", BigRIPSFocalPlane_name, &b_BigRIPSFocalPlane_name);
   fChain->SetBranchAddress("BigRIPSFocalPlane.fDataState", BigRIPSFocalPlane_fDataState, &b_BigRIPSFocalPlane_fDataState);
   fChain->SetBranchAddress("BigRIPSFocalPlane.opt_vector", BigRIPSFocalPlane_opt_vector, &b_BigRIPSFocalPlane_opt_vector);
   fChain->SetBranchAddress("BigRIPSFocalPlane.nfired_ppacx", BigRIPSFocalPlane_nfired_ppacx, &b_BigRIPSFocalPlane_nfired_ppacx);
   fChain->SetBranchAddress("BigRIPSFocalPlane.nfired_ppacy", BigRIPSFocalPlane_nfired_ppacy, &b_BigRIPSFocalPlane_nfired_ppacy);
   fChain->SetBranchAddress("BigRIPSFocalPlane.zpos", BigRIPSFocalPlane_zpos, &b_BigRIPSFocalPlane_zpos);
   fChain->SetBranchAddress("DALINaI", &DALINaI_, &b_DALINaI_);
   fChain->SetBranchAddress("DALINaI.fUniqueID", DALINaI_fUniqueID, &b_DALINaI_fUniqueID);
   fChain->SetBranchAddress("DALINaI.fBits", DALINaI_fBits, &b_DALINaI_fBits);
   fChain->SetBranchAddress("DALINaI.id", DALINaI_id, &b_DALINaI_id);
   fChain->SetBranchAddress("DALINaI.fpl", DALINaI_fpl, &b_DALINaI_fpl);
   fChain->SetBranchAddress("DALINaI.name", DALINaI_name, &b_DALINaI_name);
   fChain->SetBranchAddress("DALINaI.fDataState", DALINaI_fDataState, &b_DALINaI_fDataState);
   fChain->SetBranchAddress("DALINaI.fADC", DALINaI_fADC, &b_DALINaI_fADC);
   fChain->SetBranchAddress("DALINaI.fTDC", DALINaI_fTDC, &b_DALINaI_fTDC);
   fChain->SetBranchAddress("DALINaI.layer", DALINaI_layer, &b_DALINaI_layer);
   fChain->SetBranchAddress("DALINaI.theta", DALINaI_theta, &b_DALINaI_theta);
   fChain->SetBranchAddress("DALINaI.fXPos", DALINaI_fXPos, &b_DALINaI_fXPos);
   fChain->SetBranchAddress("DALINaI.fYPos", DALINaI_fYPos, &b_DALINaI_fYPos);
   fChain->SetBranchAddress("DALINaI.fZPos", DALINaI_fZPos, &b_DALINaI_fZPos);
   fChain->SetBranchAddress("DALINaI.costheta", DALINaI_costheta, &b_DALINaI_costheta);
   fChain->SetBranchAddress("DALINaI.fEnergy", DALINaI_fEnergy, &b_DALINaI_fEnergy);
   fChain->SetBranchAddress("DALINaI.fDoppCorEnergy", DALINaI_fDoppCorEnergy, &b_DALINaI_fDoppCorEnergy);
   fChain->SetBranchAddress("DALINaI.fEnergyWithoutT", DALINaI_fEnergyWithoutT, &b_DALINaI_fEnergyWithoutT);
   fChain->SetBranchAddress("DALINaI.fTime", DALINaI_fTime, &b_DALINaI_fTime);
   fChain->SetBranchAddress("DALINaI.fTimeOffseted", DALINaI_fTimeOffseted, &b_DALINaI_fTimeOffseted);
   fChain->SetBranchAddress("DALINaI.fTimeTrueEnergy", DALINaI_fTimeTrueEnergy, &b_DALINaI_fTimeTrueEnergy);
   fChain->SetBranchAddress("DALINaI.fTimeTrueDoppCorEnergy", DALINaI_fTimeTrueDoppCorEnergy, &b_DALINaI_fTimeTrueDoppCorEnergy);
   fChain->SetBranchAddress("DALINaI.fTimeTrueTime", DALINaI_fTimeTrueTime, &b_DALINaI_fTimeTrueTime);
   fChain->SetBranchAddress("DALINaI.fTimeTrueTimeOffseted", DALINaI_fTimeTrueTimeOffseted, &b_DALINaI_fTimeTrueTimeOffseted);
   fChain->SetBranchAddress("BigRIPSRIPS", &BigRIPSRIPS_, &b_BigRIPSRIPS_);
   fChain->SetBranchAddress("BigRIPSRIPS.fUniqueID", BigRIPSRIPS_fUniqueID, &b_BigRIPSRIPS_fUniqueID);
   fChain->SetBranchAddress("BigRIPSRIPS.fBits", BigRIPSRIPS_fBits, &b_BigRIPSRIPS_fBits);
   fChain->SetBranchAddress("BigRIPSRIPS.id", BigRIPSRIPS_id, &b_BigRIPSRIPS_id);
   fChain->SetBranchAddress("BigRIPSRIPS.fpl", BigRIPSRIPS_fpl, &b_BigRIPSRIPS_fpl);
   fChain->SetBranchAddress("BigRIPSRIPS.name", BigRIPSRIPS_name, &b_BigRIPSRIPS_name);
   fChain->SetBranchAddress("BigRIPSRIPS.fDataState", BigRIPSRIPS_fDataState, &b_BigRIPSRIPS_fDataState);
   fChain->SetBranchAddress("BigRIPSRIPS.upstream_fpl", BigRIPSRIPS_upstream_fpl, &b_BigRIPSRIPS_upstream_fpl);
   fChain->SetBranchAddress("BigRIPSRIPS.downstream_fpl", BigRIPSRIPS_downstream_fpl, &b_BigRIPSRIPS_downstream_fpl);
   fChain->SetBranchAddress("BigRIPSRIPS.center_brho", BigRIPSRIPS_center_brho, &b_BigRIPSRIPS_center_brho);
   fChain->SetBranchAddress("BigRIPSRIPS.brho", BigRIPSRIPS_brho, &b_BigRIPSRIPS_brho);
   fChain->SetBranchAddress("BigRIPSRIPS.length", BigRIPSRIPS_length, &b_BigRIPSRIPS_length);
   fChain->SetBranchAddress("BigRIPSRIPS.matrix", BigRIPSRIPS_matrix, &b_BigRIPSRIPS_matrix);
   fChain->SetBranchAddress("BigRIPSRIPS.delta", BigRIPSRIPS_delta, &b_BigRIPSRIPS_delta);
   fChain->SetBranchAddress("BigRIPSRIPS.angle", BigRIPSRIPS_angle, &b_BigRIPSRIPS_angle);
   fChain->SetBranchAddress("BigRIPSTOF", &BigRIPSTOF_, &b_BigRIPSTOF_);
   fChain->SetBranchAddress("BigRIPSTOF.fUniqueID", BigRIPSTOF_fUniqueID, &b_BigRIPSTOF_fUniqueID);
   fChain->SetBranchAddress("BigRIPSTOF.fBits", BigRIPSTOF_fBits, &b_BigRIPSTOF_fBits);
   fChain->SetBranchAddress("BigRIPSTOF.id", BigRIPSTOF_id, &b_BigRIPSTOF_id);
   fChain->SetBranchAddress("BigRIPSTOF.fpl", BigRIPSTOF_fpl, &b_BigRIPSTOF_fpl);
   fChain->SetBranchAddress("BigRIPSTOF.name", BigRIPSTOF_name, &b_BigRIPSTOF_name);
   fChain->SetBranchAddress("BigRIPSTOF.fDataState", BigRIPSTOF_fDataState, &b_BigRIPSTOF_fDataState);
   fChain->SetBranchAddress("BigRIPSTOF.tof", BigRIPSTOF_tof, &b_BigRIPSTOF_tof);
   fChain->SetBranchAddress("BigRIPSTOF.clight", BigRIPSTOF_clight, &b_BigRIPSTOF_clight);
   fChain->SetBranchAddress("BigRIPSTOF.length", BigRIPSTOF_length, &b_BigRIPSTOF_length);
   fChain->SetBranchAddress("BigRIPSTOF.upstream_plname", BigRIPSTOF_upstream_plname, &b_BigRIPSTOF_upstream_plname);
   fChain->SetBranchAddress("BigRIPSTOF.downstream_plname", BigRIPSTOF_downstream_plname, &b_BigRIPSTOF_downstream_plname);
   fChain->SetBranchAddress("BigRIPSTOF.upstream_plfpl", BigRIPSTOF_upstream_plfpl, &b_BigRIPSTOF_upstream_plfpl);
   fChain->SetBranchAddress("BigRIPSTOF.downstream_plfpl", BigRIPSTOF_downstream_plfpl, &b_BigRIPSTOF_downstream_plfpl);
   fChain->SetBranchAddress("BigRIPSTOF.time_offset", BigRIPSTOF_time_offset, &b_BigRIPSTOF_time_offset);
   fChain->SetBranchAddress("BigRIPSBeam", &BigRIPSBeam_, &b_BigRIPSBeam_);
   fChain->SetBranchAddress("BigRIPSBeam.fUniqueID", BigRIPSBeam_fUniqueID, &b_BigRIPSBeam_fUniqueID);
   fChain->SetBranchAddress("BigRIPSBeam.fBits", BigRIPSBeam_fBits, &b_BigRIPSBeam_fBits);
   fChain->SetBranchAddress("BigRIPSBeam.id", BigRIPSBeam_id, &b_BigRIPSBeam_id);
   fChain->SetBranchAddress("BigRIPSBeam.fpl", BigRIPSBeam_fpl, &b_BigRIPSBeam_fpl);
   fChain->SetBranchAddress("BigRIPSBeam.name", BigRIPSBeam_name, &b_BigRIPSBeam_name);
   fChain->SetBranchAddress("BigRIPSBeam.fDataState", BigRIPSBeam_fDataState, &b_BigRIPSBeam_fDataState);
   fChain->SetBranchAddress("BigRIPSBeam.brho", BigRIPSBeam_brho, &b_BigRIPSBeam_brho);
   fChain->SetBranchAddress("BigRIPSBeam.aoq", BigRIPSBeam_aoq, &b_BigRIPSBeam_aoq);
   fChain->SetBranchAddress("BigRIPSBeam.zet", BigRIPSBeam_zet, &b_BigRIPSBeam_zet);
   fChain->SetBranchAddress("BigRIPSBeam.beta", BigRIPSBeam_beta, &b_BigRIPSBeam_beta);
   fChain->SetBranchAddress("BigRIPSBeam.clight", BigRIPSBeam_clight, &b_BigRIPSBeam_clight);
   fChain->SetBranchAddress("BigRIPSBeam.mnucleon", BigRIPSBeam_mnucleon, &b_BigRIPSBeam_mnucleon);
   fChain->SetBranchAddress("BigRIPSBeam.nrips", BigRIPSBeam_nrips, &b_BigRIPSBeam_nrips);
   fChain->SetBranchAddress("BigRIPSBeam.ripsname[2]", BigRIPSBeam_ripsname, &b_BigRIPSBeam_ripsname);
   fChain->SetBranchAddress("BigRIPSBeam.tofname", BigRIPSBeam_tofname, &b_BigRIPSBeam_tofname);
   fChain->SetBranchAddress("BigRIPSBeam.icname", BigRIPSBeam_icname, &b_BigRIPSBeam_icname);
   fChain->SetBranchAddress("xtar", &xtar, &b_xtar);
   fChain->SetBranchAddress("ytar", &ytar, &b_ytar);
   fChain->SetBranchAddress("scattangle", &scattangle, &b_scattangle);
   fChain->SetBranchAddress("scattphi", &scattphi, &b_scattphi);
   fChain->SetBranchAddress("scatta", &scatta, &b_scatta);
   fChain->SetBranchAddress("scattb", &scattb, &b_scattb);
   fChain->SetBranchAddress("f8ppacx", f8ppacx, &b_f8ppacx);
   fChain->SetBranchAddress("f8ppacy", f8ppacy, &b_f8ppacy);
   //fChain->SetBranchAddress("f8ppacz", f8ppacz, &b_f8ppacz);
   fChain->SetBranchAddress("f8posx", f8posx, &b_f8posx);
   fChain->SetBranchAddress("f8posy", f8posy, &b_f8posy);
   //fChain->SetBranchAddress("f8posz", f8posz, &b_f8posz);
   fChain->SetBranchAddress("tof",tof,&b_tof);
   fChain->SetBranchAddress("beta",beta,&b_beta);   
   fChain->SetBranchAddress("zet",zet,&b_zet);
   fChain->SetBranchAddress("aoq",aoq,&b_aoq);
   fChain->SetBranchAddress("zetc",zetc,&b_zetc);
   fChain->SetBranchAddress("aoqc",aoqc,&b_aoqc);
   fChain->SetBranchAddress("delta",delta,&b_delta);
   fChain->SetBranchAddress("fgoodppacfocus",fgoodppacfocus,&b_fgoodppacfocus);
   fChain->SetBranchAddress("fgoodppacfocusor",fgoodppacfocusor,&b_fgoodppacfocusor);
   fChain->SetBranchAddress("dalimultwotime", &dalimultwotime, &b_dalimultwotime);
   fChain->SetBranchAddress("dalimult", &dalimult, &b_dalimult);
   fChain->SetBranchAddress("dalitimetruemult", &dalitimetruemult, &b_dalitimetruemult);
   fChain->SetBranchAddress("dalimultthres", &dalimultthres, &b_dalimultthres);
   fChain->SetBranchAddress("dalitimetruemultthres", &dalitimetruemultthres, &b_dalitimetruemultthres);

   fChain->SetBranchAddress("F3X", &F3X, &b_F3X);
   fChain->SetBranchAddress("F3A", &F3A, &b_F3A);
   fChain->SetBranchAddress("F3Y", &F3Y, &b_F3Y);
   fChain->SetBranchAddress("F3B", &F3B, &b_F3B);
   fChain->SetBranchAddress("F5X", &F5X, &b_F5X);
   fChain->SetBranchAddress("F5A", &F5A, &b_F5A);
   fChain->SetBranchAddress("F5Y", &F5Y, &b_F5Y);
   fChain->SetBranchAddress("F5B", &F5B, &b_F5B);
   fChain->SetBranchAddress("F7X", &F7X, &b_F7X);
   fChain->SetBranchAddress("F7A", &F7A, &b_F7A);
   fChain->SetBranchAddress("F7Y", &F7Y, &b_F7Y);
   fChain->SetBranchAddress("F7B", &F7B, &b_F7B);
   fChain->SetBranchAddress("F8X", &F8X, &b_F8X);
   fChain->SetBranchAddress("F8A", &F8A, &b_F8A);
   fChain->SetBranchAddress("F8Y", &F8Y, &b_F8Y);
   fChain->SetBranchAddress("F8B", &F8B, &b_F8B);
   fChain->SetBranchAddress("F9X", &F9X, &b_F9X);
   fChain->SetBranchAddress("F9A", &F9A, &b_F9A);
   fChain->SetBranchAddress("F9Y", &F9Y, &b_F9Y);
   fChain->SetBranchAddress("F9B", &F9B, &b_F9B);
   fChain->SetBranchAddress("F11X", &F11X, &b_F11X);
   fChain->SetBranchAddress("F11A", &F11A, &b_F11A);
   fChain->SetBranchAddress("F11Y", &F11Y, &b_F11Y);
   fChain->SetBranchAddress("F11B", &F11B, &b_F11B);
   fChain->SetBranchAddress("F3PLA_TL_raw", &F3PLA_TL_raw, &b_F3PLA_TL_raw);
   fChain->SetBranchAddress("F3PLA_TR_raw", &F3PLA_TR_raw, &b_F3PLA_TR_raw);
   fChain->SetBranchAddress("F3PLA_TL", &F3PLA_TL, &b_F3PLA_TL);
   fChain->SetBranchAddress("F3PLA_TR", &F3PLA_TR, &b_F3PLA_TR);
   fChain->SetBranchAddress("F3PLA_T", &F3PLA_T, &b_F3PLA_T);
   fChain->SetBranchAddress("F7PLA_TL_raw", &F7PLA_TL_raw, &b_F7PLA_TL_raw);
   fChain->SetBranchAddress("F7PLA_TR_raw", &F7PLA_TR_raw, &b_F7PLA_TR_raw);
   fChain->SetBranchAddress("F7PLA_TL", &F7PLA_TL, &b_F7PLA_TL);
   fChain->SetBranchAddress("F7PLA_TR", &F7PLA_TR, &b_F7PLA_TR);
   fChain->SetBranchAddress("F7PLA_T", &F7PLA_T, &b_F7PLA_T);
   fChain->SetBranchAddress("F8PLA_TL_raw", &F8PLA_TL_raw, &b_F8PLA_TL_raw);
   fChain->SetBranchAddress("F8PLA_TR_raw", &F8PLA_TR_raw, &b_F8PLA_TR_raw);
   fChain->SetBranchAddress("F8PLA_TL", &F8PLA_TL, &b_F8PLA_TL);
   fChain->SetBranchAddress("F8PLA_TR", &F8PLA_TR, &b_F8PLA_TR);
   fChain->SetBranchAddress("F8PLA_T", &F8PLA_T, &b_F8PLA_T);
   fChain->SetBranchAddress("F11PLA_TL_raw", &F11PLA_TL_raw, &b_F11PLA_TL_raw);
   fChain->SetBranchAddress("F11PLA_TR_raw", &F11PLA_TR_raw, &b_F11PLA_TR_raw);
   fChain->SetBranchAddress("F11PLA_TL", &F11PLA_TL, &b_F11PLA_TL);
   fChain->SetBranchAddress("F11PLA_TR", &F11PLA_TR, &b_F11PLA_TR);
   fChain->SetBranchAddress("F11PLA_T", &F11PLA_T, &b_F11PLA_T);
   fChain->SetBranchAddress("F3PLA_QL_raw", &F3PLA_QL_raw, &b_F3PLA_QL_raw);
   fChain->SetBranchAddress("F3PLA_QR_raw", &F3PLA_QR_raw, &b_F3PLA_QR_raw);
   fChain->SetBranchAddress("F3PLA_Q_ave", &F3PLA_Q_ave, &b_F3PLA_Q_ave);
   fChain->SetBranchAddress("F7PLA_QL_raw", &F7PLA_QL_raw, &b_F7PLA_QL_raw);
   fChain->SetBranchAddress("F7PLA_QR_raw", &F7PLA_QR_raw, &b_F7PLA_QR_raw);
   fChain->SetBranchAddress("F7PLA_Q_ave", &F7PLA_Q_ave, &b_F7PLA_Q_ave);
   fChain->SetBranchAddress("F8PLA_QL_raw", &F8PLA_QL_raw, &b_F8PLA_QL_raw);
   fChain->SetBranchAddress("F8PLA_QR_raw", &F8PLA_QR_raw, &b_F8PLA_QR_raw);
   fChain->SetBranchAddress("F8PLA_Q_ave", &F8PLA_Q_ave, &b_F8PLA_Q_ave);
   fChain->SetBranchAddress("F11PLA_QL_raw", &F11PLA_QL_raw, &b_F11PLA_QL_raw);
   fChain->SetBranchAddress("F11PLA_QR_raw", &F11PLA_QR_raw, &b_F11PLA_QR_raw);
   fChain->SetBranchAddress("F11PLA_Q_ave", &F11PLA_Q_ave, &b_F11PLA_Q_ave);
   fChain->SetBranchAddress("F7ICNumHit", &F7ICNumHit, &b_F7ICNumHit);
   fChain->SetBranchAddress("F7ICEnergySqSum", &F7ICEnergySqSum, &b_F7ICEnergySqSum);
   fChain->SetBranchAddress("F7ICEnergyAvSum", &F7ICEnergyAvSum, &b_F7ICEnergyAvSum);
   fChain->SetBranchAddress("F11ICNumHit", &F11ICNumHit, &b_F11ICNumHit);
   fChain->SetBranchAddress("F11ICEnergySqSum", &F11ICEnergySqSum, &b_F11ICEnergySqSum);
   fChain->SetBranchAddress("F11ICEnergyAvSum", &F11ICEnergyAvSum, &b_F11ICEnergyAvSum);
  
   /*fChain->SetBranchAddress("ratio_beta_57_brho_57",&ratio_beta_57_brho_57,&b_ratio_beta_57_brho_57);
   fChain->SetBranchAddress("ratio_beta_35_brho_35",&ratio_beta_35_brho_35,&b_ratio_beta_35_brho_35);
   fChain->SetBranchAddress("brho_35",&brho_35,&b_brho_35);
   fChain->SetBranchAddress("ratio_beta_911_brho_911",&ratio_beta_911_brho_911,&b_ratio_beta_911_brho_911);	
   fChain->SetBranchAddress("ratio_beta_89_brho_89",&ratio_beta_89_brho_89,&b_ratio_beta_89_brho_89);
   fChain->SetBranchAddress("brho_89",&brho_89,&b_brho_89);	
   fChain->SetBranchAddress("beta_57",&beta_57,&b_beta_57);	
   fChain->SetBranchAddress("beta_911",&beta_911,&b_beta_911);*/
	

   Notify();
}

Bool_t Event::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


