

{
#include "TFile.h"
#include "TTree.h"
#include <TArrow.h>
    
    
    TH1F* hist_0;
    TH1F* hist_1;
    
    TFile *f1 = new TFile("/home/shailendra/geant4/new_efficiency/Experimental_rootfiles/background_subtraction/c_137Cs_-2000_2000_15_4Feb.root");
    hist_0 = (TH1F*)f1->Get("hp_doppler_0");
   
    
    
    
    TFile *f2 = new TFile("converted_empty-target_0001.root");
    hist_1 = (TH1F*)f2->Get("hp_doppler_1");
    TFile *fout = new TFile("subtracted137Cswithoutaddback.root","RECREATE");
    
    hist_1->Scale(140);
    hist_0->Add(hist_1,-1);
    hist_0->Draw("HIST");
    
   // hist_0->Draw("HIST");
    //hist_1->Draw("HIST,same");
    hist_0->Write();
    
    fout->ls();
    fout->Close();
}
//Double_t scale = hist_0->GetYaxis()->GetBinWidth(34.81)*(hist_0->Integral());
//hist_0->Scale(scale);
