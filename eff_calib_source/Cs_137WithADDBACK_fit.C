/// 23.12.2020
#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TStyle.h"

enum ParIndex_t {
   Bkg0=0, Bkg1=1, Bkg2,
   SigScale, SigSigma, SigMean,
   N_PAR};
// Use this map to (re-)name parameters for the plot
const std::map<ParIndex_t,std::string> parNames{
   {Bkg0, "Bkg0"}, {Bkg1, "Bkg1"}, {Bkg2, "Bkg2"},
   {SigScale, "Gauss scale"}, {SigSigma, "Gauss #sigma"}, {SigMean, "Gauss #mu"}
};
 
//Background function
Double_t background(Double_t *x, Double_t *par) {
   return par[Bkg0] + par[Bkg1]*x[0] + par[Bkg2]*x[0]*x[0];
}
 
// Gauss Peak function (the signal function)
Double_t signal(Double_t *x, Double_t *par) {
   return par[SigScale]*TMath::Gaus(x[0], par[SigMean], par[SigSigma], true);
}
 
// Sum of background and peak function. We pass x and the fit parameters
// down to the signal and background functions.
Double_t fitFunction(Double_t *x, Double_t *par) {
   return background(x, par) + signal(x, par);
}

 
// Fit "fitFunction" to the histogram, and draw results on the canvas `c1`.
void FitRoutine(TCanvas* c1, TH1* histo, float fitxmin, float fitxmax, TString filename){
   c1->cd();
   // create a TF1 with the range and N_PAR parameters (six by default)
   TF1 fitFcn("fitFcn",fitFunction,450,1050,N_PAR);
   fitFcn.SetNpx(500);
   fitFcn.SetLineWidth(2);
   fitFcn.SetLineColor(kRed);
  
   // Assign the names from the map "parNames".
   for (auto& idx_name : parNames) {
     fitFcn.SetParName(idx_name.first, idx_name.second.c_str());
   }
 
   // Fit. First set ok-ish starting values for the parameters(I even dont understand it completely)
   fitFcn.SetParameters(0,0,0,0,1000,0,1000,0,2000,0);
   //fitFcn.SetParameters(0,0,0,0,200,2000);
   histo->GetXaxis()->SetRange(60,280);
   histo->Fit("fitFcn","VR+","ep");
  
 
   // improve the picture:
   // Draw signal and background functions separately
   TF1 backFcn("backFcn",background,535,775,N_PAR);
   backFcn.SetLineColor(kBlue);
   TF1 signalFcn("signalFcn",signal,500,850,N_PAR);
   signalFcn.SetLineColor(kBlue);
   signalFcn.SetNpx(500);
  
  
  //Integral of a function
    double w=histo->GetBinWidth(1);
    //int integral1 = backFcn.Integral(585,750)/w;
    int integral2 = fitFcn.Integral(480,980)/w;
    int integral3 = fitFcn.IntegralError(480,980)/w;

    //cout<<"Integral of backFcn::::::"<<integral1<<endl;
    cout<<"Integral of fitFcn ::::::"<<integral2<<endl;
    cout<<"Integralerror of fitFcn ::::::"<<integral3<<endl;
   // Retrieve fit parameters, and copy them to the signal and background functions
   Double_t par[N_PAR];
   fitFcn.GetParameters(par);
 
   backFcn.SetParameters(par);
   backFcn.DrawCopy("same");
 
   signalFcn.SetParameters(par);
   signalFcn.SetLineColor(kGreen);
   signalFcn.DrawCopy("same");
 }
   /* / draw the legend
   TLegend legend(0.55,0.9,0.90,0.75);
   legend.SetTextFont(20);
   legend.SetTextSize(0.03);
   legend.AddEntry(histo,"Exp. DataWithAddback","lpe");
   legend.AddEntry(&backFcn,"Background","l");
   legend.AddEntry(&signalFcn,"Signal","l");
   legend.AddEntry(&fitFcn,"Signal+Bkgd","l");
  // legend.AddEntry(&fitFcn,"Area of peak=1625065","l");
   legend.DrawClone();
   histo->Draw("esame");
   c1->SaveAs(filename);
*/
 

// Fitting histogram from the root file
void Cs_137WithADDBACK_fit() {

   gStyle->SetOptFit(1111);
   gStyle->SetOptStat(0);
 
   // fit range from fitxmin to fitxmax
   float fitxmin=100;
   float fitxmax=3000;

    TCanvas *c1 = new TCanvas("c1","Fitting 137Cs",10,10,700,700);
    TFile* f = new TFile("/home/shailendra/geant4/new_efficiency/fitting/137Cs/experimental/subtracted137Cswithaddback.root");
    //TH1F *h = (TH1F*) f->Get("hp_doppler_0");
    TH1F *h = (TH1F*) f->Get("hp_doppler_3");
   if (!h){
      cout << "h not found"<<endl;
      return;
   }
   h->SetMarkerStyle(2);
   h->SetMarkerSize(1);
   
   // The output of the fit. 
   FitRoutine(c1,h, fitxmin, fitxmax,"137Cs.pdf");
}
