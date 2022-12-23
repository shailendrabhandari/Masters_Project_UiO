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





////For peak 1
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


///For Peak 2


 
//Background function
Double_t background1(Double_t *x, Double_t *par) {
   return par[Bkg0] + par[Bkg1]*x[0] + par[Bkg2]*x[0]*x[0];
}
 
// Gauss Peak function (the signal function)
Double_t signal1(Double_t *x, Double_t *par) {
   return par[SigScale]*TMath::Gaus(x[0], par[SigMean], par[SigSigma], true);
}
 

// Sum of background and peak function. We pass x and the fit parameters
// down to the signal and background functions.
Double_t fitFunction1(Double_t *x, Double_t *par) {
   return background1(x, par) + signal1(x, par);
}
 
// Fit "fitFunction" to the histogram, and draw results on the canvas `c1`.
void FitRoutine(TCanvas* c1, TH1* histo, float fitxmin, float fitxmax, TString filename){
   c1->cd();


   // create a TF1 with the range and N_PAR parameters (six by default)
   TF1 fitFcn("fitFcn",fitFunction1,900,1250,N_PAR);
   fitFcn.SetNpx(500);
   fitFcn.SetLineWidth(2);
   fitFcn.SetLineColor(kRed);

   TF1 fitFcn1("fitFcn1",fitFunction1,1265,1620,N_PAR);
   fitFcn1.SetNpx(500);
   fitFcn1.SetLineWidth(2);
   fitFcn1.SetLineColor(kRed);
 
   // Assign the names from the map "parNames".
   for (auto& idx_name : parNames) {
     fitFcn.SetParName(idx_name.first, idx_name.second.c_str());
   }
 
 
   // Assign the names from the map "parNames".
   for (auto& idx_name : parNames) {
     fitFcn1.SetParName(idx_name.first, idx_name.second.c_str());
   }
 
   // Fit. First set ok-ish starting values for the parameters(I even dont understand it completely)
   fitFcn.SetParameters(0,0,0,0,1000,0,1000,0,2000,0);


   fitFcn1.SetParameters(0,0,0,0,1000,0,1000,0,2000,0);


   //fitFcn.SetParameters(0,0,0,0,200,2000);
   histo->GetXaxis()->SetRange(10,500);


   histo->Fit("fitFcn","VR+","ep");
   histo->Fit("fitFcn1","VR1+","ep1");
   
 
   // improve the picture:
   // Draw signal and background functions separately

   TF1 backFcn("backFcn",background,1060,1290,N_PAR);
   backFcn.SetLineColor(kBlue);
   TF1 signalFcn("signalFcn",signal,800,1750,N_PAR);
   signalFcn.SetLineColor(kBlue);
   signalFcn.SetNpx(500);




   TF1 backFcn1("backFcn1",background1,1020,1600,N_PAR);
   backFcn1.SetLineColor(kBlue);
   TF1 signalFcn1("signalFcn1",signal1,1100,1700,N_PAR);
   signalFcn1.SetLineColor(kBlue);
   signalFcn1.SetNpx(500);
 
   // Retrieve fit parameters, and copy them to the signal and background functions

   Double_t par[N_PAR];
   fitFcn.GetParameters(par);
 
   backFcn.SetParameters(par);
   //backFcn.DrawCopy("same");
 
   signalFcn.SetParameters(par);
   signalFcn.SetLineColor(kGreen);
   signalFcn.DrawCopy("same");



   //Double_t par[N_PAR];
   fitFcn1.GetParameters(par);
 
   backFcn1.SetParameters(par);
   backFcn1.DrawCopy("same");
 
   signalFcn1.SetParameters(par);
   signalFcn1.SetLineColor(kGreen);
   signalFcn1.DrawCopy("same");


/*
///Area of 1st peak
 

double w=histo->GetBinWidth(1);
int integral1 = backFcn.Integral(1055,1241)/w;
int integral2 = fitFcn.Integral(805,1255)/w;

cout<<"Integral of Background function::::::"<<integral1<<endl;
cout<<"Integral of Fit function::::::"<<integral2<<endl;
 
*/
//Area of second peak

double w=histo->GetBinWidth(1);
int integral1 = backFcn1.Integral(1070,1700)/w;
int integral2 = fitFcn1.Integral(1265,1770)/w;

cout<<"Integral of Background function::::::"<<integral1<<endl;
cout<<"Integral of Fit function::::::"<<integral2<<endl;


}
 
/*
   // draw the legend
   TLegend legend(0.65,0.7,0.90,0.85);
   legend.SetTextFont(72);
   legend.SetTextSize(0.03);
   legend.AddEntry(histo,"Exp. Data","lpe");
   legend.AddEntry(&backFcn,"Background","l");
   legend.AddEntry(&signalFcn1,"Signal","l");
   legend.AddEntry(&fitFcn1,"Signal+Bkgd","l");
   legend.DrawClone();
   histo->Draw("esame");
   c1->SaveAs(filename);



*/
 

// Fitting histogram from the root file
void fitCo60() {
   gStyle->SetOptFit(1111);
   gStyle->SetOptStat(0);
 
   // fit range from fitxmin to fitxmax
   float fitxmin=100;
   float fitxmax=3000;

    TCanvas *c1 = new TCanvas("c1","Fitting 137Cs",10,10,700,700);
    TFile* f = new TFile("60CoReconstructor.root");
   // TH1F *h = (TH1F*) f->Get("doppler");           //withoutaddback
   TH1F *h = (TH1F*) f->Get("doppler");           // withaddback
   if (!h){
      cout << "h not found"<<endl;
      return;
   }
   h->SetMarkerStyle(2);
   h->SetMarkerSize(2);
   // The output of the fit. 
   FitRoutine(c1,h, fitxmin, fitxmax,"137Cs.pdf");
}
