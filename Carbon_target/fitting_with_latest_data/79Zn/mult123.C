#include <TLegend.h>
#include <TF1.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TFile.h>
#include <map>
#include <TGraph.h>
#include <TStyle.h>
#include <TROOT.h>
#include "fit_spectrum.h"



void mult123(){
    gStyle->SetOptStat(kFALSE);
    gStyle->SetOptStat(kFALSE);
    gStyle->SetPadGridX(false);
    gStyle->SetPadGridY(false);
    
    int nbOfPeaks = 2;
    
    char temp[300];
    int minBin = 0;
    int maxBin = 4800;
    int binning = 20;
    int numBin = (maxBin-minBin)/binning;
    
    
    //****************************************************************************
    // The simulated  peaks
    std::map<int,TFile *>sim;

    //sim[0]   = new TFile("79Zn_GOLD_reconstructor_FFIND1.root");
    sim[0]   = new TFile("1264_79Zn_C_Reconstructor_mult123.root");
    
    //The experimental data:
    TFile *experimental = new TFile("three_Analysis_yes_dalitrigger_60dets_off.root","READ");
    TCanvas *fCanvas=new TCanvas("Canvas","79Zn",800,600);
    fCanvas->cd();
    
    
    //****************************************************************************
    //The  spectra:
    //The hpeak are only a projection of the two-dimensional spectra
    std::map<int,TH1F *>hsim;
    
    //for(int i=0;i<nbOfPeaks;i++)
    //{
       // sprintf(temp,"hsim%i",i);
        hsim[0] = new TH1F(temp,temp,numBin,minBin,maxBin);
        hsim[0] = (TH1F*)sim[0]->Get("doppler");
    //}
    
    //The experimental spectra:
    TH1F *hexp = new TH1F("hexp","hexp",numBin,minBin,maxBin);
    //  hexp = (TH1F*)exp->Get("h_mult1");   //mult1
    hexp = (TH1F*)experimental->Get("hEnergy_mult123_79Zn");
    //hp_doppler_22
   //hexp->Rebin(20);
    hexp->SetStats(0);
    hexp->SetFillColor(0);
    hexp->SetLineColor(1);
    hexp->SetLineWidth(2);
    
    hexp->GetXaxis()->SetNdivisions(110);
    hexp->GetYaxis()->SetNdivisions(110);
    
    hexp->GetXaxis()->SetRangeUser(300,5000);
    hexp->GetYaxis()->SetRangeUser(1,350);
    
    hexp->GetYaxis()->SetTitle("Counts"); //(16 keV/bin)");
    hexp->GetXaxis()->SetTitle("Energy (keV)");
    
    hexp->GetXaxis()->SetTitleOffset(0.70);
    hexp->GetXaxis()->CenterTitle();
    hexp->GetYaxis()->SetTitleOffset(0.70);
    hexp->GetYaxis()->CenterTitle();
    
    hexp->GetXaxis()->SetTitleFont(62);
    hexp->GetYaxis()->SetTitleFont(62);
    
    hexp->GetXaxis()->SetTitleSize(0.06);
    hexp->GetYaxis()->SetTitleSize(0.06);
    
    hexp->GetXaxis()->SetLabelSize(0.04);
    hexp->GetYaxis()->SetLabelSize(0.04);
    hexp->SetTitle(0);
    
    peak0g  = new TGraph(hsim[0]);
    ;
    Double_t fitmin=10.;
    Double_t fitmax=5000.;
    
    TF1 *whole = new TF1( "whole",ex_respf,fitmin,fitmax,2*nbOfPeaks+4); // *2 because 2 parameters for each peak, +4 for the exponential
    
    whole->SetParName(4,"1260 keV");
   
    
    
    // /* //01.11.2018
    Double_t norma1 = 0.08
;
    Double_t value1260 = norma1*0.040;
    /*Double_t value530 = norma1*0.49;
    Double_t value1417 = norma1*0.38;
    Double_t value1483 = norma1*1;
    
    Double_t norma2 = 0.015;
    Double_t value883 = norma2*1;
    Double_t value950 = norma2*0.36;
    Double_t value990 = value1260;*/
    // */
    
    //   Double_t p1260 = 4.5e-03;
    //   Double_t p990 = p1260;
    
    whole->FixParameter( 4, value1260); //  1260 peak 1
    whole->FixParameter( 7, 12.0);
    
   
    
    
    whole->SetLineColor(2);
    whole->SetLineWidth(2);
    whole->SetNpx(1000);
    hexp->Fit(whole,"R");
    //   (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hexp);
    
    cout << endl << "Chi2: " << whole->GetChisquare() << endl;
    
    TF1 *expon= new TF1( "expon",expf,fitmin,fitmax,4); // exp
    expon->SetParameters(whole->GetParameter(0),whole->GetParameter(1),
                         whole->GetParameter(2),whole->GetParameter(3));
    expon->SetLineColor(4);
    expon->SetLineWidth(2);
    expon->SetLineStyle(7);
    expon->Draw("same");
    
    TF1 *peak0f = new TF1("peak0f",resp0,fitmin,fitmax,2); //
    peak0f->SetParameter(0,whole->GetParameter(4));
    peak0f->SetParameter(1,whole->GetParameter(5));
    peak0f->SetLineColor(kMagenta-2);
    peak0f->SetLineWidth(2);
    //peak0f->SetLineStyle(9);
    peak0f->SetNpx(200);
    peak0f->Draw("same");
    //    (peak0f->GetHistogram())->SaveAs("75Cu_resp_1260.C");
    
   
    
    hexp->Draw("same");
    //    Double_t param[30];
    //   whole->GetParameters(param);
    Double_t chi2 = whole->GetChisquare(); ///whole->GetNDF();
    Double_t chi2ndf = whole->GetChisquare()/whole->GetNDF();
    // /*
    TLegend *leg = new TLegend(0.5,0.7,0.9,0.9);
    leg->AddEntry(hexp,"Experiment_C_Mult123_79Zn");
    leg->AddEntry(expon,"Double exponential BG");
    leg->AddEntry(peak0f,"GEANT4 Simulation 1260 keV");
    leg->AddEntry(whole," Peak + expo BG");

    leg->AddEntry((TObject*)0, "", "");
    //leg->AddEntry(whole, Form("Chi2 = %4.2f", chi2), "l");
    leg->AddEntry(whole, Form("Chi2/NDF = %4.2f", chi2ndf), "l");
    
    leg->SetFillColor(0);
    leg->Draw("Same");
    
    TLegend *leg2 = new TLegend(0.5,0.5,0.9,0.7);
    Char_t message[80];
    leg2->AddEntry((TObject*)0,"APMLITUDES:","");
    for(Int_t i=0;i<5;i++)
    {
        if (i % 2) i++;
        sprintf(message,"%s = %.5f",
                whole->GetParName(i), whole->GetParameter(i)
                );
        leg2->AddEntry((TObject*)0,message,"");
        //leg2->AddEntry((TObject*)0,whole->Chi2()," ");
    }
    // leg2->AddEntry((TObject*)0,"Chi2",whole->GetChisquare());
    leg2->SetFillColor(0);
    leg2->Draw("Same");
    //  */

    fCanvas->Print("79Zn_mult123_C_target.pdf");
}



