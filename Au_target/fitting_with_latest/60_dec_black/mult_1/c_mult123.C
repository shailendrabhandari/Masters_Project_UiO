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



void c_mult123(){
    gStyle->SetOptStat(kFALSE);
    gStyle->SetOptStat(kFALSE);
    gStyle->SetPadGridX(false);
    gStyle->SetPadGridY(false);
    
    int nbOfPeaks =4;
    
    char temp[300];
    int minBin = 0;
    int maxBin = 4800;
    int binning = 25;
    int numBin = (maxBin-minBin)/binning;
    
    
    //****************************************************************************
    // The simulated  peaks
    std::map<int,TFile *>sim;

    sim[0]   = new TFile("730_78Zn_Au_Reconstructor_mult1.root"); 
    sim[1]   = new TFile("889_78Zn_Au_Reconstructor_mult1.root");
    sim[2]   = new TFile("78Zn_Reconstructor_Au908_mult1.root");
    sim[3]   = new TFile("78Zn_Reconstructor_Au1215_mult1.root");
    
    //The experimental data:
    TFile *experimental = new TFile("three_Analysis_yes_dalitrigger_90dets_off.root","READ");
    TCanvas *fCanvas=new TCanvas("Canvas","78Zn",800,600);
    fCanvas->cd();
    
    
    //****************************************************************************
    //The  spectra:
    //The hpeak are only a projection of the two-dimensional spectra
    std::map<int,TH1F *>hsim;
    
    //for(int i=0;i<nbOfPeaks;i++)
    //{
        //sprintf(temp,"hsim%i",i);
        hsim[0] = new TH1F(temp,temp,numBin,minBin,maxBin);
        hsim[0] = (TH1F*)sim[0]->Get("doppler");

        hsim[1] = new TH1F(temp,temp,numBin,minBin,maxBin);
        hsim[1] = (TH1F*)sim[1]->Get("doppler");

        hsim[2] = new TH1F(temp,temp,numBin,minBin,maxBin);
        hsim[2] = (TH1F*)sim[2]->Get("doppler");

        hsim[3] = new TH1F(temp,temp,numBin,minBin,maxBin);
        hsim[3] = (TH1F*)sim[3]->Get("doppler");

   // }
    
    //The experimental spectra:
    TH1F *hexp = new TH1F("hexp","hexp",numBin,minBin,maxBin);
    //  hexp = (TH1F*)exp->Get("h_mult1");   //mult1
    hexp = (TH1F*)experimental->Get("hEnergy_mult1_78Zn");
    //hp_doppler_22
    //hexp->Rebin(20);   ///Double check later...
    hexp->SetStats(0);
    hexp->SetFillColor(0);
    hexp->SetLineColor(1);
    hexp->SetLineWidth(2);
    
    hexp->GetXaxis()->SetNdivisions(110);
    hexp->GetYaxis()->SetNdivisions(110);
    
    hexp->GetXaxis()->SetRangeUser(460,5000);
    hexp->GetYaxis()->SetRangeUser(0,8000);
    
    hexp->GetYaxis()->SetTitle("Counts (25 keV/bin)");
    hexp->GetXaxis()->SetTitle("Energy (keV)");
    
    hexp->GetXaxis()->SetTitleOffset(0.71);
    hexp->GetXaxis()->CenterTitle();
    hexp->GetYaxis()->SetTitleOffset(0.80);
    hexp->GetYaxis()->CenterTitle();
    
    hexp->GetXaxis()->SetTitleFont(62);
    hexp->GetYaxis()->SetTitleFont(62);
    
    hexp->GetXaxis()->SetTitleSize(0.06);
    hexp->GetYaxis()->SetTitleSize(0.06);
    
    hexp->GetXaxis()->SetLabelSize(0.03);
    hexp->GetYaxis()->SetLabelSize(0.03);
    hexp->SetTitle(0);
    
    peak0g  = new TGraph(hsim[0]);
    peak1g  = new TGraph(hsim[1]);
    peak2g  = new TGraph(hsim[2]);
    peak3g  = new TGraph(hsim[3]);
    Double_t fitmin=250.;
    Double_t fitmax=5000.;
    
    TF1 *whole = new TF1( "whole",ex_respf,fitmin,fitmax,2*nbOfPeaks+4); // *2 because 2 parameters for each peak, +4 for the exponential
    
    whole->SetParName(4,"730 keV");
    whole->SetParName(6,"889 keV");
    whole->SetParName(8,"908 keV");
    whole->SetParName(10,"1215 keV");
   
    
    

    Double_t norma1 = 0.39;
    Double_t value730 = norma1*0.24;
   
   /* Double_t norma2 = 0.15;
    Double_t value889 = norma2*0.10;

    Double_t norma3 = 0.15;
    Double_t value908 = norma3*0.10;




    Double_t norma4 = 0.18;
    Double_t value1215 = norma4*0.12;
    
    */
    
    
    whole->FixParameter( 4, value730); //  730 peak 1
    whole->FixParameter(11, 18.0);
    
 /*   whole->FixParameter( 6, value889); //  889 peak 2
    whole->FixParameter( 7, 14.0);

    whole->FixParameter( 8, value908); //  908 peak 3
    whole->FixParameter( 9, 18.0);


    whole->FixParameter( 10, value1215); //  1215 peak 4
    whole->FixParameter( 11, 22.0);
    */
    
    whole->SetLineColor(2);
    whole->SetLineWidth(2);
    whole->SetNpx(500);
    hexp->Fit(whole,"R");
    //   (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hexp);
    
    cout << endl << "Chi2: " << whole->GetChisquare() << endl;
    
    TF1 *expon= new TF1( "expon",expf,fitmin,fitmax,4); // exp
    expon->SetParameters(whole->GetParameter(0),whole->GetParameter(1),
                         whole->GetParameter(2),whole->GetParameter(3));
    expon->SetLineColor(4);
    expon->SetLineWidth(3);
    expon->SetLineStyle(7);
    expon->Draw("same");
    
    TF1 *peak0f = new TF1("peak0f",resp0,fitmin,fitmax,2); //
    peak0f->SetParameter(0,whole->GetParameter(4));
    peak0f->SetParameter(1,whole->GetParameter(5));
    peak0f->SetLineColor(kMagenta-2);
    peak0f->SetLineWidth(2);
    //peak0f->SetLineStyle(9);
    peak0f->SetNpx(500);
    peak0f->Draw("same");
    //    (peak0f->GetHistogram())->SaveAs("75Cu_resp_1260.C");
    



    TF1 *peak1f = new TF1("peak1f",resp1,fitmin,fitmax,2); //
    peak1f->SetParameter(0,whole->GetParameter(6));
    peak1f->SetParameter(1,whole->GetParameter(7));
    peak1f->SetLineColor(kGreen);
    peak1f->SetLineWidth(2);
    //peak1f->SetLineStyle(10);
    peak1f->SetNpx(500);
    peak1f->Draw("same");

    TF1 *peak2f = new TF1("peak2f",resp2,fitmin,fitmax,2); //
    peak2f->SetParameter(0,whole->GetParameter(8));
    peak2f->SetParameter(1,whole->GetParameter(9));
    peak2f->SetLineColor(kPlum);
    peak2f->SetLineWidth(2);
    //peak2f->SetLineStyle(11);
    peak2f->SetNpx(500);
    peak2f->Draw("same");




    TF1 *peak3f = new TF1("peak3f",resp3,fitmin,fitmax,2); //
    peak3f->SetParameter(0,whole->GetParameter(10));
    peak3f->SetParameter(1,whole->GetParameter(11));
    peak3f->SetLineColor(kCyan);
    peak3f->SetLineWidth(2);
    //peak3f->SetLineStyle(12);
    peak3f->SetNpx(500);
    peak3f->Draw("same");



    
    hexp->Draw("same");
    //   Double_t param[30];
    //   whole->GetParameters(param);
   // Double_t chi2 = whole->GetChisquare(); ///whole->GetNDF();
    Double_t chi2ndf = whole->GetChisquare()/whole->GetNDF();
    // /*
    TLegend *leg = new TLegend(0.5,0.7,0.9,0.9);
    leg->AddEntry(hexp,"Experiment_Au_Mult1 78Zn");
    leg->AddEntry(expon,"Double exponential BG");
    leg->AddEntry(peak0f,"GEANT4 Simulation 730 keV");
    leg->AddEntry(peak1f,"GEANT4 Simulation 889 keV");
    leg->AddEntry(peak2f,"GEANT4 Simulation 908 keV");
    leg->AddEntry(peak3f,"GEANT4 Simulation 1215 keV");
   
    leg->AddEntry(whole," Peaks + expo BG");
    leg->AddEntry((TObject*)0, "", "");
    //leg->AddEntry(whole, Form("Chi2 = %4.2f", chi2), "l");
   // leg->AddEntry(whole, Form("Chi2/NDF = %4.2f", chi2ndf), "l");
    
    leg->SetFillColor(0);
    leg->Draw("Same");
    
    TLegend *leg2 = new TLegend(0.5,0.5,0.9,0.7);
    Char_t message[60];
    leg2->AddEntry((TObject*)0,"APMLITUDES:","");
    //leg2->AddEntry((TObject*)1,"APMLITUDES:","");
    for(Int_t i=0;i<11;i++)
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

    fCanvas->Print("78Zn_Au_mult1_no_AddBack.pdf");
}



