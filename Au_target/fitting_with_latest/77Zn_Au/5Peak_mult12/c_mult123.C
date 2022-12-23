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
    
    int nbOfPeaks =5;
    
    char temp[300];
    int minBin = 0;
    int maxBin = 4800;
    int binning = 25;
    int numBin = (maxBin-minBin)/binning;
    
    
    //****************************************************************************
    // The simulated  peaks
    std::map<int,TFile *>sim;


    sim[0]   = new TFile("77Zn_Au_887_Reconstructor_MULT12.root"); 
    sim[1]   = new TFile("77Zn_Au_1363_Reconstructor_MULT12.root");
    sim[2]   = new TFile("77Zn_Au_1409_Reconstructor_MULT12.root");
    sim[3]   = new TFile("77Zn_Au_1583_Reconstructor_MULT12.root");
    sim[4]   = new TFile("77Zn_Au_2037_Reconstructor_MULT12.root");
    //sim[5]   = new TFile("77Zn_C_2626_Reconstructor_MULT123.root");
    
    //The experimental data:
    TFile *experimental = new TFile("three_Analysis_yes_dalitrigger_90dets_off.root","READ");
    TCanvas *fCanvas=new TCanvas("Canvas","77Zn",800,600);
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

        hsim[4] = new TH1F(temp,temp,numBin,minBin,maxBin);
        hsim[4] = (TH1F*)sim[4]->Get("doppler");

       // hsim[5] = new TH1F(temp,temp,numBin,minBin,maxBin);
       // hsim[5] = (TH1F*)sim[5]->Get("doppler");

   // }
    
    //The experimental spectra:
    TH1F *hexp = new TH1F("hexp","hexp",numBin,minBin,maxBin);
    //  hexp = (TH1F*)exp->Get("h_mult1");   //mult1
    hexp = (TH1F*)experimental->Get("hEnergy_mult123_77Zn");
    //hp_doppler_22
    //hexp->Rebin(20);   ///Double check later...
    hexp->SetStats(0);
    hexp->SetFillColor(0);
    hexp->SetLineColor(1);
    hexp->SetLineWidth(2);
    
    hexp->GetXaxis()->SetNdivisions(110);
    hexp->GetYaxis()->SetNdivisions(110);
    
    hexp->GetXaxis()->SetRangeUser(450,3500);
    hexp->GetYaxis()->SetRangeUser(0,600);
    
    hexp->GetYaxis()->SetTitle("Counts (25 keV/bin");
    hexp->GetXaxis()->SetTitle("Energy (keV)");
    
    hexp->GetXaxis()->SetTitleOffset(0.65);
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
    peak1g  = new TGraph(hsim[1]);
    peak2g  = new TGraph(hsim[2]);
    peak3g  = new TGraph(hsim[3]);
    peak4g  = new TGraph(hsim[4]);
    //peak5g  = new TGraph(hsim[5]);

    Double_t fitmin=00.;
    Double_t fitmax=7000.;
    
    TF1 *whole = new TF1( "whole",ex_respf,fitmin,fitmax,2*nbOfPeaks+4); // *2 because 2 parameters for each peak, +4 for the exponential
    
    whole->SetParName(4,"887 keV");
    whole->SetParName(6,"1363 keV");
    whole->SetParName(8,"1409 keV");
    whole->SetParName(10,"1583 keV");
    whole->SetParName(12,"2037 keV");
   // whole->SetParName(14,"2626 keV");
   
    
    

    Double_t norma1 = 0.18;
    Double_t value887 = norma1*0.018;
   
    Double_t norma2 = 0.12;
    Double_t value1363 = norma2*0.015;



    Double_t norma3 = 0.1;
    Double_t value1409 = norma3*0.018;




    Double_t norma4 = 0.16;
    Double_t value1583 = norma4*0.026;


    Double_t norma5 = 0.18;
    Double_t value2037 = norma5*0.024;


   // Double_t norma6 = 0.14;
   // Double_t value2626 = norma6*0.018;
    
    
    
    
    whole->FixParameter( 4, value887); //    691 peak 1
    whole->FixParameter( 5, 0.0);
    
    whole->FixParameter( 6, value1363); //    887 peak 2
    whole->FixParameter( 7, 0.0);

    whole->FixParameter( 8, value1409); //   1170 peak 3
    whole->FixParameter( 9, 0.0);


    whole->FixParameter( 10, value1583); //  1409 peak 4
    whole->FixParameter( 11, 0.0);



    whole->FixParameter( 12, value2037); //  1936 peak 5
    whole->FixParameter( 13, 0.0);
    


    //whole->FixParameter( 14, value2626); //  2626 peak 6
    //whole->FixParameter( 15, 3.0);
    
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



    TF1 *peak4f = new TF1("peak4f",resp4,fitmin,fitmax,2); //
    peak4f->SetParameter(0,whole->GetParameter(12));
    peak4f->SetParameter(1,whole->GetParameter(13));
    peak4f->SetLineColor(kOrange);
    peak4f->SetLineWidth(2);
    //peak4f->SetLineStyle(12);
    peak4f->SetNpx(500);
    peak4f->Draw("same");


   /* TF1 *peak5f = new TF1("peak5f",resp5,fitmin,fitmax,2); //
    peak5f->SetParameter(0,whole->GetParameter(14));
    peak5f->SetParameter(1,whole->GetParameter(15));
    peak5f->SetLineColor(kCyan);
    peak5f->SetLineWidth(2);
    //peak5f->SetLineStyle(12);
    peak5f->SetNpx(500);
    peak5f->Draw("same");

*/



    
    hexp->Draw("same");
    //   Double_t param[30];
    //   whole->GetParameters(param);
   // Double_t chi2 = whole->GetChisquare(); ///whole->GetNDF();
    Double_t chi2ndf = whole->GetChisquare()/whole->GetNDF();
    // /*
    TLegend *leg = new TLegend(0.5,0.7,0.9,0.9);
    leg->AddEntry(hexp,"Experiment_Au_Mult12_77Zn");
    leg->AddEntry(expon,"Double exponential BG");
    leg->AddEntry(peak0f,"GEANT4 Simulation 887 keV");
    leg->AddEntry(peak1f,"GEANT4 Simulation 1363 keV");
    leg->AddEntry(peak2f,"GEANT4 Simulation 1409 keV");
    leg->AddEntry(peak3f,"GEANT4 Simulation 1583 keV");
    leg->AddEntry(peak4f,"GEANT4 Simulation 2037 keV");
    //leg->AddEntry(peak5f,"GEANT4 Simulation 2626 keV");
   
    leg->AddEntry(whole," Peak + expo BG");
    leg->AddEntry((TObject*)0, "", "");
    //leg->AddEntry(whole, Form("Chi2 = %4.2f", chi2), "l");
   // leg->AddEntry(whole, Form("Chi2/NDF = %4.2f", chi2ndf), "l");
    
    leg->SetFillColor(0);
    leg->Draw("Same");
    
    TLegend *leg2 = new TLegend(0.5,0.5,0.9,0.7);
    Char_t message[60];
    leg2->AddEntry((TObject*)0,"APMLITUDES:","");
    //leg2->AddEntry((TObject*)1,"APMLITUDES:","");
    for(Int_t i=0;i<13;i++)
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

    fCanvas->Print("77Zn_mult123_Au_no_AddBack.pdf");
}



