{

    TFile *f1 = new TFile("/media/shailendra/New Volume/Data_UIO/carbon_analysis/Addback/gamma_Carbon_target_run.root");
    TH1F *h1 = (TH1F*) f1->Get("hEnergy_mult123_78Zn");   

    TFile *f2 = new TFile("/media/shailendra/New Volume/Data_UIO/carbon_analysis/Addback/gamma_Carbon_target_run.root");
    TH1F *h2 = (TH1F*) f2->Get("hEnergy_mult123_wa_78Zn");  

   
   h1->SetLineColor(kBlue);
   h1->SetLineWidth(2);
   h2->SetLineColor(kRed);
   h2->SetLineWidth(2);

   h1->SetTitle("Experimental spectrum of 78Zn_C_target without addback");
   h2->SetTitle("Experimental spectrum of 78Zn_C_target with addback");   

   THStack *st = new THStack();
   st->Add(h1);
   st->Add(h2);



   st->Draw("nostack");
   
   st->GetXaxis()->SetTitle("Energy(keV)");
   st->GetYaxis()->SetTitle("Counts");
   st->GetXaxis()->SetTitleOffset(0.78);
   st->GetXaxis()->CenterTitle();
   st->GetYaxis()->SetTitleOffset(0.96);
   st->GetYaxis()->CenterTitle();
    
   st->GetXaxis()->SetTitleFont(62);
   st->GetYaxis()->SetTitleFont(62);
    
   st->GetXaxis()->SetTitleSize(0.05);
   st->GetYaxis()->SetTitleSize(0.05);
    
   st->GetXaxis()->SetLabelSize(0.03);
   st->GetYaxis()->SetLabelSize(0.03);
   st->SetTitle(0);
   
   gPad->BuildLegend();
}
