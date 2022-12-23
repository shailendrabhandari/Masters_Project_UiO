 {


TFile *f = new TFile("three_Analysis_yes_dalitrigger.root");
TH2F *ggmatrix;
TH1F *p_ggmatrix;
TH1F* gamma_single_mult1;
TH1F* gamma_single_mult12;
TH1F* gamma_single_mult123;
TH1F* gamma_single_mult1234;

     

     
ggmatrix = (TH2F*)f->Get("gg_mult234_78Zn");
p_ggmatrix  = (TH1F*)ggmatrix->ProjectionY();


//gg_mult23_p_4096
//gg_multAll_p_4096
gamma_single_mult1 = (TH1F*)f->Get("hEnergy_mult1_wa_78Zn");
gamma_single_mult12 = (TH1F*)f->Get("hEnergy_mult12_wa_78Zn");
gamma_single_mult123 = (TH1F*)f->Get("hEnergy_mult123_wa_78Zn");
//gamma_single_mult1234 = (TH1F*)f->Get("hEnergy_mult1234_wa_78Zn");


     
hgammaprompt_730 = (TH1F*)ggmatrix->ProjectionY("gate on 730 keV",680,740);
hgammaleft_730 = (TH1F*)ggmatrix->ProjectionY("BG left",480,600);
hgammaright_730 = (TH1F*)ggmatrix->ProjectionY("BG right",760,780);

hgammaprompt_830 = (TH1F*)ggmatrix->ProjectionY("gate on 889 keV",800,885);
hgammaleft_830 = (TH1F*)ggmatrix->ProjectionY("BG left",530,590);
hgammaright_830 = (TH1F*)ggmatrix->ProjectionY("BG right",930,950);
     
     hgammaprompt_1230 = (TH1F*)ggmatrix->ProjectionY("gate on 1215 keV",1160,1270);
     hgammaleft_1230 = (TH1F*)ggmatrix->ProjectionY("BG left",1030,1060);
     hgammaright_1230 = (TH1F*)ggmatrix->ProjectionY("BG right",1330,1390);

TH1F* clean_730 = new TH1F(*hgammaprompt_730);
TH1F* bg_total_730 = new TH1F(*hgammaleft_730);
bg_total_730->Add(hgammaright_730,+1.);
clean_730->Add(bg_total_730, -1./2.);
//clean_730->Add(bg_total_730, -1./(p_ggmatrix->Integral(645,760)+p_ggmatrix->Integral(560,590))/p_ggmatrix->Integral(810,840));
 (TH1F*)clean_730->Clone("gate  with BG subtraction  730_keV");

TH1F* clean_830 = new TH1F(*hgammaprompt_830);
TH1F* bg_total_830 = new TH1F(*hgammaleft_830);
bg_total_830->Add(hgammaright_830,+1.);
clean_830->Add(bg_total_830, -1./2.);
(TH1F*)clean_830->Clone("gate  with BG subtraction  889_keV");

TH1F* clean_1230 = new TH1F(*hgammaprompt_1230);
TH1F* bg_total_1230 = new TH1F(*hgammaleft_1230);
bg_total_1230->Add(hgammaright_1230,+1.);
clean_1230->Add(bg_total_1230, -1./2.);
(TH1F*)clean_1230->Clone("gate  with BG subtraction  1215_keV");


/*TH1F* clean_880 = new TH1F(*hgammaprompt_880);
TH1F* bg_total_880 = new TH1F(*hgammaleft_880);
bg_total_880->Add(hgammaright_880,+1.);
clean_880->Add(bg_total_880, -1./2.);
 (TH1F*)clean_880->Clone("gate  with BG subtraction  880_keV");*/


gStyle->SetOptStat(0);
gStyle->SetTitleFontSize(0.07);
/*
TCanvas *c1 = new TCanvas("c1","matrix",700,900);
c1->SetTopMargin(10.0);
c1->SetBottomMargin(10.0);
c1->SetLeftMargin(10.0);
c1->SetRightMargin(10.0);
c1->ToggleEventStatus();
 c1->Divide(1,3);
c1->cd(1);
gamma_single_mult1->SetTitle("Gamma singles mult=1");
gamma_single_mult1->GetXaxis()->SetTitle("Energy (keV)");
gamma_single_mult1->GetXaxis()->SetTitleFont(42);
gamma_single_mult1->GetXaxis()->SetTitleSize(0.06);
gamma_single_mult1->GetXaxis()->SetTitleOffset(0.7);
gamma_single_mult1->GetXaxis()->SetLabelFont(42);
gamma_single_mult1->GetXaxis()->SetLabelSize();
///////////////////////////////////////
///////////////////////////////////////
gamma_single_mult1->GetYaxis()->SetTitle("Counts (16 keV/bin)");
gamma_single_mult1->GetYaxis()->SetTitleFont(42);
gamma_single_mult1->GetYaxis()->SetTitleSize(0.05);
gamma_single_mult1->GetYaxis()->SetTitleOffset(0.6);
gamma_single_mult1->GetYaxis()->SetLabelFont(42);
gamma_single_mult1->GetYaxis()->SetLabelSize();
gamma_single_mult1->SetLineColor(2);/////rengini degistirdim
gamma_single_mult1->SetLineWidth(3);///line thickness
//gamma_single_mult1->Rebin(16);
gamma_single_mult1->GetXaxis()->SetRange(15,120);
gamma_single_mult1->Draw();
c1->cd(2);
gamma_single_mult12->SetTitle("Gamma singles mult=12");
gamma_single_mult12->GetXaxis()->SetTitle("Energy (keV)");
gamma_single_mult12->GetXaxis()->SetTitleFont(42);
gamma_single_mult12->GetXaxis()->SetTitleSize(0.06);
gamma_single_mult12->GetXaxis()->SetTitleOffset(0.7);
gamma_single_mult12->GetXaxis()->SetLabelFont(42);
gamma_single_mult12->GetXaxis()->SetLabelSize();
///////////////////////////////////////
///////////////////////////////////////
gamma_single_mult12->GetYaxis()->SetTitle("Counts (16 keV/bin)");
gamma_single_mult12->GetYaxis()->SetTitleFont(42);
gamma_single_mult12->GetYaxis()->SetTitleSize(0.05);
gamma_single_mult12->GetYaxis()->SetTitleOffset(0.6);
gamma_single_mult12->GetYaxis()->SetLabelFont(42);
gamma_single_mult12->GetYaxis()->SetLabelSize();
gamma_single_mult12->SetLineColor(3);/////rengini degistirdim
gamma_single_mult12->SetLineWidth(3);///line thickness
//gamma_single_mult12->Rebin(16);
gamma_single_mult12->GetXaxis()->SetRange(15,120);
gamma_single_mult12->Draw();

c1->cd(3);
gamma_single_mult123->SetTitle("Gamma singles mult=123");
gamma_single_mult123->GetXaxis()->SetTitle("Energy (keV)");
gamma_single_mult123->GetXaxis()->SetTitleFont(42);
gamma_single_mult123->GetXaxis()->SetTitleSize(0.06);
gamma_single_mult123->GetXaxis()->SetTitleOffset(0.7);
gamma_single_mult123->GetXaxis()->SetLabelFont(42);
gamma_single_mult123->GetXaxis()->SetLabelSize();
///////////////////////////////////////
///////////////////////////////////////
gamma_single_mult123->GetYaxis()->SetTitle("Counts (16 keV/bin)");
gamma_single_mult123->GetYaxis()->SetTitleFont(42);
gamma_single_mult123->GetYaxis()->SetTitleSize(0.05);
gamma_single_mult123->GetYaxis()->SetTitleOffset(0.6);
gamma_single_mult123->GetYaxis()->SetLabelFont(42);
gamma_single_mult123->GetYaxis()->SetLabelSize();
gamma_single_mult123->SetLineColor(4);/////rengini degistirdim
gamma_single_mult123->SetLineWidth(3);///line thickness
//gamma_single_mult123->Rebin(16);
gamma_single_mult123->GetXaxis()->SetRange(15,120);
//gamma_single_mult123->SetMinimum(-100);

gamma_single_mult123->Draw();

     
c1->cd(4);
gamma_single_mult1234->SetTitle("Gamma singles mult=1234");
gamma_single_mult1234->GetXaxis()->SetTitle("Energy (keV)");
gamma_single_mult1234->GetXaxis()->SetTitleFont(42);
gamma_single_mult1234->GetXaxis()->SetTitleSize(0.06);
gamma_single_mult1234->GetXaxis()->SetTitleOffset(0.7);
gamma_single_mult1234->GetXaxis()->SetLabelFont(42);
gamma_single_mult1234->GetXaxis()->SetLabelSize();
///////////////////////////////////////
///////////////////////////////////////
gamma_single_mult1234->GetYaxis()->SetTitle("Counts (16 keV/bin)");
gamma_single_mult1234->GetYaxis()->SetTitleFont(42);
gamma_single_mult1234->GetYaxis()->SetTitleSize(0.05);
gamma_single_mult1234->GetYaxis()->SetTitleOffset(0.6);
gamma_single_mult1234->GetYaxis()->SetLabelFont(42);
gamma_single_mult1234->GetYaxis()->SetLabelSize();
gamma_single_mult1234->SetLineColor(5);/////rengini degistirdim
gamma_single_mult1234->SetLineWidth(3);///line thickness
//gamma_single_mult1234->Rebin(16);
gamma_single_mult1234->Draw();*/
   
     
     /*
     
      TCanvas *c2 = new TCanvas("c2","Projection",700,800);
     p_ggmatrix->SetTitle("Projection");
     p_ggmatrix->GetXaxis()->SetTitle("Energy (keV)");
     p_ggmatrix->GetXaxis()->SetTitleFont(42);
     p_ggmatrix->GetXaxis()->SetTitleSize(0.06);
     p_ggmatrix->GetXaxis()->SetTitleOffset(0.7);
     p_ggmatrix->GetXaxis()->SetLabelFont(42);
     p_ggmatrix->GetXaxis()->SetLabelSize();
     ///////////////////////////////////////
     ///////////////////////////////////////
     p_ggmatrix->GetYaxis()->SetTitle("Counts (16 keV/bin)");
     p_ggmatrix->GetYaxis()->SetTitleFont(42);
     p_ggmatrix->GetYaxis()->SetTitleSize(0.05);
     p_ggmatrix->GetYaxis()->SetTitleOffset(0.6);
     p_ggmatrix->GetYaxis()->SetLabelFont(42);
     p_ggmatrix->GetYaxis()->SetLabelSize();
     p_ggmatrix->SetLineColor(2);/////rengini degistirdim
     p_ggmatrix->SetLineWidth(3);///line thickness
     p_ggmatrix->Rebin(16);
     p_ggmatrix->GetXaxis()->SetRange(15,200);
     p_ggmatrix->Draw();

     
     TCanvas *c3 = new TCanvas("c3","gated spectra",700,800);
     c3->ToggleEventStatus();
     c3->SetTitle("with and without BG subtraction");
     c3->Divide(1,2);
     c3->cd(1);
     hgammaprompt_730->SetTitle("gate on 730 keV without BG subtraction");
     hgammaprompt_730->SetTitleSize(1);
     hgammaprompt_730->GetXaxis()->SetTitle("Energy (keV)");
     hgammaprompt_730->GetXaxis()->SetTitleFont(42);
     hgammaprompt_730->GetXaxis()->SetTitleSize(0.06);
     hgammaprompt_730->GetXaxis()->SetTitleOffset(0.7);
     hgammaprompt_730->GetXaxis()->SetLabelFont(42);
     hgammaprompt_730->GetXaxis()->SetLabelSize();
     ///////////////////////////////////////
     ///////////////////////////////////////
     hgammaprompt_730->GetYaxis()->SetTitle("Counts (16 keV/bin)");
     hgammaprompt_730->GetYaxis()->SetTitleFont(42);
     hgammaprompt_730->GetYaxis()->SetTitleSize(0.06);
     hgammaprompt_730->GetYaxis()->SetTitleOffset(0.6);
     hgammaprompt_730->GetYaxis()->SetLabelFont(42);
     hgammaprompt_730->GetYaxis()->SetLabelSize();
     hgammaprompt_730->SetLineColor(1);/////rengini degistirdim
     hgammaprompt_730->SetLineWidth(3);///line thickness
     hgammaprompt_730->Rebin(16);
     hgammaprompt_730->GetXaxis()->SetRange(10,140);
     hgammaprompt_730->SetMinimum(-1);
     hgammaprompt_730->Draw();
     cout << "done 1/1 " << endl;
     c3->cd(2);
     clean_730->SetTitle("gate on 730 keV with BG subtraction");
     clean_730->SetTitleSize(1);
     clean_730->GetXaxis()->SetTitle("Energy (keV)");
     clean_730->GetXaxis()->SetTitleFont(42);
     clean_730->GetXaxis()->SetTitleSize(0.06);
     clean_730->GetXaxis()->SetTitleOffset(0.7);
     clean_730->GetXaxis()->SetLabelFont(42);
     clean_730->GetXaxis()->SetLabelSize();
     ///////////////////////////////////////
     ///////////////////////////////////////
     clean_730->GetYaxis()->SetTitle("Counts (16 keV/bin)");
     clean_730->GetYaxis()->SetTitleFont(42);
     clean_730->GetYaxis()->SetTitleSize(0.06);
     clean_730->GetYaxis()->SetTitleOffset(0.6);
     clean_730->GetYaxis()->SetLabelFont(42);
     clean_730->GetYaxis()->SetLabelSize();
     clean_730->SetLineColor(1);/////rengini degistirdim
     clean_730->SetLineWidth(3);///line thickness
     clean_730->Rebin(16);
     clean_730->GetXaxis()->SetRange(10,140);
     clean_730->SetMinimum(-1);
     clean_730->Draw();
     cout << "done 2/2 " << endl;
     
     
     TCanvas *c4 = new TCanvas("c4","gated spectra_2",700,800);
     c4->ToggleEventStatus();
     c4->SetTitle("without and without BG subtraction");
     c4->Divide(1,2);
     c4->cd(1);
     hgammaprompt_830->SetTitle("gate on 830 keV without BG subtraction");
     hgammaprompt_830->SetTitleSize(1);
     hgammaprompt_830->GetXaxis()->SetTitle("Energy (keV)");
     hgammaprompt_830->GetXaxis()->SetTitleFont(42);
     hgammaprompt_830->GetXaxis()->SetTitleSize(0.06);
     hgammaprompt_830->GetXaxis()->SetTitleOffset(0.7);
     hgammaprompt_830->GetXaxis()->SetLabelFont(42);
     hgammaprompt_830->GetXaxis()->SetLabelSize();
     ///////////////////////////////////////
     ///////////////////////////////////////
     hgammaprompt_830->GetYaxis()->SetTitle("Counts (16 keV/bin)");
     hgammaprompt_830->GetYaxis()->SetTitleFont(42);
     hgammaprompt_830->GetYaxis()->SetTitleSize(0.06);
     hgammaprompt_830->GetYaxis()->SetTitleOffset(0.6);
     hgammaprompt_830->GetYaxis()->SetLabelFont(42);
     hgammaprompt_830->GetYaxis()->SetLabelSize();
     hgammaprompt_830->SetLineColor(1);/////rengini degistirdim
     hgammaprompt_830->SetLineWidth(3);///line thickness
     hgammaprompt_830->Rebin(16);
     hgammaprompt_830->GetXaxis()->SetRange(10,140);
     hgammaprompt_830->SetMinimum(-1);
     hgammaprompt_830->Draw();
     cout << "done 1/1 " << endl;
     c4->cd(2);
     clean_830->SetTitle("gate on 830 keV with BG subtraction");
     clean_830->SetTitleSize(1);
     clean_830->GetXaxis()->SetTitle("Energy (keV)");
     clean_830->GetXaxis()->SetTitleFont(42);
     clean_830->GetXaxis()->SetTitleSize(0.06);
     clean_830->GetXaxis()->SetTitleOffset(0.7);
     clean_830->GetXaxis()->SetLabelFont(42);
     clean_830->GetXaxis()->SetLabelSize();
     ///////////////////////////////////////
     ///////////////////////////////////////
     clean_830->GetYaxis()->SetTitle("Counts (16 keV/bin)");
     clean_830->GetYaxis()->SetTitleFont(42);
     clean_830->GetYaxis()->SetTitleSize(0.06);
     clean_830->GetYaxis()->SetTitleOffset(0.6);
     clean_830->GetYaxis()->SetLabelFont(42);
     clean_830->GetYaxis()->SetLabelSize();
     clean_830->SetLineColor(1);/////rengini degistirdim
     clean_830->SetLineWidth(3);///line thickness
     clean_830->Rebin(16);
     clean_830->GetXaxis()->SetRange(10,140);
     clean_830->SetMinimum(-1);
     clean_830->Draw();
     cout << "done 2/2 " << endl;


     TCanvas *c5 = new TCanvas("c5","gated spectra_3",700,800);
     c5->ToggleEventStatus();
     c5->SetTitle("without and without BG subtraction");
     c5->Divide(1,2);
     c5->cd(1);
     hgammaprompt_1230->SetTitle("gate on 1230 keV without BG subtraction");
     hgammaprompt_1230->SetTitleSize(1);
     hgammaprompt_1230->GetXaxis()->SetTitle("Energy (keV)");
     hgammaprompt_1230->GetXaxis()->SetTitleFont(42);
     hgammaprompt_1230->GetXaxis()->SetTitleSize(0.06);
     hgammaprompt_1230->GetXaxis()->SetTitleOffset(0.7);
     hgammaprompt_1230->GetXaxis()->SetLabelFont(42);
     hgammaprompt_1230->GetXaxis()->SetLabelSize();
     ///////////////////////////////////////
     ///////////////////////////////////////
     hgammaprompt_1230->GetYaxis()->SetTitle("Counts (16 keV/bin)");
     hgammaprompt_1230->GetYaxis()->SetTitleFont(42);
     hgammaprompt_1230->GetYaxis()->SetTitleSize(0.06);
     hgammaprompt_1230->GetYaxis()->SetTitleOffset(0.6);
     hgammaprompt_1230->GetYaxis()->SetLabelFont(42);
     hgammaprompt_1230->GetYaxis()->SetLabelSize();
     hgammaprompt_1230->SetLineColor(1);/////rengini degistirdim
     hgammaprompt_1230->SetLineWidth(3);///line thickness
     hgammaprompt_1230->Rebin(16);
     hgammaprompt_1230->GetXaxis()->SetRange(10,140);
     hgammaprompt_1230->SetMinimum(-1);
     hgammaprompt_1230->Draw();
     cout << "done 1/1 " << endl;
     c5->cd(2);
     clean_1230->SetTitle("gate on 1230 keV with BG subtraction");
     clean_1230->SetTitleSize(1);
     clean_1230->GetXaxis()->SetTitle("Energy (keV)");
     clean_1230->GetXaxis()->SetTitleFont(42);
     clean_1230->GetXaxis()->SetTitleSize(0.06);
     clean_1230->GetXaxis()->SetTitleOffset(0.7);
     clean_1230->GetXaxis()->SetLabelFont(42);
     clean_1230->GetXaxis()->SetLabelSize();
     ///////////////////////////////////////
     ///////////////////////////////////////
     clean_1230->GetYaxis()->SetTitle("Counts (16 keV/bin)");
     clean_1230->GetYaxis()->SetTitleFont(42);
     clean_1230->GetYaxis()->SetTitleSize(0.06);
     clean_1230->GetYaxis()->SetTitleOffset(0.6);
     clean_1230->GetYaxis()->SetLabelFont(42);
     clean_1230->GetYaxis()->SetLabelSize();
     clean_1230->SetLineColor(1);/////rengini degistirdim
     clean_1230->SetLineWidth(3);///line thickness
     clean_1230->Rebin(16);
     clean_1230->GetXaxis()->SetRange(10,140);
     clean_1230->SetMinimum(-1);
     clean_1230->Draw();
     cout << "done 2/2 " << endl;
    */

     
     TCanvas *c6 = new TCanvas("c6","gated spectra_$",700,800);
     c6->ToggleEventStatus();
     c6->SetTitle("with BG subtraction");
     c6->Divide(1,3);
     c6->cd(1);
     clean_730->SetTitle("Gated on 730 keV");
     clean_730->SetTitleSize(1);
     clean_730->GetXaxis()->SetTitle("Energy (keV)");
     clean_730->GetXaxis()->SetTitleFont(62);
     clean_730->GetXaxis()->SetTitleSize(0.06);
     clean_730->GetXaxis()->SetTitleOffset(0.7);
     clean_730->GetXaxis()->SetLabelFont(42);
     clean_730->GetXaxis()->SetLabelSize();
     ///////////////////////////////////////
     ///////////////////////////////////////
     clean_730->GetYaxis()->SetTitle("Counts (16 keV/bin)");
     clean_730->GetYaxis()->SetTitleFont(62);
     clean_730->GetYaxis()->SetTitleSize(0.06);
     clean_730->GetYaxis()->SetTitleOffset(0.45);
     clean_730->GetYaxis()->SetLabelFont(42);
     clean_730->GetYaxis()->SetLabelSize();
     clean_730->SetLineColor(1);/////rengini degistirdim
     clean_730->SetLineWidth(3);///line thickness
     clean_730->Rebin(16);
     clean_730->GetXaxis()->SetRange(10,140);
     clean_730->SetMinimum(-1);
     clean_730->Draw();
     cout << "done 1/3 " << endl;
     c6->cd(2);
     clean_830->SetTitle("Gated on 889 keV");
     clean_830->SetTitleSize(1);
     clean_830->GetXaxis()->SetTitle("Energy (keV)");
     clean_830->GetXaxis()->SetTitleFont(62);
     clean_830->GetXaxis()->SetTitleSize(0.06);
     clean_830->GetXaxis()->SetTitleOffset(0.7);
     clean_830->GetXaxis()->SetLabelFont(42);
     clean_830->GetXaxis()->SetLabelSize();
     ///////////////////////////////////////
     ///////////////////////////////////////
     clean_830->GetYaxis()->SetTitle("Counts (16 keV/bin)");
     clean_830->GetYaxis()->SetTitleFont(62);
     clean_830->GetYaxis()->SetTitleSize(0.06);
     clean_830->GetYaxis()->SetTitleOffset(0.45);
     clean_830->GetYaxis()->SetLabelFont(42);
     clean_830->GetYaxis()->SetLabelSize();
     clean_830->SetLineColor(2);/////rengini degistirdim
     clean_830->SetLineWidth(3);///line thickness
     clean_830->Rebin(16);
     clean_830->GetXaxis()->SetRange(10,140);
     clean_830->SetMinimum(-1);
     clean_830->Draw();
     cout << "done 2/3 " << endl;
     c6->cd(3);
     clean_1230->SetTitle("Gated on 1215 keV");
     clean_1230->SetTitleSize(1);
     clean_1230->GetXaxis()->SetTitle("Energy (keV)");
     clean_1230->GetXaxis()->SetTitleFont(62);
     clean_1230->GetXaxis()->SetTitleSize(0.06);
     clean_1230->GetXaxis()->SetTitleOffset(0.7);
     clean_1230->GetXaxis()->SetLabelFont(42);
     clean_1230->GetXaxis()->SetLabelSize();
     ///////////////////////////////////////
     ///////////////////////////////////////
     clean_1230->GetYaxis()->SetTitle("Counts (16 keV/bin)");
     clean_1230->GetYaxis()->SetTitleFont(62);
     clean_1230->GetYaxis()->SetTitleSize(0.06);
     clean_1230->GetYaxis()->SetTitleOffset(0.45);
     clean_1230->GetYaxis()->SetLabelFont(62);
     clean_1230->GetYaxis()->SetLabelSize();
     clean_1230->SetLineColor(4);/////rengini degistirdim
     clean_1230->SetLineWidth(3);///line thickness
     clean_1230->Rebin(16);
     clean_1230->GetXaxis()->SetRange(10,140);
     clean_1230->SetMinimum(-1);
     clean_1230->Draw();
     cout << "done 3/3 " << endl;
     
     
     
/*TFile *fout = new TFile("78Zn_spectra.root","RECREATE");
  
    fout->cd();
    
   p_ggmatrix ->Write();
   ggmatrix->Write();
   gamma_single_mult1 ->Write();
   
   hgammaprompt_730->Write();
	
	
  clean_730->Write();
  fout->ls();
  fout->Close();
*/

}

