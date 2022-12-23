 {


TFile *f = new TFile("three_Analysis_yes_dalitrigger_All_DALI_on_11_may_2021.root");
TH2F *ggmatrix;
TH1F *p_ggmatrix;
TH1F* gamma_single_mult1;
TH1F* gamma_single_mult2;
TH1F* gamma_single_mult3;

TH1F* gamma_single_mult12;
TH1F* gamma_single_mult123;
TH1F* gamma_single_mult1234;

    
     
ggmatrix = (TH2F*)f->Get("gg_mult234_79Zn");
p_ggmatrix  = (TH1F*)ggmatrix->ProjectionY();


//gg_mult23_p_4096
//gg_multAll_p_4096
gamma_single_mult1 = (TH1F*)f->Get("hEnergy_mult1_wa_79Zn");
gamma_single_mult2 = (TH1F*)f->Get("hEnergy_mult2_wa_79Zn");
gamma_single_mult3 = (TH1F*)f->Get("hEnergy_mult3_wa_79Zn");
//gamma_single_mult12 = (TH1F*)f->Get("hEnergy_mult4_wa_79Zn");
gamma_single_mult12 = (TH1F*)f->Get("hEnergy_mult12_wa_79Zn");
gamma_single_mult123 = (TH1F*)f->Get("hEnergy_mult123_wa_79Zn");
//gamma_single_mult1234 = (TH1F*)f->Get("hEnergy_mult1234_wa_79Zn");


     
hgammaprompt_980 = (TH1F*)ggmatrix->ProjectionY("gate on 980 keV",890,1000);
hgammaleft_980 = (TH1F*)ggmatrix->ProjectionY("BG left",510,560);
hgammaright_980 = (TH1F*)ggmatrix->ProjectionY("BG right",1050,1080);

hgammaprompt_1380 = (TH1F*)ggmatrix->ProjectionY("gate on 1380 keV",1340,1470);
hgammaleft_1380 = (TH1F*)ggmatrix->ProjectionY("BG left",1030,1090);
hgammaright_1380 = (TH1F*)ggmatrix->ProjectionY("BG right",1500,1600);
     
     hgammaprompt_1260 = (TH1F*)ggmatrix->ProjectionY("gate on 1260 keV",1170,1290);
     hgammaleft_1260 = (TH1F*)ggmatrix->ProjectionY("BG left",1030,1090);
     hgammaright_1260 = (TH1F*)ggmatrix->ProjectionY("BG right",1500,1600);

TH1F* clean_980 = new TH1F(*hgammaprompt_980);
TH1F* bg_total_980 = new TH1F(*hgammaleft_980);
bg_total_980->Add(hgammaright_980,+1.);
clean_980->Add(bg_total_980, -1./2.);
//clean_980->Add(bg_total_980, -1./(p_ggmatrix->Integral(645,760)+p_ggmatrix->Integral(560,590))/p_ggmatrix->Integral(810,840));
 (TH1F*)clean_980->Clone("gate  with BG subtraction  980_keV");

TH1F* clean_1380 = new TH1F(*hgammaprompt_1380);
TH1F* bg_total_1380 = new TH1F(*hgammaleft_1380);
bg_total_1380->Add(hgammaright_1380,+1.);
clean_1380->Add(bg_total_1380, -1./2.);
(TH1F*)clean_1380->Clone("gate  with BG subtraction  1380_keV");

TH1F* clean_1260 = new TH1F(*hgammaprompt_1260);
TH1F* bg_total_1260 = new TH1F(*hgammaleft_1260);
bg_total_1260->Add(hgammaright_1260,+1.);
clean_1260->Add(bg_total_1260, -1./2.);
(TH1F*)clean_1260->Clone("gate  with BG subtraction  1260_keV");


/*TH1F* clean_880 = new TH1F(*hgammaprompt_880);
TH1F* bg_total_880 = new TH1F(*hgammaleft_880);
bg_total_880->Add(hgammaright_880,+1.);
clean_880->Add(bg_total_880, -1./2.);
 (TH1F*)clean_880->Clone("gate  with BG subtraction  880_keV");*/


gStyle->SetOptStat(0);
gStyle->SetTitleFontSize(0.07);

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
gamma_single_mult2->SetTitle("Gamma singles mult=2");   //("Gamma singles mult=12");
gamma_single_mult2->GetXaxis()->SetTitle("Energy (keV)");
gamma_single_mult2->GetXaxis()->SetTitleFont(42);
gamma_single_mult2->GetXaxis()->SetTitleSize(0.06);
gamma_single_mult2->GetXaxis()->SetTitleOffset(0.7);
gamma_single_mult2->GetXaxis()->SetLabelFont(42);
gamma_single_mult2->GetXaxis()->SetLabelSize();
///////////////////////////////////////
///////////////////////////////////////
gamma_single_mult2->GetYaxis()->SetTitle("Counts (16 keV/bin)");
gamma_single_mult2->GetYaxis()->SetTitleFont(42);
gamma_single_mult2->GetYaxis()->SetTitleSize(0.05);
gamma_single_mult2->GetYaxis()->SetTitleOffset(0.6);
gamma_single_mult2->GetYaxis()->SetLabelFont(42);
gamma_single_mult2->GetYaxis()->SetLabelSize();
gamma_single_mult2->SetLineColor(3);/////rengini degistirdim
gamma_single_mult2->SetLineWidth(3);///line thickness
//gamma_single_mult2->Rebin(16);
gamma_single_mult2->GetXaxis()->SetRange(15,120);
gamma_single_mult2->Draw();

c1->cd(3);
gamma_single_mult3->SetTitle("Gamma singles mult=3");   //("Gamma singles mult=123");
gamma_single_mult3->GetXaxis()->SetTitle("Energy (keV)");
gamma_single_mult3->GetXaxis()->SetTitleFont(42);
gamma_single_mult3->GetXaxis()->SetTitleSize(0.06);
gamma_single_mult3->GetXaxis()->SetTitleOffset(0.7);
gamma_single_mult3->GetXaxis()->SetLabelFont(42);
gamma_single_mult3->GetXaxis()->SetLabelSize();
///////////////////////////////////////
///////////////////////////////////////
gamma_single_mult3->GetYaxis()->SetTitle("Counts (16 keV/bin)");
gamma_single_mult3->GetYaxis()->SetTitleFont(42);
gamma_single_mult3->GetYaxis()->SetTitleSize(0.05);
gamma_single_mult3->GetYaxis()->SetTitleOffset(0.6);
gamma_single_mult3->GetYaxis()->SetLabelFont(42);
gamma_single_mult3->GetYaxis()->SetLabelSize();
gamma_single_mult3->SetLineColor(4);/////rengini degistirdim
gamma_single_mult3->SetLineWidth(3);///line thickness
//gamma_single_mult3->Rebin(16);
gamma_single_mult3->GetXaxis()->SetRange(15,120);
//gamma_single_mult3->SetMinimum(-100);

gamma_single_mult3->Draw();

     
/*c1->cd(4);
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
     hgammaprompt_980->SetTitle("Gate on 983 keV without BG subtraction");
     hgammaprompt_980->SetTitleSize(1);
     hgammaprompt_980->GetXaxis()->SetTitle("Energy (keV)");
     hgammaprompt_980->GetXaxis()->SetTitleFont(42);
     hgammaprompt_980->GetXaxis()->SetTitleSize(0.06);
     hgammaprompt_980->GetXaxis()->SetTitleOffset(0.7);
     hgammaprompt_980->GetXaxis()->SetLabelFont(42);
     hgammaprompt_980->GetXaxis()->SetLabelSize();
     ///////////////////////////////////////
     ///////////////////////////////////////
     hgammaprompt_980->GetYaxis()->SetTitle("Counts (16 keV/bin)");
     hgammaprompt_980->GetYaxis()->SetTitleFont(42);
     hgammaprompt_980->GetYaxis()->SetTitleSize(0.06);
     hgammaprompt_980->GetYaxis()->SetTitleOffset(0.6);
     hgammaprompt_980->GetYaxis()->SetLabelFont(42);
     hgammaprompt_980->GetYaxis()->SetLabelSize();
     hgammaprompt_980->SetLineColor(1);/////rengini degistirdim
     hgammaprompt_980->SetLineWidth(3);///line thickness
     hgammaprompt_980->Rebin(16);
     hgammaprompt_980->GetXaxis()->SetRange(10,140);
     hgammaprompt_980->SetMinimum(-1);
     hgammaprompt_980->Draw();
     cout << "done 1/1 " << endl;
     c3->cd(2);
     clean_980->SetTitle("Gate on 983 keV with BG subtraction");
     clean_980->SetTitleSize(1);
     clean_980->GetXaxis()->SetTitle("Energy (keV)");
     clean_980->GetXaxis()->SetTitleFont(42);
     clean_980->GetXaxis()->SetTitleSize(0.06);
     clean_980->GetXaxis()->SetTitleOffset(0.7);
     clean_980->GetXaxis()->SetLabelFont(42);
     clean_980->GetXaxis()->SetLabelSize();
     ///////////////////////////////////////
     ///////////////////////////////////////
     clean_980->GetYaxis()->SetTitle("Counts (16 keV/bin)");
     clean_980->GetYaxis()->SetTitleFont(42);
     clean_980->GetYaxis()->SetTitleSize(0.06);
     clean_980->GetYaxis()->SetTitleOffset(0.6);
     clean_980->GetYaxis()->SetLabelFont(42);
     clean_980->GetYaxis()->SetLabelSize();
     clean_980->SetLineColor(1);/////rengini degistirdim
     clean_980->SetLineWidth(3);///line thickness
     clean_980->Rebin(16);
     clean_980->GetXaxis()->SetRange(10,140);
     clean_980->SetMinimum(-1);
     clean_980->Draw();
     cout << "done 2/2 " << endl;
     
     
     TCanvas *c4 = new TCanvas("c4","gated spectra_2",700,800);
     c4->ToggleEventStatus();
     c4->SetTitle("without and without BG subtraction");
     c4->Divide(1,2);
     c4->cd(1);
     hgammaprompt_1380->SetTitle("Gate on 1380 keV without BG subtraction");
     hgammaprompt_1380->SetTitleSize(1);
     hgammaprompt_1380->GetXaxis()->SetTitle("Energy (keV)");
     hgammaprompt_1380->GetXaxis()->SetTitleFont(42);
     hgammaprompt_1380->GetXaxis()->SetTitleSize(0.06);
     hgammaprompt_1380->GetXaxis()->SetTitleOffset(0.7);
     hgammaprompt_1380->GetXaxis()->SetLabelFont(42);
     hgammaprompt_1380->GetXaxis()->SetLabelSize();
     ///////////////////////////////////////
     ///////////////////////////////////////
     hgammaprompt_1380->GetYaxis()->SetTitle("Counts (16 keV/bin)");
     hgammaprompt_1380->GetYaxis()->SetTitleFont(42);
     hgammaprompt_1380->GetYaxis()->SetTitleSize(0.06);
     hgammaprompt_1380->GetYaxis()->SetTitleOffset(0.6);
     hgammaprompt_1380->GetYaxis()->SetLabelFont(42);
     hgammaprompt_1380->GetYaxis()->SetLabelSize();
     hgammaprompt_1380->SetLineColor(1);/////rengini degistirdim
     hgammaprompt_1380->SetLineWidth(3);///line thickness
     hgammaprompt_1380->Rebin(16);
     hgammaprompt_1380->GetXaxis()->SetRange(10,140);
     hgammaprompt_1380->SetMinimum(-1);
     hgammaprompt_1380->Draw();
     cout << "done 1/1 " << endl;
     c4->cd(2);
     clean_1380->SetTitle("Gate on 1380 keV with BG subtraction");
     clean_1380->SetTitleSize(1);
     clean_1380->GetXaxis()->SetTitle("Energy (keV)");
     clean_1380->GetXaxis()->SetTitleFont(42);
     clean_1380->GetXaxis()->SetTitleSize(0.06);
     clean_1380->GetXaxis()->SetTitleOffset(0.7);
     clean_1380->GetXaxis()->SetLabelFont(42);
     clean_1380->GetXaxis()->SetLabelSize();
     ///////////////////////////////////////
     ///////////////////////////////////////
     clean_1380->GetYaxis()->SetTitle("Counts (16 keV/bin)");
     clean_1380->GetYaxis()->SetTitleFont(42);
     clean_1380->GetYaxis()->SetTitleSize(0.06);
     clean_1380->GetYaxis()->SetTitleOffset(0.6);
     clean_1380->GetYaxis()->SetLabelFont(42);
     clean_1380->GetYaxis()->SetLabelSize();
     clean_1380->SetLineColor(1);/////rengini degistirdim
     clean_1380->SetLineWidth(3);///line thickness
     clean_1380->Rebin(16);
     clean_1380->GetXaxis()->SetRange(10,140);
     clean_1380->SetMinimum(-1);
     clean_1380->Draw();
     cout << "done 2/2 " << endl;


     TCanvas *c5 = new TCanvas("c5","gated spectra_3",700,800);
     c5->ToggleEventStatus();
     c5->SetTitle("without and without BG subtraction");
     c5->Divide(1,2);
     c5->cd(1);
     hgammaprompt_1260->SetTitle("Gate on 1260 keV without BG subtraction");
     hgammaprompt_1260->SetTitleSize(1);
     hgammaprompt_1260->GetXaxis()->SetTitle("Energy (keV)");
     hgammaprompt_1260->GetXaxis()->SetTitleFont(42);
     hgammaprompt_1260->GetXaxis()->SetTitleSize(0.06);
     hgammaprompt_1260->GetXaxis()->SetTitleOffset(0.7);
     hgammaprompt_1260->GetXaxis()->SetLabelFont(42);
     hgammaprompt_1260->GetXaxis()->SetLabelSize();
     ///////////////////////////////////////
     ///////////////////////////////////////
     hgammaprompt_1260->GetYaxis()->SetTitle("Counts (16 keV/bin)");
     hgammaprompt_1260->GetYaxis()->SetTitleFont(42);
     hgammaprompt_1260->GetYaxis()->SetTitleSize(0.06);
     hgammaprompt_1260->GetYaxis()->SetTitleOffset(0.6);
     hgammaprompt_1260->GetYaxis()->SetLabelFont(42);
     hgammaprompt_1260->GetYaxis()->SetLabelSize();
     hgammaprompt_1260->SetLineColor(1);/////rengini degistirdim
     hgammaprompt_1260->SetLineWidth(3);///line thickness
     hgammaprompt_1260->Rebin(16);
     hgammaprompt_1260->GetXaxis()->SetRange(10,140);
     hgammaprompt_1260->SetMinimum(-1);
     hgammaprompt_1260->Draw();
     cout << "done 1/1 " << endl;
     c5->cd(2);
     clean_1260->SetTitle("Gate on 1260 keV with BG subtraction");
     clean_1260->SetTitleSize(1);
     clean_1260->GetXaxis()->SetTitle("Energy (keV)");
     clean_1260->GetXaxis()->SetTitleFont(42);
     clean_1260->GetXaxis()->SetTitleSize(0.06);
     clean_1260->GetXaxis()->SetTitleOffset(0.7);
     clean_1260->GetXaxis()->SetLabelFont(42);
     clean_1260->GetXaxis()->SetLabelSize();
     ///////////////////////////////////////
     ///////////////////////////////////////
     clean_1260->GetYaxis()->SetTitle("Counts (16 keV/bin)");
     clean_1260->GetYaxis()->SetTitleFont(42);
     clean_1260->GetYaxis()->SetTitleSize(0.06);
     clean_1260->GetYaxis()->SetTitleOffset(0.6);
     clean_1260->GetYaxis()->SetLabelFont(42);
     clean_1260->GetYaxis()->SetLabelSize();
     clean_1260->SetLineColor(1);/////rengini degistirdim
     clean_1260->SetLineWidth(3);///line thickness
     clean_1260->Rebin(16);
     clean_1260->GetXaxis()->SetRange(10,140);
     clean_1260->SetMinimum(-1);
     clean_1260->Draw();
     cout << "done 2/2 " << endl;
    

     /*
     TCanvas *c6 = new TCanvas("c6","gated spectra_$",700,800);
     c6->ToggleEventStatus();
     c6->SetTitle("with BG subtraction");
     c6->Divide(1,3);
     c6->cd(1);
     clean_980->SetTitle("gate on 980 keV");
     clean_980->SetTitleSize(1);
     clean_980->GetXaxis()->SetTitle("Energy (keV)");
     clean_980->GetXaxis()->SetTitleFont(42);
     clean_980->GetXaxis()->SetTitleSize(0.06);
     clean_980->GetXaxis()->SetTitleOffset(0.7);
     clean_980->GetXaxis()->SetLabelFont(42);
     clean_980->GetXaxis()->SetLabelSize();
     ///////////////////////////////////////
     ///////////////////////////////////////
     clean_980->GetYaxis()->SetTitle("Counts (16 keV/bin)");
     clean_980->GetYaxis()->SetTitleFont(42);
     clean_980->GetYaxis()->SetTitleSize(0.06);
     clean_980->GetYaxis()->SetTitleOffset(0.6);
     clean_980->GetYaxis()->SetLabelFont(42);
     clean_980->GetYaxis()->SetLabelSize();
     clean_980->SetLineColor(1);/////rengini degistirdim
     clean_980->SetLineWidth(3);///line thickness
     clean_980->Rebin(16);
     clean_980->GetXaxis()->SetRange(10,140);
     clean_980->SetMinimum(-1);
     clean_980->Draw();
     cout << "done 1/3 " << endl;
     c6->cd(2);
     clean_1380->SetTitle("gate on 1380 keV");
     clean_1380->SetTitleSize(1);
     clean_1380->GetXaxis()->SetTitle("Energy (keV)");
     clean_1380->GetXaxis()->SetTitleFont(42);
     clean_1380->GetXaxis()->SetTitleSize(0.06);
     clean_1380->GetXaxis()->SetTitleOffset(0.7);
     clean_1380->GetXaxis()->SetLabelFont(42);
     clean_1380->GetXaxis()->SetLabelSize();
     ///////////////////////////////////////
     ///////////////////////////////////////
     clean_1380->GetYaxis()->SetTitle("Counts (16 keV/bin)");
     clean_1380->GetYaxis()->SetTitleFont(42);
     clean_1380->GetYaxis()->SetTitleSize(0.06);
     clean_1380->GetYaxis()->SetTitleOffset(0.6);
     clean_1380->GetYaxis()->SetLabelFont(42);
     clean_1380->GetYaxis()->SetLabelSize();
     clean_1380->SetLineColor(1);/////rengini degistirdim
     clean_1380->SetLineWidth(3);///line thickness
     clean_1380->Rebin(16);
     clean_1380->GetXaxis()->SetRange(10,140);
     clean_1380->SetMinimum(-1);
     clean_1380->Draw();
     cout << "done 2/3 " << endl;
     c6->cd(3);
     clean_1260->SetTitle("gate on 1260 keV");
     clean_1260->SetTitleSize(1);
     clean_1260->GetXaxis()->SetTitle("Energy (keV)");
     clean_1260->GetXaxis()->SetTitleFont(42);
     clean_1260->GetXaxis()->SetTitleSize(0.06);
     clean_1260->GetXaxis()->SetTitleOffset(0.7);
     clean_1260->GetXaxis()->SetLabelFont(42);
     clean_1260->GetXaxis()->SetLabelSize();
     ///////////////////////////////////////
     ///////////////////////////////////////
     clean_1260->GetYaxis()->SetTitle("Counts (16 keV/bin)");
     clean_1260->GetYaxis()->SetTitleFont(42);
     clean_1260->GetYaxis()->SetTitleSize(0.06);
     clean_1260->GetYaxis()->SetTitleOffset(0.6);
     clean_1260->GetYaxis()->SetLabelFont(42);
     clean_1260->GetYaxis()->SetLabelSize();
     clean_1260->SetLineColor(1);/////rengini degistirdim
     clean_1260->SetLineWidth(3);///line thickness
     clean_1260->Rebin(16);
     clean_1260->GetXaxis()->SetRange(10,140);
     clean_1260->SetMinimum(-1);
     clean_1260->Draw();
     cout << "done 3/3 " << endl;
   
*/
     
     
/*TFile *fout = new TFile("79Zn_spectra.root","RECREATE");
  
    fout->cd();
    
   p_ggmatrix ->Write();
   ggmatrix->Write();
   gamma_single_mult1 ->Write();
   
   hgammaprompt_980->Write();
	
	
  clean_980->Write();
  fout->ls();
  fout->Close();
*/

}

