{
  TFile *f = new TFile("137Cs_-2000_2000_15_4Feb.root");
  TH2F* h[10];
  TH1F* hp[10];
  TFile *fout = new TFile("converted_137cs_calibrated_histograms.root","RECREATE");
  for(int i=0;i<10;i++){
    h[i] = (TH2F*)f->Get(Form("h_doppler[%d]",i));
    h[i]->SetName(Form("h_doppler_%d",i));
    fout->cd();
    //h[i]->Write();
    hp[i] = (TH1F*)h[i]->ProjectionY(Form("hp_doppler_%d",i),0,186);
    hp[i]->Write();
  }
  fout->ls();
  fout->Close();
}
