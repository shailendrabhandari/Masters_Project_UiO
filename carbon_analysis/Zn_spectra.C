{
gStyle->SetPalette();

//// to save banana gates /////
//TFile *f=new TFile("gates_dali_spectra.root","UPDATE");
//ge->Write(); 
f->Close(); 

//TFile *f=new TFile("gates_dali_spectra1.root","UPDATE");
TCut  *br_79Zn, *zd_79Zn; *br_78Zn, *zd_78Zn , *br_77Zn, *zd_77Zn ;
 
br_79Zn = (TCut *)gDirectory->Get("br_79Zn");
zd_79Zn = (TCut *)gDirectory->Get("zd_79Zn");
br_78Zn = (TCut *)gDirectory->Get("br_78Zn");
zd_78Zn = (TCut *)gDirectory->Get("zd_78Zn");
br_77Zn = (TCut *)gDirectory->Get("br_77Zn");
zd_77Zn = (TCut *)gDirectory->Get("zd_77Zn");

 f->Close();


TChain chain("tree");

    //chain.Add("/Volumes/Seagate2/rootfiles_79Zn/Gamma150034.root");// 79Zn + 197Au

    //chain.Add("../te.root");
  

////For Carbon Target
    chain.Add("/media/shailendra/Backup Plus/New_root_files_may/rootfiles_77Cu+carbon/Gamma150093.root");
    chain.Add("/media/shailendra/Backup Plus/New_root_files_may/rootfiles_77Cu+carbon/Gamma150094.root");
    chain.Add("/media/shailendra/Backup Plus/New_root_files_may/rootfiles_77Cu+carbon/Gamma150095.root");
    /*chain.Add("/media/shailendra/Backup Plus/New_root_files_may/rootfiles_77Cu+carbon/Gamma150096.root");
    chain.Add("/media/shailendra/Backup Plus/New_root_files_may/rootfiles_77Cu+carbon/Gamma150097.root");
    chain.Add("/media/shailendra/Backup Plus/New_root_files_may/rootfiles_77Cu+carbon/Gamma150098.root");
    chain.Add("/media/shailendra/Backup Plus/New_root_files_may/rootfiles_77Cu+carbon/Gamma150099.root");
    chain.Add("/media/shailendra/Backup Plus/New_root_files_may/rootfiles_77Cu+carbon/Gamma150100.root");
    chain.Add("/media/shailendra/Backup Plus/New_root_files_may/rootfiles_77Cu+carbon/Gamma150101.root");
    chain.Add("/media/shailendra/Backup Plus/New_root_files_may/rootfiles_77Cu+carbon/Gamma150102.root");
    chain.Add("/media/shailendra/Backup Plus/New_root_files_may/rootfiles_77Cu+carbon/Gamma150103.root");
    chain.Add("/media/shailendra/Backup Plus/New_root_files_may/rootfiles_77Cu+carbon/Gamma150104.root");
    chain.Add("/media/shailendra/Backup Plus/New_root_files_may/rootfiles_77Cu+carbon/Gamma150105.root");
    chain.Add("/media/shailendra/Backup Plus/New_root_files_may/rootfiles_77Cu+carbon/Gamma150106.root");
    chain.Add("/media/shailendra/Backup Plus/New_root_files_may/rootfiles_77Cu+carbon/Gamma150107.root");
    chain.Add("/media/shailendra/Backup Plus/New_root_files_may/rootfiles_77Cu+carbon/Gamma150108.root");
    chain.Add("/media/shailendra/Backup Plus/New_root_files_may/rootfiles_77Cu+carbon/Gamma150109.root");
    chain.Add("/media/shailendra/Backup Plus/New_root_files_may/rootfiles_77Cu+carbon/Gamma150110.root");
    chain.Add("/media/shailendra/Backup Plus/New_root_files_may/rootfiles_77Cu+carbon/Gamma150111.root");
    chain.Add("/media/shailendra/Backup Plus/New_root_files_may/rootfiles_77Cu+carbon/Gamma150112.root");
    chain.Add("/media/shailendra/Backup Plus/New_root_files_may/rootfiles_77Cu+carbon/Gamma150113.root");
    chain.Add("/media/shailendra/Backup Plus/New_root_files_may/rootfiles_77Cu+carbon/Gamma150114.root");
    chain.Add("/media/shailendra/Backup Plus/New_root_files_may/rootfiles_77Cu+carbon/Gamma150115.root");
    chain.Add("/media/shailendra/Backup Plus/New_root_files_may/rootfiles_77Cu+carbon/Gamma150116.root");*/







////For Gold Target

   /* chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150034.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150035.root");
     chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150036.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150037.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150038.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150039.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150040.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150041.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150042.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150043.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150044.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150045.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150046.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150047.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150048.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150049.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150050.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150051.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150052.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150053.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150054.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150055.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150056.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150057.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150058.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150059.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150060.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150061.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150062.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150063.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150064.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150065.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150066.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150067.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150068.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150069.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150070.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150071.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150072.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150073.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150074.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150075.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150076.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150077.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150078.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150079.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150080.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150081.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150082.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150083.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150084.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150085.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150086.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150087.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150088.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150089.root");
    chain.Add("/media/shailendra/Backup Plus/rootfiles_79Zn_2/Gamma150090.root");*/




/*
//////Creating PID Gates
TCanvas *pid = new TCanvas("pid","pid",700,800);
pid->Divide(1,2);
pid->cd(1);
pid_1->SetLogz();
chain.Draw("zet[2]:aoq[2]>>h1(1000,2.3,2.8,1000,26,34)","","col");
pid->cd(2);
pid_2->SetLogz();
chain.Draw("zet[5]:aoq[5]>>h2(1000,2.3,2.8,1000,26,34)","","col");*/

// tu put pid gates
    
//chain.Draw("zetc[2]:aoqc[2]>>h1(1000,2.3,2.8,1000,26,34)","","col"); // BR pid
//chain.Draw("zetc[5]:aoqc[5]>>h2(1000,2.3,2.8,1000,26,34)","","col"); //ZD pid

    /*
///Time alignmet /TimeOffseted
TCanvas *dali_time = new TCanvas("dali_time","dali_time",1200,800);
dali_time->Divide(1,3);
dali_time->cd(1);
dali_time_1->SetLogy();
chain.Draw("DALINaI.fTDC>>h1(1000,3000,14000)","DALINaI.id==1","");
dali_time->cd(2);
//dali_time_2->Setlogy();
chain.Draw("DALINaI.fTime>>h2(1000,-1500,-900)","DALINaI.id==1","");
dali_time->cd(3);
chain.Draw("DALINaI.fTimeOffseted>>h3(1000,-100,2000)","DALINaI.id==1","");*/

 //hw   
////TimeOfsetted selection for 79Zn
chain.Draw("DALINaI.fTimeOffseted>>h1(1000,-100,100)","dalimult<2&&fTimeOffseted>-100&&fTimeOffseted<100&&br_78Zn&&zd_78Zn",""); 
//("histogram","cuts&conditions", "options")
 //-10,7 //79Zn


///TimeOfsetted selection for 77Zn
//chain.Draw("DALINaI.fTimeOffseted:DALINaI.id>>(200,0,200,5000,-100,100)","dalimult<2&&fTimeOffseted>-10&&fTimeOffseted<10&&br_77Zn&&zd_77Zn","col");
//-10,8  77Zn



//79Zn  /////
//chain.Draw("DALINaI.fTimeOffseted:DALINaI.fDoppCorEnergy>>h1(256,0,4095,200,-100,100)","DALINaI.id>40&&dalimult<8&&fTimeOffseted>-80&&fTimeOffseted<80&&br_79Zn&&zd_79Zn","colz");
//chain.Draw("DALINaI.fTimeOffseted:DALINaI.fDoppCorEnergy>>Zn(256,0,4095,200,-100,100)","DALINaI.id>40&&dalimult<8&&zd_79Zn","colz");
// -10,8 //79Zn
//DALINaI.id==1

//chain.Draw("DALINaI.fDoppCorEnergy:DALINaI.id>>Dali_ID_78Zn(186,1,187,600,0,6000)","dalimult<3&&fTimeOffseted>-13&&fTimeOffseted<4&&br_78Zn&&zd_78Zn","colz");
//77Zn

//chain.Draw("DALINaI.fTimeOffseted:DALINaI.fDoppCorEnergy>>h1(256,0,4095,200,-100,100)","DALINaI.id>40&&dalimult<8&&fTimeOffseted>-100&&fTimeOffseted<100&&br_77Zn&&zd_77Zn","colz");
// -10,8 //77Zn


///79Zn Gamma-ray spectrum
//chain.Draw("DALINaI.fDoppCorEnergy>>79Zn(300,300,3000)","DALINaI.id>40&&dalimult<2&&fTimeOffseted>-10&&fTimeOffseted<6&&br_79Zn&&zd_79Zn","");

////77Zn
//chain.Draw("DALINaI.fDoppCorEnergy>>77Zn(1000,0,4000)","DALINaI.id>60&&dalimult<2&&fTimeOffseted>-11&&fTimeOffseted<6&&br_77Zn&&zd_77Zn","");

/*
TCanvas *dali = new TCanvas("dali","dali",1200,800);
dali->Divide(2,3);
dali->cd(1);
chain.Draw("beta[0]","br_79Zn&&zd_79Zn","");
dali->cd(2);
chain.Draw("beta[1]","br_79Zn&&zd_79Zn","");
dali->cd(3);
chain.Draw("beta[2]","br_79Zn&&zd_79Zn","");
dali->cd(4);
chain.Draw("beta[3]","br_79Zn&&zd_79Zn","");
dali->cd(5);
chain.Draw("beta[4]","br_79Zn&&zd_79Zn","");
dali->cd(6);
chain.Draw("beta[5]","br_79Zn&&zd_79Zn","");

*/


/*chain.Draw("DALINaI.fDoppCorEnergy>>Zn78_1(400,400,1200)","DALINaI.id>60&&dalimult<2&&fTimeOffseted>-13&&fTimeOffseted<0&&br_78Zn&&zd_78Zn","");
   

TFile *myfile = TFile::Open("78zn1.root", "recreate");
chain.Write();
delete myfile;*/
















//chain.Draw("DALINaI.fDoppCorEnergy:DALINaI.id>>Dali_ID_77Zn(186,1,187,600,0,6000)","dalimult<2&&fTimeOffseted>-15&&fTimeOffseted<40&&br_77Zn&&zd_77Zn","col");
//chain.Draw("DALINaI.fDoppCorEnergy>>77Zn(600,0,6000)","DALINaI.id>40&&dalimult<2&&fTimeOffseted>-10&&fTimeOffseted<10&&br_77Zn&&zd_77Zn","");
//chain.Draw("DALINaI.fDoppCorEnergy:DALINaI.id>>Dali_ID_78Zn(186,1,187,600,0,6000)","dalimult<2&&fTimeOffseted>-13&&fTimeOffseted<4&&br_78Zn&&zd_78Zn","colz");


//TCanvas *dali = new TCanvas("dali","dali",1200,800);

//dali->Divide(2,2);
/*
//eda version
dali->cd(1);
chain.Draw("DALINaI.fDoppCorEnergy:DALINaI.id>>Dali_ID_77Cu(186,1,187,600,0,6000)","dalimult<2&&fTimeOffseted>-10&&fTimeOffseted<10&&br_77Cu&&zd_77Cu","col");
cout << "done 1/4 " << endl;
dali->cd(2);
chain.Draw("DALINaI.fDoppCorEnergy>>77Cu(600,0,6000)","DALINaI.id>40&&dalimult<2&&fTimeOffseted>-10&&fTimeOffseted<10&&br_77Cu&&zd_77Cu","");
cout << "done 2/4 " << endl;
dali->cd(3);
chain.Draw("DALINaI.fDoppCorEnergy:DALINaI.id>>Dali_ID_78Zn(186,1,187,600,0,6000)","dalimult<2&&fTimeOffseted>-10&&fTimeOffseted<10&&br_78Zn&&zd_78Zn","col");
cout << "done 3/4 " << endl;
dali->cd(4);
chain.Draw("DALINaI.fDoppCorEnergy>>78Zn(600,0,6000)","DALINaI.id>40&&dalimult<2&&fTimeOffseted>-10&&fTimeOffseted<10&&br_78Zn&&zd_78Zn","");
cout << "done 4/4 " << endl;
*/
/*
//kwimmer
dali->cd(1);
chain.Draw("DALINaI.fDoppCorEnergy:DALINaI.id>>Dali_ID_77Cu(186,1,187,600,0,6000)","dalimult<4&&fTimeOffseted>-10&&fTimeOffseted<10&&br_77Cu&&zd_77Cu","col");
cout << "done 1/2 " << endl;
dali->cd(2);
TH1F *cu = (TH1F*)Dali_ID_77Cu->ProjectionY("cup",41,187);
cu->Draw();
dali->cd(3);
chain.Draw("DALINaI.fDoppCorEnergy:DALINaI.id>>Dali_ID_78Zn(186,1,187,600,0,6000)","dalimult<2&&fTimeOffseted>-10&&fTimeOffseted<10&&br_78Zn&&zd_78Zn 	","col");
cout << "done 2/2 " << endl;
dali->cd(4);
TH1F *zn = (TH1F*)Dali_ID_78Zn->ProjectionY("znp",41,187);
zn->Draw();
*/





/* //////beta_BR and beta_ZD values for 78Zn
    TCanvas *beta78Zn = new TCanvas("beta78Zn","beta78Zn",1000,1000);
    beta78Zn->Divide(2,1);
    beta78Zn->cd(1);
    chain.Draw("beta_37>>beta37(125,0.4,0.7)","br_78Zn&&zd_78Zn","");
    //beta37->Fit("gaus","","",0.545,0.567);
    beta78Zn->cd(2);
    chain.Draw("beta_811>>beta811(125,0.4,0.7)","br_78Zn&&zd_78Zn","");
    //beta811->Fit("gaus","","",0.499,0.521);*/







}




