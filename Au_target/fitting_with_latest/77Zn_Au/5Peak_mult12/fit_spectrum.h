// Shailendra_04Feb 2021_____________________________________________________________________
//For 77,78,79Zn Mult1-6

TGraph *peak0g  ;
TGraph *peak1g  ;
TGraph *peak2g ;
TGraph *peak3g  ;
TGraph *peak4g ;
TGraph *peak5g ;
/*TGraph *peak6g ;
TGraph *peak7g ;
TGraph *peak8g ;
TGraph *peak9g ;
TGraph *peak10g ;
TGraph *peak11g ;
TGraph *peak12g ;
TGraph *peak13g ;*/



Double_t expf(Double_t *x, Double_t *par) //triple exp
{
    return TMath::Exp(par[0]+par[1]*x[0]) + TMath::Exp(par[2]+par[3]*x[0]) ;
}



Double_t resp0(Double_t *x, Double_t *par)
{
    return peak0g->Eval(x[0]-par[1])*par[0];
}



Double_t resp1(Double_t *x, Double_t *par)
{
    return peak1g->Eval(x[0]-par[1])*par[0];

}




Double_t resp2(Double_t *x, Double_t *par)
{
    return peak2g->Eval(x[0]-par[1])*par[0];
}





Double_t resp3(Double_t *x, Double_t *par)
{
    return peak3g->Eval(x[0]-par[1])*par[0];
}



Double_t resp4(Double_t *x, Double_t *par)
{
    return peak4g->Eval(x[0]-par[1])*par[0];
}


/*Double_t resp5(Double_t *x, Double_t *par)
{
    return peak5g->Eval(x[0]-par[1])*par[0];
}

*/






Double_t ex_respf(Double_t *x, Double_t *par)
{
    return expf(x,par) + resp0(x,par+4) +resp1(x,par+6) +resp2(x,par+8)+resp3(x,par+10)+resp4(x,par+12);//+resp5(x,par+14);
}


