#include <iostream>
#include <TApplication.h>
#include <TMatrix.h>
#include <TROOT.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TTree.h>
#include <TDirectoryFile.h>
#include <TLorentzVector.h>
#include <TString.h>
#include <TStopwatch.h>
#include <TPaveText.h>
#include <cstdlib>
#include <iomanip>

#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#include "TVector2.h"
#include "TFile.h"
#include "TString.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMath.h"
#include "TText.h"
#include "TRandom3.h"
#include "TArray.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TMinuit.h"

#include "TCanvas.h"
#include "TLatex.h"
#include "TImage.h"
#include "TGaxis.h"
#include "TGraphErrors.h"
#include "TPaveStats.h"
using namespace std;



double PolFunction(double *x, double *par);
double BWFunction(double *x, double *par);
double BWplusPol(double *x, double *par);
double GausFunction(double *x, double *par);
double GausplusPol(double *x, double *par);
double Voigtian(Double_t *x, Double_t *par);
double VoigtplusPol(double *x, double *par);
double Rapidity(double pt, double pz, double m);

//////  ****  Setup : start ****  //////
int SelectRap= 1; // 0:|y|<1.4  , 1:|y|<1.2   2:|y|<1.0    3:|y|<0.8   4:|y|<0.6    5:|y|<0.4
Double_t rapidity[6]={2.0, 1.6, 1.4, 1.2, 1.0, 0.8}; 
Double_t  dely=    rapidity[SelectRap]  ;  
double raprange[2]={-dely/2, dely/2};


const int ptbin = 8;    // # of pt bins
const int ptbinloop = 8;    // # of pt bins


//-STANDARD   
// 
bool UseMix =     kFALSE; 
bool TypeCount =  kFALSE ; 
bool iMixPol =    kTRUE ;  //   kFALSE; 
int POLdegree= 3 ; // 3; // pol-gauss Pol3 ,S
bool   UseTree =   kTRUE ;// use input from TTRee or histos // v45-->kTRUE ONLY

  Double_t    mix1min = 1.170 ; // 1.175 ; - A lot depends on that values - MAIN 
  Double_t    mix1max = 1.185 ; // 1.180 ;
  Double_t    mix2min = 1.205 ;
  Double_t    mix2max = 1.220;

 double FitRangeSigma0[2]={1.182,1.205}; 

 double FitRangeLSigma0[2]={1.170,1.215};  //-- data Bins0

 double FitRangeL2Sigma0[2]={1.16,1.230};

int BinPtRange = 1;

bool iMixBGMC =  kFALSE ;
//OFFI  
double SigmaCountBound[2]={1.186,1.200}; 
double SigmaMassCount[2]={1.182, 1.202}; 

//////  ****  Setup : end ****  //////


bool reject;
double yieldRange[2]={0};
double FitRange[2]={0};

double ptedges[ptbin+1+1]={0};

Int_t iBinPt[40]=               { 12,  17,  21,  25,  29,  35,  41, 51,  81, 241,281,321,361,401 };
double ptedges_Sigma[ptbin+1] = { 1.1, 1.6, 2.0, 2.4, 2.8, 3.4, 4., 5.0, 8.0 };
double pt_points[ptbin] =       { 1.35,1.8, 2.2, 2.6, 3.0, 3.8,    4.5, 6.5 };
double pt_points_e[ptbin] =     {  .25, .2,  .2, .2,   .3,  .3,    0.5, 1.5 };

double Sigma0Mass[ptbin]={0};
double Sigma0Width[ptbin]={0};
double Sigma0Mass_e[ptbin]={0};
double Sigma0Width_e[ptbin]={0};
double MissedYield[ptbin]={0};
double MissedYield_e[ptbin]={0};
double yield[ptbin]={0};
double yield_e[ptbin]={0};
double spectrum[ptbin]={0};
double spectrum_e[ptbin]={0};
double spectrum_esys[ptbin]={0};
double spectrum_etot[ptbin]={0};

double ratioev[ptbin]={0};
double ratioev_e[ptbin]={0};

double ratioTrueSub[ptbin]={0};
double ratioTrueSub_e[ptbin]={0};


double signifi[ptbin]={0};
double signifi_e[ptbin]={0};
double GenEv[ptbin]={0};
double GenEv_e[ptbin]={0};
double spectrum_mc[ptbin]={0};
double spectrum_mce[ptbin]={0};
double ratiomc[ptbin]={0};
double ratiomc_e[ptbin]={0};
double ratiomc_etot[ptbin]={0};

double EffA[ptbin]={0}; double EffA_e[ptbin]={0};
double EffP[ptbin]={0}; double EffP_e[ptbin]={0};
double EffAP[ptbin]={0}; double EffAP_e[ptbin]={0};
double EffAS[ptbin]={0}; double EffAS_e[ptbin]={0};
double EffPS[ptbin]={0}; double EffPS_e[ptbin]={0};
double EffAPS[ptbin]={0}; double EffAPS_e[ptbin]={0};
  // **  Calculate effieciency :start  ** //
    Double_t rec_num[40]; Double_t input_num[40];
    Double_t recA_num[40]; Double_t inputA_num[40];
    Double_t recP_num[40]; Double_t inputP_num[40];
    Double_t recAP_num[40]; Double_t inputAP_num[40];
    Double_t recAS_num[40]; Double_t inputAS_num[40];
    Double_t recPS_num[40]; Double_t inputPS_num[40];
    Double_t recAPS_num[40]; Double_t inputAPS_num[40];

 
double temp_yield=0, temp_yield_e=0, temp_yield1=0;
Double_t feloss[8] ={1,1,1,1,1,1,1,1};
Double_t fASigSig = 2.;
 Double_t flukageantAbase[8] ={ 0.935, 0.95, 0.96, 0.966, 0.970, 0.975, 0.982, 1.0 };   // USE !!!
Double_t flukageantA[8] ; 
Double_t flukageantP[8] ;
Double_t etacorr =   dely/ rapidity[4] ; // is not used

Double_t rTrueRec[8] ={ 0.79, 0.82, 0.85, 0.87, 0.91, 0.957, 1.1, 1.18 };

Double_t  Eff[ptbin]  ={9.33415e-05, 0.000400704, 0.000873701, 0.00160246, 0.00249329,  0.00327019,  0.00408349,  0.00489897};
                          
Double_t  Eff_e[ptbin]={5.70201e-06, 1.91066e-05, 3.89734e-05, 7.21558e-05,0.000104931,0.000179498,0.000241532,0.000339698};

 Double_t fsystotal[ptbin]= {0};


Double_t fsyserr = 0.085;  // Note 1052 -- Mat. budget pp @ 7 TeV for TWO photons!!!

Double_t fsysINEL = 0.; //  0.062;

Double_t fsystfluct[ptbin]= {0};
 Double_t fsystLam[ptbin]= {0.0707,  0.0650,   0.0606,  0.0565,  0.0530,  0.0460,  0.0405,  0.0285};

 Double_t fsystGam[ptbin]= {0.044,  0.044,  0.044,  0.044,  0.044,  0.044, 0.044,  0.044 };


  Double_t fsystMass[ptbin]= {0.052, 0.044, 0.039, 0.0325, 0.028, 0.023, 0.018, 0.009 } ;   // v 442 - 5jun18
 
 Double_t npi0ev[100], enpi0ev[100],  nSumpi0ev[100],  enSumpi0ev[100]  ;

Double_t mSigmax=1.230;  Double_t mSigmin=1.130;                                                                          
                                    
Double_t mLammax=1.122;  Double_t mLammin=1.110;                                                                          
                                     
Int_t nchMgg = 100;                                                                                                       
   
Int_t nmassbin = 0;

Double_t ptMaxxM = 0; 
Double_t ptMinxM = 0; 


void sitest(){

  for(int i=0; i<ptbin+1; i++){    // WHY ?

    fsystfluct[i] = sqrt( fsystLam[i]*fsystLam[i] +  fsystGam[i]*fsystGam[i] +  fsystMass[i]*fsystMass[i] );

    Double_t fsysttotal =  sqrt(  fsystfluct[i]* fsystfluct[i] +  fsyserr* fsyserr + fsysINEL*fsysINEL );

    if(i==0) ptedges[i]=0;
    else {
      ptedges[i] = ptedges_Sigma[i-1];
      ptedges[ptbin+1] = ptedges_Sigma[ptbin];
    }
  } 
  


  int PARBINS_temp;
  PARBINS_temp = POLdegree+1+4;
  const int PARBINS=PARBINS_temp;
  int offset=0;
  offset=4;
  double bkg_params[POLdegree+1];
    
  TString *polname=new TString();
  if(POLdegree==1) polname->Append("pol1(0)");
  if(POLdegree==2) polname->Append("pol2(0)");
  if(POLdegree==3) polname->Append("pol3(0)");
  if(POLdegree==4) polname->Append("pol4(0)");
  if(POLdegree==5) polname->Append("pol5(0)");
    
  TF1 *myPol[ptbin];

  FitRange[0] = FitRangeSigma0[0]; 
  FitRange[1] = FitRangeSigma0[1];

  for(int i=0; i<ptbin; i++){
    TString *name = new TString("myPol");
    *name +=i+1;
    myPol[i] = new TF1(name->Data(),polname->Data(),0,2);
    myPol[i]->SetLineColor(3);
    if( i >= BinPtRange ) {   
      FitRange[0] = FitRangeLSigma0[0];  FitRange[1] = FitRangeLSigma0[1];    
    }
    myPol[i]->SetRange(FitRange[0],FitRange[1]);
  }

 Int_t iFlukaGeant = 0;
  for(int i = 0; i< 8;  i++){
    if( iFlukaGeant == 1 ){      flukageantA[i] =  flukageantAbase[i] ;           flukageantP[i] = 0.99 ;    }
    else if( iFlukaGeant == 0 ){      flukageantA[i] = 1 ;      flukageantP[i] = 1 ;    }
  }

  TFile *fmix  = TFile::Open("./himc10bcdefp4-mc-woPID.root");
  TFile *file0 = TFile::Open("./HISTOS-MC-4nov22.root");
  TDirectoryFile *MyTask = (TDirectoryFile*)file0->Get("MyTask");
  TList *Histos = (TList*)MyTask->Get("Histos");
  TH2F *hMassPt2D =   (TH2F*)Histos->FindObject("hLamGv0") ;
  TH2F *hMixPt2D =   (TH2F*)Histos->FindObject("hmLamGv0") ;

  TFile *fmc   = TFile::Open("./HISTOS-MC-4nov22.root");
  TDirectoryFile *MyTaskMC = (TDirectoryFile*)fmc->Get("MyTask");
  TList *HistosMC = (TList*)MyTaskMC->Get("Histos");
  TH1F * sig0gen =   (TH1F*)HistosMC->FindObject("hmc4piSig0Eta") ;
 TFile * file3 = TFile::Open("./mc_pythia8_sum.root");

  Double_t nPP = 326159711.   ;  // hev->GetBinContent(15) ; 
  Double_t nPPmc =  326159711. ; //  hevmc->GetBinContent(15); printf("npp mc = %f \n",nPPmc);  printf(" nPPmc =   %f \n", nPPmc ); 
  printf(" nPP=nPPmc = %f \n", nPP  );   
  cout << "open histogram : Rapidity : " << SelectRap << endl;
  cout << "  rapidity range 0 : " << raprange[0] << "  rapidity range 1 : "<< raprange[1] << endl;
    
    
  Int_t TypeSigma = 2;
  char hname[55];
  // This variable are for Mixed event background //

    
  fmix->Close();
  TH1D *hMassPt[8]; TH1D *hMixPt[8];  TH1D *hMixBGSub[8];   TH1D *hMixBGSub2[8];
  TH1D *hMass0Pt[8];
  TH1D *hMass1Pt[8];
  TH1D *hMass2Pt[8];

  TH1D *hMassPtTrueRec[8]; 
  TH1D *hMass1PtTrueRec[8]; 
  TH1D *hReEv[8]; 
  TH1D *hMiEv[8];   
  TH1D *hReEvSig[8];
  TH1D *hTrueBG[8];
  TH1D *hMiEvScale[8];
  TH1D *hSigTree[ptbin];   

  for(int i=0; i<ptbin; i++){
    hReEv[i]= new TH1D(Form("hReEv_%d",i),Form("Signal+Bg_%d",i),nchMgg,mSigmin,mSigmax);
    hMiEv[i]= new TH1D(Form("hMiEv_%d",i),Form("Mixed Bg_%d",i),nchMgg,mSigmin,mSigmax);
    hMiEvScale[i]= new TH1D(Form("hMiEvScale_%d",i),Form("Mixed Scaled Bg_%d",i),nchMgg,mSigmin,mSigmax);
    hReEvSig[i]= new TH1D(Form("hReEvSig_%d",i),Form("Signal after subtr mix Bg_%d",i),nchMgg,mSigmin,mSigmax);
    hSigTree[i]= new TH1D(Form("hSigTree_%d",i),Form("SigmaTTree_%d",i),nchMgg,mSigmin,mSigmax);
    hSigTree[i]->SetMarkerStyle(21);
    hSigTree[i]->SetMarkerSize(0.5);
    hSigTree[i]->SetLineColor(kRed-3);
    hSigTree[i]->SetMarkerColor(kRed-3);

    hMassPt[i]= new TH1D(Form("hMassPt_%d",i),Form("MassPt_%d",i),nchMgg,mSigmin,mSigmax);
  }
  TH1D *hSigMass = new TH1D("hSigMass","hSigMass", nchMgg,mSigmin,mSigmax );
  
  TString nameBinPt,  nameBinPtTrue,  titleBinPt, nameBinPtmix, titleBinPtmix, nameBinPtmix0,titleBinPtmix0 ,  nameBinPtmix1,titleBinPtmix1 ;


// Part 2. Signal + Background distributions from TTree output for Systematic Studies +++

 printf(" >>>>>>>>>>>>>>>>>>>>  Loop over extended pT bins, kk  \n \n");

  const int ptbinALL = 14; 
  double ptedges_SigmaALL[ptbinALL+1] = {0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.6, 2.0, 2.4, 2.8, 3.4, 4., 5.0, 8.0, 11. };
  double spectrum_mcALL[ptbinALL+1] = {0.};
  double spectrum_mcALLe[ptbinALL+1] = {0.};
// Main loop over pT bins, kk 
  for(int kk = 0; kk<ptbinloop; kk++ ){

    printf("\n pT bin %d \n", kk );

    nameBinPt   = Form("pt%d",iBinPt[kk]);
    hMassPt[kk] = hMassPt2D->ProjectionX(nameBinPt, iBinPt[kk], iBinPt[kk+1]-1);
    hMassPt[kk] ->SetTitle(titleBinPt);
    hMassPt[kk] ->SetLineColor(kBlue);
  }

  printf("...start can \n");
  TCanvas *can11 = new TCanvas("can11","can11",13,34,1200,600);
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat(1111);

  can11->Range(-1.25,-0.2625,11.25,2.3625);
  can11->SetFillColor(10);
  can11->SetBorderMode(0);
  can11->SetBorderSize(2);
  can11->SetFrameFillColor(0);
  can11->SetFrameBorderMode(0);
  can11->SetFrameBorderMode(0);
  can11->Divide(4,2);    
  for(int ll=0; ll<ptbinloop; ll++){
    can11->cd(ll+1);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    hMassPt[ll]->GetXaxis()->SetTitle("p_{T} (GeV/c)");      
      
    hMassPt[ll]->Draw("E");
  } 

 for(int kk = 0; kk<ptbinloop; kk++ ){

   if(kk == 2 ) { 
     FitRange[0] = 1.173 ;  FitRange[1] = 1.209; 
    SigmaCountBound[0]= 1.188 ; SigmaCountBound[0]= 1.197 ;
    POLdegree = 3 ;
   }

   if ( kk >= 4 && kk < 6) {
     FitRange[0] = FitRangeL2Sigma0[0]; 
     FitRange[1] = FitRangeL2Sigma0[1];
   }

   if ( kk == 4 ) {    SigmaCountBound[0]= 1.182 ; SigmaCountBound[0]= 1.120 ; }
   if ( kk == 4 ) {  FitRange[0] = 1.16 ;   FitRange[1] = 1.220 ; }
   if ( kk == 3 ) {  FitRange[0] = 1.164 ;   FitRange[1] = 1.2100 ; }

    printf(" >>>>>>>>>>>>>>>>>>>>  Part 4.  Background subtraction with different methods +++ \n \n" );

    if(UseMix ==kFALSE  &&  TypeCount == kFALSE ){
      // 1st Iteration : fit BG  
  
      
      reject=kTRUE;
            
      TF1 *myBkg =  myBkg =new TF1("myBkg",PolFunction,FitRange[0],FitRange[1],POLdegree+1);
      myBkg->SetParameters(6,0);
      myBkg->SetLineColor(2);
      hMassPt[kk]->Fit(myBkg,"NQ","",FitRange[0],FitRange[1]);// 1st fit
      //hMassPtP[kk]->Fit(myBkg,"NQ","",FitRange[0],FitRange[1]);// 1st fit
    
      if ( kk == 0 ){
	printf(" f1 \n");
      }

        
      TF1 *FullFit_t2;
      TF1 *FullFit;
            
      FullFit_t2=new TF1("FullFit_t2",VoigtplusPol,FitRange[0],FitRange[1],PARBINS+1);
      FullFit_t2->SetParameter(0,1.0);// mean
      FullFit_t2->SetParLimits(0,0.01,10.);// 0.0001 to 0.01 for BW only; 0.01 to 10 for true Voigtian
      FullFit_t2->SetParameter(1,1.19);// mean
      FullFit_t2->SetParLimits(1,1.175,1.21);
      FullFit_t2->SetParameter(2,.002);// abb Gaussian sigma
      FullFit_t2->SetParLimits(2,0.0001,.015);
      FullFit_t2->FixParameter(3,.0001);// BW width
      FullFit_t2->SetParName(0,"Norm");
      FullFit_t2->SetParName(1,"Mass");
      FullFit_t2->SetParName(2,"#sigma");
      FullFit_t2->SetParName(3,"#Gamma");
      FullFit_t2->SetParameter(4,0);
      FullFit_t2->SetParameter(5,0);
      FullFit_t2->SetParameter(6,0);
          FullFit_t2->SetParameter(7,0);
      FullFit_t2->SetParameter(8,0);
      FullFit_t2->SetParameter(9,0);
      
 
  if ( kk == 0 ){
     printf(" f2 \n");

    }
      // 2nd interation : fit Peak
            
      reject=kFALSE;
      double pars[PARBINS];
      double pars_e[PARBINS];
            
      for(int polbin=0; polbin<PARBINS-offset; polbin++) bkg_params[polbin] = myBkg->GetParameter(polbin);
      for(int polbin=0; polbin<PARBINS-offset; polbin++)
	FullFit_t2->FixParameter(offset+polbin, myBkg->GetParameter(polbin));
            
      hMassPt[kk]->Fit(FullFit_t2,"IMQ","",FitRange[0],FitRange[1]);// 2nd fit //IMQN+      

  if ( kk == 0 ){
     printf(" f3 \n");
    }



      for(int polbin=0; polbin<PARBINS; polbin++) pars[polbin] = FullFit_t2->GetParameter(polbin);
            
      ///////////////////////////////////////
      // 3rd iteration : fit Total
     
      FullFit=new TF1("FullFit",VoigtplusPol,FitRange[0],FitRange[1],PARBINS);
      FullFit->SetParLimits(0,0.01,10.);// 0.0001 to 0.01 for BW only; 0.01 to 10 for true Voigtian
      FullFit->SetParLimits(1,1.175,1.21); 
      if ( kk == 0) FullFit->SetParLimits(1,1.190,1.198); 

      FullFit->SetParLimits(2,0.001,.005);
        FullFit->FixParameter(3,0.0);
      FullFit->SetParName(0,"Norm");
      FullFit->SetParName(1,"Mass");
      FullFit->SetParName(2,"#sigma");
      FullFit->SetParName(3,"#Gamma");
      FullFit->SetParameter(4,0);
      FullFit->SetParameter(5,0);
      FullFit->SetParameter(6,0);
      FullFit->SetParameter(7,0);
      FullFit->SetParameter(8,0);
      FullFit->SetParameter(9,0);
            
      for(int polbin=0; polbin<PARBINS; polbin++) FullFit->SetParameter(polbin, pars[polbin]);
            
      hMassPt[kk]->Fit(FullFit,"IMEQ","",FitRange[0],FitRange[1]);// 3rd fit //IMEQ+

      if ( kk == 0 ){
	printf(" f4 \n");
      }
      
      for(int polbin=0; polbin<PARBINS-offset; polbin++) bkg_params[polbin] = FullFit->GetParameter(polbin+offset);
      for(int polbin=0; polbin<PARBINS; polbin++) pars[polbin] = FullFit->GetParameter(polbin);
      for(int polbin=0; polbin<PARBINS; polbin++) pars_e[polbin] = FullFit->GetParError(polbin);
      
      for(int polbin=0; polbin<POLdegree+1; polbin++) myBkg->SetParameter(polbin, bkg_params[polbin]);
      for(int polbin=0; polbin<POLdegree+1; polbin++) myPol[kk]->FixParameter(polbin, bkg_params[polbin]);
            
      // **  Get mass and width  ** //
      Sigma0Mass[kk] = FullFit->GetParameter(1);
      Sigma0Mass_e[kk] = FullFit->GetParError(1);
      Sigma0Width[kk] = FullFit->GetParameter(2);
      Sigma0Width_e[kk] = FullFit->GetParError(2);
            
      // **  Get raw yield and spectrum  ** //

      MissedYield[kk] = 1000.*(FullFit->Integral(0.,SigmaCountBound[0]) + FullFit->Integral(SigmaCountBound[1],100.)
			       - myBkg->Integral(0.,SigmaCountBound[0]) - myBkg->Integral(SigmaCountBound[1],100.));
            
      MissedYield_e[kk] = MissedYield[kk]*FullFit->GetParError(0)/FullFit->GetParameter(0);

      temp_yield = hMassPt[kk]->Integral(hMassPt[kk]->GetXaxis()->FindBin(SigmaCountBound[0]),
      					 hMassPt[kk]->GetXaxis()->FindBin(SigmaCountBound[1]-.0001));
      temp_yield1 = hMassPt[kk]->Integral(hMassPt[kk]->GetXaxis()->FindBin(SigmaCountBound[0]),
					  hMassPt[kk]->GetXaxis()->FindBin(SigmaCountBound[1]-.0001)); 

      temp_yield -= 1000.*(myBkg->Integral(SigmaCountBound[0],SigmaCountBound[1]));
      // this line is for calculating significance
      Double_t lowPt =  Sigma0Mass[kk] -  3.0 * Sigma0Width[kk];
      Double_t highPt =  Sigma0Mass[kk] +  3.0 * Sigma0Width[kk];

      temp_yield += MissedYield[kk];
      for(int massBin=hMassPt[kk]->GetXaxis()->FindBin(SigmaCountBound[0]);
	  massBin<=hMassPt[kk]->GetXaxis()->FindBin(SigmaCountBound[1]); massBin++){
	temp_yield_e += pow(hMassPt[kk]->GetBinError(massBin),2);
      }
            
      temp_yield_e += pow(MissedYield_e[kk],2);
      temp_yield_e = sqrt(temp_yield_e);
            
      temp_yield1 = hMassPt[kk]->Integral(hMassPt[kk]->GetXaxis()->FindBin(lowPt),
      			  hMassPt[kk]->GetXaxis()->FindBin(highPt-.0001));

      signifi[kk]=temp_yield/sqrt(temp_yield1);
      signifi_e[kk]=temp_yield/sqrt(temp_yield1)*sqrt((temp_yield_e/temp_yield)*
		    (temp_yield_e/temp_yield)+(sqrt(temp_yield1)/temp_yield1)*(sqrt(temp_yield1)/temp_yield1));

            
      printf("Low pT %f bin %f High pT %f bin %f  \n", lowPt, SigmaCountBound[0], highPt, SigmaCountBound[1]-.0001 );          
      cout<<"Signal/Bkg = "<<temp_yield/(1000.*myBkg->Integral(SigmaCountBound[0],SigmaCountBound[1]))<<endl;
      cout<<"Significance 1 = "<<temp_yield/sqrt(temp_yield1)<<endl;
      printf("TempYeild %f TempYield1 %f MissedY %f  \n",temp_yield, temp_yield1, MissedYield[kk] ); 
     cout<<"Included/Total = "<<(temp_yield-MissedYield[kk])/temp_yield<<endl;
      cout<<"yield = "<<temp_yield<<" +- "<<temp_yield_e<<  "Rel. Error Yield = " << temp_yield_e/temp_yield << endl;
      cout<<"MC input = "<<input_num[kk]<< endl ;
      cout<<"MC true rec = "<< rec_num[kk]  << endl;
      Double_t data_mc_ratio =      temp_yield/rec_num[kk];       
      cout<<"Data/rec_eff = "<< temp_yield/rec_num[kk]  << endl;  

      yield[kk] = temp_yield;
      yield_e[kk] = temp_yield_e;
      if ( kk == 7 ){
	printf(" f5 \n");
	hMassPt[kk] ->  Draw("E");

      }
    }  // end of UseMix FALSE and TypeCount False ==  Pol-Gauss fit
    printf(" end of Pol-Gauss fit with  UseMix FALSE and TypeCount False \n \n");

    // START Pol-Count PolBG and BinCounting !!!
  
  } // end of pt loop with kk variable == // END of SIG+BG fits
        
  printf("...start can \n");
  TCanvas *can = new TCanvas("can","can",13,34,1200,600);
  //    gStyle->SetOptFit(1111);
  gStyle->SetOptStat(1111);

  can->Range(-1.25,-0.2625,11.25,2.3625);
  can->SetFillColor(10);
  can->SetBorderMode(0);
  can->SetBorderSize(2);
  can->SetFrameFillColor(0);
  can->SetFrameBorderMode(0);
  can->SetFrameBorderMode(0);
  can->Divide(4,2);    
  for(int ll=0; ll<ptbinloop; ll++){
    can->cd(ll+1);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    hMassPt[ll]->GetXaxis()->SetTitle("p_{T} (GeV/c)");      
      
    hMassPt[ll]->Draw("E");

    myPol[ll]->DrawCopy("same");
  }
   
  can->Draw("");

     //+++
    // Part 3. Calculation of acceptance*efficiency ++++

   
   TCanvas *can1 = new TCanvas("can1","can1",13,34,800,600);
   sig0gen->Draw("E");

   TH1F * sig0genrec = (TH1F*)sig0gen->Clone("sig0gen");
  for(int kk = 0; kk<ptbin; kk++ ){

    Double_t ptMin = hMassPt2D->GetYaxis()->GetBinLowEdge(iBinPt[kk]);
    Double_t ptMax = hMassPt2D->GetYaxis()->GetBinUpEdge (iBinPt[kk+1]-1 );
      
    nameBinPt   = Form("pt%d",iBinPt[kk]);
    nameBinPtTrue   = Form("TrueRec pt%d",iBinPt[kk]);
    nameBinPtmix0   = Form("Mix pt0%d",iBinPt[kk]);
    nameBinPtmix1   = Form("Mix pt1%d",iBinPt[kk]);
    nameBinPtmix   = Form("Mix pt%d",iBinPt[kk]);
    titleBinPt  = Form("M_{#Lambda #gamma}, %.1f<p_{T}<%.1f GeV/c",ptMin,ptMax);
    titleBinPtmix0  = Form("Mix M_{#Lambda #gamma}, %.1f<p_{T}<%.1f GeV/c",ptMin,ptMax);
    titleBinPtmix1  = Form("Mix M_{#Lambda #gamma}, %.1f<p_{T}<%.1f GeV/c",ptMin,ptMax);
    titleBinPtmix  = Form("Mix M_{#Lambda #gamma}, %.1f<p_{T}<%.1f GeV/c",ptMin,ptMax);

   // input data in Pt bins //
    if( UseTree == kFALSE )  hMassPt[kk] = hMassPt2D->ProjectionX(nameBinPt, iBinPt[kk], iBinPt[kk+1]-1);
    hMassPt[kk] ->SetTitle(titleBinPt);
    hMassPt[kk] ->SetLineColor(kBlue);
    Double_t nevgenrec = 100000; // CHECK abb-8nov22

    int   nPtgen     = sig0gen->GetNbinsX();
    double nevgen = 0;  
    for(Int_t i = 1; i<= nPtgen ; i++){
      Double_t Ptmingen=sig0gen->GetXaxis()->GetBinLowEdge(i) ;
      Double_t Ptmaxgen=sig0gen->GetXaxis()->GetBinUpEdge(i) ;
      if( Ptmingen >= ptMin && Ptmaxgen <= ptMax ){
	Double_t  ngeni = sig0gen->GetBinContent(i) ;
	Double_t  ngenreci = sig0genrec->GetBinContent(i) ;
       	nevgen = nevgen + ngeni;
      }
     }
       
    input_num[kk] = nevgen;
    rec_num[kk] = nevgenrec;      
    Eff[kk] = rec_num[kk]/input_num[kk];
    Eff_e[kk] = sqrt( pow(sqrt(rec_num[kk])/input_num[kk],2) + pow(sqrt(input_num[kk])*rec_num[kk]/pow(input_num[kk],2),2));        

    double base = Eff[kk]*(raprange[1]-raprange[0])*(2*pt_points_e[kk])*nPP*feloss[kk]*fASigSig  ; //Events passing PhysSelectionCuts
    double base_mc =(raprange[1]-raprange[0])*(2*pt_points_e[kk])*nPPmc*fASigSig;


    printf("Efficiency = %f +- %f   MC_rec_count = %f  MC_input_count = %f \n", Eff[kk], Eff_e[kk], rec_num[kk], input_num[kk]);

    cout<<"Efficiency = "<<Eff[kk]<<"  +- "<<Eff_e[kk]<<"   MC_rec_count = "<<rec_num[kk]<< "   MC_input_count =  "<< input_num[kk]<< endl;

    cout<<"Rel Error of Efficiency = "<< Eff_e[kk] / Eff[kk]<< endl; 

    printf(" base %f  base_mc %f \n", base, base_mc );

    printf(" **  Calculate effieciency : end **<<<<<<<<<<< \n");
      
    // Part 3. Calculation of acceptance*efficiency ++++

   
   TH1F * sig0genrec = (TH1F*)sig0gen->Clone("sig0gen");
  
    nameBinPt   = Form("pt%d",iBinPt[kk]);
    nameBinPtTrue   = Form("TrueRec pt%d",iBinPt[kk]);
    nameBinPtmix0   = Form("Mix pt0%d",iBinPt[kk]);
    nameBinPtmix1   = Form("Mix pt1%d",iBinPt[kk]);
    nameBinPtmix   = Form("Mix pt%d",iBinPt[kk]);
    titleBinPt  = Form("M_{#Lambda #gamma}, %.1f<p_{T}<%.1f GeV/c",ptMin,ptMax);
    titleBinPtmix0  = Form("Mix M_{#Lambda #gamma}, %.1f<p_{T}<%.1f GeV/c",ptMin,ptMax);
    titleBinPtmix1  = Form("Mix M_{#Lambda #gamma}, %.1f<p_{T}<%.1f GeV/c",ptMin,ptMax);
    titleBinPtmix  = Form("Mix M_{#Lambda #gamma}, %.1f<p_{T}<%.1f GeV/c",ptMin,ptMax);

   // input data in Pt bins //
    if( UseTree == kFALSE )  hMassPt[kk] = hMassPt2D->ProjectionX(nameBinPt, iBinPt[kk], iBinPt[kk+1]-1);
    hMassPt[kk] ->SetTitle(titleBinPt);
    hMassPt[kk] ->SetLineColor(kBlue);
    for(Int_t i = 1; i<= nPtgen ; i++){
      Double_t Ptmingen=sig0gen->GetXaxis()->GetBinLowEdge(i) ;
      Double_t Ptmaxgen=sig0gen->GetXaxis()->GetBinUpEdge(i) ;
      if( Ptmingen >= ptMin && Ptmaxgen <= ptMax ){
	Double_t  ngeni = sig0gen->GetBinContent(i) ;
	Double_t  ngenreci = sig0genrec->GetBinContent(i) ;
       	nevgen = nevgen + ngeni;
      }
     }
 
       
    input_num[kk] = nevgen;
    rec_num[kk] = nevgenrec;      
    Eff[kk] = rec_num[kk]/input_num[kk];
    Eff_e[kk] = sqrt( pow(sqrt(rec_num[kk])/input_num[kk],2) + pow(sqrt(input_num[kk])*rec_num[kk]/pow(input_num[kk],2),2));        
 printf("Efficiency = %f +- %f   MC_rec_count = %f  MC_input_count = %f \n", Eff[kk], Eff_e[kk], rec_num[kk], input_num[kk]);

    cout<<"Efficiency = "<<Eff[kk]<<"  +- "<<Eff_e[kk]<<"   MC_rec_count = "<<rec_num[kk]<< "   MC_input_count =  "<< input_num[kk]<< endl;

    cout<<"Rel Error of Efficiency = "<< Eff_e[kk] / Eff[kk]<< endl; 

    printf(" base %f  base_mc %f \n", base, base_mc );

    printf(" **  Calculate effieciency : end **<<<<<<<<<<< \n");
      spectrum[kk] = temp_yield;
      spectrum[kk] /= base;
            
      spectrum_e[kk] = pow(temp_yield_e/base,2);
      spectrum_e[kk] += pow(temp_yield/base*Eff_e[kk]/Eff[kk],2);
      spectrum_e[kk] = sqrt(spectrum_e[kk]);
      printf(" end of Pol-Gauss fit with  UseMix FALSE and TypeCount False \n \n");
 



    
  

  printf(" ============================== output ============================== \n \n \n " ); 

  gStyle->SetOptStat(0);
  TFile *fm68= new TFile("./mcPyt68Rebinned.root");
    TGraphErrors *mcSigma0Pyt6Rebin = (TGraphErrors*)fm68 ->Get("mcSigma0Pyt6Rebin");
    TGraphErrors *mcSigma0Pyt8Rebin = (TGraphErrors*)fm68 ->Get("mcSigma0Pyt8Rebin");
  mcSigma0Pyt8Rebin->SetLineColor(kMagenta);
  mcSigma0Pyt8Rebin->SetMarkerColor(kMagenta);

  printf(" +222222222222222++++++++++++++++++++++++++++++++ add ratio Sigma/Lambda(p_T) +++++++++++++++++++++ \n") ;
   
    if(1) return;
   
    TGraphErrors *gr_Eff = new TGraphErrors(ptbin,pt_points, Eff, pt_points_e, Eff_e);
    gr_Eff->SetMarkerStyle(20);
    gr_Eff->SetMinimum(1.e-6);
    gr_Eff->SetMaximum(.05);
    gr_Eff->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    gr_Eff->GetYaxis()->SetTitle("Efficiency");
    gr_Eff->SetTitle("#Sigma^{0}+cc Efficiency"); 

    TH1D   *h_Eff = new TH1D("h_Eff", "AcceptanceEfficiency", ptbin+1,ptedges);
    h_Eff->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_Eff->GetYaxis()->SetTitle("Acceptance*Efficiency*B.R.");


    TH1D   *h_GenTrue = new TH1D("h_GenTrue", "GenTrue", ptbin+1,ptedges);
    TH1D   *h_RecTrue = new TH1D("h_RecTrue", "RecTrue", ptbin+1,ptedges);
    h_RecTrue->SetMarkerStyle(21);
    h_RecTrue->SetMarkerSize(.7);
    h_RecTrue->SetTitle("Rec/True");
    //    h_RecTrue->SetTitle("pp #sqrt{s}= 7 TeV");
    h_RecTrue->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_RecTrue->GetYaxis()->SetTitle("counts");
    h_RecTrue->SetLineColor(kBlue-3);
    h_RecTrue->SetMarkerColor(kBlue-3);


    TH1D   *h_Ratio = new TH1D("h_Ratio", "Ratio", ptbin+1,ptedges);

    h_Ratio->SetMarkerStyle(21);
    h_Ratio->SetMarkerSize(.5);
    h_Ratio->SetTitle("Rec/True");
    //    h_Ratio->SetTitle("pp #sqrt{s}= 7 TeV");
    h_Ratio->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_Ratio->GetYaxis()->SetTitle("ratio");
    h_Ratio->SetLineColor(kBlue-3);
    h_Ratio->SetMarkerColor(kBlue-3);

    TH1D   *h_rls = new TH1D("h_Ratio", "Ratio", ptbin+1,ptedges);

    h_rls->SetMarkerStyle(21);
    h_rls->SetMarkerSize(.5);
    h_rls->SetTitle("Rec/True");
    //    h_Ratio->SetTitle("pp #sqrt{s}= 7 TeV");
    h_rls->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_rls->GetYaxis()->SetTitle("ratio");
    h_rls->SetLineColor(kBlue-3);
    h_rls->SetMarkerColor(kBlue-3);



    TH1D   *hInput = new TH1D("hInput", "hInput", ptbin+1,ptedges);
    for(int i=0; i<ptbin; i++) {
        h_Eff->SetBinContent(i+2, Eff[i]);
        h_Eff->SetBinError(i+2, Eff_e[i]);

        h_GenTrue->SetBinContent(i+2,  input_num[i]  );
        h_GenTrue->SetBinError(i+2,   sqrt( input_num[i] ) );

        h_RecTrue->SetBinContent(i+2,  rec_num[i]  );
        h_RecTrue->SetBinError(i+2,   sqrt( rec_num[i] ) );
 
	hInput->SetBinContent(i+2,  input_num[i]  );
        hInput->SetBinError(i+2,   sqrt( input_num[i] ) );
   }

    TH1D *h_Rawspectrum = new TH1D("h_Rawspectrum","Sigma0 raw spectrum",ptbin+1,ptedges);
    h_Rawspectrum->SetMarkerStyle(20);
    h_Rawspectrum->SetMarkerSize(1);
    //    h_Rawspectrum->SetTitle("#Sigma0 raw spectrum");
    h_Rawspectrum->SetTitle("pp #sqrt{s}= 7 TeV");
    h_Rawspectrum->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_Rawspectrum->GetYaxis()->SetTitle("counts");
    // h_Rawspectrum->GetYaxis()->SetTitle("1/dydp_{T} (GeV/c)^{-1}");
    h_Rawspectrum->SetMinimum(3.0e-6);
    h_Rawspectrum->GetXaxis()->SetLimits(0,5.6);

    TH1D *h_NormRawspectrum = new TH1D("h_NormRawspectrum","Sigma0 norm raw spectrum",ptbin+1,ptedges);
    h_NormRawspectrum->SetMarkerStyle(20);
    h_NormRawspectrum->SetMarkerSize(1);
    h_NormRawspectrum->SetTitle("pp #sqrt{s}= 7 TeV");
    h_NormRawspectrum->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_NormRawspectrum->GetYaxis()->SetTitle("counts");

    Double_t  xPyt6, xPyt8,  yPyt6, yPyt8, rat68;
    for(int i=0; i<ptbin; i++) {
 
	h_Rawspectrum->SetBinContent(i+2, yield[i]);
	h_Rawspectrum->SetBinError(i+2, yield_e[i]);

	h_NormRawspectrum->SetBinContent(i+2, yield[i]/nPP );
        h_NormRawspectrum->SetBinError(i+2, yield_e[i]/nPP );

	printf("yeild/npp %f  bin %d",  yield[i]/nPP, i );
	cout<<"    yeild/npp = "<< yield[i]/nPP << endl;
	
	ratioev[i] =  yield[i]/rec_num[i] ;
	ratioev_e[i] =	ratioev[i] * sqrt( (yield_e[i]/yield[i])*(yield_e[i]/yield[i]) + 1/rec_num[i] ) ;
        h_Ratio->SetBinContent(i+2, ratioev[i]);
        h_Ratio->SetBinError(i+2, ratioev_e[i]);
    }

    TCanvas *Rec = new TCanvas("Rec","Reconstructed",0,0,800,800);
    Rec->Divide (1,2) ;
    Rec->cd(1);
    h_RecTrue->Draw("E");
    TLegend *legl = new TLegend(0.6,0.5,0.9,0.9); // TLegend( X0, Y0, X1, Y1 )
    legl->AddEntry(h_RecTrue,"True","lp");
    legl->AddEntry( h_Rawspectrum,"After BG subtr.","lp");
    legl->Draw("same");
    h_Rawspectrum->Draw("sameE");
    Rec->cd(2);
    gStyle->SetOptFit(111);
    h_Ratio->Fit("pol0","","",1.1,8.);
    h_Ratio->Draw("E");

    TH1D *h_spectrum = new TH1D("h_spectrum","Sigma0 raw spectrum",ptbin+1,ptedges);
    h_spectrum->SetMarkerStyle(21);
    h_spectrum->SetMarkerSize(.5);
    h_spectrum->SetTitle("pp #sqrt{s}= 7 TeV");
    h_spectrum->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_spectrum->GetYaxis()->SetTitle("1/N_{E}d^{2}N/dydp_{T} (GeV/c)^{-1}");
    h_spectrum->SetLineColor(kRed);
    h_spectrum->SetMarkerColor(kRed);
    
    for(int i=0; i<ptbin; i++) {
        h_spectrum->SetBinContent(i+2, spectrum[i]);
        h_spectrum->SetBinError(i+2, spectrum_e[i]);
    }

    TCanvas *CStat = new TCanvas("CSTAT","StatErrors",0,0,800,800);
    CStat->Divide (2,2) ;
    CStat->cd (1) ;
    h_Rawspectrum->Draw("E");
     CStat->cd (2) ;  
    TH1D *h_staterr = new TH1D("h_staterr","Sigma0 Stat Errors",ptbin+1,ptedges);
    h_staterr->SetMarkerStyle(21);
    h_staterr->SetMarkerSize(.5);
    h_staterr->SetTitle("pp #sqrt{s}= 7 TeV");
    h_staterr->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_staterr->GetYaxis()->SetTitle("1/N_{E}d^{2}N/dydp_{T} (GeV/c)^{-1}");
    h_staterr->SetLineColor(kRed);
    
    for(int i=0; i<ptbin; i++) {
        h_staterr->SetBinContent(i+2, yield_e[i]);
        h_staterr->SetBinError(i+2, sqrt(yield_e[i]) );
    }
    h_staterr->Draw("E");
    
     CStat->cd (3) ;  
    TH1D *h_relstaterr = new TH1D("h_relstaterr","Sigma0 Stat Errors",ptbin+1,ptedges);
    h_relstaterr->SetMarkerStyle(21);
    h_relstaterr->SetMarkerSize(.5);
    h_relstaterr->SetTitle("pp #sqrt{s}= 7 TeV");
    h_relstaterr->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_relstaterr->GetYaxis()->SetTitle("1/N_{E}d^{2}N/dydp_{T} (GeV/c)^{-1}");
    h_relstaterr->SetLineColor(kBlue);
    
    for(int i=0; i<ptbin; i++) {
        h_relstaterr->SetBinContent(i+2, yield_e[i]/yield[i] );
        h_relstaterr->SetBinError(i+2, sqrt(yield_e[i])/yield[i] );
    }
    gStyle->SetOptFit(1111);
    h_relstaterr->Fit("pol0","","",1.,10.);
    h_relstaterr->Draw("E");

    TH1D *h_spectrum_mc = new TH1D("h_spectrum_mc","Sigma0 mc raw spectrum",ptbin+1,ptedges);
    h_spectrum_mc->SetMarkerStyle(20);
    h_spectrum_mc->SetMarkerSize(.5);
    // h_spectrum_mc->SetTitle("#Sigma0 spectrum");
    h_spectrum_mc->SetTitle("pp #sqrt{s}= 7 TeV");
    h_spectrum_mc->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_spectrum_mc->GetYaxis()->SetTitle("1/N_{E}d^{2}N/dydp_{T} (GeV/c)^{-1}");
    h_spectrum_mc->SetLineColor(kRed);
    
    for(int i=0; i<ptbin; i++) {
        h_spectrum_mc->SetBinContent(i+2, spectrum_mc[i]);
        h_spectrum_mc->SetBinError(i+2, spectrum_mce[i]);
    }
    // h_spectrum_mc->Draw();
        
    TGraphErrors *SigMass = new TGraphErrors(ptbin,pt_points, Sigma0Mass, pt_points_e, Sigma0Mass_e);
    SigMass->SetMarkerStyle(22);
    SigMass->SetMinimum(1.19);
    SigMass->SetMaximum(1.198);
    SigMass->GetXaxis()->SetTitle("p_{T} ");
    SigMass->GetYaxis()->SetTitle("(GeV/c^{2})");
    SigMass->SetTitle("#Sigma^{0} mass vs p_{T}");
    //    SigMass->Draw("");
    
    TGraphErrors *SigWidth = new TGraphErrors(ptbin,pt_points, Sigma0Width, pt_points_e, Sigma0Width_e);
    SigWidth->SetMarkerStyle(23);
    SigWidth->SetMinimum(0.);
    SigWidth->SetMaximum(.005);
    SigWidth->GetXaxis()->SetTitle("p_{T} ");
    SigWidth->GetYaxis()->SetTitle("(GeV/c)");
    SigWidth->SetTitle("#Sigma^{0} width vs p_{T}");
    //    SigWidth->Draw("");
    
    TGraphErrors *SigSignifi = new TGraphErrors(ptbin,pt_points, signifi, pt_points_e, signifi_e);
    SigSignifi->SetMarkerStyle(20);
    SigSignifi->SetMinimum(0.);
    SigSignifi->SetMaximum(15.);
    SigSignifi->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    SigSignifi->GetYaxis()->SetTitle("Sig/#sqrt{Sig+BG}");
    SigSignifi->SetTitle("#Sigma^{0} significance vs p_{T}");
    gStyle->SetOptStat(0);
    //      SigSignifi->Draw("");
    //      h_Rawspectrum->Draw("E");

    TH1D *h_SigMass  = new TH1D("h_SigMass","Mass in Pt bins",ptbin+1,ptedges);
    TH1D *h_SigWidth = new TH1D("h_SigWidth","Width of mass in Pt bins",ptbin+1,ptedges);
    TH1D *h_SigSignifi = new TH1D("h_SigSignifi","Significance",ptbin+1,ptedges);
    h_SigSignifi->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_SigSignifi->GetYaxis()->SetTitle("counts");

    TH1D *h_SigYields = new TH1D("h_SigYields","Yields",ptbin+1,ptedges);
    h_SigYields->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_SigYields->GetYaxis()->SetTitle("counts");


    for(int i=0; i<ptbin; i++) {
      h_SigMass->SetBinContent(i+2, Sigma0Mass[i]);
      h_SigMass->SetBinError(i+2,  Sigma0Mass_e[i]);

      h_SigWidth->SetBinContent(i+2, Sigma0Width[i]);
      h_SigWidth->SetBinError(i+2,  Sigma0Width_e[i]);

      h_SigSignifi->SetBinContent(i+2, signifi[i]);
      h_SigSignifi->SetBinError(i+2,   signifi_e[i]);

      h_SigYields->SetBinContent(i+2, yield[i]);
      h_SigYields->SetBinError(i+2,  yield_e[i]);

    }
    h_SigMass->GetXaxis()->SetTitle("p_{T} (GeV/c)  ");
    h_SigMass->GetYaxis()->SetTitle("M (GeV/c^{2})  ");

    h_SigWidth->GetXaxis()->SetTitle("p_{T},(GeV/c)  ");
    h_SigWidth->GetYaxis()->SetTitle("#sigma_{M},(GeV/c^{2})  ");

    TCanvas *InfoHisto = new TCanvas("InfoHisto","Raw yield, Mass, Width, Significance",0,0,800,800);
    InfoHisto->Divide(2,2);
    InfoHisto->cd(1);
    h_Rawspectrum->Draw("E");
    InfoHisto->cd(2);
    SigSignifi->Draw("");   
    InfoHisto->cd(3);
    SigMass->Draw("");
    InfoHisto->cd(4);
    SigWidth->Draw("");
    
    TH1D *h_specsyst = new TH1D("h_specsyst","Sigma0 spectrum, systematic errors",ptbin+1,ptedges);
    h_specsyst->SetMarkerStyle(21);
    h_specsyst->SetMarkerSize(.5);
    //    h_specsyst ->SetTitle("#Sigma0 spectrum");
    h_specsyst ->SetTitle("pp #sqrt{s}= 7 TeV");
    h_specsyst->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_specsyst->SetTitle("1/N_{E}d^{2}N/dydp_{T} (GeV/c)^{-1}");
    h_specsyst->SetLineColor(kRed);
    h_specsyst->SetMarkerColor(kRed);
    
    for(int i=0; i<ptbin; i++) {
        h_spectrum->SetBinContent(i+2, spectrum[i]);
        h_spectrum->SetBinError(i+2, spectrum_e[i]);
    }

    
    TH1D *h_spectotal = new TH1D("h_spectotal","Sigma0 spectrum, total erros",ptbin+1,ptedges);
    h_spectotal->SetMarkerStyle(21);
    h_spectotal->SetMarkerSize(1.1);
    //    h_spectotal ->SetTitle("#Sigma0 spectrum");
    h_spectotal->SetTitle("pp #sqrt{s}= 7 TeV");
    h_spectotal ->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_spectotal->GetYaxis()->SetTitle("1/N_{inel}d^{2}N/dydp_{T} (GeV/c)^{-1}");
    h_spectotal->SetLineColor(kRed);
    h_spectotal->SetMarkerColor(kRed);
    

    TH1D *ratio = new TH1D("ratio","Rec./Sim.",ptbin+1,ptedges);
    ratio->SetLineColor(1);
    
    TH1D *hratio68 = new TH1D("hratio68","Sim-Pyt6/Sim-Pyt8",ptbin+1,ptedges);
    hratio68->SetLineColor(1);
    
    TH1D *hratio68dtmc = new TH1D("hratio68dtmc","Sig0 DTPyt8.",ptbin+1,ptedges);
    hratio68dtmc->SetLineColor(kMagenta);
    
    TH1D *ratiototal = new TH1D("ratiototal","DT/Pyt6",ptbin+1,ptedges);
    ratio->SetLineColor(1);

    TH1D *ratioStatDtMC = new TH1D("ratioStatDtMC","Rec-stat-err./Sim.",ptbin+1,ptedges);
    ratio->SetLineColor(1);

    Double_t ratiomc68[8], ratiomc68_e[8]; 
    Double_t ratiomc_e[8], ratiomc_esyst[8];
    for(int i=0; i<ptbin; i++) {

      mcSigma0Pyt6Rebin->GetPoint(i,xPyt6,yPyt6);
      mcSigma0Pyt8Rebin->GetPoint(i,xPyt8,yPyt8);
      rat68 = yPyt6/yPyt8 ;
     
      ratiomc68[i] = ratiomc[i] * rat68;
      ratiomc68_e[i] =  ratiomc_e[i]* rat68;
    
      ratio->SetBinContent(i+2, ratiomc[i]);
      ratio->SetBinError(i+2, ratiomc_e[i]);

      fsystotal[i] =  sqrt( fsyserr*fsyserr + fsystfluct[i]*fsystfluct[i] ); 
      printf(" Fsyst Err %f \n ",fsystotal[i] );
      //      fsystotal[i] = 0.116;

      spectrum_esys[i] =   spectrum[i] *  fsystotal[i];
      spectrum_etot[i] = sqrt( spectrum_esys[i]*spectrum_esys[i]+spectrum_e[i]*spectrum_e[i]);  
      printf("Sigma0  %f  Total  Err %f Syst %f Stat %f Bin %d \n ", 
	     spectrum[i],  spectrum_etot[i], spectrum_esys[i],spectrum_e[i], i  );


      h_specsyst->SetBinContent(i+2, spectrum[i]);      
      h_specsyst->SetBinError(i+2, spectrum_esys[i]);

      h_spectotal->SetBinContent(i+2, spectrum[i]);  
      h_spectotal->SetBinError(i+2, spectrum_etot[i]);

      ratiomc_etot[i] = ratiomc[i]*sqrt( spectrum_etot[i]*spectrum_etot[i]/spectrum[i]/spectrum[i]
                                       +spectrum_mce[i]*spectrum_mce[i]/spectrum_mc[i]/spectrum_mc[i] );
  

      ratiomc_esyst[i] = ratiomc[i]*sqrt( spectrum_esys[i]*spectrum_esys[i]/spectrum[i]/spectrum[i]
                                       +spectrum_mce[i]*spectrum_mce[i]/spectrum_mc[i]/spectrum_mc[i] );
  

      ratiomc_e[i] = ratiomc[i]*sqrt( spectrum_e[i]*spectrum_e[i]/spectrum[i]/spectrum[i]
                                       +spectrum_mce[i]*spectrum_mce[i]/spectrum_mc[i]/spectrum_mc[i] );


      ratiototal->SetBinContent(i+2, ratiomc[i]);
      //      ratiototal->SetBinError(i+2, ratiomc_etot[i]);
      ratiototal->SetBinError(i+2, ratiomc_esyst[i]);

      hratio68-> SetBinContent(i+2, rat68);
      hratio68-> SetBinError(i+2,rat68*0.01 );

      hratio68dtmc-> SetBinContent(i+2, ratiomc[i]*rat68);
      hratio68dtmc->SetBinError(i+2,  ratiomc_e[i]*rat68 );


      ratioStatDtMC->SetBinContent(i+2, ratiomc[i]);
      ratioStatDtMC->SetBinError(i+2, ratiomc_e[i]);

      printf("ratioStatDtMC %f +- %f bin %d \n \n \n ", ratiomc[i],   ratiomc_e[i], i );
 
    }
   double systerr_e[ptbin];
   double lamsysterr_e[ptbin];
   
   double systerrMat[ptbin];
   double systerrINEL[ptbin];
    for(int i=0; i<ptbin; i++) {
      lamsysterr_e[i] = 2.e-2;
      systerr_e[i] = 1.e-6;
      systerrMat[i] = fsyserr;
      systerrINEL[i] = fsysINEL ; 
   }

   TGraphErrors *h_SystLam = new TGraphErrors(ptbin,pt_points, fsystLam, pt_points_e, systerr_e);
   TGraphErrors *h_SystLam2 = new TGraphErrors(ptbin,pt_points, fsystLam, pt_points_e, lamsysterr_e);

   TGraphErrors *h_SystGam = new TGraphErrors(ptbin,pt_points, fsystGam, pt_points_e, systerr_e);
   TGraphErrors *h_SystGam2 = new TGraphErrors(ptbin,pt_points, fsystGam, pt_points_e, lamsysterr_e);


   TGraphErrors *h_SystMass = new TGraphErrors(ptbin,pt_points, fsystMass, pt_points_e, systerr_e);
   TGraphErrors *h_SystMat = new TGraphErrors(ptbin,pt_points, systerrMat, pt_points_e, systerr_e);
   TGraphErrors *h_SystINEL = new TGraphErrors(ptbin,pt_points, systerrINEL, pt_points_e, systerr_e);
   TGraphErrors *h_SystTot = new TGraphErrors(ptbin,pt_points, fsystotal, pt_points_e, systerr_e);

   h_SystTot->SetMarkerStyle(21);
   h_SystTot->SetMarkerSize(1.0);    
   h_SystTot->SetLineColor(1);
   h_SystTot->SetMarkerColor(1);
   h_SystTot ->GetXaxis()->SetTitle("p_{T} (GeV/c)");

   h_SystMass->SetMarkerStyle(23);
   h_SystMass->SetMarkerSize(0.8);
   h_SystMass->SetLineColor(kRed+3);
   h_SystMass->SetMarkerColor(kRed+3);
 

   h_SystGam->SetMarkerStyle(22);
   h_SystGam->SetMarkerSize(0.8);
   h_SystGam->SetLineColor(2);
   h_SystGam->SetMarkerColor(2);


   h_SystLam->SetMarkerStyle(26);
   h_SystLam->SetMarkerSize(0.8);
   h_SystLam->SetLineColor(3);
   h_SystLam->SetMarkerColor(3);

   h_SystMat->SetMarkerStyle(24);
   h_SystMat->SetMarkerSize(0.8);
   h_SystMat->SetLineColor(4);
   h_SystMat->SetMarkerColor(4);

   h_SystINEL->SetMarkerStyle(27);
   h_SystINEL->SetMarkerSize(0.8);
   h_SystINEL->SetLineColor(5);
   h_SystINEL->SetMarkerColor(5);
    TGraphErrors *gr_stat = new TGraphErrors(ptbin,pt_points,  spectrum , pt_points_e, spectrum_e );
    gr_stat->SetMarkerStyle(20);
    gr_stat->SetMinimum(1.e-6);
    gr_stat->SetMaximum(.05);
    gr_stat->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    gr_stat->GetYaxis()->SetTitle(" ");
    gr_stat->SetTitle(" ");

    ratiototal->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    ratiototal->GetYaxis()->SetTitle("ratio=rec./sim.");


    TCanvas *spectrum_tot = new TCanvas("spectrum_tot","Data/MC Sigma0 spectrum",0,0,800,1000);
    gStyle->SetOptStat(111);
    gStyle->SetOptFit(1111);
    spectrum_tot->Divide (1,2) ;
    spectrum_tot->cd(1);
    h_spectrum->Draw("E");
    spectrum_tot->cd(2);
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(11111);

    ratiototal->Draw("E1");
    ratioStatDtMC->Draw("E1,same");

    if (UseMix ==kTRUE  &&  TypeCount == kFALSE  ){        // Mix-Gauss  MIxed BG and Gauss Fit   
    printf(" +++++++++++++++++++++++++++++++++ add ratio Sigma/Lambda(p_T) +++++++++++++++++++++ \n") ;

  printf("++ writing to SigMac started   \n"  );

  TFile *fout0 = TFile::Open("SigMacDT.root","RECREATE");
 
  fout0->Close();


  TFile *fout = TFile::Open("SigMac2v15DTtail.root","RECREATE");

  hMassPt[0]->Write("pt0");  hMassPt[1]->Write("pt1");  hMassPt[2]->Write("pt2");  hMassPt[3]->Write("pt3"); 
  hMassPt[4]->Write("pt4");  hMassPt[5]->Write("pt5");  hMassPt[6]->Write("pt6");  hMassPt[7]->Write("pt7"); 
 printf(" write done \n");

  fout->Close();
    
 
}   // end of main

//________________________________________________________________________
double PolFunction(double *x, double *par){
    
    if (reject && x[0] > yieldRange[0] && x[0] < yieldRange[1]) {
        TF1::RejectPoint();
        return 0;
    }
    
    if(POLdegree==1) return par[0] + par[1]*x[0];
    else if(POLdegree==2) return par[0] + par[1]*x[0] + par[2]*pow(x[0],2);
    else if(POLdegree==3) return par[0] + par[1]*x[0] + par[2]*pow(x[0],2) + par[3]*pow(x[0],3);
    else if(POLdegree==4) return par[0] + par[1]*x[0] + par[2]*pow(x[0],2) + par[3]*pow(x[0],3) + par[4]*pow(x[0],4);
    else if(POLdegree==5) return par[0] + par[1]*x[0] + par[2]*pow(x[0],2) + par[3]*pow(x[0],3) + par[4]*pow(x[0],4) + par[5]*pow(x[0],5);
    else return 0;
    
}
//________________________________________________________________________

double BWFunction(double *x, double *par){
    return (par[0]*par[2]/(1000.*2.*3.1415926))/( pow(x[0]-par[1],2) + par[2]*par[2]/4.);
}
//________________________________________________________________________
double BWplusPol(double *x, double *par){
    return BWFunction(x,par) + PolFunction(x,&par[3]);
}
//________________________________________________________________________
double GausFunction(double *x, double *par){
    return (par[0]*.001)/(sqrt(2*3.14159*par[2]*par[2]))*exp(-pow((x[0]-par[1])/(sqrt(2)*par[2]),2));
}
//________________________________________________________________________
double GausplusPol(double *x, double *par){
    return GausFunction(x,par) + PolFunction(x,&par[3]);
}
//________________________________________________________________________
double  Voigtian(Double_t *x, Double_t *par)
{// code taken directly from RooVoigtian in ROOT
    Bool_t _doFast;
    _doFast = kFALSE;
    Double_t Norm = par[0];
    Double_t mean = par[1];
    Double_t s = fabs(par[2]);// sigma; Gaussian sigma
    Double_t w = fabs(par[3]);// Width; BW width
    
    Double_t arg = x[0] - mean;
    Double_t _invRootPi = 1./sqrt(atan2(0.,-1.));
    
    // return constant for zero width and sigma
    if (s==0. && w==0.) return 1.;
    
    // Breit-Wigner for zero sigma
    if (s==0.) return (Norm*1./(arg*arg+0.25*w*w));
    Double_t coef= -0.5/(s*s);
    
    // Gauss for zero width
    if (w==0.) return Norm*exp(coef*arg*arg);
    
    // actual Voigtian for non-trivial width and sigma
    Double_t c = 1./(sqrt(2.)*s);
    Double_t a = 0.5*c*w;
    Double_t u = c*arg;
    
    std::complex<Double_t> z(u,a) ;
    std::complex<Double_t> v(0.) ;
    
    v = RooMath::faddeeva(z);
    
    return Norm*c*_invRootPi*v.real();
    
    
}
//________________________________________________________________________
double VoigtplusPol(double *x, double *par){
    return  Voigtian(x,par) + PolFunction(x,&par[4]);
}

//=============================================================================
void PPRstyle()
{
  //////////////////////////////////////////////////////////////////////
  //
  // ROOT style macro for the TRD TDR
  //
  //////////////////////////////////////////////////////////////////////

  gStyle->SetPalette(1);
  gStyle->SetCanvasBorderMode(-1);
  gStyle->SetCanvasBorderSize(1);
  gStyle->SetCanvasColor(10);

  gStyle->SetFrameFillColor(10);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetFrameBorderMode(-1);
  //  gStyle->SetFrameLineWidth(1.2);
  gStyle->SetFrameLineColor(1);

  gStyle->SetHistFillColor(kBlack);
  gStyle->SetHistLineWidth(2);
  gStyle->SetHistLineColor(1);

  gStyle->SetPadColor(10);
  gStyle->SetPadBorderSize(1);
  gStyle->SetPadBorderMode(-1);

  gStyle->SetStatColor(10);
  gStyle->SetTitleColor(kBlack,"X");
  gStyle->SetTitleColor(kBlack,"Y");

  gStyle->SetLabelSize(0.024,"X");
  gStyle->SetLabelSize(0.024,"Y");
  gStyle->SetLabelSize(0.024,"Z");
  gStyle->SetTitleSize(0.024,"X");
  gStyle->SetTitleSize(0.024,"Y");
  gStyle->SetTitleSize(0.024,"Z");
  gStyle->SetTitleFont(42,"X");
  gStyle->SetTitleFont(42,"Y");
  gStyle->SetTitleFont(42,"X");
  gStyle->SetLabelFont(42,"X");
  gStyle->SetLabelFont(42,"Y");
  gStyle->SetLabelFont(42,"Z");
  gStyle->SetStatFont(42);

  gStyle->SetTitleOffset(1.0,"X");
  gStyle->SetTitleOffset(1.4,"Y");

//   gStyle->SetFillColor(kWhite);
  gStyle->SetTitleFillColor(kWhite);

  gStyle->SetOptDate(0);
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

}
 //<<<<<<<<<<<<<<<<<<<<<<<< 
Double_t Rapidity(Double_t pt, Double_t pz, Double_t m)
{
    //
    // calculates rapidity keeping the sign in case E == pz
    //
    
    Double_t energy = TMath::Sqrt(pt*pt+pz*pz+m*m);
    if (energy != TMath::Abs(pz))
    return 0.5*TMath::Log((energy+pz)/(energy-pz));
    
    //    Printf("W- mt=0");
    return TMath::Sign(1.e30,pz);
}
