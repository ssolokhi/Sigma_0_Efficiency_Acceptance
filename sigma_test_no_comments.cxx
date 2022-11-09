#include <iostream>
#include <TApplication.h>
#include <TMatrix.h>
#include <TROOT.h>
#include <TH1F.h>
double PolFunction(double *x, double *par);
double GausFunction(double *x, double *par);
double GausplusPol(double *x, double *par);
double Voigtian(Double_t *x, Double_t *par);
double VoigtplusPol(double *x, double *par);
double Rapidity(double pt, double pz, double m);

//////  ****  Setup : start ****  //////
int SelectRap = 0;
Double_t rapidity[6] = {2.0, 1.6, 1.4, 1.2, 1.0, 0.8};
Double_t dely = rapidity[SelectRap];
double raprange[2] = {-dely / 2, dely / 2};

const int ptbin = 8; // # of pt bins
bool UseMix = kFALSE;
bool TypeCount = kFALSE;
bool iMixPol = kTRUE; //   kFALSE;
int POLdegree = 3;    // 3; // pol-gauss Pol3 ,
bool UseTree = kTRUE; // kFALSE; // use input from TTRee or histos // v45-->kTRUE ONLY

Double_t mix1min = 1.170; // 1.175 ; - A lot depends on that values - MAIN
Double_t mix1max = 1.185; // 1.180 ;
Double_t mix2min = 1.205;
Double_t mix2max = 1.220;

double FitRangeSigma0[2] = {1.182, 1.205};

double FitRangeLSigma0[2] = {1.170, 1.215}; //-- data Bins0

double FitRangeL2Sigma0[2] = {1.16, 1.230};

int BinPtRange = 1;

bool iMixBGMC = kFALSE;
// OFFI
double SigmaCountBound[2] = {1.186, 1.200};

double SigmaMassCount[2] = {1.182, 1.202};

//////  ****  Setup : end ****  //////
bool reject;
double yieldRange[2] = {0};
double FitRange[2] = {0};
double ptedges[ptbin + 1 + 1] = {0};

Int_t iBinPt[40] = {12, 17, 21, 25, 29, 35, 41, 51, 81, 241, 281, 321, 361, 401};
double ptedges_Sigma[ptbin + 1] = {1.1, 1.6, 2.0, 2.4, 2.8, 3.4, 4., 5.0, 8.0};
double pt_points[ptbin] = {1.35, 1.8, 2.2, 2.6, 3.0, 3.8, 4.5, 6.5};
double pt_points_e[ptbin] = {.25, .2, .2, .2, .3, .3, 0.5, 1.5};
double Sigma0Mass[ptbin] = {0};
double Sigma0Width[ptbin] = {0};
double Sigma0Mass_e[ptbin] = {0};
double Sigma0Width_e[ptbin] = {0};
double MissedYield[ptbin] = {0};
double MissedYield_e[ptbin] = {0};
double yield[ptbin] = {0};
double yield_e[ptbin] = {0};
double spectrum[ptbin] = {0};
double spectrum_e[ptbin] = {0};
double spectrum_esys[ptbin] = {0};
double spectrum_etot[ptbin] = {0};

double ratioev[ptbin] = {0};
double ratioev_e[ptbin] = {0};

double ratioTrueSub[ptbin] = {0};
double ratioTrueSub_e[ptbin] = {0};

double signifi[ptbin] = {0};
double signifi_e[ptbin] = {0};
double GenEv[ptbin] = {0};
double GenEv_e[ptbin] = {0};
double spectrum_mc[ptbin] = {0};
double spectrum_mce[ptbin] = {0};
double ratiomc[ptbin] = {0};
double ratiomc_e[ptbin] = {0};
double ratiomc_etot[ptbin] = {0};

double EffA[ptbin] = {0};
double EffA_e[ptbin] = {0};
double EffP[ptbin] = {0};
double EffP_e[ptbin] = {0};
double EffAP[ptbin] = {0};
double EffAP_e[ptbin] = {0};
double EffAS[ptbin] = {0};
double EffAS_e[ptbin] = {0};
double EffPS[ptbin] = {0};
double EffPS_e[ptbin] = {0};
double EffAPS[ptbin] = {0};
double EffAPS_e[ptbin] = {0};

// **  Calculate effieciency :start  ** //
Double_t rec_num[40];
Double_t input_num[40];
Double_t recA_num[40];
Double_t inputA_num[40];
Double_t recP_num[40];
Double_t inputP_num[40];
Double_t recAP_num[40];
Double_t inputAP_num[40];
Double_t recAS_num[40];
Double_t inputAS_num[40];
Double_t recPS_num[40];
Double_t inputPS_num[40];
Double_t recAPS_num[40];
Double_t inputAPS_num[40];

double temp_yield = 0, temp_yield_e = 0, temp_yield1 = 0;

Double_t feloss[8] = {1, 1, 1, 1, 1, 1, 1, 1};
Double_t fASigSig = 2.;
Double_t flukageantAbase[8] = {0.935, 0.95, 0.96, 0.966, 0.970, 0.975, 0.982, 1.0}; // USE !!!
Double_t flukageantA[8];
Double_t flukageantP[8];
Double_t etacorr = dely / rapidity[4]; // is not used

Double_t rTrueRec[8] = {0.79, 0.82, 0.85, 0.87, 0.91, 0.957, 1.1, 1.18};
Double_t Eff[ptbin] = {9.33415e-05, 0.000400704, 0.000873701, 0.00160246, 0.00249329, 0.00327019, 0.00408349, 0.00489897};
Double_t Eff_e[ptbin] = {5.70201e-06, 1.91066e-05, 3.89734e-05, 7.21558e-05, 0.000104931, 0.000179498, 0.000241532, 0.000339698};
Double_t fsystotal[ptbin] = {0};

Double_t fsyserr = 0.085; // Mat. budget pp @ 7 TeV for TWO photons!!!

Double_t fsysINEL = 0.; //  0.062;

Double_t fsystfluct[ptbin] = {0};
//  expo smoothed lambda
Double_t fsystLam[ptbin] = {0.0707, 0.0650, 0.0606, 0.0565, 0.0530, 0.0460, 0.0405, 0.0285};
Double_t fsystGam[ptbin] = {0.044, 0.044, 0.044, 0.044, 0.044, 0.044, 0.044, 0.044};
Double_t fsystMass[ptbin] = {0.052, 0.044, 0.039, 0.0325, 0.028, 0.023, 0.018, 0.009}; // v 442 - 5jun18
//=============================
Double_t npi0ev[100], enpi0ev[100], nSumpi0ev[100], enSumpi0ev[100];
Double_t mSigmax = 1.230;
Double_t mSigmin = 1.130;
Double_t mLammax = 1.122;
Double_t mLammin = 1.110;

Int_t nchMgg = 100;
Int_t nmassbin = 0;

Double_t ptMaxxM = 0;
Double_t ptMinxM = 0;

void sigmatest()
{
  TCanvas *canPyt6Pyt801 = new TCanvas("canPyt6Pyt801", "canPyt6Pyt801", 13, 34, 800, 800);

  for (int i = 0; i < ptbin + 1; i++)
  { // WHY ?
    fsystfluct[i] = sqrt(fsystLam[i] * fsystLam[i] + fsystGam[i] * fsystGam[i] + fsystMass[i] * fsystMass[i]);
    Double_t fsysttotal = sqrt(fsystfluct[i] * fsystfluct[i] + fsyserr * fsyserr + fsysINEL * fsysINEL);
    printf("totFluctSyst %f Lam %f Gam %f Mass %f  INEL %f i %d \n", fsystfluct[i], fsystLam[i], fsystGam[i], fsystMass[i], fsysINEL, i);
    printf("Tot   %f \n", fsysttotal);

    if (i == 0)
      ptedges[i] = 0;
    else
    {
      ptedges[i] = ptedges_Sigma[i - 1];
      ptedges[ptbin + 1] = ptedges_Sigma[ptbin];
    }
  }

  int PARBINS_temp;
  PARBINS_temp = POLdegree + 1 + 4;
  const int PARBINS = PARBINS_temp;
  int offset = 0;
  offset = 4;
  double bkg_params[POLdegree + 1];

  TString *polname = new TString();
  if (POLdegree == 1)
    polname->Append("pol1(0)");
  if (POLdegree == 2)
    polname->Append("pol2(0)");
  if (POLdegree == 3)
    polname->Append("pol3(0)");
  if (POLdegree == 4)
    polname->Append("pol4(0)");
  if (POLdegree == 5)
    polname->Append("pol5(0)");

  TF1 *myPol[ptbin];

  FitRange[0] = FitRangeSigma0[0];
  FitRange[1] = FitRangeSigma0[1];

  for (int i = 0; i < ptbin; i++)
  {
    TString *name = new TString("myPol");
    *name += i + 1;
    myPol[i] = new TF1(name->Data(), polname->Data(), 0, 2);
    myPol[i]->SetLineColor(3);
    if (i >= BinPtRange)
    {
      FitRange[0] = FitRangeLSigma0[0];
      FitRange[1] = FitRangeLSigma0[1];
    }
    myPol[i]->SetRange(FitRange[0], FitRange[1]);
  }

  Int_t iFlukaGeant = 0;
  for (int i = 0; i < 8; i++)
  {
    if (iFlukaGeant == 1)
    {
      flukageantA[i] = flukageantAbase[i];
      flukageantP[i] = 0.99;
    }
    else if (iFlukaGeant == 0)
    {
      flukageantA[i] = 1;
      flukageantP[i] = 1;
    }
  }

  TFile *file0 = TFile::Open("./higood_10bcdefp4-data-woPID.root");
  TFile *fmix = TFile::Open("./higood_10bcdefp4-data-woPID.root");
  TFile *fmc = TFile::Open("./himc10bcdefp4-mc-woPID.root");

  TFile *file3 = TFile::Open("./mc_pythia8_sum.root");

  if (file0 == NULL)
    return; // this is for empty directory
  else
    cout << "open data root file ---> DONE" << endl;

  if (fmc == NULL)
    return; // this is for empty directorya
  else
    cout << "open mc root file " << endl;

  TObjArray *histESD = (TObjArray *)file0->Get("PWGGA_CaloConv/CaloConv");
  TObjArray *histESDmc = (TObjArray *)fmc->Get("PWGGA_CaloConv/CaloConv");
  TObjArray *histESDmix = (TObjArray *)fmix->Get("PWGGA_CaloConv/CaloConv");

  cout << "read histogram from Data" << endl;

  TDirectoryFile *td = (TDirectoryFile *)file0->Get("PWGGA_CaloConv");
  TList *l = (TList *)td->Get("CaloConv");
  TList *list = (TList *)l->FindObject("ConvEventCuts_00000103");
  TH1F *h1 = (TH1F *)list->FindObject("ESD_EventCuts 00000103");

  Double_t nInput = h1->GetBinContent(1);
  Double_t nOffTrig = h1->GetBinContent(2);
  Double_t nnovtxcontr = h1->GetBinContent(3);
  Double_t nvtxZgt10 = h1->GetBinContent(4);
  Double_t npileup = h1->GetBinContent(5);
  Double_t nvtxZlt10 = h1->GetBinContent(7);

  Double_t nMinBias = nvtxZlt10 + nvtxZgt10 + npileup + nnovtxcontr;
  Double_t nSumm = nMinBias + nOffTrig;

  printf("\n novtx %f Zgt10 %f  pileup %f Zlt10 %f MInBias %f \n", nnovtxcontr, nvtxZgt10, nvtxZlt10, npileup, nMinBias);
  printf("nInput %f nSumm %f \n \n", nInput, nSumm);

  Double_t nNormEv = 1;

  // come 29may18 factor 1/0.852 is due to the normalization on INEL events Error ~+- 0.06
  nNormEv = (nvtxZlt10 + nvtxZlt10 * nnovtxcontr / (nvtxZlt10 + nvtxZgt10)) / 0.852; // / 0.852 ;
  Double_t nNormEvAdd = nvtxZlt10 * nnovtxcontr / (nvtxZlt10 + nvtxZgt10);

  cout << "Ev.Min.Bias= " << nMinBias << "  Norm= " << nNormEv << "  Zv<10= " << nvtxZlt10 / nMinBias << "  Zv>10= " << nvtxZgt10 / nMinBias << "  NoVtx= " << nnovtxcontr / nMinBias << "  PileUp= " << npileup / nMinBias << endl;
  printf("nInput %f nSumm %f \n \n", nInput, nSumm);

  printf("\n \n");

  printf("nvtxZlt10 %f nNormEvAdd %f   nNormEv = %f r1 %f r2 %f  \n", nvtxZlt10, nNormEvAdd, nNormEv, nvtxZlt10 / nNormEv, nNormEvAdd / nNormEv);
  printf(" done 1 \n");

  Double_t nPP = 326159711.;   // hev->GetBinContent(15) ;
  Double_t nPPmc = 326159711.; //  hevmc->GetBinContent(15); printf("npp mc = %f \n",nPPmc);  printf(" nPPmc =   %f \n", nPPmc );
  printf(" nPP=nPPmc = %f \n", nPP);

  TH2F *hMassPt2D_A = (TH2F *)histESD->FindObject(Form("hDtPSig0MvsPtRec%d", SelectRap));
  TH2F *hMassPt2D_P = (TH2F *)histESD->FindObject(Form("hDtASig0MvsPtRec%d", SelectRap));

  TH2F *hMassPt2D = (TH2F *)hMassPt2D_A->Clone("hMassPt2D");
  TH2F *hMassPt2Dp = (TH2F *)hMassPt2D_P->Clone("hMassPt2Dp");
  hMassPt2D->Add(hMassPt2Dp, 1.);

  TH2F *hMassPt2DTrueRec_A = (TH2F *)histESDmc->FindObject(Form("hMCASig0MvsPtRec%d", SelectRap));
  TH2F *hMassPt2DTrueRec_P = (TH2F *)histESDmc->FindObject(Form("hMCPSig0MvsPtRec%d", SelectRap));
  TH2F *hMassPt2DTrueRec = (TH2F *)hMassPt2DTrueRec_A->Clone("hMassPt2DTrueRec");
  TH2F *hMassPt2DpTrueRec = (TH2F *)hMassPt2DTrueRec_P->Clone("hMassPt2DpTrueRec");
  hMassPt2DTrueRec->Add(hMassPt2DpTrueRec, 1.);

  cout << "open histogram : Rapidity : " << SelectRap << endl;
  cout << "  rapidity range 0 : " << raprange[0] << "  rapidity range 1 : " << raprange[1] << endl;

  Int_t TypeSigma = 2;
  char hname[55];
  // This variable are for Mixed event background //
  sprintf(hname, "hDtASig0MvsPtMix%d", SelectRap);
  TH2F *hMixPt2D_A = (TH2F *)histESDmix->FindObject(hname);
  sprintf(hname, "hDtPSig0MvsPtMix%d", SelectRap);
  TH2F *hMixPt2D_P = (TH2F *)histESDmix->FindObject(hname);
  TH2F *hMixPt2D = (TH2F *)hMixPt2D_A->Clone("hMixPt2D");
  TH2F *hMixPt2p = (TH2F *)hMixPt2D_P->Clone("hMixPt2p");
  hMixPt2D->Add(hMixPt2p, 1.);

  fmix->Close();
  TH1D *hMassPt[8];
  TH1D *hMixPt[8];
  TH1D *hMixBGSub[8];
  TH1D *hMixBGSub2[8];
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

  for (int i = 0; i < ptbin; i++)
  {
    hReEv[i] = new TH1D(Form("hReEv_%d", i), Form("Signal+Bg_%d", i), nchMgg, mSigmin, mSigmax);
    hMiEv[i] = new TH1D(Form("hMiEv_%d", i), Form("Mixed Bg_%d", i), nchMgg, mSigmin, mSigmax);
    hMiEvScale[i] = new TH1D(Form("hMiEvScale_%d", i), Form("Mixed Scaled Bg_%d", i), nchMgg, mSigmin, mSigmax);
    hReEvSig[i] = new TH1D(Form("hReEvSig_%d", i), Form("Signal after subtr mix Bg_%d", i), nchMgg, mSigmin, mSigmax);
    hSigTree[i] = new TH1D(Form("hSigTree_%d", i), Form("SigmaTTree_%d", i), nchMgg, mSigmin, mSigmax);
    hSigTree[i]->SetMarkerStyle(21);
    hSigTree[i]->SetMarkerSize(0.5);
    hSigTree[i]->SetLineColor(kRed - 3);
    hSigTree[i]->SetMarkerColor(kRed - 3);

    hMassPt[i] = new TH1D(Form("hMassPt_%d", i), Form("MassPt_%d", i), nchMgg, mSigmin, mSigmax);
    hMass0Pt[i] = new TH1D(Form("hMass0Pt_%d", i), Form("Mass0Pt_%d", i), nchMgg, mSigmin, mSigmax);
    hMass1Pt[i] = new TH1D(Form("hMass1Pt_%d", i), Form("Mass1Pt_%d", i), nchMgg, mSigmin, mSigmax);
    hMass2Pt[i] = new TH1D(Form("hMass2Pt_%d", i), Form("Mass2Pt_%d", i), nchMgg, mSigmin, mSigmax);

    hMassPtTrueRec[i] = new TH1D(Form("hMassPtTrueRec_%d", i), Form("MassPtTrueRec_%d", i), nchMgg, mSigmin, mSigmax);
    hMass1PtTrueRec[i] = new TH1D(Form("hMass1PtTrueRec_%d", i), Form("Mass1PtTrueRec_%d", i), nchMgg, mSigmin, mSigmax);
  }
  TH1D *hSigMass = new TH1D("hSigMass", "hSigMass", nchMgg, mSigmin, mSigmax);

  TString nameBinPt, nameBinPtTrue, titleBinPt, nameBinPtmix, titleBinPtmix, nameBinPtmix0, titleBinPtmix0, nameBinPtmix1, titleBinPtmix1;

  // Part 2. Signal + Background distributions from TTree output for Systematic Studies +++
  // DATA from TTRee
  Double_t fCutLambdaIPNeg = -1;
  Double_t fCutLambdaIPPos = -1;
  Double_t fCutLambdaDCADaughters = 1000;
  Double_t fCutLambdaCosPA = -1;
  Double_t fCutTPCPIDNSigmas = 1e+6;
  Double_t fCutLeastNumberOfCrossedRows = 70;
  Double_t fCutLambdaMass[2] = {1, 2};
  Double_t fCutLambdaArmPt[2] = {0, 0};
  Double_t fCutLambdaArmAlpha[2] = {0, 0};
  Double_t fCutLambdaRadius = 0;
  Double_t fCutLambdaEta = 1000;

  Double_t fCutGammaCosPA = -1;
  Double_t fCutGammaRadius[2] = {0, 0};
  Double_t fCutGammaArmPt[2] = {0, 0};
  Double_t fCutGammaArmAlpha[2] = {0, 0};
  Double_t fCutGammaZConv = 0;
  Double_t fCutGammaChi2 = 0;
  Double_t fCutGammaPt = 0;
  Double_t fCutGammaEta = 1000;

  Double_t fCutSigmaArmPt = 1000;
  Double_t fCutSigmaArmAlpha = 1000;
  Double_t fCutSigmaMass = 1.40;
  Double_t fCutSigma0Rap[2] = {0.5, 0.80};

  float fLambdaMod, fLambdaMass, fLambdaOnFly, fLambdaPx, fLambdaPt, fLambdaPz, fLambdaArmPt, fLambdaArmAlpha, fLambdaEta, fLambdaCosPointingAngle, fLambdaDCAtoPVPos, fLambdaDCAtoPVNeg, fLambdaDCADaughters, fLambdaRadius, fLambdacTau, fLambdaK0sDiff, fLambdaPosClust, fLambdaNegClust, fLambdaRatRows, fLambdaChi2V0, fLambdaSigmasProton,
      fGammaMass, fGammaPx, fGammaPt, fGammaPz, fGammaCosPointingAngle, fGammaDCADaughters, fGammaDCAtoPVNeg, fGammaDCAtoPVPos, fGammaEta, fGammaArmPt, fGammaArmAlpha, fGammaZConv, fGammaChi2, fGammaRadius,
      fGammaDCAzToPrimVtx,
      fGammaPtPos, fGammaPtNeg, fGammaNsigePos, fGammaNsigeNeg, fGammaNsigPiPos, fGammaNsigPiNeg,
      fGammaNsigKPos, fGammaNsigKNeg,
      fGammaPsi, fGammaRatRows,
      fSigmaMass, fSigmaPt, fSigmaArmPt, fSigmaArmAlpha,
      fLambdaTPx, fLambdaTPt, fLambdaTPz,
      fGammaTPx, fGammaTPt, fGammaTPz,
      fSigmaTPx, fSigmaTPt, fSigmaTPz,
      fLambdaTPxRec, fLambdaTPtRec, fLambdaTPzRec,
      fGammaTPxRec, fGammaTPtRec, fGammaTPzRec,
      fSigmaTPxRec, fSigmaTPtRec, fSigmaTPzRec,

      fLambdaTMassRec, fGammaTMassRec, fSigmaTMassRec,

      fCentr;
  float fSigmaRapidity;

  TDirectoryFile *dir = (TDirectoryFile *)file0->Get("PWGGA_CaloConv/CaloConv");
  TList *tDataList = (TList *)dir->FindObject("TList");
  TTree *tTreeEvent = (TTree *)tDataList->FindObject("tTreeEvent");

  tTreeEvent->SetBranchAddress("fLambdaMod", &fLambdaMod);
  tTreeEvent->SetBranchAddress("fLambdaMass", &fLambdaMass);
  tTreeEvent->SetBranchAddress("fLambdaOnFly", &fLambdaOnFly);
  tTreeEvent->SetBranchAddress("fLambdaPx", &fLambdaPx);
  tTreeEvent->SetBranchAddress("fLambdaPt", &fLambdaPt);
  tTreeEvent->SetBranchAddress("fLambdaPz", &fLambdaPz);
  tTreeEvent->SetBranchAddress("fLambdaArmPt", &fLambdaArmPt);
  tTreeEvent->SetBranchAddress("fLambdaArmAlpha", &fLambdaArmAlpha);
  tTreeEvent->SetBranchAddress("fLambdaEta", &fLambdaEta);
  tTreeEvent->SetBranchAddress("fLambdaCosPointingAngle", &fLambdaCosPointingAngle);
  tTreeEvent->SetBranchAddress("fLambdaDCAtoPVPos", &fLambdaDCAtoPVPos);
  tTreeEvent->SetBranchAddress("fLambdaDCAtoPVNeg", &fLambdaDCAtoPVNeg);
  tTreeEvent->SetBranchAddress("fLambdaDCADaughters", &fLambdaDCADaughters);
  tTreeEvent->SetBranchAddress("fLambdaRadius", &fLambdaRadius);
  tTreeEvent->SetBranchAddress("fLambdacTau", &fLambdacTau);
  tTreeEvent->SetBranchAddress("fLambdaK0sDiff", &fLambdaK0sDiff);
  tTreeEvent->SetBranchAddress("fLambdaPosClust", &fLambdaPosClust);
  tTreeEvent->SetBranchAddress("fLambdaNegClust", &fLambdaPosClust);
  tTreeEvent->SetBranchAddress("fLambdaRatRows", &fLambdaRatRows);
  tTreeEvent->SetBranchAddress("fLambdaChi2V0", &fLambdaChi2V0);
  tTreeEvent->SetBranchAddress("fLambdaSigmasProton", &fLambdaSigmasProton);
  tTreeEvent->SetBranchAddress("fGammaMass", &fGammaMass);
  tTreeEvent->SetBranchAddress("fGammaPx", &fGammaPx);
  tTreeEvent->SetBranchAddress("fGammaPt", &fGammaPt);
  tTreeEvent->SetBranchAddress("fGammaPz", &fGammaPz);
  tTreeEvent->SetBranchAddress("fGammaCosPointingAngle", &fGammaCosPointingAngle);
  tTreeEvent->SetBranchAddress("fGammaDCADaughters", &fGammaDCADaughters);
  tTreeEvent->SetBranchAddress("fGammaDCAtoPVNeg", &fGammaDCAtoPVNeg);
  tTreeEvent->SetBranchAddress("fGammaDCAtoPVPos", &fGammaDCAtoPVPos);
  tTreeEvent->SetBranchAddress("fGammaEta", &fGammaEta);
  tTreeEvent->SetBranchAddress("fGammaArmPt", &fGammaArmPt);
  tTreeEvent->SetBranchAddress("fGammaArmAlpha", &fGammaArmAlpha);
  tTreeEvent->SetBranchAddress("fGammaZConv", &fGammaZConv);
  tTreeEvent->SetBranchAddress("fGammaChi2", &fGammaChi2);
  tTreeEvent->SetBranchAddress("fGammaRadius", &fGammaRadius);
  tTreeEvent->SetBranchAddress("fGammaDCAzToPrimVtx", &fGammaDCAzToPrimVtx);
  tTreeEvent->SetBranchAddress("fGammaPtPos", &fGammaPtPos);
  tTreeEvent->SetBranchAddress("fGammaPtNeg", &fGammaPtNeg);
  tTreeEvent->SetBranchAddress("fGammaNsigePos", &fGammaNsigePos);
  tTreeEvent->SetBranchAddress("fGammaNsigeNeg", &fGammaNsigeNeg);
  tTreeEvent->SetBranchAddress("fGammaNsigPiPos", &fGammaNsigPiPos);
  tTreeEvent->SetBranchAddress("fGammaNsigPiNeg", &fGammaNsigPiNeg);
  tTreeEvent->SetBranchAddress("fGammaNsigKPos", &fGammaNsigKPos);
  tTreeEvent->SetBranchAddress("fGammaNsigKNeg", &fGammaNsigKNeg);
  tTreeEvent->SetBranchAddress("fGammaPsi", &fGammaPsi);
  tTreeEvent->SetBranchAddress("fGammaRatRows", &fGammaRatRows);
  tTreeEvent->SetBranchAddress("fSigmaMass", &fSigmaMass);
  tTreeEvent->SetBranchAddress("fSigmaPt", &fSigmaPt);
  tTreeEvent->SetBranchAddress("fSigmaArmPt", &fSigmaArmPt);
  tTreeEvent->SetBranchAddress("fSigmaArmAlpha", &fSigmaArmAlpha);
  tTreeEvent->SetBranchAddress("fLambdaTPx", &fLambdaTPx);
  tTreeEvent->SetBranchAddress("fLambdaTPt", &fLambdaTPt);
  tTreeEvent->SetBranchAddress("fLambdaTPz", &fLambdaTPz);
  tTreeEvent->SetBranchAddress("fGammaTPx", &fGammaTPx);
  tTreeEvent->SetBranchAddress("fGammaTPt", &fGammaTPt);
  tTreeEvent->SetBranchAddress("fGammaTPz", &fGammaTPz);
  tTreeEvent->SetBranchAddress("fSigmaTPx", &fSigmaTPx);
  tTreeEvent->SetBranchAddress("fSigmaTPt", &fSigmaTPt);
  tTreeEvent->SetBranchAddress("fSigmaTPz", &fSigmaTPz);
  tTreeEvent->SetBranchAddress("fLambdaTPxRec", &fLambdaTPxRec);
  tTreeEvent->SetBranchAddress("fLambdaTPtRec", &fLambdaTPtRec);
  tTreeEvent->SetBranchAddress("fLambdaTPzRec", &fLambdaTPzRec);
  tTreeEvent->SetBranchAddress("fGammaTPxRec", &fGammaTPxRec);
  tTreeEvent->SetBranchAddress("fGammaTPtRec", &fGammaTPtRec);
  tTreeEvent->SetBranchAddress("fGammaTPzRec", &fGammaTPzRec);
  tTreeEvent->SetBranchAddress("fSigmaTPxRec", &fSigmaTPxRec);
  tTreeEvent->SetBranchAddress("fSigmaTPtRec", &fSigmaTPtRec);
  tTreeEvent->SetBranchAddress("fSigmaTPzRec", &fSigmaTPzRec);
  tTreeEvent->SetBranchAddress("fLambdaTMassRec", &fLambdaTMassRec);
  tTreeEvent->SetBranchAddress("fGammaTMassRec", &fGammaTMassRec);
  tTreeEvent->SetBranchAddress("fSigmaTMassRec", &fLambdaTMassRec);

  // CutsApplication
  Double_t cutmin = -20; // -20; // -20; // -20; // -30 ;
  Double_t cutmax = 30;  // 5; //  10; // 20;  //  30 ;

  Double_t cutmin1 = 0.;
  Double_t cutmax1 = 200.0;
  Double_t cutmin2 = -1.e-8;
  Double_t cutmax2 = 1.e-8;
  Double_t cutmin3 = 0.0;
  Double_t cutmax3 = 1.20;
  Double_t cutmin4 = 0.0;
  Double_t cutmax4 = 30.;

  Double_t cutmin5 = 0.99;
  Double_t cutmax5 = 1.001;
  Double_t cutmin6 = 0.;
  Double_t cutmax6 = 140;
  Double_t cutmin7 = 0.;
  Double_t cutmax7 = 140;
  Double_t cutmin8 = 0.;
  Double_t cutmax8 = 2.;

  Double_t cutmin9 = 0.;
  Double_t cutmax9 = 160.;

  Double_t cutmin10 = 0.;
  Double_t cutmax10 = 2.0;
  Double_t cutmin11 = 0.;
  Double_t cutmax11 = 2.0;

  Double_t cutmin12 = -8.0;
  Double_t cutmax12 = 8.0;
  Double_t cutmin13 = -8.0;
  Double_t cutmax13 = 8.0;
  Double_t cutmin14 = -10.0;
  Double_t cutmax14 = 10.0;
  Double_t cutmin15 = -10.0;
  Double_t cutmax15 = 10.0;

  TH1D *hCutDT = new TH1D("hCutDT", "; ", 100, cutmin, cutmax);
  TH1D *hCutMC = new TH1D("hCutMC", "; ", 100, cutmin, cutmax);

  TH1D *hCutDT0 = new TH1D("hCutDT0", "; ", 100, cutmin, cutmax);
  TH1D *hCutMC0 = new TH1D("hCutMC0", "; ", 100, cutmin, cutmax);

  TH1D *hCutDT1 = new TH1D("hCutDT1", "; Lambda Radius, cm", 1000, cutmin1, cutmax1);
  TH1D *hCutMC1 = new TH1D("hCutMC1", "; Lambda Radius, cm", 1000, cutmin1, cutmax1);

  TH1D *hCutDT2 = new TH1D("hCutDT2", "; Lambda Mass, GeV/c^{2}", 100, cutmin2, cutmax2);
  TH1D *hCutMC2 = new TH1D("hCutMC2", "; Lambda Mass, GeV/c^{2}", 100, cutmin2, cutmax2);

  TH1D *hCutDT3 = new TH1D("hCutDT3", "; Lambda Ratio Rows", 100, cutmin3, cutmax3);
  TH1D *hCutMC3 = new TH1D("hCutMC3", "; Lambda Ratio Rows", 100, cutmin3, cutmax3);

  TH1D *hCutDT4 = new TH1D("hCutDT4", "; Lambda K0sDiff", 100, cutmin4, cutmax4);
  TH1D *hCutMC4 = new TH1D("hCutMC4", "; Lambda K0sDiff", 100, cutmin4, cutmax4);

  TH1D *hCutDT5 = new TH1D("hCutDT5", "; Lambda Cos(PointingAngle)", 100, cutmin5, cutmax5);
  TH1D *hCutMC5 = new TH1D("hCutMC5", "; Lambda Cos(PointingAngle)", 100, cutmin5, cutmax5);

  TH1D *hCutDT6 = new TH1D("hCutDT6", "; Lambda DCA to PV pos. track, cm", 100, cutmin6, cutmax6);
  TH1D *hCutMC6 = new TH1D("hCutMC6", "; Lambda DCA to PV pos. track, cm", 100, cutmin6, cutmax6);

  TH1D *hCutDT7 = new TH1D("hCutDT7", "; Lambda DCA to PV neg. track, cm", 100, cutmin7, cutmax7);
  TH1D *hCutMC7 = new TH1D("hCutMC7", "; Lambda DCA to PV neg. track, cm", 100, cutmin7, cutmax7);

  TH1D *hCutDT8 = new TH1D("hCutDT8", "; Lambda DCA of daughters, cm", 1000, cutmin8, cutmax8);
  TH1D *hCutMC8 = new TH1D("hCutMC8", "; Lambda DCA of daughters, cm", 1000, cutmin8, cutmax8);

  TH1D *hCutDT9 = new TH1D("hCutDT9", "; Lambda SigmasProton", 100, cutmin9, cutmax9);
  TH1D *hCutMC9 = new TH1D("hCutMC9", "; Lambda SigmasProton", 100, cutmin9, cutmax9);

  TH1D *hCutDT10 = new TH1D("hCutDT10", "; Lambda cTau ", 100, cutmin10, cutmax10);
  TH1D *hCutMC10 = new TH1D("hCutMC10", "; Lambda cTau", 100, cutmin10, cutmax10);

  TH1D *hCutDT11 = new TH1D("hCutDT11", "; Lambda Arm Alpha", 100, cutmin11, cutmax11);
  TH1D *hCutMC11 = new TH1D("hCutMC11", "; Lambda Arm Alpha ", 100, cutmin11, cutmax11);

  TH1D *hCutDT12 = new TH1D("hCutDT12", "; Lambda Arm pT", 100, cutmin12, cutmax12);
  TH1D *hCutMC12 = new TH1D("hCutMC12", "; Lambda Arm pT", 100, cutmin12, cutmax12);

  TH1D *hCutDT13 = new TH1D("hCutDT13", "; Lambda Mod", 20, -200, 200);
  TH1D *hCutMC13 = new TH1D("hCutMC13", "; Lambda Mod", 20, -200, 200);

  TH1D *hCutDT14 = new TH1D("hCutDT14", "; gamma N sigmas PiPos", 100, cutmin14, cutmax14);
  TH1D *hCutMC14 = new TH1D("hCutMC14", "; gamma N sigmas PiPos", 100, cutmin14, cutmax14);

  TH1D *hCutDT15 = new TH1D("hCutDT15", "; gamma N sigmas PiNeg", 100, cutmin15, cutmax15);
  TH1D *hCutMC15 = new TH1D("hCutMC15", "; gamma N sigmas PiNeg", 100, cutmin15, cutmax15);

  TH2D *h2CutDT14 = new TH2D("h2CutDT14", "gamma N sigmas; nSigma PiPos;  PiNeg", 40, cutmin14, cutmax14, 40, cutmin14, cutmax14);
  TH2D *h2CutMC14 = new TH2D("h2CutMC14", "gamma N sigmas; nSigma PiPos; PiNeg", 40, cutmin14, cutmax14, 40, cutmin14, cutmax14);

  TH2D *h2CutDT14mix = new TH2D("h2CutDT14mix", "gamma N sigmas; nSigma ePos;  eNeg", 40, cutmin14, cutmax14, 40, cutmin14, cutmax14);

  Double_t CutNsigPi = 0.0;
  Double_t CutNsigPiMAX = 0.;
  Double_t CutNsigPiDIFF = 0.;

  Double_t varvalue14old = 0.;

  for (Long64_t iEntry = 0; iEntry < tTreeEvent->GetEntries(); iEntry++)
  {
    tTreeEvent->GetEntry(iEntry);
    Double_t SigmaPz = fGammaPz + fLambdaPz;
    Double_t SigmaRap = Rapidity(fSigmaPt, SigmaPz, fSigmaMass);
    if (TMath::Abs(SigmaRap) > 0.5)
      continue;

    Double_t varvalue = fLambdaSigmasProton; // fGammaCosPointingAngle;
    Double_t varvalue1 = fLambdaRadius;
    Double_t varvalue2 = fLambdaMass;
    Double_t varvalue3 = fLambdaRatRows;
    Double_t varvalue4 = fLambdaK0sDiff;

    Double_t varvalue5 = fLambdaCosPointingAngle;
    Double_t varvalue6 = fLambdaDCAtoPVPos;
    Double_t varvalue7 = fLambdaDCAtoPVNeg;
    Double_t varvalue8 = fLambdaDCADaughters;

    Double_t varvalue9 = fLambdaSigmasProton;

    Double_t varvalue10 = fLambdacTau;
    Double_t varvalue11 = fLambdaArmAlpha;
    Double_t varvalue12 = fLambdaArmPt;
    Double_t varvalue13 = fLambdaMod;
    Double_t varvalue14 = fGammaNsigPiPos;
    Double_t varvalue15 = fGammaNsigPiNeg;
    Double_t varvalue16 = fGammaNsigePos;
    Double_t varvalue17 = fGammaNsigeNeg;

    if (!(varvalue > cutmin && varvalue <= cutmax))
      continue;
    for (int kk = 0; kk < ptbin; kk++)
    {
      if (fSigmaPt >= ptedges_Sigma[kk] && fSigmaPt < ptedges_Sigma[kk + 1])
      {
        hSigTree[kk]->Fill(fSigmaMass);
        if (UseTree == kTRUE)
        {
          hMassPt[kk]->Fill(fGammaPsi);
          hCutDT14->Fill(varvalue14);
          hCutDT15->Fill(varvalue15);
          h2CutDT14->Fill(varvalue14, varvalue15);
          h2CutDT14mix->Fill(varvalue16, varvalue17);
          varvalue14old = varvalue14;

          if (fLambdaMod == 100)
            hCutDT->Fill(varvalue);
          if (fLambdaMod == -100)
            hCutDT0->Fill(varvalue);

          hCutDT1->Fill(varvalue1);
          hCutDT2->Fill(varvalue2);
          hCutDT3->Fill(varvalue3);
          hCutDT4->Fill(varvalue4);
          hCutDT5->Fill(varvalue5);
          hCutDT6->Fill(varvalue6);
          hCutDT7->Fill(varvalue7);
          hCutDT8->Fill(varvalue8);
          hCutDT9->Fill(varvalue9);
          hCutDT10->Fill(varvalue10);
          hCutDT11->Fill(varvalue11);
          hCutDT12->Fill(varvalue12);
          hCutDT13->Fill(varvalue13);
        }
      }
    }
    hSigMass->Fill(fSigmaMass);
  }
  printf("TTree DATA done \n \n");

  // MC acceptance*effi corrections from true reconstructed Sigma0
  TDirectoryFile *dirMC = (TDirectoryFile *)fmc->Get("PWGGA_CaloConv/CaloConv");
  TList *tDataListMC = (TList *)dirMC->FindObject("TList");
  TTree *tTreeEventMC = (TTree *)tDataListMC->FindObject("tTreeEvent");
  tTreeEventMC->SetBranchAddress("fLambdaMod", &fLambdaMod);
  tTreeEventMC->SetBranchAddress("fLambdaMass", &fLambdaMass);
  tTreeEventMC->SetBranchAddress("fLambdaOnFly", &fLambdaOnFly);
  tTreeEventMC->SetBranchAddress("fLambdaPx", &fLambdaPx);
  tTreeEventMC->SetBranchAddress("fLambdaPt", &fLambdaPt);
  tTreeEventMC->SetBranchAddress("fLambdaPz", &fLambdaPz);
  tTreeEventMC->SetBranchAddress("fLambdaEta", &fLambdaEta);
  tTreeEventMC->SetBranchAddress("fLambdaCosPointingAngle", &fLambdaCosPointingAngle);
  tTreeEventMC->SetBranchAddress("fLambdaDCAtoPVPos", &fLambdaDCAtoPVPos);
  tTreeEventMC->SetBranchAddress("fLambdaDCAtoPVNeg", &fLambdaDCAtoPVNeg);
  tTreeEventMC->SetBranchAddress("fLambdaDCADaughters", &fLambdaDCADaughters);
  tTreeEventMC->SetBranchAddress("fLambdaRadius", &fLambdaRadius);
  tTreeEventMC->SetBranchAddress("fLambdacTau", &fLambdacTau);
  tTreeEventMC->SetBranchAddress("fLambdaK0sDiff", &fLambdaK0sDiff);
  tTreeEventMC->SetBranchAddress("fLambdaPosClust", &fLambdaPosClust);
  tTreeEventMC->SetBranchAddress("fLambdaNegClust", &fLambdaPosClust);
  tTreeEventMC->SetBranchAddress("fLambdaRatRows", &fLambdaRatRows);
  tTreeEventMC->SetBranchAddress("fLambdaChi2V0", &fLambdaChi2V0);
  tTreeEventMC->SetBranchAddress("fLambdaSigmasProton", &fLambdaSigmasProton);
  tTreeEventMC->SetBranchAddress("fGammaMass", &fGammaMass);
  tTreeEventMC->SetBranchAddress("fGammaPx", &fGammaPx);
  tTreeEventMC->SetBranchAddress("fGammaPt", &fGammaPt);
  tTreeEventMC->SetBranchAddress("fGammaPz", &fGammaPz);
  tTreeEventMC->SetBranchAddress("fGammaCosPointingAngle", &fGammaCosPointingAngle);
  tTreeEventMC->SetBranchAddress("fGammaDCADaughters", &fGammaDCADaughters);
  tTreeEventMC->SetBranchAddress("fGammaDCAtoPVNeg", &fGammaDCAtoPVNeg);
  tTreeEventMC->SetBranchAddress("fGammaDCAtoPVPos", &fGammaDCAtoPVPos);
  tTreeEventMC->SetBranchAddress("fGammaEta", &fGammaEta);
  tTreeEventMC->SetBranchAddress("fGammaArmPt", &fGammaArmPt);
  tTreeEventMC->SetBranchAddress("fGammaArmAlpha", &fGammaArmAlpha);
  tTreeEventMC->SetBranchAddress("fGammaZConv", &fGammaZConv);
  tTreeEventMC->SetBranchAddress("fGammaChi2", &fGammaChi2);
  tTreeEventMC->SetBranchAddress("fGammaRadius", &fGammaRadius);
  tTreeEventMC->SetBranchAddress("fGammaDCAzToPrimVtx", &fGammaDCAzToPrimVtx);
  tTreeEventMC->SetBranchAddress("fGammaPtPos", &fGammaPtPos);
  tTreeEventMC->SetBranchAddress("fGammaPtNeg", &fGammaPtNeg);
  tTreeEventMC->SetBranchAddress("fGammaNsigePos", &fGammaNsigePos);
  tTreeEventMC->SetBranchAddress("fGammaNsigeNeg", &fGammaNsigeNeg);
  tTreeEventMC->SetBranchAddress("fGammaNsigPiPos", &fGammaNsigPiPos);
  tTreeEventMC->SetBranchAddress("fGammaNsigPiNeg", &fGammaNsigPiNeg);
  tTreeEventMC->SetBranchAddress("fGammaNsigKPos", &fGammaNsigKPos);
  tTreeEventMC->SetBranchAddress("fGammaNsigKNeg", &fGammaNsigKNeg);
  tTreeEventMC->SetBranchAddress("fGammaPsi", &fGammaPsi);
  tTreeEventMC->SetBranchAddress("fGammaRatRows", &fGammaRatRows);
  tTreeEventMC->SetBranchAddress("fSigmaMass", &fSigmaMass);
  tTreeEventMC->SetBranchAddress("fSigmaRapidity", &fSigmaRapidity);
  tTreeEventMC->SetBranchAddress("fSigmaPt", &fSigmaPt);
  tTreeEventMC->SetBranchAddress("fSigmaArmPt", &fSigmaArmPt);
  tTreeEventMC->SetBranchAddress("fSigmaArmAlpha", &fSigmaArmAlpha);
  tTreeEventMC->SetBranchAddress("fLambdaTPx", &fLambdaTPx);
  tTreeEventMC->SetBranchAddress("fLambdaTPt", &fLambdaTPt);
  tTreeEventMC->SetBranchAddress("fLambdaTPz", &fLambdaTPz);
  tTreeEventMC->SetBranchAddress("fGammaTPx", &fGammaTPx);
  tTreeEventMC->SetBranchAddress("fGammaTPt", &fGammaTPt);
  tTreeEventMC->SetBranchAddress("fGammaTPz", &fGammaTPz);
  tTreeEventMC->SetBranchAddress("fSigmaTPx", &fSigmaTPx);
  tTreeEventMC->SetBranchAddress("fSigmaTPt", &fSigmaTPt);
  tTreeEventMC->SetBranchAddress("fSigmaTPz", &fSigmaTPz);
  tTreeEventMC->SetBranchAddress("fLambdaTPxRec", &fLambdaTPxRec);
  tTreeEventMC->SetBranchAddress("fLambdaTPtRec", &fLambdaTPtRec);
  tTreeEventMC->SetBranchAddress("fLambdaTPzRec", &fLambdaTPzRec);
  tTreeEventMC->SetBranchAddress("fGammaTPxRec", &fGammaTPxRec);
  tTreeEventMC->SetBranchAddress("fGammaTPtRec", &fGammaTPtRec);
  tTreeEventMC->SetBranchAddress("fGammaTPzRec", &fGammaTPzRec);
  tTreeEventMC->SetBranchAddress("fSigmaTPxRec", &fSigmaTPxRec);
  tTreeEventMC->SetBranchAddress("fSigmaTPtRec", &fSigmaTPtRec);
  tTreeEventMC->SetBranchAddress("fSigmaTPzRec", &fSigmaTPzRec);
  tTreeEventMC->SetBranchAddress("fLambdaTMassRec", &fLambdaTMassRec);
  tTreeEventMC->SetBranchAddress("fGammaTMassRec", &fGammaTMassRec);
  tTreeEventMC->SetBranchAddress("fSigmaTMassRec", &fSigmaTMassRec);

  for (Long64_t iEntry = 0; iEntry < tTreeEventMC->GetEntries(); iEntry++)
  {
    tTreeEventMC->GetEntry(iEntry);
    Double_t SigmaPz = fGammaPz + fLambdaPz;
    if (fabs(fSigmaRapidity) > 0.5)
      continue;
    // CutsApplication
    if (fLambdaMod == -100)
      continue;

    Double_t varvalue = fLambdaSigmasProton; // fGammaCosPointingAngle;

    Double_t varvalue1 = fLambdaRadius;  // fGammaRadius ; // fLambdaRadius;
    Double_t varvalue2 = fLambdaMass;    // fGammaMass; // fLambdaMass;
    Double_t varvalue3 = fLambdaRatRows; // fGammaRatRows; // fLambdaRatRows;
    Double_t varvalue4 = fLambdaK0sDiff; // fGammaChi2; // fLambdaK0sDiff ;

    Double_t varvalue5 = fLambdaCosPointingAngle; // fGammaCosPointingAngle; // fLambdaCosPointingAngle;
    Double_t varvalue6 = fLambdaDCAtoPVPos;       // fGammaDCAtoPVPos; // fLambdaDCAtoPVPos;
    Double_t varvalue7 = fLambdaDCAtoPVNeg;       // fGammaDCAtoPVNeg; // fLambdaDCAtoPVNeg;
    Double_t varvalue8 = fLambdaDCADaughters;     // fGammaDCADaughters; //fLambdaDCADaughters;

    Double_t varvalue9 = fLambdaSigmasProton; // fGammaDCAzToPrimVtx; // fLambdaSigmasProton;
    Double_t varvalue10 = fLambdacTau;        // fGammaPtPos; // fLambdacTau;
    Double_t varvalue11 = fLambdaArmAlpha;    // fGammaPtNeg; //fLambdaArmAlpha;
    Double_t varvalue12 = fLambdaArmPt;       // fGammaNsigePos; // fLambdaArmPt;
    Double_t varvalue13 = fLambdaArmPt;       // fGammaNsigeNeg; // fLambdaArmPt;
    Double_t varvalue14 = fGammaNsigPiPos;    // fLambdaArmPt;
    Double_t varvalue15 = fGammaNsigPiNeg;    // fLambdaArmPt;
    Double_t varvalue16 = fGammaNsigKPos;     // fLambdaArmPt;
    Double_t varvalue17 = fGammaNsigKNeg;     // fLambdaArmPt;

    if (fGammaTPt < 0)
      continue;
    if (fLambdaMod == 100)
      hCutMC->Fill(varvalue);
    if (fLambdaMod == -100)
      hCutMC0->Fill(varvalue);

    hCutMC1->Fill(varvalue1);
    hCutMC2->Fill(varvalue2);
    hCutMC3->Fill(varvalue3);
    hCutMC4->Fill(varvalue4);
    hCutMC5->Fill(varvalue5);
    hCutMC6->Fill(varvalue6);
    hCutMC7->Fill(varvalue7);
    hCutMC8->Fill(varvalue8);
    hCutMC9->Fill(varvalue9);
    hCutMC10->Fill(varvalue10);
    hCutMC11->Fill(varvalue11);

    hCutMC12->Fill(varvalue12);
    hCutMC13->Fill(varvalue13);

    for (int kk = 0; kk < ptbin; kk++)
    {
      if (fSigmaPt >= ptedges_Sigma[kk] && fSigmaPt < ptedges_Sigma[kk + 1])
      {
        if (UseTree == kTRUE)
        {
          hMassPtTrueRec[kk]->Fill(fGammaPsi);
          hMass1PtTrueRec[kk]->Fill(fGammaPsi);
          hCutMC14->Fill(varvalue14);
        }
      }
    }
  }
  printf("  <<<<<<<<<<<<<<<<<<End of MC acceptance from TTree \n \n");

  // Loop over extended pT bins, kk
  //""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  const int ptbinALL = 14;
  double ptedges_SigmaALL[ptbinALL + 1] = {0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.6, 2.0, 2.4, 2.8, 3.4, 4., 5.0, 8.0, 11.};
  double spectrum_mcALL[ptbinALL + 1] = {0.};
  double spectrum_mcALLe[ptbinALL + 1] = {0.};

  TH1F *sig0genP = (TH1F *)histESDmc->FindObject(Form("hMCPSig0PtGen%d", SelectRap));
  TH1F *sig0genA = (TH1F *)histESDmc->FindObject(Form("hMCASig0PtGen%d", SelectRap));
  TH1F *sig0gen = (TH1F *)sig0genP->Clone("sig0gen");
  sig0gen->Add(sig0genA, 1.);

  for (int kk = 0; kk < ptbinALL; kk++)
  {
    double base_mc = (raprange[1] - raprange[0]) * (2 * pt_points_e[kk]) * nPPmc * fASigSig;
    int nPtgen = sig0gen->GetNbinsX();
    double nevgener = 0;
    double ptMin = ptedges_SigmaALL[kk];
    double ptMax = ptedges_SigmaALL[kk + 1];
    for (Int_t i = 1; i <= nPtgen; i++)
    {
      Double_t Ptmingen = sig0gen->GetXaxis()->GetBinLowEdge(i);
      Double_t Ptmaxgen = sig0gen->GetXaxis()->GetBinUpEdge(i);
      if (Ptmingen >= ptMin && Ptmaxgen <= ptMax)
      {
        Double_t ngenir = sig0gen->GetBinContent(i);
        nevgener = nevgener + ngenir;
      }
    }
    spectrum_mcALL[kk] = nevgener / base_mc;
    spectrum_mcALLe[kk] = sqrt(nevgener) / base_mc;
    printf(" sp %f i %d base %f pT %f -- %f \n", spectrum_mcALL[kk], kk, base_mc, ptedges_SigmaALL[kk], ptedges_SigmaALL[kk + 1]);
  }
  TH1D *h_spectrum_mcALL = new TH1D("h_spectrum", "Sigma0 raw spectrum", ptbinALL + 1, ptedges_SigmaALL);
  h_spectrum_mcALL->SetMarkerStyle(21);
  h_spectrum_mcALL->SetMarkerSize(.5);
  h_spectrum_mcALL->SetTitle("pp #sqrt{s}= 7 TeV");
  h_spectrum_mcALL->GetYaxis()->SetTitle("1/N_{E}d^{2}N/dydp_{T} (GeV/c)^{-1}");
  h_spectrum_mcALL->SetLineColor(kRed);
  for (int i = 0; i < ptbinALL; i++)
  {
    h_spectrum_mcALL->SetBinContent(i + 2, spectrum_mcALL[i]);
    h_spectrum_mcALL->SetBinError(i + 2, spectrum_mcALLe[i]);
  }
  //""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  // Main loop over pT bins, kk
  for (int kk = 0; kk < ptbin; kk++)
  {
    printf("\n pT bin %d \n", kk);
    if (kk >= BinPtRange)
    {
      FitRange[0] = FitRangeLSigma0[0];
      FitRange[1] = FitRangeLSigma0[1];
    }

    if (kk >= 4 && kk < 6)
    {
      FitRange[0] = FitRangeL2Sigma0[0];
      FitRange[1] = FitRangeL2Sigma0[1];
    }

    if (kk == 3)
    {
      FitRange[0] = 1.170;
      FitRange[1] = 1.200;
    }

    if (kk == 7)
    {
      FitRange[0] = 1.172;
      FitRange[1] = 1.212;
    }
    // Part 3. Calculation of acceptance*efficiency ++++
    // read histogram from MC. NOTE here is sum Sigma0 + ASigma0 in data and MC
    TH2F *hMassPt2D_Amc = (TH2F *)histESDmc->FindObject(Form("hDtASig0MvsPtRec%d", SelectRap));
    TH2F *hMassPt2D_Pmc = (TH2F *)histESDmc->FindObject(Form("hDtPSig0MvsPtRec%d", SelectRap));
    TH2F *hMassPt2D2mc = (TH2F *)hMassPt2D_Amc->Clone("hMassPt2D2mc");
    TH2F *hMassPt2pmc = (TH2F *)hMassPt2D_Pmc->Clone("hMassPt2pmc");
    hMassPt2D2mc->Add(hMassPt2pmc, 1.);

    TH1F *sig0genP = (TH1F *)histESDmc->FindObject(Form("hMCPSig0PtGen%d", SelectRap));
    TH1F *sig0genPS = (TH1F *)sig0genP->Clone("sig0genPS");
    sig0genPS->Scale(flukageantP[kk]);
    TH1F *sig0genrecP = (TH1F *)histESDmc->FindObject(Form("hMCPSig0PtRec%d", SelectRap));

    TH1F *sig0genA = (TH1F *)histESDmc->FindObject(Form("hMCASig0PtGen%d", SelectRap));
    TH1F *sig0genAS = (TH1F *)sig0genA->Clone("sig0genAS");
    sig0genAS->Scale(flukageantA[kk]);
    TH1F *sig0genrecA = (TH1F *)histESDmc->FindObject(Form("hMCASig0PtRec%d", SelectRap));

    TH1F *sig0gen = (TH1F *)sig0genPS->Clone("sig0gen");
    TH1F *sig0genrec = (TH1F *)sig0genrecP->Clone("sig0genrec");
    sig0gen->Add(sig0genAS, 1.);
    sig0genrec->Add(sig0genrecA, 1.);
    Double_t ptMin = hMassPt2D->GetYaxis()->GetBinLowEdge(iBinPt[kk]);
    Double_t ptMax = hMassPt2D->GetYaxis()->GetBinUpEdge(iBinPt[kk + 1] - 1);

    nameBinPt = Form("pt%d", iBinPt[kk]);
    nameBinPtTrue = Form("TrueRec pt%d", iBinPt[kk]);
    nameBinPtmix0 = Form("Mix pt0%d", iBinPt[kk]);
    nameBinPtmix1 = Form("Mix pt1%d", iBinPt[kk]);
    nameBinPtmix = Form("Mix pt%d", iBinPt[kk]);
    titleBinPt = Form("M_{#Lambda #gamma}, %.1f<p_{T}<%.1f GeV/c", ptMin, ptMax);
    titleBinPtmix0 = Form("Mix M_{#Lambda #gamma}, %.1f<p_{T}<%.1f GeV/c", ptMin, ptMax);
    titleBinPtmix1 = Form("Mix M_{#Lambda #gamma}, %.1f<p_{T}<%.1f GeV/c", ptMin, ptMax);
    titleBinPtmix = Form("Mix M_{#Lambda #gamma}, %.1f<p_{T}<%.1f GeV/c", ptMin, ptMax);

    // input data in Pt bins //
    if (UseTree == kFALSE)
      hMassPt[kk] = hMassPt2D->ProjectionX(nameBinPt, iBinPt[kk], iBinPt[kk + 1] - 1);
    hMassPt[kk]->SetTitle(titleBinPt);
    hMassPt[kk]->SetLineColor(kBlue);

    // with reconstructed pT for effi calculation
    if (UseTree == kFALSE)
      hMassPtTrueRec[kk] = hMassPt2DTrueRec->ProjectionX(nameBinPtTrue, iBinPt[kk], iBinPt[kk + 1] - 1);
    hMassPtTrueRec[kk]->SetTitle(titleBinPt);
    hMassPtTrueRec[kk]->SetLineColor(kMagenta);
    hMassPtTrueRec[kk]->SetMarkerColor(kMagenta);
    hMassPtTrueRec[kk]->SetMarkerStyle(20);
    hMassPtTrueRec[kk]->SetMarkerSize(0.9);
    int nptb = hMassPtTrueRec[kk]->GetNbinsX();

    int iminbin = hMassPtTrueRec[kk]->GetXaxis()->FindBin(SigmaMassCount[0]);
    int imaxbin = hMassPtTrueRec[kk]->GetXaxis()->FindBin(SigmaMassCount[1]);

    Double_t nevgenrec = 0;
    Double_t nevgenrecMassCount = 0;
    for (Int_t ib = 1; ib <= nptb + 1; ib++)
    { // note nptb+1 to get full statistics, including last bin, "peak-effi"
      Double_t ngenreci = hMassPtTrueRec[kk]->GetBinContent(ib);
      nevgenrec = nevgenrec + ngenreci;
      if (ib >= iminbin && ib <= imaxbin)
      {
        nevgenrecMassCount = nevgenrec + ngenreci;
      }
    }
    int nPtgen = sig0gen->GetNbinsX();
    double nevgen = 0;
    for (Int_t i = 1; i <= nPtgen; i++)
    {
      Double_t Ptmingen = sig0gen->GetXaxis()->GetBinLowEdge(i);
      Double_t Ptmaxgen = sig0gen->GetXaxis()->GetBinUpEdge(i);
      if (Ptmingen >= ptMin && Ptmaxgen <= ptMax)
      {
        Double_t ngeni = sig0gen->GetBinContent(i);
        Double_t ngenreci = sig0genrec->GetBinContent(i);
        nevgen = nevgen + ngeni;
      }
    }
    input_num[kk] = nevgen;
    rec_num[kk] = nevgenrec;
    printf("Efficiency = %f +- %f   MC_rec_count = %f  MC_input_count = %f \n", Eff[kk], Eff_e[kk], rec_num[kk], input_num[kk]);

    cout << "Efficiency = " << Eff[kk] << "  +- " << Eff_e[kk] << "   MC_rec_count = " << rec_num[kk] << "   MC_input_count =  " << input_num[kk] << endl;

    cout << "Rel Error of Efficiency = " << Eff_e[kk] / Eff[kk] << endl;

    double base = Eff[kk] * (raprange[1] - raprange[0]) * (2 * pt_points_e[kk]) * nPP * feloss[kk] * fASigSig; // Events passing PhysSelectionCuts
    double base_mc = (raprange[1] - raprange[0]) * (2 * pt_points_e[kk]) * nPPmc * fASigSig;
    // **  Calculate effieciency : end **<<<<<<<<<<< //

    // Part 4.  Background subtraction with different methods +++
    //        >>>>>>>>>>>>>>>>>>>>  Pol-Gauss START
    if (UseMix == kFALSE && TypeCount == kFALSE)
    {
      /////////////////////////////////////
      // 1st Iteration : fit BG

      reject = kTRUE;

      TF1 *myBkg = myBkg = new TF1("myBkg", PolFunction, FitRange[0], FitRange[1], POLdegree + 1);
      myBkg->SetParameters(6, 0);
      myBkg->SetLineColor(2);
      hMassPt[kk]->Fit(myBkg, "NQ", "", FitRange[0], FitRange[1]); // 1st fit
      TF1 *FullFit_t2;
      TF1 *FullFit;

      FullFit_t2 = new TF1("FullFit_t2", VoigtplusPol, FitRange[0], FitRange[1], PARBINS + 1);
      FullFit_t2->SetParameter(0, 1.0);       // mean
      FullFit_t2->SetParLimits(0, 0.01, 10.); // 0.0001 to 0.01 for BW only; 0.01 to 10 for true Voigtian
      FullFit_t2->SetParameter(1, 1.19);      // mean
      FullFit_t2->SetParLimits(1, 1.175, 1.21);
      FullFit_t2->SetParameter(2, .002); // abb Gaussian sigma
      FullFit_t2->SetParLimits(2, 0.0001, .015);
      FullFit_t2->FixParameter(3, .0001); // BW width
      FullFit_t2->SetParName(0, "Norm");
      FullFit_t2->SetParName(1, "Mass");
      FullFit_t2->SetParName(2, "#sigma");
      FullFit_t2->SetParName(3, "#Gamma");
      FullFit_t2->SetParameter(4, 0);
      FullFit_t2->SetParameter(5, 0);
      FullFit_t2->SetParameter(6, 0);
      FullFit_t2->SetParameter(7, 0);
      FullFit_t2->SetParameter(8, 0);
      FullFit_t2->SetParameter(9, 0);

      /////////////////////////////////////////////
      // 2nd interation : fit Peak
      reject = kFALSE;
      double pars[PARBINS];
      double pars_e[PARBINS];

      for (int polbin = 0; polbin < PARBINS - offset; polbin++)
        bkg_params[polbin] = myBkg->GetParameter(polbin);
      for (int polbin = 0; polbin < PARBINS - offset; polbin++)
        FullFit_t2->FixParameter(offset + polbin, myBkg->GetParameter(polbin));

      hMassPt[kk]->Fit(FullFit_t2, "IMQ", "", FitRange[0], FitRange[1]); // 2nd fit //IMQN+
      for (int polbin = 0; polbin < PARBINS; polbin++)
        pars[polbin] = FullFit_t2->GetParameter(polbin);

      ///////////////////////////////////////
      // 3rd iteration : fit Total
      FullFit = new TF1("FullFit", VoigtplusPol, FitRange[0], FitRange[1], PARBINS);
      FullFit->SetParLimits(0, 0.01, 10.); // 0.0001 to 0.01 for BW only; 0.01 to 10 for true Voigtian
      FullFit->SetParLimits(1, 1.175, 1.21);
      if (kk == 0)
        FullFit->SetParLimits(1, 1.190, 1.198);
      FullFit->SetParLimits(2, 0.001, .005);
      FullFit->FixParameter(3, 0.0);
      FullFit->SetParName(0, "Norm");
      FullFit->SetParName(1, "Mass");
      FullFit->SetParName(2, "#sigma");
      FullFit->SetParName(3, "#Gamma");
      FullFit->SetParameter(4, 0);
      FullFit->SetParameter(5, 0);
      FullFit->SetParameter(6, 0);
      FullFit->SetParameter(7, 0);
      FullFit->SetParameter(8, 0);
      FullFit->SetParameter(9, 0);

      for (int polbin = 0; polbin < PARBINS; polbin++)
        FullFit->SetParameter(polbin, pars[polbin]);

      hMassPt[kk]->Fit(FullFit, "IMEQ", "", FitRange[0], FitRange[1]); // 3rd fit //IMEQ+
      for (int polbin = 0; polbin < PARBINS - offset; polbin++)
        bkg_params[polbin] = FullFit->GetParameter(polbin + offset);
      for (int polbin = 0; polbin < PARBINS; polbin++)
        pars[polbin] = FullFit->GetParameter(polbin);
      for (int polbin = 0; polbin < PARBINS; polbin++)
        pars_e[polbin] = FullFit->GetParError(polbin);

      for (int polbin = 0; polbin < POLdegree + 1; polbin++)
        myBkg->SetParameter(polbin, bkg_params[polbin]);
      for (int polbin = 0; polbin < POLdegree + 1; polbin++)
        myPol[kk]->FixParameter(polbin, bkg_params[polbin]);

      // **  Get mass and width  ** //
      Sigma0Mass[kk] = FullFit->GetParameter(1);
      Sigma0Mass_e[kk] = FullFit->GetParError(1);
      Sigma0Width[kk] = FullFit->GetParameter(2);
      Sigma0Width_e[kk] = FullFit->GetParError(2);

      // **  Get raw yield and spectrum  ** //
      MissedYield[kk] = 1000. * (FullFit->Integral(0., SigmaCountBound[0]) + FullFit->Integral(SigmaCountBound[1], 100.) - myBkg->Integral(0., SigmaCountBound[0]) - myBkg->Integral(SigmaCountBound[1], 100.));

      MissedYield_e[kk] = MissedYield[kk] * FullFit->GetParError(0) / FullFit->GetParameter(0);

      temp_yield = hMassPt[kk]->Integral(hMassPt[kk]->GetXaxis()->FindBin(SigmaCountBound[0]),
                                         hMassPt[kk]->GetXaxis()->FindBin(SigmaCountBound[1] - .0001));
      temp_yield1 = hMassPt[kk]->Integral(hMassPt[kk]->GetXaxis()->FindBin(SigmaCountBound[0]),
                                          hMassPt[kk]->GetXaxis()->FindBin(SigmaCountBound[1] - .0001));

      temp_yield -= 1000. * (myBkg->Integral(SigmaCountBound[0], SigmaCountBound[1]));
      // this line is for calculating significance

      Double_t lowPt = Sigma0Mass[kk] - 3.0 * Sigma0Width[kk];
      Double_t highPt = Sigma0Mass[kk] + 3.0 * Sigma0Width[kk];

      temp_yield += MissedYield[kk];
      for (int massBin = hMassPt[kk]->GetXaxis()->FindBin(SigmaCountBound[0]);
           massBin <= hMassPt[kk]->GetXaxis()->FindBin(SigmaCountBound[1]); massBin++)
      {
        temp_yield_e += pow(hMassPt[kk]->GetBinError(massBin), 2);
      }

      temp_yield_e += pow(MissedYield_e[kk], 2);
      temp_yield_e = sqrt(temp_yield_e);

      temp_yield1 = hMassPt[kk]->Integral(hMassPt[kk]->GetXaxis()->FindBin(lowPt),
                                          hMassPt[kk]->GetXaxis()->FindBin(highPt - .0001));

      signifi[kk] = temp_yield / sqrt(temp_yield1);
      signifi_e[kk] = temp_yield / sqrt(temp_yield1) * sqrt((temp_yield_e / temp_yield) * (temp_yield_e / temp_yield) + (sqrt(temp_yield1) / temp_yield1) * (sqrt(temp_yield1) / temp_yield1));

      printf("Low pT %f bin %f High pT %f bin %f  \n", lowPt, SigmaCountBound[0], highPt, SigmaCountBound[1] - .0001);
      cout << "Signal/Bkg = " << temp_yield / (1000. * myBkg->Integral(SigmaCountBound[0], SigmaCountBound[1])) << endl;
      cout << "Significance 1 = " << temp_yield / sqrt(temp_yield1) << endl;
      printf("TempYeild %f TempYield1 %f MissedY %f  \n", temp_yield, temp_yield1, MissedYield[kk]);
      cout << "Included/Total = " << (temp_yield - MissedYield[kk]) / temp_yield << endl;
      cout << "yield = " << temp_yield << " +- " << temp_yield_e << "Rel. Error Yield = " << temp_yield_e / temp_yield << endl;
      cout << "MC input = " << input_num[kk] << endl;
      cout << "MC true rec = " << rec_num[kk] << endl;
      Double_t data_mc_ratio = temp_yield / rec_num[kk];
      cout << "Data/rec_eff = " << temp_yield / rec_num[kk] << endl;

      yield[kk] = temp_yield;
      yield_e[kk] = temp_yield_e;

      spectrum[kk] = temp_yield;
      spectrum[kk] /= base;

      spectrum_e[kk] = pow(temp_yield_e / base, 2);
      spectrum_e[kk] += pow(temp_yield / base * Eff_e[kk] / Eff[kk], 2);
      spectrum_e[kk] = sqrt(spectrum_e[kk]);
      printf(" end of Pol-Gauss fit with  UseMix FALSE and TypeCount False \n \n");

    } // end of UseMix FALSE and TypeCount False ==  Pol-Gauss fit

    // START Pol-Count PolBG and BinCounting !!!
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if (UseMix == kFALSE && TypeCount == kTRUE)
    { // START PolBG and BinCounting !!!

      printf("POL-COUNT  \n Pt bin %d PolCount  1st Iteration : fit BG>>>>>>>>>>>>>>>> \n", kk);

      //      if( kk == 0 ){	FitRange[0] = FitRangeSigma0[0];	FitRange[1] = FitRangeSigma0[1];      }

      if (kk == 0)
      {
        FitRange[0] = 1.170;
        FitRange[1] = 1.220;
      }
      if (kk == 4)
      {
        FitRange[0] = 1.160;
        FitRange[1] = 1.230;
      }

      printf("PolBGBinCoun rangeMin %f rangeMax %f \n", FitRange[0], FitRange[1]);

      reject = kTRUE;

      TF1 *myBkg = new TF1("myBkg", PolFunction, FitRange[0], FitRange[1], POLdegree + 1);
      myBkg->SetParameters(6, 0);
      myBkg->SetLineColor(kBlue);
      hMassPt[kk]->Fit(myBkg, "NQ", "", FitRange[0], FitRange[1]); // 1st fit
      TF1 *FullFit_t2;
      TF1 *FullFit;
      TF1 *FullFit_t6;

      FullFit_t2 = new TF1("FullFit_t2", VoigtplusPol, FitRange[0], FitRange[1], PARBINS + 1);
      FullFit_t2->SetParameter(0, 1.0);       // mean
      FullFit_t2->SetParLimits(0, 0.01, 10.); // 0.0001 to 0.01 for BW only; 0.01 to 10 for true Voigtian
      FullFit_t2->SetParameter(1, 1.19);      // mean
      FullFit_t2->SetParLimits(1, 1.175, 1.21);
      FullFit_t2->SetParameter(2, 0.002); // abb Gaussian sigma
      FullFit_t2->SetParLimits(2, 0.0001, .015);
      FullFit_t2->FixParameter(3, .0001); // BW width
      FullFit_t2->SetParName(0, "Norm");
      FullFit_t2->SetParName(1, "Mass");
      FullFit_t2->SetParName(2, "#sigma");
      FullFit_t2->SetParName(3, "#Gamma");
      FullFit_t2->SetParameter(4, 0);
      FullFit_t2->SetParameter(5, 0);
      FullFit_t2->SetParameter(6, 0);
      FullFit_t2->SetParameter(7, 0);
      FullFit_t2->SetParameter(8, 0);
      FullFit_t2->SetParameter(9, 0);
      printf(" PolCoun 2nd interation :  fit Peak \n");

      reject = kFALSE;
      double pars[PARBINS];
      double pars_e[PARBINS];

      for (int polbin = 0; polbin < PARBINS - offset; polbin++)
        bkg_params[polbin] = myBkg->GetParameter(polbin);
      for (int polbin = 0; polbin < PARBINS - offset; polbin++)
        FullFit_t2->FixParameter(offset + polbin, myBkg->GetParameter(polbin));

      hMassPt[kk]->Fit(FullFit_t2, "IMQ", "", FitRange[0], FitRange[1]); // 2nd fit //IMQN+
      for (int polbin = 0; polbin < PARBINS; polbin++)
        pars[polbin] = FullFit_t2->GetParameter(polbin);

      printf(" PolCoun  3rd iteration : fit Total \n");
      FullFit_t6 = new TF1("FullFit_t6", VoigtplusPol, FitRange[0], FitRange[1], PARBINS);
      FullFit_t6->SetParLimits(0, 0.01, 10.); // 0.0001 to 0.01 for BW only; 0.01 to 10 for true Voigtian
      FullFit_t6->SetParLimits(1, 1.175, 1.21);
      FullFit_t6->SetParLimits(2, 0.001, .005);
      if (kk == 0)
      {
        FullFit_t6->SetParameter(2, 0.003); // abb Gaussian sigma
        FullFit_t6->SetParLimits(2, 0.0025, .0035);
      }
      FullFit_t6->FixParameter(3, 0.0);
      FullFit_t6->SetParName(0, "Norm");
      FullFit_t6->SetParName(1, "Mass");
      FullFit_t6->SetParName(2, "#sigma");
      FullFit_t6->SetParName(3, "#Gamma");
      FullFit_t6->SetParameter(4, 0);
      FullFit_t6->SetParameter(5, 0);
      FullFit_t6->SetParameter(6, 0);
      FullFit_t6->SetParameter(7, 0);
      FullFit_t6->SetParameter(8, 0);
      FullFit_t6->SetParameter(9, 0);

      for (int polbin = 0; polbin < PARBINS; polbin++)
        FullFit_t6->SetParameter(polbin, pars[polbin]);

      hMassPt[kk]->Fit(FullFit_t6, "IMEQ", "", FitRange[0], FitRange[1]); // 3rd fit //IMEQ+  // line 897 JSz

      for (int polbin = 0; polbin < PARBINS - offset; polbin++)
        bkg_params[polbin] = FullFit_t6->GetParameter(polbin + offset);
      for (int polbin = 0; polbin < PARBINS; polbin++)
        pars[polbin] = FullFit_t6->GetParameter(polbin);
      for (int polbin = 0; polbin < PARBINS; polbin++)
        pars_e[polbin] = FullFit_t6->GetParError(polbin);

      for (int polbin = 0; polbin < POLdegree + 1; polbin++)
        myBkg->SetParameter(polbin, bkg_params[polbin]);
      for (int polbin = 0; polbin < POLdegree + 1; polbin++)
        myPol[kk]->FixParameter(polbin, bkg_params[polbin]);

      // **  Get mass and width  ** //
      Sigma0Mass[kk] = FullFit_t6->GetParameter(1);
      Sigma0Mass_e[kk] = FullFit_t6->GetParError(1);
      Sigma0Width[kk] = FullFit_t6->GetParameter(2);
      Sigma0Width_e[kk] = FullFit_t6->GetParError(2);

      // **  Get raw yield and spectrum  ** //
      MissedYield[kk] = 1000. * (FullFit_t6->Integral(0., SigmaCountBound[0]) + FullFit_t6->Integral(SigmaCountBound[1], 100.) - myBkg->Integral(0., SigmaCountBound[0]) - myBkg->Integral(SigmaCountBound[1], 100.));

      MissedYield_e[kk] = MissedYield[kk] * FullFit_t6->GetParError(0) / FullFit_t6->GetParameter(0);
      temp_yield = hMassPt[kk]->Integral(hMassPt[kk]->GetXaxis()->FindBin(SigmaCountBound[0]),
                                         hMassPt[kk]->GetXaxis()->FindBin(SigmaCountBound[1] - .0001));
      temp_yield1 = hMassPt[kk]->Integral(hMassPt[kk]->GetXaxis()->FindBin(SigmaCountBound[0]),
                                          hMassPt[kk]->GetXaxis()->FindBin(SigmaCountBound[1] - .0001));
      // this line is for calculating significance

      temp_yield -= 1000. * (myBkg->Integral(SigmaCountBound[0], SigmaCountBound[1]));
      temp_yield += MissedYield[kk];
      for (int massBin = hMassPt[kk]->GetXaxis()->FindBin(SigmaCountBound[0]);
           massBin <= hMassPt[kk]->GetXaxis()->FindBin(SigmaCountBound[1]); massBin++)
      {
        temp_yield_e += pow(hMassPt[kk]->GetBinError(massBin), 2);
      }
      // ABB scheme: subtract Pol BG from the data inside the region of mass resolution:
      Int_t iMin = 0;
      Int_t iMax = 0;
      Double_t nSigEv = 0;
      if (TypeCount)
      {

        printf("POL-COUNT now \n");

        iMin = hMassPt[kk]->FindBin(SigmaCountBound[0]); // was 5 or 2 sigma == WHY???
        iMax = hMassPt[kk]->FindBin(SigmaCountBound[1]);

        if (kk == 0)
        {
          iMin = hMassPt[kk]->FindBin(1.189);
          iMax = hMassPt[kk]->FindBin(1.200);
        }

        printf("Pol3-BinCounts iMin %d iMax %d  \n", iMin, iMax);
        // have to subtract pol BG from here ! -> Use JS codes from l.975

        Double_t nTotEv = 0;
        Double_t nBgEv = 0;
        Int_t iMassBin = 0;
        for (Int_t iBin = iMin; iBin < iMax; iBin++)
        {

          nTotEv += hMassPt[kk]->GetBinContent(iBin);
          printf(" events %f ibin %d \n", hMassPt[kk]->GetBinContent(iBin), iBin);
        }
        nBgEv = 1000. * myBkg->Integral(SigmaCountBound[0], SigmaCountBound[1]);

        if (kk == 0)
          nBgEv = 1000. * myBkg->Integral(1.189, 1.200);
        nSigEv = nTotEv - nBgEv;

        printf("---SUM sig %f ev tot %f Bg %f kk %d  \n", nSigEv, nTotEv, nBgEv, kk);
      }
      printf("---temp_yield sig %f  tot yield1 %f BG %f  ratio-sig %f \n",
             temp_yield, temp_yield1, 1000. * myBkg->Integral(SigmaCountBound[0], SigmaCountBound[1]), nSigEv / temp_yield);

      temp_yield_e += pow(MissedYield_e[kk], 2);
      temp_yield_e = sqrt(temp_yield_e);

      signifi[kk] = temp_yield / sqrt(temp_yield1);
      signifi_e[kk] = temp_yield / sqrt(temp_yield1) * sqrt((temp_yield_e / temp_yield) * (temp_yield_e / temp_yield) + (sqrt(temp_yield1) / temp_yield1) * (sqrt(temp_yield1) / temp_yield1));

      cout << "PolCount Signal/Bkg = " << temp_yield / (1000. * myBkg->Integral(SigmaCountBound[0], SigmaCountBound[1])) << endl;
      cout << "PolCount Significance = " << temp_yield / sqrt(temp_yield1) << endl;
      cout << "PolCount Included/Total = " << (temp_yield - MissedYield[kk]) / temp_yield << endl;
      cout << "PolCount yield = " << temp_yield << " +- " << temp_yield_e << " Rel Error Yield= " << temp_yield_e / temp_yield << endl;
      cout << "PolCopunt MC input = " << input_num[kk] << endl;

      yield[kk] = temp_yield;
      yield_e[kk] = temp_yield_e;

      spectrum[kk] = temp_yield;
      spectrum[kk] /= base;

      spectrum_e[kk] = pow(temp_yield_e / base, 2);
      spectrum_e[kk] += pow(temp_yield / base * Eff_e[kk] / Eff[kk], 2);
      spectrum_e[kk] = sqrt(spectrum_e[kk]);
      printf("  END PolBG and BinCounting \n");
    } // END Pol-Count:  PolBG and BinCounting

    // START MixBG and Gaus Peak Mix-Gauss
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if (UseMix == kTRUE && TypeCount == kFALSE)
    { // START MixBG and Gaus Peak
      printf(" \n  *** Mixed event BG ***  BinpT %d \n", kk);

      hMassPt[kk]->SetLineColor(kBlue);

      hMassPt[kk] = hMassPt2D->ProjectionX(nameBinPtmix0, iBinPt[kk], iBinPt[kk + 1] - 1);
      hMassPt[kk]->SetTitle(titleBinPtmix0);
      hMassPt[kk]->SetLineColor(kBlue + 2);
      hMassPt[kk]->GetXaxis()->SetTitle("M (GeV/c^{2})");

      hMixBGSub[kk] = hMassPt2D->ProjectionX(nameBinPtmix1, iBinPt[kk], iBinPt[kk + 1] - 1);
      hMixBGSub[kk]->SetTitle(titleBinPtmix1);
      hMixBGSub[kk]->SetLineColor(kMagenta + 5);

      hMixBGSub2[kk] = hMassPt2D->ProjectionX(nameBinPtmix1, iBinPt[kk], iBinPt[kk + 1] - 1);
      hMixBGSub2[kk]->SetTitle(titleBinPtmix1);
      hMixBGSub2[kk]->SetLineColor(kMagenta + 5);

      hMixPt[kk] = hMixPt2D->ProjectionX(nameBinPtmix, iBinPt[kk], iBinPt[kk + 1] - 1);
      hMixPt[kk]->SetTitle(titleBinPtmix);
      hMixPt[kk]->SetLineColor(kRed - 3);
      hMixPt[kk]->GetXaxis()->SetTitle("M (GeV/c^{2})");

      double sumRe = 0;
      double sumMi = 0;
      double sumReL = 0;
      double sumMiL = 0;
      double sumReR = 0;
      double sumMiR = 0;

      int nptb = hMassPt[kk]->GetNbinsX();
      double ReEv[nptb], MiEv[nptb], MiEvScale[nptb], ReEvSig[nptb];
      double ReEv_e[nptb], MiEv_e[nptb], MiEvScale_e[nptb], ReEvSig_e[nptb];

      printf(" mix1min %f mix1max %f mix2min %f mix2max %f \n", mix1min, mix1max, mix2min, mix2max);
      for (Int_t ib = 1; ib <= nptb; ib++)
      {

        if (hMassPt[kk]->GetXaxis()->GetBinUpEdge(ib) >= mix1min && hMassPt[kk]->GetXaxis()->GetBinUpEdge(ib) <= mix1max)
          sumReL = sumReL + hMassPt[kk]->GetBinContent(ib);
        if (hMassPt[kk]->GetXaxis()->GetBinUpEdge(ib) >= mix2min && hMassPt[kk]->GetXaxis()->GetBinUpEdge(ib) <= mix2max)
          sumReR = sumReR + hMassPt[kk]->GetBinContent(ib);

        if (hMixPt[kk]->GetXaxis()->GetBinUpEdge(ib) >= mix1min && hMixPt[kk]->GetXaxis()->GetBinUpEdge(ib) <= mix1max)
          sumMiL = sumMiL + hMixPt[kk]->GetBinContent(ib);
        if (hMixPt[kk]->GetXaxis()->GetBinUpEdge(ib) >= mix2min && hMixPt[kk]->GetXaxis()->GetBinUpEdge(ib) <= mix2max)
          sumMiR = sumMiR + hMixPt[kk]->GetBinContent(ib);

        sumRe = sumReL + sumReR;
        sumMi = sumMiL + sumMiR;

        ReEv[ib] = hMassPt[kk]->GetBinContent(ib);
        MiEv[ib] = hMixPt[kk]->GetBinContent(ib);

        ReEv_e[ib] = sqrt(ReEv[ib]);
        MiEv_e[ib] = sqrt(MiEv[ib]);
      }
      printf(" sumReL %f sumMiL %f r= %f  \n", sumReL, sumMiL, sumMiL / sumReL);
      printf(" sumReR %f sumMiR %f r= %f \n", sumReR, sumMiR, sumMiR / sumReR);
      sumRe = sumReL + sumReR;
      sumMi = sumMiL + sumMiR;
      printf(" sumRe %f sumMi %f r= %f \n", sumRe, sumMi, sumMi / sumRe);

      if (sumMi != 0)
        hMixPt[kk]->Scale(sumRe / sumMi);
      else
        hMixPt[kk]->Scale(1.);

      Double_t fscale = sumRe / sumMi;
      for (Int_t ib = 1; ib <= nptb; ib++)
      {
        MiEv[ib] = MiEv[ib] * fscale;
        MiEv_e[ib] = MiEv[ib] * fscale;

        ReEvSig[ib] = ReEv[ib] - MiEv[ib];
        ReEvSig_e[ib] = sqrt(ReEv_e[ib] * ReEv_e[ib] + MiEv_e[ib] * MiEv_e[ib]);

        if (ReEv[ib] == 0)
          ReEvSig[ib] = 0.;
      }

      for (int i = 1; i <= nptb; i++)
      {

        hReEvSig[kk]->SetBinContent(i, ReEvSig[i]);
        hReEvSig[kk]->SetBinError(i, ReEvSig_e[i]);

        hReEv[kk]->SetBinContent(i, ReEv[i]);
        hReEv[kk]->SetBinError(i, ReEv_e[i]);

        hMiEv[kk]->SetBinContent(i, MiEv[i]);
        hMiEv[kk]->SetBinError(i, MiEv_e[i]);

        hMiEvScale[kk]->SetBinContent(i, MiEvScale[i]);
        hMiEvScale[kk]->SetBinError(i, MiEvScale_e[i]);
      }

      if (iMixPol == kFALSE)
      { // for counting of events inside mix1max and mix1min

        hMixBGSub[kk]->Add(hMixPt[kk], -1);
        hMixBGSub2[kk]->Add(hMixBGSub[kk], 1);

        temp_yield = hReEvSig[kk]->Integral(hReEvSig[kk]->GetXaxis()->FindBin(mix1max),
                                            hReEvSig[kk]->GetXaxis()->FindBin(mix2min - .0001));
        temp_yield1 = hReEv[kk]->Integral(hReEvSig[kk]->GetXaxis()->FindBin(mix1max),
                                          hReEv[kk]->GetXaxis()->FindBin(mix2min - .0001));

        for (int massBin = hReEvSig[kk]->GetXaxis()->FindBin(SigmaCountBound[0]);
             massBin <= hReEvSig[kk]->GetXaxis()->FindBin(SigmaCountBound[1]); massBin++)
        {
          temp_yield_e += pow(hReEvSig[kk]->GetBinError(massBin), 2);
        }

        temp_yield = 0;
        int iminbin = hReEvSig[kk]->GetXaxis()->FindBin(SigmaMassCount[0]);
        int imaxbin = hReEvSig[kk]->GetXaxis()->FindBin(SigmaMassCount[1]);
        for (int i = 1; i <= nptb; i++)
        {
          if (i >= iminbin && i <= imaxbin)
          {
            Double_t cont = hReEvSig[kk]->GetBinContent(i);
            temp_yield += cont;
            printf("cont %f Yield %f ---- Massbin %d \n", cont, temp_yield, i);
          }
        }
        if (temp_yield <= 0)
          temp_yield_e = 0;
        else
          temp_yield_e = sqrt(fabs(temp_yield));
        printf("true-rec %f in mass-range %f ratio %f \n", nevgenrec, nevgenrecMassCount, nevgenrecMassCount / nevgenrec);

        TF1 *GaussFit;
        GaussFit = new TF1("GaussFit", GausFunction, FitRange[0], FitRange[1], 3);
        GaussFit->SetParameter(1, 1.1926); // mean
        GaussFit->SetParLimits(1, 1.1900, 1.1970);
        GaussFit->SetParameter(2, 0.002); // Gaussian sigma
        GaussFit->SetParLimits(2, 0.001, 0.010);
        GaussFit->SetLineColor(kBlack);

        printf("mixBG \n");

        hReEvSig[kk]->Fit(GaussFit, "IMQ", "", FitRange[0], FitRange[1]);

        Sigma0Mass[kk] = GaussFit->GetParameter(1);
        Sigma0Mass_e[kk] = GaussFit->GetParError(1);
        Sigma0Width[kk] = GaussFit->GetParameter(2);
        Sigma0Width_e[kk] = GaussFit->GetParError(2);

        temp_yield_e = sqrt(temp_yield_e);

        yield[kk] = temp_yield;
        yield_e[kk] = temp_yield_e;

        signifi[kk] = temp_yield / sqrt(temp_yield1);
        signifi_e[kk] = temp_yield / sqrt(temp_yield1) * sqrt((temp_yield_e / temp_yield) * (temp_yield_e / temp_yield) + (sqrt(temp_yield1) / temp_yield1) * (sqrt(temp_yield1) / temp_yield1));

        spectrum[kk] = temp_yield;
        spectrum[kk] /= base;

        spectrum_e[kk] = pow(temp_yield_e / base, 2);
        spectrum_e[kk] += pow(temp_yield / base * Eff_e[kk] / Eff[kk], 2);
        spectrum_e[kk] = sqrt(spectrum_e[kk]);

        ratiomc[kk] = spectrum[kk] / spectrum_mc[kk];
        ratiomc_e[kk] = ratiomc[kk] * sqrt(spectrum_e[kk] * spectrum_e[kk] / spectrum[kk] /
                                               spectrum[kk] +
                                           spectrum_mce[kk] * spectrum_mce[kk] / spectrum_mc[kk] / spectrum_mc[kk]);

      } // end of iMIxPol FALSE -- bin counting
      if (iMixPol == kTRUE)
      {
        /////////////////////////////////////
        // 1st Iteration : fit BG

        reject = kTRUE;

        TF1 *myBkg = myBkg = new TF1("myBkg", PolFunction, FitRange[0], FitRange[1], POLdegree + 1);
        myBkg->SetParameters(6, 0);
        myBkg->SetLineColor(3);
        hReEvSig[kk]->Fit(myBkg, "NQ", "", FitRange[0], FitRange[1]); // 1st fit

        TF1 *FullFit_t2;
        TF1 *FullFit;

        FullFit_t2 = new TF1("FullFit_t2", VoigtplusPol, FitRange[0], FitRange[1], PARBINS + 1);
        FullFit_t2->SetParameter(0, 1.0);       // mean
        FullFit_t2->SetParLimits(0, 0.01, 10.); // 0.0001 to 0.01 for BW only; 0.01 to 10 for true Voigtian
        FullFit_t2->SetParameter(1, 1.19);      // mean
        FullFit_t2->SetParLimits(1, 1.190, 1.199);
        FullFit_t2->SetParameter(2, .002); // abb Gaussian sigma
        FullFit_t2->SetParLimits(2, 0.0001, .015);
        FullFit_t2->FixParameter(3, .0001); // BW width
        FullFit_t2->SetParName(0, "Norm");
        FullFit_t2->SetParName(1, "Mass");
        FullFit_t2->SetParName(2, "#sigma");
        FullFit_t2->SetParName(3, "#Gamma");
        FullFit_t2->SetParameter(4, 0);
        FullFit_t2->SetParameter(5, 0);
        FullFit_t2->SetParameter(6, 0);
        FullFit_t2->SetParameter(7, 0);
        FullFit_t2->SetParameter(8, 0);
        FullFit_t2->SetParameter(9, 0);

        /////////////////////////////////////////////
        // 2nd interation : fit Peak
        reject = kFALSE;
        double pars[PARBINS];
        double pars_e[PARBINS];

        for (int polbin = 0; polbin < PARBINS - offset; polbin++)
          bkg_params[polbin] = myBkg->GetParameter(polbin);
        for (int polbin = 0; polbin < PARBINS - offset; polbin++)
          FullFit_t2->FixParameter(offset + polbin, myBkg->GetParameter(polbin));

        hReEvSig[kk]->Fit(FullFit_t2, "IMQ", "", FitRange[0], FitRange[1]); // 2nd fit //IMQN+
        for (int polbin = 0; polbin < PARBINS; polbin++)
          pars[polbin] = FullFit_t2->GetParameter(polbin);

        ///////////////////////////////////////
        // 3rd iteration : fit Total
        FullFit = new TF1("FullFit", VoigtplusPol, FitRange[0], FitRange[1], PARBINS);
        FullFit->SetParLimits(0, 0.01, 10.); // 0.0001 to 0.01 for BW only; 0.01 to 10 for true Voigtian
        FullFit->SetParLimits(1, 1.175, 1.21);
        FullFit->SetParLimits(2, 0.001, .005);
        FullFit->FixParameter(3, 0.0);
        FullFit->SetParName(0, "Norm");
        FullFit->SetParName(1, "Mass");
        FullFit->SetParName(2, "#sigma");
        FullFit->SetParName(3, "#Gamma");
        FullFit->SetParameter(4, 0);
        FullFit->SetParameter(5, 0);
        FullFit->SetParameter(6, 0);
        FullFit->SetParameter(7, 0);
        FullFit->SetParameter(8, 0);
        FullFit->SetParameter(9, 0);
        FullFit->SetLineColor(kBlack);

        for (int polbin = 0; polbin < PARBINS; polbin++)
          FullFit->SetParameter(polbin, pars[polbin]);

        hReEvSig[kk]->Fit(FullFit, "IMEQ", "", FitRange[0], FitRange[1]); // 3rd fit //IMEQ+

        for (int polbin = 0; polbin < PARBINS - offset; polbin++)
          bkg_params[polbin] = FullFit->GetParameter(polbin + offset);
        for (int polbin = 0; polbin < PARBINS; polbin++)
          pars[polbin] = FullFit->GetParameter(polbin);
        for (int polbin = 0; polbin < PARBINS; polbin++)
          pars_e[polbin] = FullFit->GetParError(polbin);

        for (int polbin = 0; polbin < POLdegree + 1; polbin++)
          myBkg->SetParameter(polbin, bkg_params[polbin]);
        for (int polbin = 0; polbin < POLdegree + 1; polbin++)
          myPol[kk]->FixParameter(polbin, bkg_params[polbin]);

        // **  Get mass and width  ** //
        Sigma0Mass[kk] = FullFit->GetParameter(1);
        Sigma0Mass_e[kk] = FullFit->GetParError(1);
        Sigma0Width[kk] = FullFit->GetParameter(2);
        Sigma0Width_e[kk] = FullFit->GetParError(2);

        // **  Get raw yield and spectrum  ** //
        MissedYield[kk] = 1000. * (FullFit->Integral(0., SigmaCountBound[0]) + FullFit->Integral(SigmaCountBound[1], 100.) - myBkg->Integral(0., SigmaCountBound[0]) - myBkg->Integral(SigmaCountBound[1], 100.));

        MissedYield_e[kk] = MissedYield[kk] * FullFit->GetParError(0) / FullFit->GetParameter(0);

        temp_yield = hReEvSig[kk]->Integral(hMassPt[kk]->GetXaxis()->FindBin(SigmaCountBound[0]),
                                            hReEvSig[kk]->GetXaxis()->FindBin(SigmaCountBound[1] - .0001));
        temp_yield1 = hReEv[kk]->Integral(hMassPt[kk]->GetXaxis()->FindBin(SigmaCountBound[0]),
                                          hReEvSig[kk]->GetXaxis()->FindBin(SigmaCountBound[1] - .0001));
        // this line is for calculating significance
        cout << "0 Signal= " << temp_yield << "  BG= " << temp_yield1 << endl;

        temp_yield -= 1000. * (myBkg->Integral(SigmaCountBound[0], SigmaCountBound[1]));
        temp_yield += MissedYield[kk];

        cout << "1 Signal= " << temp_yield << "  BG= " << temp_yield1 << endl;

        for (int massBin = hReEvSig[kk]->GetXaxis()->FindBin(SigmaCountBound[0]);
             massBin <= hReEvSig[kk]->GetXaxis()->FindBin(SigmaCountBound[1]); massBin++)
        {
          temp_yield_e += pow(hReEvSig[kk]->GetBinError(massBin), 2);
        }

        temp_yield_e += pow(MissedYield_e[kk], 2);
        temp_yield_e = sqrt(temp_yield_e);

        signifi[kk] = temp_yield / sqrt(temp_yield1);
        signifi_e[kk] = temp_yield / sqrt(temp_yield1) * sqrt((temp_yield_e / temp_yield) * (temp_yield_e / temp_yield) + (sqrt(temp_yield1) / temp_yield1) * (sqrt(temp_yield1) / temp_yield1));

        cout << "Signal= " << temp_yield << "  BG= " << 1000. * myBkg->Integral(SigmaCountBound[0], SigmaCountBound[1]) << "  Signal/Bkg = " << temp_yield / (1000. * myBkg->Integral(SigmaCountBound[0], SigmaCountBound[1])) << endl;

        cout << "Significance = " << temp_yield / sqrt(temp_yield1) << endl;
        cout << "Included/Total = " << (temp_yield - MissedYield[kk]) / temp_yield << endl;
        cout << "yield = " << temp_yield << " +- " << temp_yield_e << " Rel Error Yield= " << temp_yield_e / temp_yield << endl;
        cout << "MC input = " << input_num[kk] << endl;
        cout << "MC true rec = " << rec_num[kk] << "ratio = rec/true " << temp_yield / rec_num[kk] << endl;
        Double_t mc_needs = temp_yield / Eff[kk];
        cout << "MCneeds = Data/rec_eff = " << temp_yield / Eff[kk] << endl;
        Double_t factor_mc = mc_needs / input_num[kk];
        cout << "factor_mc = " << factor_mc << endl;

        yield[kk] = temp_yield;
        yield_e[kk] = temp_yield_e;

        spectrum[kk] = temp_yield;
        spectrum[kk] /= base;

        spectrum_e[kk] = pow(temp_yield_e / base, 2);
        spectrum_e[kk] += pow(temp_yield / base * Eff_e[kk] / Eff[kk], 2);

        spectrum_e[kk] = sqrt(spectrum_e[kk]);
      } // end of iMixPol

    } // End of Use Mixed BG <<<<<<<<<<<<<

    spectrum_mc[kk] = nevgen / base_mc;
    spectrum_mce[kk] = sqrt(nevgen) / base_mc;

    ratiomc[kk] = spectrum[kk] / spectrum_mc[kk];
    ratiomc_e[kk] = ratiomc[kk] * sqrt(spectrum_e[kk] * spectrum_e[kk] / spectrum[kk] /
                                           spectrum[kk] +
                                       spectrum_mce[kk] * spectrum_mce[kk] / spectrum_mc[kk] / spectrum_mc[kk]);
  } // end of pt loop with kk variable
  //==============================// output //==============================//

  gStyle->SetOptStat(0);

  TFile *fm68 = new TFile("./mcPyt68Rebinned.root");
  TGraphErrors *mcSigma0Pyt6Rebin = (TGraphErrors *)fm68->Get("mcSigma0Pyt6Rebin");
  TGraphErrors *mcSigma0Pyt8Rebin = (TGraphErrors *)fm68->Get("mcSigma0Pyt8Rebin");

  TCanvas *canPyt6Pyt80 = new TCanvas("canPyt6Pyt80", "canPyt6Pyt80", 13, 34, 800, 800);
  mcSigma0Pyt6Rebin->Draw("");

  mcSigma0Pyt8Rebin->SetLineColor(kMagenta);
  mcSigma0Pyt8Rebin->SetMarkerColor(kMagenta);
  mcSigma0Pyt8Rebin->Draw("same");

  printf(" +222222222222222++++++++++++++++++++++++++++++++ add ratio Sigma/Lambda(p_T) +++++++++++++++++++++ \n");

  TGraphErrors *gr_Eff = new TGraphErrors(ptbin, pt_points, Eff, pt_points_e, Eff_e);
  gr_Eff->SetMarkerStyle(20);
  gr_Eff->SetMinimum(1.e-6);
  gr_Eff->SetMaximum(.05);
  gr_Eff->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  gr_Eff->GetYaxis()->SetTitle("Efficiency");
  gr_Eff->SetTitle("#Sigma^{0}+cc Efficiency");

  TH1D *h_Eff = new TH1D("h_Eff", "AcceptanceEfficiency", ptbin + 1, ptedges);
  h_Eff->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_Eff->GetYaxis()->SetTitle("Acceptance*Efficiency*B.R.");

  TH1D *h_GenTrue = new TH1D("h_GenTrue", "GenTrue", ptbin + 1, ptedges);
  TH1D *h_RecTrue = new TH1D("h_RecTrue", "RecTrue", ptbin + 1, ptedges);
  h_RecTrue->SetMarkerStyle(21);
  h_RecTrue->SetMarkerSize(.7);
  h_RecTrue->SetTitle("Rec/True");
  h_RecTrue->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_RecTrue->GetYaxis()->SetTitle("counts");
  h_RecTrue->SetLineColor(kBlue - 3);
  h_RecTrue->SetMarkerColor(kBlue - 3);

  TH1D *h_Ratio = new TH1D("h_Ratio", "Ratio", ptbin + 1, ptedges);
  h_Ratio->SetMarkerStyle(21);
  h_Ratio->SetMarkerSize(.5);
  h_Ratio->SetTitle("Rec/True");
  h_Ratio->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_Ratio->GetYaxis()->SetTitle("ratio");
  h_Ratio->SetLineColor(kBlue - 3);
  h_Ratio->SetMarkerColor(kBlue - 3);

  TH1D *h_rls = new TH1D("h_Ratio", "Ratio", ptbin + 1, ptedges);

  h_rls->SetMarkerStyle(21);
  h_rls->SetMarkerSize(.5);
  h_rls->SetTitle("Rec/True");
  h_rls->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_rls->GetYaxis()->SetTitle("ratio");
  h_rls->SetLineColor(kBlue - 3);
  h_rls->SetMarkerColor(kBlue - 3);

  TH1D *hInput = new TH1D("hInput", "hInput", ptbin + 1, ptedges);
  for (int i = 0; i < ptbin; i++)
  {
    h_Eff->SetBinContent(i + 2, Eff[i]);
    h_Eff->SetBinError(i + 2, Eff_e[i]);

    h_GenTrue->SetBinContent(i + 2, input_num[i]);
    h_GenTrue->SetBinError(i + 2, sqrt(input_num[i]));

    h_RecTrue->SetBinContent(i + 2, rec_num[i]);
    h_RecTrue->SetBinError(i + 2, sqrt(rec_num[i]));

    hInput->SetBinContent(i + 2, input_num[i]);
    hInput->SetBinError(i + 2, sqrt(input_num[i]));
  }

  TH1D *h_Rawspectrum = new TH1D("h_Rawspectrum", "Sigma0 raw spectrum", ptbin + 1, ptedges);
  h_Rawspectrum->SetMarkerStyle(20);
  h_Rawspectrum->SetMarkerSize(1);
  h_Rawspectrum->SetTitle("pp #sqrt{s}= 7 TeV");
  h_Rawspectrum->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_Rawspectrum->GetYaxis()->SetTitle("counts");
  h_Rawspectrum->SetMinimum(3.0e-6);
  h_Rawspectrum->GetXaxis()->SetLimits(0, 5.6);

  TH1D *h_NormRawspectrum = new TH1D("h_NormRawspectrum", "Sigma0 norm raw spectrum", ptbin + 1, ptedges);
  h_NormRawspectrum->SetMarkerStyle(20);
  h_NormRawspectrum->SetMarkerSize(1);
  h_NormRawspectrum->SetTitle("pp #sqrt{s}= 7 TeV");
  h_NormRawspectrum->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_NormRawspectrum->GetYaxis()->SetTitle("counts");

  Double_t xPyt6, xPyt8, yPyt6, yPyt8, rat68;
  for (int i = 0; i < ptbin; i++)
  {

    h_Rawspectrum->SetBinContent(i + 2, yield[i]);
    h_Rawspectrum->SetBinError(i + 2, yield_e[i]);

    h_NormRawspectrum->SetBinContent(i + 2, yield[i] / nPP);
    h_NormRawspectrum->SetBinError(i + 2, yield_e[i] / nPP);

    printf("yeild/npp %f  bin %d", yield[i] / nPP, i);
    cout << "    yeild/npp = " << yield[i] / nPP << endl;

    ratioev[i] = yield[i] / rec_num[i];
    ratioev_e[i] = ratioev[i] * sqrt((yield_e[i] / yield[i]) * (yield_e[i] / yield[i]) + 1 / rec_num[i]);

    h_Ratio->SetBinContent(i + 2, ratioev[i]);
    h_Ratio->SetBinError(i + 2, ratioev_e[i]);
  }

  TCanvas *Rec = new TCanvas("Rec", "Reconstructed", 0, 0, 800, 800);
  Rec->Divide(1, 2);
  Rec->cd(1);
  h_RecTrue->Draw("E");
  TLegend *legl = new TLegend(0.6, 0.5, 0.9, 0.9); // TLegend( X0, Y0, X1, Y1 )
  legl->AddEntry(h_RecTrue, "True", "lp");
  legl->AddEntry(h_Rawspectrum, "After BG subtr.", "lp");
  legl->Draw("same");
  h_Rawspectrum->Draw("sameE");
  Rec->cd(2);
  gStyle->SetOptFit(111);
  h_Ratio->Fit("pol0", "", "", 1.1, 8.);
  h_Ratio->Draw("E");

  printf(" -RRR-AAA---TTT---III---OOO \n ");

  TH1D *h_spectrum = new TH1D("h_spectrum", "Sigma0 raw spectrum", ptbin + 1, ptedges);
  h_spectrum->SetMarkerStyle(21);
  h_spectrum->SetMarkerSize(.5);
  h_spectrum->SetTitle("pp #sqrt{s}= 7 TeV");
  h_spectrum->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_spectrum->GetYaxis()->SetTitle("1/N_{E}d^{2}N/dydp_{T} (GeV/c)^{-1}");
  h_spectrum->SetLineColor(kRed);
  h_spectrum->SetMarkerColor(kRed);

  for (int i = 0; i < ptbin; i++)
  {
    h_spectrum->SetBinContent(i + 2, spectrum[i]);
    h_spectrum->SetBinError(i + 2, spectrum_e[i]);
  }

  TCanvas *CStat = new TCanvas("CSTAT", "StatErrors", 0, 0, 800, 800);
  CStat->Divide(2, 2);
  CStat->cd(1);
  h_Rawspectrum->Draw("E");
  CStat->cd(2);
  TH1D *h_staterr = new TH1D("h_staterr", "Sigma0 Stat Errors", ptbin + 1, ptedges);
  h_staterr->SetMarkerStyle(21);
  h_staterr->SetMarkerSize(.5);
  h_staterr->SetTitle("pp #sqrt{s}= 7 TeV");
  h_staterr->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_staterr->GetYaxis()->SetTitle("1/N_{E}d^{2}N/dydp_{T} (GeV/c)^{-1}");
  h_staterr->SetLineColor(kRed);

  for (int i = 0; i < ptbin; i++)
  {
    h_staterr->SetBinContent(i + 2, yield_e[i]);
    h_staterr->SetBinError(i + 2, sqrt(yield_e[i]));
  }
  h_staterr->Draw("E");

  CStat->cd(3);
  TH1D *h_relstaterr = new TH1D("h_relstaterr", "Sigma0 Stat Errors", ptbin + 1, ptedges);
  h_relstaterr->SetMarkerStyle(21);
  h_relstaterr->SetMarkerSize(.5);
  h_relstaterr->SetTitle("pp #sqrt{s}= 7 TeV");
  h_relstaterr->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_relstaterr->GetYaxis()->SetTitle("1/N_{E}d^{2}N/dydp_{T} (GeV/c)^{-1}");
  h_relstaterr->SetLineColor(kBlue);

  for (int i = 0; i < ptbin; i++)
  {
    h_relstaterr->SetBinContent(i + 2, yield_e[i] / yield[i]);
    h_relstaterr->SetBinError(i + 2, sqrt(yield_e[i]) / yield[i]);
  }
  gStyle->SetOptFit(1111);
  h_relstaterr->Fit("pol0", "", "", 1., 10.);
  h_relstaterr->Draw("E");

  TH1D *h_spectrum_mc = new TH1D("h_spectrum_mc", "Sigma0 mc raw spectrum", ptbin + 1, ptedges);
  h_spectrum_mc->SetMarkerStyle(20);
  h_spectrum_mc->SetMarkerSize(.5);
  h_spectrum_mc->SetTitle("pp #sqrt{s}= 7 TeV");
  h_spectrum_mc->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_spectrum_mc->GetYaxis()->SetTitle("1/N_{E}d^{2}N/dydp_{T} (GeV/c)^{-1}");
  h_spectrum_mc->SetLineColor(kRed);

  for (int i = 0; i < ptbin; i++)
  {
    h_spectrum_mc->SetBinContent(i + 2, spectrum_mc[i]);
    h_spectrum_mc->SetBinError(i + 2, spectrum_mce[i]);
  }

  TGraphErrors *SigMass = new TGraphErrors(ptbin, pt_points, Sigma0Mass, pt_points_e, Sigma0Mass_e);
  SigMass->SetMarkerStyle(22);
  SigMass->SetMinimum(1.19);
  SigMass->SetMaximum(1.198);
  SigMass->GetXaxis()->SetTitle("p_{T} ");
  SigMass->GetYaxis()->SetTitle("(GeV/c^{2})");
  SigMass->SetTitle("#Sigma^{0} mass vs p_{T}");

  TGraphErrors *SigWidth = new TGraphErrors(ptbin, pt_points, Sigma0Width, pt_points_e, Sigma0Width_e);
  SigWidth->SetMarkerStyle(23);
  SigWidth->SetMinimum(0.);
  SigWidth->SetMaximum(.005);
  SigWidth->GetXaxis()->SetTitle("p_{T} ");
  SigWidth->GetYaxis()->SetTitle("(GeV/c)");
  SigWidth->SetTitle("#Sigma^{0} width vs p_{T}");

  TGraphErrors *SigSignifi = new TGraphErrors(ptbin, pt_points, signifi, pt_points_e, signifi_e);
  SigSignifi->SetMarkerStyle(20);
  SigSignifi->SetMinimum(0.);
  SigSignifi->SetMaximum(15.);
  SigSignifi->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  SigSignifi->GetYaxis()->SetTitle("Sig/#sqrt{Sig+BG}");
  SigSignifi->SetTitle("#Sigma^{0} significance vs p_{T}");
  gStyle->SetOptStat(0);

  TH1D *h_SigMass = new TH1D("h_SigMass", "Mass in Pt bins", ptbin + 1, ptedges);
  TH1D *h_SigWidth = new TH1D("h_SigWidth", "Width of mass in Pt bins", ptbin + 1, ptedges);
  TH1D *h_SigSignifi = new TH1D("h_SigSignifi", "Significance", ptbin + 1, ptedges);
  h_SigSignifi->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_SigSignifi->GetYaxis()->SetTitle("counts");

  TH1D *h_SigYields = new TH1D("h_SigYields", "Yields", ptbin + 1, ptedges);
  h_SigYields->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_SigYields->GetYaxis()->SetTitle("counts");

  for (int i = 0; i < ptbin; i++)
  {
    h_SigMass->SetBinContent(i + 2, Sigma0Mass[i]);
    h_SigMass->SetBinError(i + 2, Sigma0Mass_e[i]);

    h_SigWidth->SetBinContent(i + 2, Sigma0Width[i]);
    h_SigWidth->SetBinError(i + 2, Sigma0Width_e[i]);

    h_SigSignifi->SetBinContent(i + 2, signifi[i]);
    h_SigSignifi->SetBinError(i + 2, signifi_e[i]);

    h_SigYields->SetBinContent(i + 2, yield[i]);
    h_SigYields->SetBinError(i + 2, yield_e[i]);
  }
  h_SigMass->GetXaxis()->SetTitle("p_{T} (GeV/c)  ");
  h_SigMass->GetYaxis()->SetTitle("M (GeV/c^{2})  ");

  h_SigWidth->GetXaxis()->SetTitle("p_{T},(GeV/c)  ");
  h_SigWidth->GetYaxis()->SetTitle("#sigma_{M},(GeV/c^{2})  ");

  gStyle->SetOptStat(0);
  TCanvas *cM = new TCanvas("cM", "", 0, 0, 800, 800);
  gStyle->SetOptStat(0);
  h_SigMass->Draw("E");

  gStyle->SetOptStat(0);
  TCanvas *cW = new TCanvas("cW", "", 0, 0, 800, 800);
  gStyle->SetOptStat(0);
  h_SigWidth->Draw("E");

  TCanvas *InfoHisto = new TCanvas("InfoHisto", "Raw yield, Mass, Width, Significance", 0, 0, 800, 800);
  InfoHisto->Divide(2, 2);
  InfoHisto->cd(1);
  h_Rawspectrum->Draw("E");
  InfoHisto->cd(2);
  SigSignifi->Draw("");
  InfoHisto->cd(3);
  SigMass->Draw("");
  InfoHisto->cd(4);
  SigWidth->Draw("");

  TH1D *h_specsyst = new TH1D("h_specsyst", "Sigma0 spectrum, systematic errors", ptbin + 1, ptedges);
  h_specsyst->SetMarkerStyle(21);
  h_specsyst->SetMarkerSize(.5);
  h_specsyst->SetTitle("pp #sqrt{s}= 7 TeV");
  h_specsyst->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_specsyst->SetTitle("1/N_{E}d^{2}N/dydp_{T} (GeV/c)^{-1}");
  h_specsyst->SetLineColor(kRed);
  h_specsyst->SetMarkerColor(kRed);

  for (int i = 0; i < ptbin; i++)
  {
    h_spectrum->SetBinContent(i + 2, spectrum[i]);
    h_spectrum->SetBinError(i + 2, spectrum_e[i]);
  }

  TH1D *h_spectotal = new TH1D("h_spectotal", "Sigma0 spectrum, total erros", ptbin + 1, ptedges);
  h_spectotal->SetMarkerStyle(21);
  h_spectotal->SetMarkerSize(1.1);
  h_spectotal->SetTitle("pp #sqrt{s}= 7 TeV");
  h_spectotal->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_spectotal->GetYaxis()->SetTitle("1/N_{inel}d^{2}N/dydp_{T} (GeV/c)^{-1}");
  h_spectotal->SetLineColor(kRed);
  h_spectotal->SetMarkerColor(kRed);

  TH1D *ratio = new TH1D("ratio", "Rec./Sim.", ptbin + 1, ptedges);
  ratio->SetLineColor(1);

  TH1D *hratio68 = new TH1D("hratio68", "Sim-Pyt6/Sim-Pyt8", ptbin + 1, ptedges);
  hratio68->SetLineColor(1);

  TH1D *hratio68dtmc = new TH1D("hratio68dtmc", "Sig0 DTPyt8.", ptbin + 1, ptedges);
  hratio68dtmc->SetLineColor(kMagenta);

  TH1D *ratiototal = new TH1D("ratiototal", "DT/Pyt6", ptbin + 1, ptedges);
  ratio->SetLineColor(1);

  TH1D *ratioStatDtMC = new TH1D("ratioStatDtMC", "Rec-stat-err./Sim.", ptbin + 1, ptedges);
  ratio->SetLineColor(1);

  Double_t ratiomc68[8], ratiomc68_e[8];
  Double_t ratiomc_e[8], ratiomc_esyst[8];
  for (int i = 0; i < ptbin; i++)
  {

    mcSigma0Pyt6Rebin->GetPoint(i, xPyt6, yPyt6);
    mcSigma0Pyt8Rebin->GetPoint(i, xPyt8, yPyt8);
    rat68 = yPyt6 / yPyt8;

    ratiomc68[i] = ratiomc[i] * rat68;
    ratiomc68_e[i] = ratiomc_e[i] * rat68;

    ratio->SetBinContent(i + 2, ratiomc[i]);
    ratio->SetBinError(i + 2, ratiomc_e[i]);

    fsystotal[i] = sqrt(fsyserr * fsyserr + fsystfluct[i] * fsystfluct[i]);
    printf(" Fsyst Err %f \n ", fsystotal[i]);

    spectrum_esys[i] = spectrum[i] * fsystotal[i];
    spectrum_etot[i] = sqrt(spectrum_esys[i] * spectrum_esys[i] + spectrum_e[i] * spectrum_e[i]);
    printf("Sigma0  %f  Total  Err %f Syst %f Stat %f Bin %d \n ",
           spectrum[i], spectrum_etot[i], spectrum_esys[i], spectrum_e[i], i);

    h_specsyst->SetBinContent(i + 2, spectrum[i]);
    h_specsyst->SetBinError(i + 2, spectrum_esys[i]);

    h_spectotal->SetBinContent(i + 2, spectrum[i]);
    h_spectotal->SetBinError(i + 2, spectrum_etot[i]);

    ratiomc_etot[i] = ratiomc[i] * sqrt(spectrum_etot[i] * spectrum_etot[i] / spectrum[i] / spectrum[i] + spectrum_mce[i] * spectrum_mce[i] / spectrum_mc[i] / spectrum_mc[i]);

    ratiomc_esyst[i] = ratiomc[i] * sqrt(spectrum_esys[i] * spectrum_esys[i] / spectrum[i] / spectrum[i] + spectrum_mce[i] * spectrum_mce[i] / spectrum_mc[i] / spectrum_mc[i]);

    ratiomc_e[i] = ratiomc[i] * sqrt(spectrum_e[i] * spectrum_e[i] / spectrum[i] / spectrum[i] + spectrum_mce[i] * spectrum_mce[i] / spectrum_mc[i] / spectrum_mc[i]);

    ratiototal->SetBinContent(i + 2, ratiomc[i]);
    ratiototal->SetBinError(i + 2, ratiomc_esyst[i]);

    hratio68->SetBinContent(i + 2, rat68);
    hratio68->SetBinError(i + 2, rat68 * 0.01);

    hratio68dtmc->SetBinContent(i + 2, ratiomc[i] * rat68);
    hratio68dtmc->SetBinError(i + 2, ratiomc_e[i] * rat68);

    ratioStatDtMC->SetBinContent(i + 2, ratiomc[i]);
    ratioStatDtMC->SetBinError(i + 2, ratiomc_e[i]);

    printf("ratioStatDtMC %f +- %f bin %d \n \n \n ", ratiomc[i], ratiomc_e[i], i);
  }

  double systerr_e[ptbin];
  double lamsysterr_e[ptbin];

  double systerrMat[ptbin];
  double systerrINEL[ptbin];
  for (int i = 0; i < ptbin; i++)
  {
    lamsysterr_e[i] = 2.e-2;
    systerr_e[i] = 1.e-6;
    systerrMat[i] = fsyserr;
    systerrINEL[i] = fsysINEL;
  }

  TGraphErrors *h_SystLam = new TGraphErrors(ptbin, pt_points, fsystLam, pt_points_e, systerr_e);
  TGraphErrors *h_SystLam2 = new TGraphErrors(ptbin, pt_points, fsystLam, pt_points_e, lamsysterr_e);

  TGraphErrors *h_SystGam = new TGraphErrors(ptbin, pt_points, fsystGam, pt_points_e, systerr_e);
  TGraphErrors *h_SystGam2 = new TGraphErrors(ptbin, pt_points, fsystGam, pt_points_e, lamsysterr_e);

  TGraphErrors *h_SystMass = new TGraphErrors(ptbin, pt_points, fsystMass, pt_points_e, systerr_e);
  TGraphErrors *h_SystMat = new TGraphErrors(ptbin, pt_points, systerrMat, pt_points_e, systerr_e);
  TGraphErrors *h_SystINEL = new TGraphErrors(ptbin, pt_points, systerrINEL, pt_points_e, systerr_e);
  TGraphErrors *h_SystTot = new TGraphErrors(ptbin, pt_points, fsystotal, pt_points_e, systerr_e);

  h_SystTot->SetMarkerStyle(21);
  h_SystTot->SetMarkerSize(1.0);
  h_SystTot->SetLineColor(1);
  h_SystTot->SetMarkerColor(1);
  h_SystTot->GetXaxis()->SetTitle("p_{T} (GeV/c)");

  h_SystMass->SetMarkerStyle(23);
  h_SystMass->SetMarkerSize(0.8);
  h_SystMass->SetLineColor(kRed + 3);
  h_SystMass->SetMarkerColor(kRed + 3);

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

  TCanvas *cSyst = new TCanvas("cSyst", "Systematics", 0, 0, 800, 800);
  h_SystTot->Draw("");
  h_SystMass->Draw("same");
  h_SystGam->Draw("same");
  h_SystLam->Draw("same");
  h_SystMat->Draw("same");

  TLegend *legl1 = new TLegend(0.6, 0.5, 0.9, 0.9); // TLegend( X0, Y0, X1, Y1 )
  legl1->AddEntry(h_SystTot, "Total", "lp");
  legl1->AddEntry(h_SystMass, "signal extraction", "lp");
  legl1->AddEntry(h_SystGam, "gamma selection", "lp");
  legl1->AddEntry(h_SystLam, "lambda selection", "lp");
  legl1->AddEntry(h_SystMat, "material budget", "lp");

  legl1->Draw("same");

  TCanvas *cSyst2 = new TCanvas("cSyst2", "Systematics", 0, 0, 800, 800);

  h_SystLam2->Draw("");

  TCanvas *cSyst3 = new TCanvas("cSyst3", "Systematics", 0, 0, 800, 800);

  h_SystGam2->Draw("");

  TGraphErrors *gr_stat = new TGraphErrors(ptbin, pt_points, spectrum, pt_points_e, spectrum_e);
  gr_stat->SetMarkerStyle(20);
  gr_stat->SetMinimum(1.e-6);
  gr_stat->SetMaximum(.05);
  gr_stat->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  gr_stat->GetYaxis()->SetTitle(" ");
  gr_stat->SetTitle(" ");

  ratiototal->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  ratiototal->GetYaxis()->SetTitle("ratio=rec./sim.");

  TCanvas *spectrum_tot = new TCanvas("spectrum_tot", "Data/MC Sigma0 spectrum", 0, 0, 800, 1000);
  gStyle->SetOptStat(111);
  gStyle->SetOptFit(1111);
  spectrum_tot->Divide(1, 2);
  spectrum_tot->cd(1);
  h_spectrum->Draw("E");
  h_spectrum_mc->Draw("Esame");
  spectrum_tot->cd(2);
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(11111);

  ratiototal->Draw("E1");
  ratioStatDtMC->Draw("E1,same");

  TCanvas *spSt = new TCanvas("spSt", "Data/MC Sigma0 ratio", 0, 0, 800, 1000);
  ratiototal->Draw("E1");
  ratioStatDtMC->Draw("E1,same");

  printf("...start can \n");

  TCanvas *can = new TCanvas("can", "can", 13, 34, 1200, 600);
  gStyle->SetOptStat(1111);

  can->Range(-1.25, -0.2625, 11.25, 2.3625);
  can->SetFillColor(10);
  can->SetBorderMode(0);
  can->SetBorderSize(2);
  can->SetFrameFillColor(0);
  can->SetFrameBorderMode(0);
  can->SetFrameBorderMode(0);
  can->Divide(4, 2);
  for (int ll = 0; ll < ptbin; ll++)
  {
    can->cd(ll + 1);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    hMassPt[ll]->Draw("E");

    if (UseMix == kFALSE)
    {
      myPol[ll]->DrawCopy("same");
    }
    else
    {
      gStyle->SetOptStat(0);
      gStyle->SetOptFit(0);

      hMixPt[ll]->Draw("same");

      gStyle->SetOptStat(0);
      gStyle->SetOptFit(0);
      hReEvSig[ll]->Draw("sameE");
    }
  }

  printf(" +++++++++++++++++++++++++++++++++ add ratio Sigma/Lambda(p_T) +++++++++++++++++++++ \n");

  Double_t xPtd, yPtd, exPtd, eyPtd, pt[8], dpt[8], ratiodmc[8], ept[8], eratiodmc[8];
  Double_t ratioSigLam[8], eratioSigLam[8], eStatSig, eTotSig;
  Double_t xPtmc, yPtmc, exPtmc, eyPtmc;
  Double_t xPtSig, yPtSig, exPtSig, eyPtSig;

  TFile *f1 = new TFile("LambdaYieldOffi-17jan17.root");
  TGraphErrors *hLamSpec = (TGraphErrors *)f1->Get("LamYieldRebin");

  Double_t ratioLamSig[8];

  double lamspectr[ptbin];
  double lamspectr_e[ptbin], lamspectr_estat[ptbin];

  double eratioSystLamSig[ptbin], eratioStatLamSig[ptbin], eratioTotLamSig[ptbin];
  double eratioSystSigLam[ptbin], eratioStatSigLam[ptbin], eratioTotSigLam[ptbin];

  for (Int_t iPt = 0; iPt < hLamSpec->GetN(); iPt++)
  {
    hLamSpec->GetPoint(iPt, xPtd, yPtd);

    exPtd = hLamSpec->GetErrorX(iPt);
    eyPtd = hLamSpec->GetErrorY(iPt);

    if (iPt <= 6)
      eyPtd = 0.08 * yPtd;
    else
      eyPtd = 0.10 * yPtd;

    lamspectr[iPt] = yPtd;
    lamspectr_e[iPt] = eyPtd;
    lamspectr_estat[iPt] = 0.02 * yPtd;
    yPtSig = spectrum[iPt];       //  h_spectotal->GetBinContent(iPt);
    eyPtSig = spectrum_esys[iPt]; // h_spectotal->GetBinError(iPt) ;
    eStatSig = spectrum_e[iPt];
    eTotSig = sqrt(spectrum_esys[iPt] * spectrum_esys[iPt] + spectrum_e[iPt] * spectrum_e[iPt]);

    Double_t ptMinx = hLamSpec->GetXaxis()->GetBinLowEdge(iPt);
    Double_t ptMaxx = hLamSpec->GetXaxis()->GetBinUpEdge(iPt);
    Double_t LamVal = hLamSpec->GetXaxis()->GetBinCenter(iPt);
    printf("\n \n Lam ipt-ibin = %d  ptMin %f ptMax %f LamVal %f Lam %f +- %f  \n", iPt, ptMinx, ptMaxx, LamVal, yPtd, eyPtd);

    Double_t ptMin = h_spectotal->GetXaxis()->GetBinLowEdge(iPt);
    Double_t ptMax = h_spectotal->GetXaxis()->GetBinUpEdge(iPt);
    printf("Sig0 ipt-ibin = %d ptMin %f ptMax %f Sig0 %f +- sys %f +- stat %f \n", iPt, ptMin, ptMax, yPtSig, eyPtSig, eStatSig);

    pt[iPt] = xPtSig;
    ept[iPt] = exPtSig;
    dpt[iPt] = exPtSig;

    ratiodmc[iPt] = yPtSig / yPtd; // Sigma/Lambda  2. IF PSig or ASig, As Lambda == Const

    ratioSigLam[iPt] = yPtSig / yPtd;

    ratioLamSig[iPt] = 1. / ratioSigLam[iPt];

    eratiodmc[iPt] = ratiodmc[iPt] * sqrt(eyPtd * eyPtd / yPtd / yPtd + eyPtSig * eyPtSig / yPtSig / yPtSig); // correct

    eratioStatSigLam[iPt] = ratiodmc[iPt] * sqrt(eyPtd * eyPtd / yPtd / yPtd + eStatSig * eStatSig / yPtSig / yPtSig); // correct

    eratioTotSigLam[iPt] = ratiodmc[iPt] * sqrt(eyPtd * eyPtd / yPtd / yPtd + eTotSig * eTotSig / yPtSig / yPtSig); // correct

    // Error:
    eratioSystLamSig[iPt] = ratioLamSig[iPt] * sqrt(eyPtd * eyPtd / yPtd / yPtd + eyPtSig * eyPtSig / yPtSig / yPtSig);   // correct
    eratioStatLamSig[iPt] = ratioLamSig[iPt] * sqrt(eyPtd * eyPtd / yPtd / yPtd + eStatSig * eStatSig / yPtSig / yPtSig); // correct
    eratioTotLamSig[iPt] = ratioLamSig[iPt] * sqrt(eyPtd * eyPtd / yPtd / yPtd + eTotSig * eTotSig / yPtSig / yPtSig);    // correct

    printf("lam-x %f +- %f sig-x %f +- %f \n", xPtd, exPtd, xPtmc, exPtmc);
    printf("Sig-y %f +- stat %f  +- tot %f \n", yPtSig, eStatSig, eyPtSig);
    printf("++ratio %f L/S = %f er %f corr %f  \n", ratiodmc[iPt], ratioLamSig[iPt], eratiodmc[iPt], etacorr);
  }

  // here is ratio with Syst Errors only:
  TGraphErrors *hRatioSigLam = new TGraphErrors(ptbin, pt_points, ratiodmc, pt_points_e, eratiodmc);
  TGraphErrors *hRatioStatSigLam = new TGraphErrors(ptbin, pt_points, ratiodmc, pt_points_e, eratioStatSigLam);
  TGraphErrors *hRatioTotSigLam = new TGraphErrors(ptbin, pt_points, ratioSigLam, pt_points_e, eratioTotSigLam);
  TGraphErrors *hRatioSystSigLam = new TGraphErrors(ptbin, pt_points, ratioSigLam, pt_points_e, eratioSystSigLam);
  TGraphErrors *hRatioLamSig = new TGraphErrors(ptbin, pt_points, ratioLamSig, pt_points_e, eratiodmc);
  TGraphErrors *hRatioStatLamSig = new TGraphErrors(ptbin, pt_points, ratioLamSig, pt_points_e, eratioStatSigLam);
  TGraphErrors *hRatioTotLamSig = new TGraphErrors(ptbin, pt_points, ratioLamSig, pt_points_e, eratioTotSigLam);
  TGraphErrors *hRatioSystLamSig = new TGraphErrors(ptbin, pt_points, ratioLamSig, pt_points_e, eratioTotSigLam);

  hRatioSigLam->SetTitle("hRatio");
  gStyle->SetOptFit(111);

  TCanvas *cSig = new TCanvas("cSig", "cSigtot", 10, 10, 600, 600);
  h_spectotal->Draw("");

  TCanvas *cLam = new TCanvas("cLam", "cLamtot", 10, 10, 600, 600);

  TH1D *h_RatioSystSigLam = new TH1D("h_RatioSystSigLam", "ratio Sig/Lambda spectrum, syst erros", ptbin + 1, ptedges);
  h_RatioSystSigLam->SetMarkerStyle(21);
  h_RatioSystSigLam->SetMarkerSize(.5);
  h_RatioSystSigLam->SetTitle("pp #sqrt{s}= 7 TeV");
  h_RatioSystSigLam->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_RatioSystSigLam->GetYaxis()->SetTitle("ratio = #Sigma^{0}/#Lambda");
  h_RatioSystSigLam->SetLineColor(kRed);

  TH1D *h_RatioStatSigLam = new TH1D("h_RatioStatSigLam", "ratio Sig/Lambda spectrum, stat erros", ptbin + 1, ptedges);
  h_RatioStatSigLam->SetMarkerStyle(21);
  h_RatioStatSigLam->SetMarkerSize(.5);
  h_RatioStatSigLam->SetTitle("pp #sqrt{s}= 7 TeV");
  h_RatioStatSigLam->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_RatioStatSigLam->GetYaxis()->SetTitle("ratio = #Sigma^{0}/#Lambda");
  h_RatioStatSigLam->SetLineColor(kBlue + 2);

  TH1D *h_RatioTotSigLam = new TH1D("h_RatioTotSigLam", "ratio Sig/Lambda spectrum, total erros", ptbin + 1, ptedges);
  h_RatioTotSigLam->SetMarkerStyle(21);
  h_RatioTotSigLam->SetMarkerSize(.5);
  h_RatioTotSigLam->SetTitle("pp #sqrt{s}= 7 TeV");
  h_RatioTotSigLam->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_RatioTotSigLam->GetYaxis()->SetTitle("ratio = #Sigma^{0}/#Lambda");
  h_RatioTotSigLam->SetLineColor(kBlack);

  TH1D *h_RatioStatLamSig = new TH1D("h_RatioStatLamSig", "ratio Lambda/Sigma spectrum, stat erros", ptbin + 1, ptedges);
  h_RatioStatLamSig->SetMarkerStyle(21);
  h_RatioStatLamSig->SetMarkerSize(.5);
  h_RatioStatLamSig->SetTitle("pp #sqrt{s}= 7 TeV");
  h_RatioStatLamSig->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_RatioStatLamSig->GetYaxis()->SetTitle("ratio = #Lambda/#Sigma^{0}");
  h_RatioStatLamSig->SetLineColor(kBlue + 2);

  TH1D *h_RatioSysLamSig = new TH1D("h_RatioSysLamSig", "ratio Lambda/Sigma spectrum, sys erros", ptbin + 1, ptedges);
  h_RatioSysLamSig->SetMarkerStyle(21);
  h_RatioSysLamSig->SetMarkerSize(.5);
  h_RatioSysLamSig->SetTitle("pp #sqrt{s}= 7 TeV");
  h_RatioSysLamSig->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_RatioSysLamSig->GetYaxis()->SetTitle("ratio = #Lambda/#Sigma^{0}");
  h_RatioSysLamSig->SetLineColor(kRed);

  TH1D *h_RatioTotLamSig = new TH1D("h_RatioTotLamSig", "ratio Lambda/Sigma spectrum, total erros", ptbin + 1, ptedges);
  h_RatioTotLamSig->SetMarkerStyle(21);
  h_RatioTotLamSig->SetMarkerSize(.5);
  h_RatioTotLamSig->SetTitle("pp #sqrt{s}= 7 TeV");
  h_RatioTotLamSig->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_RatioTotLamSig->GetYaxis()->SetTitle("ratio = #Lambda/#Sigma^{0}");
  h_RatioTotLamSig->SetLineColor(kBlack);

  TH1D *h_LamSpec = new TH1D("h_LamSpec", "Lambda spectrum, total erros", ptbin + 1, ptedges);
  h_LamSpec->SetMarkerStyle(21);
  h_LamSpec->SetMarkerSize(.5);
  h_LamSpec->SetTitle("pp #sqrt{s}= 7 TeV");
  h_LamSpec->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_LamSpec->GetYaxis()->SetTitle("1/N_{inel}d^{2}N/dydp_{T} (GeV/c)^{-1}");
  h_LamSpec->SetLineColor(kBlue);

  TH1D *h_LamSpecStat = new TH1D("h_LamSpecStat", "Lambda spectrum, stat erros", ptbin + 1, ptedges);
  h_LamSpecStat->SetMarkerStyle(21);
  h_LamSpecStat->SetMarkerSize(.5);
  h_LamSpecStat->SetTitle("pp #sqrt{s}= 7 TeV");
  h_LamSpecStat->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_LamSpecStat->GetYaxis()->SetTitle("1/N_{inel}d^{2}N/dydp_{T} (GeV/c)^{-1}");
  h_LamSpecStat->SetLineColor(kBlue);

  TH1D *h_SigSpec = new TH1D("h_SigSpec", "Sigma0 spectrum, total erros", ptbin + 1, ptedges);
  h_SigSpec->SetMarkerStyle(21);
  h_SigSpec->SetMarkerSize(.5);
  h_SigSpec->SetTitle("pp #sqrt{s}= 7 TeV");
  h_SigSpec->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_SigSpec->GetYaxis()->SetTitle("1/N_{inel}d^{2}N/dydp_{T} (GeV/c)^{-1}");
  h_SigSpec->SetLineColor(kRed);
  h_SigSpec->SetMarkerColor(kRed);
  h_SigSpec->Draw("");

  for (int i = 0; i < ptbin; i++)
  {

    h_SigSpec->SetBinContent(i + 2, spectrum[i]);
    h_SigSpec->SetBinError(i + 2, spectrum_etot[i]);

    h_LamSpec->SetBinContent(i + 2, lamspectr[i]);
    h_LamSpec->SetBinError(i + 2, lamspectr_e[i]);

    h_LamSpecStat->SetBinContent(i + 2, lamspectr[i]);
    h_LamSpecStat->SetBinError(i + 2, lamspectr_estat[i]);

    h_RatioSystSigLam->SetBinContent(i + 2, ratiodmc[i]);
    h_RatioSystSigLam->SetBinError(i + 2, eratiodmc[i]);

    h_RatioStatSigLam->SetBinContent(i + 2, ratiodmc[i]);
    h_RatioStatSigLam->SetBinError(i + 2, eratioStatSigLam[i]);

    h_RatioTotSigLam->SetBinContent(i + 2, ratiodmc[i]);
    h_RatioTotSigLam->SetBinError(i + 2, eratioTotSigLam[i]);

    h_RatioTotLamSig->SetBinContent(i + 2, ratioLamSig[i]);
    h_RatioTotLamSig->SetBinError(i + 2, eratioTotLamSig[i]);

    h_RatioStatLamSig->SetBinContent(i + 2, ratioLamSig[i]);
    h_RatioStatLamSig->SetBinError(i + 2, eratioStatLamSig[i]);

    h_RatioSysLamSig->SetBinContent(i + 2, ratioLamSig[i]);
    h_RatioSysLamSig->SetBinError(i + 2, eratioSystLamSig[i]);
  }

  TCanvas *crsl = new TCanvas("crsl", "cSigLamtot", 10, 10, 600, 600);
  h_RatioSystSigLam->Draw("");
  h_RatioStatSigLam->Draw("sameE");
  h_RatioTotSigLam->Draw("sameE");

  TCanvas *clamsigtot = new TCanvas("clamsigtot", "clamsigtot", 10, 10, 600, 600);
  h_RatioTotLamSig->Draw("");
  h_RatioStatLamSig->Draw("sameE");
  h_RatioSysLamSig->Draw("sameE");

  TCanvas *cSigLam = new TCanvas("cSigLam", "cSigLam", 10, 10, 600, 600);
  cSigLam->Divide(2, 2);
  cSigLam->cd(1);
  h_RatioTotSigLam->Draw("");

  cSigLam->cd(2);
  h_RatioTotSigLam->Draw("");

  cSigLam->cd(3);
  h_RatioStatSigLam->Draw("");

  cSigLam->cd(4);
  h_RatioSystSigLam->Draw("");

  TCanvas *cLamSig = new TCanvas("cLamSig", "cLamSig", 10, 10, 600, 600);
  cLamSig->Divide(2, 2);

  cLamSig->cd(1);
  h_RatioTotLamSig->Draw("");

  cLamSig->cd(2);
  h_RatioTotLamSig->Draw("");
  h_RatioStatLamSig->Draw("sameE");
  h_RatioSysLamSig->Draw("sameE");

  cLamSig->cd(3);
  h_RatioStatLamSig->Draw("");
  cLamSig->cd(4);
  h_RatioSysLamSig->Draw("");

  TCanvas *cTotLamSig = new TCanvas("cTotLamSig", "cTotLamSig", 10, 10, 600, 600);

  h_RatioTotLamSig->Draw("E1");

  h_RatioSysLamSig->Draw("sameE1");

  TCanvas *cTotSigLam = new TCanvas("cTotSigLam", "cTotSigLam", 10, 10, 600, 600);

  h_RatioTotSigLam->Draw("E1");

  h_RatioSystSigLam->Draw("sameE1");

  //  end of ratio   */

  printf("++ writing to SigMac started   \n");

  TFile *fout0 = TFile::Open("SigMacDT.root", "RECREATE");
  h_spectrum->Write("stat");
  h_specsyst->Write("syst");
  fout0->Close();

  TFile *fout = TFile::Open("SigMac2v15DTtail.root", "RECREATE");

  h_spectrum->Write("stat");
  h_specsyst->Write("syst");

  hMassPt[0]->Write("pt0");
  hMassPt[1]->Write("pt1");
  hMassPt[2]->Write("pt2");
  hMassPt[3]->Write("pt3");
  hMassPt[4]->Write("pt4");
  hMassPt[5]->Write("pt5");
  hMassPt[6]->Write("pt6");
  hMassPt[7]->Write("pt7");

  can->Write("can");
  cSyst->Write("cSyst");

  h_Rawspectrum->Write("h_Rawspectrum");

  h_NormRawspectrum->Write("h_NormRawData");

  h_spectrum->Write("h_spectrum");
  h_spectotal->Write(" h_spectotal");
  h_spectrum_mc->Write("h_spectrum_mc");
  ratio->Write("ratio");
  ratiototal->Write("ratiototal");
  ratioStatDtMC->Write("ratioStatDtMC");

  hratio68dtmc->Write("hratio68dtmc");

  h_relstaterr->Write("h_relstaterr");

  gr_stat->Write("gr_stat");

  h_Eff->Write("h_Eff");
  h_GenTrue->Write("h_GenTrue");
  h_RecTrue->Write("h_RecTrue");

  spectrum_tot->Write("spectrum_tot");
  h_SigWidth->Write("h_SigWidth");
  h_SigMass->Write("h_SigMass");
  h_SigSignifi->Write("h_SigSignifi");
  h_SigYields->Write("h_SigYields");
  InfoHisto->Write("InfoHisto");

  hMassPtTrueRec[0]->Write("h_RecTrueSigmaPt0");
  hMassPtTrueRec[1]->Write("h_RecTrueSigmaPt1");
  hMassPtTrueRec[2]->Write("h_RecTrueSigmaPt2");
  hMassPtTrueRec[3]->Write("h_RecTrueSigmaPt3");
  hMassPtTrueRec[4]->Write("h_RecTrueSigmaPt4");
  hMassPtTrueRec[5]->Write("h_RecTrueSigmaPt5");
  hMassPtTrueRec[6]->Write("h_RecTrueSigmaPt6");
  hMassPtTrueRec[7]->Write("h_RecTrueSigmaPt7");
  hCutDT->Write("hCutDT");
  hCutMC->Write("hCutMC");
  hCutDT1->Write("hCutDT1");
  hCutMC1->Write("hCutMC1");
  hCutDT2->Write("hCutDT2");
  hCutMC2->Write("hCutMC2");
  hCutDT3->Write("hCutDT3");
  hCutMC3->Write("hCutMC3");
  hCutDT4->Write("hCutDT4");
  hCutMC4->Write("hCutMC4");
  hCutDT5->Write("hCutDT5");
  hCutMC5->Write("hCutMC5");
  hCutDT6->Write("hCutDT6");
  hCutMC6->Write("hCutMC6");
  hCutDT7->Write("hCutDT7");
  hCutMC7->Write("hCutMC7");
  hCutDT8->Write("hCutDT8");
  hCutMC8->Write("hCutMC8");
  hCutDT9->Write("hCutDT9");
  hCutMC9->Write("hCutMC9");
  hCutDT10->Write("hCutDT10");
  hCutMC10->Write("hCutMC10");
  hCutDT11->Write("hCutDT11");
  hCutMC11->Write("hCutMC11");
  hCutDT12->Write("hCutDT12");
  hCutMC12->Write("hCutMC12");
  hCutDT13->Write("hCutDT13");
  hCutMC13->Write("hCutMC13");
  hCutDT14->Write("hCutDT14");
  hCutMC14->Write("hCutMC14");
  hCutDT15->Write("hCutDT15");
  hCutMC15->Write("hCutMC15");
  h_SystLam2->Write("h_SystLam2");
  h_SystGam2->Write("h_SystGam2");
  h_SigSpec->Write("h_SigSpec");
  h_LamSpec->Write("h_LamSpec");
  h_LamSpecStat->Write("h_LamSpecStat");
  h_RatioSystSigLam->Write("h_RatioSystSigLam");
  h_RatioStatSigLam->Write("h_RatioStatSigLam");
  h_RatioTotSigLam->Write("h_RatioTotSigLam");
  h_RatioTotLamSig->Write("h_RatioTotLamSig");
  h_RatioStatLamSig->Write("h_RatioStatLamSig");
  h_RatioSysLamSig->Write("h_RatioSysLamSig");
  cSigLam->Write("cSigLam");
  cLamSig->Write("cLamSig");
  h2CutDT14->Write("h2CutDT14");
  h2CutDT14mix->Write("h2CutDT14mix");
  h2CutMC14->Write("h2CutMC14");

  printf(" write done \n");

  fout->Close();

} // end of main

//________________________________________________________________________
double PolFunction(double *x, double *par)
{

  if (reject && x[0] > yieldRange[0] && x[0] < yieldRange[1])
  {
    TF1::RejectPoint();
    return 0;
  }

  if (POLdegree == 1)
    return par[0] + par[1] * x[0];
  else if (POLdegree == 2)
    return par[0] + par[1] * x[0] + par[2] * pow(x[0], 2);
  else if (POLdegree == 3)
    return par[0] + par[1] * x[0] + par[2] * pow(x[0], 2) + par[3] * pow(x[0], 3);
  else if (POLdegree == 4)
    return par[0] + par[1] * x[0] + par[2] * pow(x[0], 2) + par[3] * pow(x[0], 3) + par[4] * pow(x[0], 4);
  else if (POLdegree == 5)
    return par[0] + par[1] * x[0] + par[2] * pow(x[0], 2) + par[3] * pow(x[0], 3) + par[4] * pow(x[0], 4) + par[5] * pow(x[0], 5);
  else
    return 0;
}
//________________________________________________________________________
double GausFunction(double *x, double *par)
{
  return (par[0] * .001) / (sqrt(2 * 3.14159 * par[2] * par[2])) * exp(-pow((x[0] - par[1]) / (sqrt(2) * par[2]), 2));
}
//________________________________________________________________________
double GausplusPol(double *x, double *par)
{
  return GausFunction(x, par) + PolFunction(x, &par[3]);
}
//________________________________________________________________________
double Voigtian(Double_t *x, Double_t *par)
{ // code taken directly from RooVoigtian in ROOT
  Bool_t _doFast;
  _doFast = kFALSE;
  Double_t Norm = par[0];
  Double_t mean = par[1];
  Double_t s = fabs(par[2]); // sigma; Gaussian sigma
  Double_t w = fabs(par[3]); // Width; BW width

  Double_t arg = x[0] - mean;
  Double_t _invRootPi = 1. / sqrt(atan2(0., -1.));

  // return constant for zero width and sigma
  if (s == 0. && w == 0.)
    return 1.;

  // Breit-Wigner for zero sigma
  if (s == 0.)
    return (Norm * 1. / (arg * arg + 0.25 * w * w));
  Double_t coef = -0.5 / (s * s);

  // Gauss for zero width
  if (w == 0.)
    return Norm * exp(coef * arg * arg);

  // actual Voigtian for non-trivial width and sigma
  Double_t c = 1. / (sqrt(2.) * s);
  Double_t a = 0.5 * c * w;
  Double_t u = c * arg;

  std::complex<Double_t> z(u, a);
  std::complex<Double_t> v(0.);

  v = RooMath::faddeeva(z);

  return Norm * c * _invRootPi * v.real();
}
//________________________________________________________________________
double VoigtplusPol(double *x, double *par)
{
  return Voigtian(x, par) + PolFunction(x, &par[4]);
}

Double_t Rapidity(Double_t pt, Double_t pz, Double_t m)
{
  Double_t energy = TMath::Sqrt(pt * pt + pz * pz + m * m);
  if (energy != TMath::Abs(pz))
    return 0.5 * TMath::Log((energy + pz) / (energy - pz));
  return TMath::Sign(1.e30, pz);
}