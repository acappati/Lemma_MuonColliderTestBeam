// **************************
// 
// usage:
//    root -l -b -q plotProfileTrend_noTargetRuns.C
//
// **************************


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <iomanip>

#include "TString.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TAxis.h"


using namespace std;



// do the fit function 
//void doTheFits()
void plotProfileTrend_noTargetRuns()
{

  // number of runs to be analyzed
  // current of the magnet +580A
  int runNumberPlus580A[6]  = {500332,500327,500325,500320,500319,500314};
  // current of the magnet -580A
  int runNumberMinus580A[6] = {500330,500329,500323,500322,500317,500316};

  float beamEnergy[6] = {16.,18.,20.,22.,24.,26.};
  string sbeamEnergy[6] = {"16GeV","18GeV","20GeV","22GeV","24GeV","26GeV"};


  // define output path and make output directory
  TString plotPathdet37 = "plotFits_det37";
  gSystem->Exec(("mkdir -p "+plotPathdet37));  
  TString plotPathdet36 = "plotFits_det36";
  gSystem->Exec(("mkdir -p "+plotPathdet36));  



  // current of the magnet +580A
  for(int i=0; i<6; i++){
  
    TString inputFileName;
    inputFileName.Form("plotsSiOccupancy_Run%d.root",runNumberPlus580A[i]);

    TFile* inputFile = TFile::Open(inputFileName);
    TH1F* hist_det37 = (TH1F*)inputFile->Get("hist_xh_det37");

    TF1* gfit = new TF1("gfit","gaus");
    gfit->SetParameter(1,hist_det37->GetMean());
    gfit->SetParameter(2,4.);
    gfit->SetLineColor(kRed);
  
    hist_det37->Fit("gfit");

    TCanvas* c = new TCanvas(); 
    c->cd();
    hist_det37->Draw("hist");
    //gfit->Draw("samel");

    //gStyle->SetOptStat(0);
    //gStyle->SetOptFit(1111);

    TPaveText* pvtext = new TPaveText(0.85,0.55,0.96,0.65,"brNDC");
    pvtext->AddText("+580A");
    pvtext->AddText(sbeamEnergy[i].c_str());
    pvtext->SetTextSize(0.03);
    pvtext->SetFillColor(kWhite);
    pvtext->SetBorderSize(1);
    pvtext->SetTextFont(40);
    pvtext->SetTextSize(0.037);
    pvtext->SetTextFont(42);
    pvtext->SetTextAlign(12);
    pvtext->Draw(); 

    c->Update();

    c->SaveAs(Form(plotPathdet37 + "/" + hist_det37->GetName() + "_run%d.png",runNumberPlus580A[i]));

  }



  // current of the magnet -580A
  for(int i=0; i<6; i++){
  
    TString inputFileName;
    inputFileName.Form("plotsSiOccupancy_Run%d.root",runNumberMinus580A[i]);

    TFile* inputFile = TFile::Open(inputFileName);
    TH1F* hist_det36 = (TH1F*)inputFile->Get("hist_xh_det36");

    TF1* gfit = new TF1("gfit","gaus");
    gfit->SetParameter(1,hist_det36->GetMean());
    gfit->SetParameter(2,4.);
    gfit->SetLineColor(kRed);
  
    hist_det36->Fit("gfit");

    TCanvas* c = new TCanvas(); 
    c->cd();
    hist_det36->Draw("hist");
    //gfit->Draw("samel");

    //gStyle->SetOptStat(0);
    //gStyle->SetOptFit(1111);

    TPaveText* pvtext = new TPaveText(0.85,0.55,0.96,0.65,"brNDC");
    pvtext->AddText("-580A");
    pvtext->AddText(sbeamEnergy[i].c_str());
    pvtext->SetTextSize(0.03);
    pvtext->SetFillColor(kWhite);
    pvtext->SetBorderSize(1);
    pvtext->SetTextFont(40);
    pvtext->SetTextSize(0.037);
    pvtext->SetTextFont(42);
    pvtext->SetTextAlign(12);
    pvtext->Draw(); 

    c->Update();

    c->SaveAs(Form(plotPathdet36 + "/" + hist_det36->GetName() + "_run%d.png",runNumberMinus580A[i]));

  }


}




// main function
// void plotProfileTrend_noTargetRuns()
// {

//   //doTheFits();


//   // current of the magnet +580A
//   int runNumberPlus580A[6]  = {500332,500327,500325,500320,500319,500314};
//   // current of the magnet -580A
//   int runNumberMinus580A[6] = {500330,500329,500323,500322,500317,500316};

//   float beamEnergy[6] = {16.,18.,20.,22.,24.,26.};

//   float meanPlus[6]      = {4.57,};
//   float meanPlus_err[6]  = {1.16,};
//   float sigmaPlus[6]     = {2.37,};
//   float sigmaPlus_err[6] = {1.96,};

//   float meanMinus[6] = {};
//   float meanMinus_err[6] = {};
//   float sigmaMinus[6] = {};
//   float sigmaMinus_err[6] = {};

// }
