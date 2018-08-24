// **************************
// 
// usage:
//    root -l -b -q plotVariables.C
//
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



// main function
void plotVariables()
{

  // define input file name
  TString inputFileName = "/afs/cern.ch/user/a/abertoli/public/lemma/LemmaMC2018_MuMu_T_FieldMf_N10000_calo-converted.root";
  //inputFileName.Form("/afs/cern.ch/user/a/abertoli/public/lemma/si-%d.root",runNumber);

  // define output path and make output directory
  TString plotPath = "plots_variables";
  //plotPath.Form("plotsSiOccupancy_Run%d",runNumber);
  gSystem->Exec(("mkdir -p "+plotPath));
  
  
   Double_t chi2m;
   Double_t x_pos_mum[12];
   Double_t x_pos_mum_err[12];
   Double_t z_x_pos_mum[12];
   Double_t x_pos_DT_mum[8];
   Double_t z_pos_DT_mum[8];
   Double_t p_mum;
   Double_t p_mup;
   Double_t chi2p;
   Int_t    good_mup_trk;
   Double_t x_pos_mup[12];
   Double_t x_pos_mup_err[12];
   Double_t z_x_pos_mup[12];
   Double_t x_pos_DT_mup[8];
   Double_t z_pos_DT_mup[8];
   Int_t    subdet[100];
   Int_t    itrack[100];
   Double_t xh[100];
   Double_t yh[100];
   Double_t zh[100];
   Int_t    nhits;
   Double_t Calo_EnDep[100];
   Int_t    event_type;
   Double_t gen_pos_mum[7];
   Double_t gen_pos_mup[7];

   
  TFile* inputFile = TFile::Open(inputFileName);
  TTree* inputTree = (TTree*)inputFile->Get("tb_output");

  inputTree->SetBranchAddress("chi2m",	        &chi2m);	     
  inputTree->SetBranchAddress("x_pos_mum",      &x_pos_mum[0]); 
  inputTree->SetBranchAddress("x_pos_mum_err",  &x_pos_mum_err[0]);
  inputTree->SetBranchAddress("z_x_pos_mum",    &z_x_pos_mum[0]);
  inputTree->SetBranchAddress("x_pos_DT_mum",   &x_pos_DT_mum[0]);
  inputTree->SetBranchAddress("z_pos_DT_mum",   &z_pos_DT_mum[0]);
  inputTree->SetBranchAddress("p_mum",          &p_mum);	     
  inputTree->SetBranchAddress("p_mup",          &p_mup);	     
  inputTree->SetBranchAddress("chi2p",          &chi2p);	     
  inputTree->SetBranchAddress("good_mup_trk",   &good_mup_trk);  
  inputTree->SetBranchAddress("x_pos_mup",      &x_pos_mup[0]); 
  inputTree->SetBranchAddress("x_pos_mup_err",  &x_pos_mup_err[0]);
  inputTree->SetBranchAddress("z_x_pos_mup",    &z_x_pos_mup[0]);
  inputTree->SetBranchAddress("x_pos_DT_mup",   &x_pos_DT_mup[0]);
  inputTree->SetBranchAddress("z_pos_DT_mup",   &z_pos_DT_mup[0]);
  inputTree->SetBranchAddress("subdet",         &subdet[0]);   
  inputTree->SetBranchAddress("itrack",         &itrack[0]);   
  inputTree->SetBranchAddress("xh",             &xh[0]);	     
  inputTree->SetBranchAddress("yh",             &yh[0]);	     
  inputTree->SetBranchAddress("zh",             &zh[0]);	     
  inputTree->SetBranchAddress("nhits",          &nhits);	     
  inputTree->SetBranchAddress("Calo_EnDep",     &Calo_EnDep[0]);
  inputTree->SetBranchAddress("event_type",     &event_type);       
  inputTree->SetBranchAddress("gen_pos_mum",    &gen_pos_mum[0]);
  inputTree->SetBranchAddress("gen_pos_mup",    &gen_pos_mup[0]);


  // def histos 
  TH1F* hist_pTot = new TH1F("hist_pTot","hist_pTot",200,40000.,47000.);
  

  // loop over tree entries 
  Long64_t entries = inputTree->GetEntries();
  for(Long64_t z=0; z<entries; ++z){

    inputTree->GetEntry(z);

    if( p_mum > 0. && p_mup > 0.) {hist_pTot->Fill(p_mum + p_mup);}


  }//end over tree entries 



  //plot histos
  TCanvas* c_pTot = new TCanvas("c_pTot","c_pTot");
  c_pTot->cd();
  hist_pTot->SetTitle("");
  hist_pTot->GetXaxis()->SetTitle("p_{#mu^{+}} + p_{#mu^{-}}");
  hist_pTot->Draw("hist");
  c_pTot->SaveAs((plotPath + "/" + c_pTot->GetName() + ".png"));
   
     

  cout<<" Plots done! =) "<<endl; 

}
