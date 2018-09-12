// **************************
// 
// usage:
//    root -l -b -q 'detectorSiOccupancy_Lemma.C(runNumber)'
//
// e.g. root -l -b -q 'detectorSiOccupancy_Lemma.C(500275)'
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
void detectorSiOccupancy_Lemma(int runNumber)
{

  // define input file name
  TString inputFileName;
  inputFileName.Form("/afs/cern.ch/user/a/abertoli/public/lemma/si-%d.root",runNumber);

  // define output path and make output directory
  TString plotPath;
  plotPath.Form("plotsSiOccupancy_Run%d",runNumber);
  gSystem->Exec(("mkdir -p "+plotPath));
  
  
  Int_t    iev;
  Int_t    nhits;
  Int_t    subdet[500];   //[nhits]
  Float_t  xh[500];       //[nhits]
  Float_t  yh[500];       //[nhits]
  Float_t  zh[500];       //[nhits]
  Int_t    itrack[500];   //[nhits]
  Double_t Calo_EnDep[25];
  Int_t    Calo_Time[25];

 
  TFile* inputFile = TFile::Open(inputFileName);
  TTree* inputTree = (TTree*)inputFile->Get("lemma");

  inputTree->SetBranchAddress("iev", &iev);
  inputTree->SetBranchAddress("nhits", &nhits);
  inputTree->SetBranchAddress("subdet", &subdet[0]);
  inputTree->SetBranchAddress("xh", &xh[0]);
  inputTree->SetBranchAddress("yh", &yh[0]);
  inputTree->SetBranchAddress("zh", &zh[0]);
  inputTree->SetBranchAddress("itrack", &itrack[0]);
  inputTree->SetBranchAddress("Calo_EnDep", &Calo_EnDep[0]);
  inputTree->SetBranchAddress("Calo_Time", &Calo_Time[0]);

  
  // histos for x and y occupancy
  TH1F* hist_xh_det10 = new TH1F("hist_xh_det10","hist_xh_det10",100,0.,2.);
  TH1F* hist_xh_det20 = new TH1F("hist_xh_det20","hist_xh_det20",100,0.,2.);
  TH1F* hist_xh_det30 = new TH1F("hist_xh_det30","hist_xh_det30",100,0.,10.);
  TH1F* hist_xh_det31 = new TH1F("hist_xh_det31","hist_xh_det31",100,0.,10.);
  TH1F* hist_xh_det32 = new TH1F("hist_xh_det32","hist_xh_det32",100,0.,10.);
  TH1F* hist_xh_det33 = new TH1F("hist_xh_det33","hist_xh_det33",100,0.,10.);
  TH1F* hist_xh_det34 = new TH1F("hist_xh_det34","hist_xh_det34",100,0.,10.);
  TH1F* hist_xh_det35 = new TH1F("hist_xh_det35","hist_xh_det35",100,0.,10.);
  TH1F* hist_xh_det36 = new TH1F("hist_xh_det36","hist_xh_det36",100,0.,20.);
  TH1F* hist_xh_det37 = new TH1F("hist_xh_det37","hist_xh_det37",100,0.,20.);

  TH1F* hist_yh_det10 = new TH1F("hist_yh_det10","hist_yh_det10",100,0.,2.);
  TH1F* hist_yh_det20 = new TH1F("hist_yh_det20","hist_yh_det20",100,0.,2.);
  TH1F* hist_yh_det30 = new TH1F("hist_yh_det30","hist_yh_det30",100,0.,10.);
  TH1F* hist_yh_det31 = new TH1F("hist_yh_det31","hist_yh_det31",100,0.,10.);
  TH1F* hist_yh_det32 = new TH1F("hist_yh_det32","hist_yh_det32",100,0.,10.);
  TH1F* hist_yh_det33 = new TH1F("hist_yh_det33","hist_yh_det33",100,0.,10.);
  TH1F* hist_yh_det34 = new TH1F("hist_yh_det34","hist_yh_det34",100,0.,10.);
  TH1F* hist_yh_det35 = new TH1F("hist_yh_det35","hist_yh_det35",100,0.,10.);
  TH1F* hist_yh_det36 = new TH1F("hist_yh_det36","hist_yh_det36",100,0.,20.);
  TH1F* hist_yh_det37 = new TH1F("hist_yh_det37","hist_yh_det37",100,0.,20.);

  // histos for number of hits per event
  TH1F* hist_xh_nhit_det10 = new TH1F("hist_xh_nhit_det10","hist_xh_nhit_det10",11,-0.5,10.5);
  TH1F* hist_xh_nhit_det20 = new TH1F("hist_xh_nhit_det20","hist_xh_nhit_det20",11,-0.5,10.5);
  TH1F* hist_xh_nhit_det30 = new TH1F("hist_xh_nhit_det30","hist_xh_nhit_det30",11,-0.5,10.5);
  TH1F* hist_xh_nhit_det31 = new TH1F("hist_xh_nhit_det31","hist_xh_nhit_det31",11,-0.5,10.5);
  TH1F* hist_xh_nhit_det32 = new TH1F("hist_xh_nhit_det32","hist_xh_nhit_det32",11,-0.5,10.5);
  TH1F* hist_xh_nhit_det33 = new TH1F("hist_xh_nhit_det33","hist_xh_nhit_det33",11,-0.5,10.5);
  TH1F* hist_xh_nhit_det34 = new TH1F("hist_xh_nhit_det34","hist_xh_nhit_det34",11,-0.5,10.5);
  TH1F* hist_xh_nhit_det35 = new TH1F("hist_xh_nhit_det35","hist_xh_nhit_det35",11,-0.5,10.5);
  TH1F* hist_xh_nhit_det36 = new TH1F("hist_xh_nhit_det36","hist_xh_nhit_det36",11,-0.5,10.5);
  TH1F* hist_xh_nhit_det37 = new TH1F("hist_xh_nhit_det37","hist_xh_nhit_det37",11,-0.5,10.5);

  TH1F* hist_yh_nhit_det10 = new TH1F("hist_yh_nhit_det10","hist_yh_nhit_det10",11,-0.5,10.5);
  TH1F* hist_yh_nhit_det20 = new TH1F("hist_yh_nhit_det20","hist_yh_nhit_det20",11,-0.5,10.5);
  TH1F* hist_yh_nhit_det30 = new TH1F("hist_yh_nhit_det30","hist_yh_nhit_det30",11,-0.5,10.5);
  TH1F* hist_yh_nhit_det31 = new TH1F("hist_yh_nhit_det31","hist_yh_nhit_det31",11,-0.5,10.5);
  TH1F* hist_yh_nhit_det32 = new TH1F("hist_yh_nhit_det32","hist_yh_nhit_det32",11,-0.5,10.5);
  TH1F* hist_yh_nhit_det33 = new TH1F("hist_yh_nhit_det33","hist_yh_nhit_det33",11,-0.5,10.5);
  TH1F* hist_yh_nhit_det34 = new TH1F("hist_yh_nhit_det34","hist_yh_nhit_det34",11,-0.5,10.5);
  TH1F* hist_yh_nhit_det35 = new TH1F("hist_yh_nhit_det35","hist_yh_nhit_det35",11,-0.5,10.5);
  TH1F* hist_yh_nhit_det36 = new TH1F("hist_yh_nhit_det36","hist_yh_nhit_det36",11,-0.5,10.5);
  TH1F* hist_yh_nhit_det37 = new TH1F("hist_yh_nhit_det37","hist_yh_nhit_det37",11,-0.5,10.5);


  // loop over tree entries 
  Long64_t entries = inputTree->GetEntries();
  for(Long64_t z=0; z<entries; ++z){

    inputTree->GetEntry(z);

    int nhits_xh_det10 = 0;
    int nhits_xh_det20 = 0;
    int nhits_xh_det30 = 0;
    int nhits_xh_det31 = 0;
    int nhits_xh_det32 = 0;
    int nhits_xh_det33 = 0;
    int nhits_xh_det34 = 0;
    int nhits_xh_det35 = 0;
    int nhits_xh_det36 = 0;
    int nhits_xh_det37 = 0;

    int nhits_yh_det10 = 0;
    int nhits_yh_det20 = 0;
    int nhits_yh_det30 = 0;
    int nhits_yh_det31 = 0;
    int nhits_yh_det32 = 0;
    int nhits_yh_det33 = 0;
    int nhits_yh_det34 = 0;
    int nhits_yh_det35 = 0;
    int nhits_yh_det36 = 0;
    int nhits_yh_det37 = 0;

    for(int i=0; i<nhits;i++){

      if(xh[i]>0.){
        if(subdet[i]==10) {hist_xh_det10->Fill(xh[i]); nhits_xh_det10++; }
        if(subdet[i]==20) {hist_xh_det20->Fill(xh[i]); nhits_xh_det20++; }
        if(subdet[i]==30) {hist_xh_det30->Fill(xh[i]); nhits_xh_det30++; }
        if(subdet[i]==31) {hist_xh_det31->Fill(xh[i]); nhits_xh_det31++; }
        if(subdet[i]==32) {hist_xh_det32->Fill(xh[i]); nhits_xh_det32++; }
        if(subdet[i]==33) {hist_xh_det33->Fill(xh[i]); nhits_xh_det33++; }
        if(subdet[i]==34) {hist_xh_det34->Fill(xh[i]); nhits_xh_det34++; }
        if(subdet[i]==35) {hist_xh_det35->Fill(xh[i]); nhits_xh_det35++; }
        if(subdet[i]==36) {hist_xh_det36->Fill(xh[i]); nhits_xh_det36++; }
        if(subdet[i]==37) {hist_xh_det37->Fill(xh[i]); nhits_xh_det37++; }
      }

      if(yh[i]>0.){
        if(subdet[i]==10) {hist_yh_det10->Fill(yh[i]); nhits_yh_det10++; }
        if(subdet[i]==20) {hist_yh_det20->Fill(yh[i]); nhits_yh_det20++; }
        if(subdet[i]==30) {hist_yh_det30->Fill(yh[i]); nhits_yh_det30++; }
        if(subdet[i]==31) {hist_yh_det31->Fill(yh[i]); nhits_yh_det31++; }
        if(subdet[i]==32) {hist_yh_det32->Fill(yh[i]); nhits_yh_det32++; }
        if(subdet[i]==33) {hist_yh_det33->Fill(yh[i]); nhits_yh_det33++; }
        if(subdet[i]==34) {hist_yh_det34->Fill(yh[i]); nhits_yh_det34++; }
        if(subdet[i]==35) {hist_yh_det35->Fill(yh[i]); nhits_yh_det35++; }
        if(subdet[i]==36) {hist_yh_det36->Fill(yh[i]); nhits_yh_det36++; }
        if(subdet[i]==37) {hist_yh_det37->Fill(yh[i]); nhits_yh_det37++; }
      } 
    }// end loop over number of hits 

    // fill histos with number of hits 
    hist_xh_nhit_det10->Fill(nhits_xh_det10);
    hist_xh_nhit_det20->Fill(nhits_xh_det20);
    hist_xh_nhit_det30->Fill(nhits_xh_det30);
    hist_xh_nhit_det31->Fill(nhits_xh_det31);
    hist_xh_nhit_det32->Fill(nhits_xh_det32);
    hist_xh_nhit_det33->Fill(nhits_xh_det33);
    hist_xh_nhit_det34->Fill(nhits_xh_det34);
    hist_xh_nhit_det35->Fill(nhits_xh_det35);
    hist_xh_nhit_det36->Fill(nhits_xh_det36);
    hist_xh_nhit_det37->Fill(nhits_xh_det37);

    hist_yh_nhit_det10->Fill(nhits_yh_det10);
    hist_yh_nhit_det20->Fill(nhits_yh_det20);
    hist_yh_nhit_det30->Fill(nhits_yh_det30);
    hist_yh_nhit_det31->Fill(nhits_yh_det31);
    hist_yh_nhit_det32->Fill(nhits_yh_det32);
    hist_yh_nhit_det33->Fill(nhits_yh_det33);
    hist_yh_nhit_det34->Fill(nhits_yh_det34);
    hist_yh_nhit_det35->Fill(nhits_yh_det35);
    hist_yh_nhit_det36->Fill(nhits_yh_det36);
    hist_yh_nhit_det37->Fill(nhits_yh_det37);


  }//end over tree entries 



  //plot histos for x and y occupancy 
  TCanvas* c_det10 = new TCanvas("c_det10","c_det10");
  c_det10->Divide(2,1);
  c_det10->cd(1);
  hist_xh_det10->SetTitle("x view");
  hist_xh_det10->Draw("hist");
  c_det10->cd(2);
  hist_yh_det10->SetTitle("y view");
  hist_yh_det10->Draw("hist");
  c_det10->SaveAs((plotPath + "/" + c_det10->GetName() + ".png"));
   
  TCanvas* c_det20 = new TCanvas("c_det20","c_det20");
  c_det20->Divide(2,1);
  c_det20->cd(1);
  hist_xh_det20->SetTitle("x view");
  hist_xh_det20->Draw("hist");
  c_det20->cd(2);
  hist_yh_det20->SetTitle("y view");
  hist_yh_det20->Draw("hist");
  c_det20->SaveAs((plotPath + "/" + c_det20->GetName() + ".png"));
   
  TCanvas* c_det30 = new TCanvas("c_det30","c_det30");
  c_det30->Divide(2,1);
  c_det30->cd(1);
  hist_xh_det30->SetTitle("x view");
  hist_xh_det30->Draw("hist");
  c_det30->cd(2);
  hist_yh_det30->SetTitle("y view");
  hist_yh_det30->Draw("hist");
  c_det30->SaveAs((plotPath + "/" + c_det30->GetName() + ".png"));
  
  TCanvas* c_det31 = new TCanvas("c_det31","c_det31");
  c_det31->Divide(2,1);
  c_det31->cd(1);
  hist_xh_det31->SetTitle("x view");
  hist_xh_det31->Draw("hist");
  c_det31->cd(2);
  hist_yh_det31->SetTitle("y view");
  hist_yh_det31->Draw("hist");
  c_det31->SaveAs((plotPath + "/" + c_det31->GetName() + ".png"));

  TCanvas* c_det32 = new TCanvas("c_det32","c_det32");
  c_det32->Divide(2,1);
  c_det32->cd(1);
  hist_xh_det32->SetTitle("x view");
  hist_xh_det32->Draw("hist");
  c_det32->cd(2);
  hist_yh_det32->SetTitle("y view");
  hist_yh_det32->Draw("hist");
  c_det32->SaveAs((plotPath + "/" + c_det32->GetName() + ".png"));
  
  TCanvas* c_det33 = new TCanvas("c_det33","c_det33");
  c_det33->Divide(2,1);
  c_det33->cd(1);
  hist_xh_det33->SetTitle("x view");
  hist_xh_det33->Draw("hist");
  c_det33->cd(2);
  hist_yh_det33->SetTitle("y view");
  hist_yh_det33->Draw("hist");
  c_det33->SaveAs((plotPath + "/" + c_det33->GetName() + ".png"));
  
  TCanvas* c_det34 = new TCanvas("c_det34","c_det34");
  c_det34->Divide(2,1);
  c_det34->cd(1);
  hist_xh_det34->SetTitle("x view");
  hist_xh_det34->Draw("hist");
  c_det34->cd(2);
  hist_yh_det34->SetTitle("y view");
  hist_yh_det34->Draw("hist");
  c_det34->SaveAs((plotPath + "/" + c_det34->GetName() + ".png"));
  
  TCanvas* c_det35 = new TCanvas("c_det35","c_det35");
  c_det35->Divide(2,1);
  c_det35->cd(1);
  hist_xh_det35->SetTitle("x view");
  hist_xh_det35->Draw("hist");
  c_det35->cd(2);
  hist_yh_det35->SetTitle("y view");
  hist_yh_det35->Draw("hist");
  c_det35->SaveAs((plotPath + "/" + c_det35->GetName() + ".png"));
  
  TCanvas* c_det36 = new TCanvas("c_det36","c_det36");
  c_det36->Divide(2,1);
  c_det36->cd(1);
  hist_xh_det36->SetTitle("x view");
  hist_xh_det36->Draw("hist");
  c_det36->cd(2);
  hist_yh_det36->SetTitle("y view");
  hist_yh_det36->Draw("hist");
  c_det36->SaveAs((plotPath + "/" + c_det36->GetName() + ".png"));
  
  TCanvas* c_det37 = new TCanvas("c_det37","c_det37");
  c_det37->Divide(2,1);
  c_det37->cd(1);
  hist_xh_det37->SetTitle("x view");
  hist_xh_det37->Draw("hist");
  c_det37->cd(2);
  hist_yh_det37->SetTitle("y view");
  hist_yh_det37->Draw("hist");
  c_det37->SaveAs((plotPath + "/" + c_det37->GetName() + ".png"));
  
    
  // plot  histos for number of hits per event
  TCanvas* c_nhits_det10 = new TCanvas("c_nhits_det10","c_nhits_det10");
  c_nhits_det10->Divide(2,1);
  c_nhits_det10->cd(1);
  hist_xh_nhit_det10->SetTitle("x view");
  hist_xh_nhit_det10->GetXaxis()->SetTitle("n hits");
  hist_xh_nhit_det10->Draw("hist");
  c_nhits_det10->cd(2);
  hist_yh_nhit_det10->SetTitle("y view");
  hist_yh_nhit_det10->GetXaxis()->SetTitle("n hits");
  hist_yh_nhit_det10->Draw("hist");
  c_nhits_det10->SaveAs((plotPath + "/" + c_nhits_det10->GetName() + ".png"));


  TCanvas* c_nhits_det20 = new TCanvas("c_nhits_det20","c_nhits_det20");
  c_nhits_det20->Divide(2,1);
  c_nhits_det20->cd(1);
  hist_xh_nhit_det20->SetTitle("x view");
  hist_xh_nhit_det20->GetXaxis()->SetTitle("n hits");
  hist_xh_nhit_det20->Draw("hist");
  c_nhits_det20->cd(2);
  hist_yh_nhit_det20->SetTitle("y view");
  hist_yh_nhit_det20->GetXaxis()->SetTitle("n hits");
  hist_yh_nhit_det20->Draw("hist");
  c_nhits_det20->SaveAs((plotPath + "/" + c_nhits_det20->GetName() + ".png"));


  TCanvas* c_nhits_det30 = new TCanvas("c_nhits_det30","c_nhits_det30");
  c_nhits_det30->Divide(2,1);
  c_nhits_det30->cd(1);
  hist_xh_nhit_det30->SetTitle("x view");
  hist_xh_nhit_det30->GetXaxis()->SetTitle("n hits");
  hist_xh_nhit_det30->Draw("hist");
  c_nhits_det30->cd(2);
  hist_yh_nhit_det30->SetTitle("y view");
  hist_yh_nhit_det30->GetXaxis()->SetTitle("n hits");
  hist_yh_nhit_det30->Draw("hist");
  c_nhits_det30->SaveAs((plotPath + "/" + c_nhits_det30->GetName() + ".png"));


  TCanvas* c_nhits_det31 = new TCanvas("c_nhits_det31","c_nhits_det31");
  c_nhits_det31->Divide(2,1);
  c_nhits_det31->cd(1);
  hist_xh_nhit_det31->SetTitle("x view");
  hist_xh_nhit_det31->GetXaxis()->SetTitle("n hits");
  hist_xh_nhit_det31->Draw("hist");
  c_nhits_det31->cd(2);
  hist_yh_nhit_det31->SetTitle("y view");
  hist_yh_nhit_det31->GetXaxis()->SetTitle("n hits");
  hist_yh_nhit_det31->Draw("hist");
  c_nhits_det31->SaveAs((plotPath + "/" + c_nhits_det31->GetName() + ".png"));


  TCanvas* c_nhits_det32 = new TCanvas("c_nhits_det32","c_nhits_det32");
  c_nhits_det32->Divide(2,1);
  c_nhits_det32->cd(1);
  hist_xh_nhit_det32->SetTitle("x view");
  hist_xh_nhit_det32->GetXaxis()->SetTitle("n hits");
  hist_xh_nhit_det32->Draw("hist");
  c_nhits_det32->cd(2);
  hist_yh_nhit_det32->SetTitle("y view");
  hist_yh_nhit_det32->GetXaxis()->SetTitle("n hits");
  hist_yh_nhit_det32->Draw("hist");
  c_nhits_det32->SaveAs((plotPath + "/" + c_nhits_det32->GetName() + ".png"));


  TCanvas* c_nhits_det33 = new TCanvas("c_nhits_det33","c_nhits_det33");
  c_nhits_det33->Divide(2,1);
  c_nhits_det33->cd(1);
  hist_xh_nhit_det33->SetTitle("x view");
  hist_xh_nhit_det33->GetXaxis()->SetTitle("n hits");
  hist_xh_nhit_det33->Draw("hist");
  c_nhits_det33->cd(2);
  hist_yh_nhit_det33->SetTitle("y view");
  hist_yh_nhit_det33->GetXaxis()->SetTitle("n hits");
  hist_yh_nhit_det33->Draw("hist");
  c_nhits_det33->SaveAs((plotPath + "/" + c_nhits_det33->GetName() + ".png"));


  TCanvas* c_nhits_det34 = new TCanvas("c_nhits_det34","c_nhits_det34");
  c_nhits_det34->Divide(2,1);
  c_nhits_det34->cd(1);
  hist_xh_nhit_det34->SetTitle("x view");
  hist_xh_nhit_det34->GetXaxis()->SetTitle("n hits");
  hist_xh_nhit_det34->Draw("hist");
  c_nhits_det34->cd(2);
  hist_yh_nhit_det34->SetTitle("y view");
  hist_yh_nhit_det34->GetXaxis()->SetTitle("n hits");
  hist_yh_nhit_det34->Draw("hist");
  c_nhits_det34->SaveAs((plotPath + "/" + c_nhits_det34->GetName() + ".png"));


  TCanvas* c_nhits_det35 = new TCanvas("c_nhits_det35","c_nhits_det35");
  c_nhits_det35->Divide(2,1);
  c_nhits_det35->cd(1);
  hist_xh_nhit_det35->SetTitle("x view");
  hist_xh_nhit_det35->GetXaxis()->SetTitle("n hits");
  hist_xh_nhit_det35->Draw("hist");
  c_nhits_det35->cd(2);
  hist_yh_nhit_det35->SetTitle("y view");
  hist_yh_nhit_det35->GetXaxis()->SetTitle("n hits");
  hist_yh_nhit_det35->Draw("hist");
  c_nhits_det35->SaveAs((plotPath + "/" + c_nhits_det35->GetName() + ".png"));


  TCanvas* c_nhits_det36 = new TCanvas("c_nhits_det36","c_nhits_det36");
  c_nhits_det36->Divide(2,1);
  c_nhits_det36->cd(1);
  hist_xh_nhit_det36->SetTitle("x view");
  hist_xh_nhit_det36->GetXaxis()->SetTitle("n hits");
  hist_xh_nhit_det36->Draw("hist");
  c_nhits_det36->cd(2);
  hist_yh_nhit_det36->SetTitle("y view");
  hist_yh_nhit_det36->GetXaxis()->SetTitle("n hits");
  hist_yh_nhit_det36->Draw("hist");
  c_nhits_det36->SaveAs((plotPath + "/" + c_nhits_det36->GetName() + ".png"));


  TCanvas* c_nhits_det37 = new TCanvas("c_nhits_det37","c_nhits_det37");
  c_nhits_det37->Divide(2,1);
  c_nhits_det37->cd(1);
  hist_xh_nhit_det37->SetTitle("x view");
  hist_xh_nhit_det37->GetXaxis()->SetTitle("n hits");
  hist_xh_nhit_det37->Draw("hist");
  c_nhits_det37->cd(2);
  hist_yh_nhit_det37->SetTitle("y view");
  hist_yh_nhit_det37->GetXaxis()->SetTitle("n hits");
  hist_yh_nhit_det37->Draw("hist");
  c_nhits_det37->SaveAs((plotPath + "/" + c_nhits_det37->GetName() + ".png"));

  

  cout<<" Plots done! =) "<<endl; 

}
