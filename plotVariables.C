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
#include "TLegend.h"
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
  TString plotPath = "LemmaVariables_MC";
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
  TH1F* hist_pMuPlus     = new TH1F("hist_pMuPlus","hist_pMuPlus",200,16000.,30000.);
  TH1F* hist_pMuMinus    = new TH1F("hist_pMuMinus","hist_pMuMinus",200,16000.,30000.);
  TH1F* hist_pTot        = new TH1F("hist_pTot","hist_pTot",200,40000.,47000.);
  TH1F* hist_chi2MuPlus  = new TH1F("hist_chi2MuPlus","hist_chi2MuPlus",50,0.,10.);
  TH1F* hist_chi2MuMinus = new TH1F("hist_chi2MuMinus","hist_chi2MuMinus",50,0.,10.);
  // TH1F* hist_ThetaMuPlus  = new TH1F("hist_ThetaMuPlus","hist_ThetaMuPlus",10,0.,10.);    //angle in bending plane
  // TH1F* hist_ThetaMuMinus = new TH1F("hist_ThetaMuMinus","hist_ThetaMuMinus",10,0.,10.);  //angle in bending plane
  

  TH1F* hist_xh_det10_MuPlus = new TH1F("hist_xh_det10_MuPlus","hist_xh_det10_MuPlus",100,   -9.5,   9.5);
  TH1F* hist_xh_det20_MuPlus = new TH1F("hist_xh_det20_MuPlus","hist_xh_det20_MuPlus",100,   -9.5,   9.5);
  TH1F* hist_xh_det30_MuPlus = new TH1F("hist_xh_det30_MuPlus","hist_xh_det30_MuPlus",100,  -46.5,  46.5);
  TH1F* hist_xh_det31_MuPlus = new TH1F("hist_xh_det31_MuPlus","hist_xh_det31_MuPlus",100,  -46.5,  46.5);
  TH1F* hist_xh_det32_MuPlus = new TH1F("hist_xh_det32_MuPlus","hist_xh_det32_MuPlus",100,   40.0, 120.0);
  TH1F* hist_xh_det33_MuPlus = new TH1F("hist_xh_det33_MuPlus","hist_xh_det33_MuPlus",100, -120.0, -40.0);
  TH1F* hist_xh_det34_MuPlus = new TH1F("hist_xh_det34_MuPlus","hist_xh_det34_MuPlus",100,   93.5, 186.5);
  TH1F* hist_xh_det35_MuPlus = new TH1F("hist_xh_det35_MuPlus","hist_xh_det35_MuPlus",100, -186.5, -93.5);
  TH1F* hist_xh_det36_MuPlus = new TH1F("hist_xh_det36_MuPlus","hist_xh_det36_MuPlus",100,  150.0, 330.0);
  TH1F* hist_xh_det37_MuPlus = new TH1F("hist_xh_det37_MuPlus","hist_xh_det37_MuPlus",100, -330.0,-150.0);

  TH1F* hist_xh_det10_MuMinus = new TH1F("hist_xh_det10_MuMinus","hist_xh_det10_MuMinus",100,   -9.5,   9.5);
  TH1F* hist_xh_det20_MuMinus = new TH1F("hist_xh_det20_MuMinus","hist_xh_det20_MuMinus",100,   -9.5,   9.5);
  TH1F* hist_xh_det30_MuMinus = new TH1F("hist_xh_det30_MuMinus","hist_xh_det30_MuMinus",100,  -46.5,  46.5);
  TH1F* hist_xh_det31_MuMinus = new TH1F("hist_xh_det31_MuMinus","hist_xh_det31_MuMinus",100,  -46.5,  46.5);
  TH1F* hist_xh_det32_MuMinus = new TH1F("hist_xh_det32_MuMinus","hist_xh_det32_MuMinus",100,   40.0, 120.0);
  TH1F* hist_xh_det33_MuMinus = new TH1F("hist_xh_det33_MuMinus","hist_xh_det33_MuMinus",100, -120.0, -40.0);
  TH1F* hist_xh_det34_MuMinus = new TH1F("hist_xh_det34_MuMinus","hist_xh_det34_MuMinus",100,   93.5, 186.5);
  TH1F* hist_xh_det35_MuMinus = new TH1F("hist_xh_det35_MuMinus","hist_xh_det35_MuMinus",100, -186.5, -93.5);
  TH1F* hist_xh_det36_MuMinus = new TH1F("hist_xh_det36_MuMinus","hist_xh_det36_MuMinus",100,  150.0, 330.0);
  TH1F* hist_xh_det37_MuMinus = new TH1F("hist_xh_det37_MuMinus","hist_xh_det37_MuMinus",100, -330.0,-150.0);


  // loop over tree entries 
  Long64_t entries = inputTree->GetEntries();
  for(Long64_t z=0; z<entries; ++z){

    inputTree->GetEntry(z);

    if( p_mup > 0.) {
      hist_pMuPlus->Fill(p_mup);      //momentum for mu plus
      hist_chi2MuPlus->Fill(chi2p);   //chi2 for mu plus tracks
    }
    if( p_mum > 0.) {
      hist_pMuMinus->Fill(p_mum);     //momentum for mu minus
      hist_chi2MuMinus->Fill(chi2m);  //chi2 for mu minus tracks
    }
    if( p_mum > 0. && p_mup > 0.) {hist_pTot->Fill(p_mum + p_mup);}  //total momentum 

    // histos for hits in Si det
    for(int i=0;i<100;i++){
    
      if(xh[i] == -9990 ){continue;}   // skip empty events

      // mu plus events
      if(itrack[i] == 1){             
	if(subdet[i] == 10) {hist_xh_det10_MuPlus->Fill(xh[i]);}
	if(subdet[i] == 20) {hist_xh_det20_MuPlus->Fill(xh[i]);}
	if(subdet[i] == 30) {hist_xh_det30_MuPlus->Fill(xh[i]);}
	if(subdet[i] == 31) {hist_xh_det31_MuPlus->Fill(xh[i]);}
	if(subdet[i] == 32) {hist_xh_det32_MuPlus->Fill(xh[i]);}
	if(subdet[i] == 33) {hist_xh_det33_MuPlus->Fill(xh[i]);}
	if(subdet[i] == 34) {hist_xh_det34_MuPlus->Fill(xh[i]);}
	if(subdet[i] == 35) {hist_xh_det35_MuPlus->Fill(xh[i]);}
	if(subdet[i] == 36) {hist_xh_det36_MuPlus->Fill(xh[i]);}
	if(subdet[i] == 37) {hist_xh_det37_MuPlus->Fill(xh[i]);}
      }
      // mu minus events
      else if(itrack[i] == -1){             
	if(subdet[i] == 10) {hist_xh_det10_MuMinus->Fill(xh[i]);}
	if(subdet[i] == 20) {hist_xh_det20_MuMinus->Fill(xh[i]);}
	if(subdet[i] == 30) {hist_xh_det30_MuMinus->Fill(xh[i]);}
	if(subdet[i] == 31) {hist_xh_det31_MuMinus->Fill(xh[i]);}
	if(subdet[i] == 32) {hist_xh_det32_MuMinus->Fill(xh[i]);}
	if(subdet[i] == 33) {hist_xh_det33_MuMinus->Fill(xh[i]);}
	if(subdet[i] == 34) {hist_xh_det34_MuMinus->Fill(xh[i]);}
	if(subdet[i] == 35) {hist_xh_det35_MuMinus->Fill(xh[i]);}
	if(subdet[i] == 36) {hist_xh_det36_MuMinus->Fill(xh[i]);}
	if(subdet[i] == 37) {hist_xh_det37_MuMinus->Fill(xh[i]);}
      }
      else continue;
    }
    
  }//end over tree entries 



  //plot histos
  TCanvas* c_pMuPlus = new TCanvas("c_pMuPlus","c_pMuPlus");
  c_pMuPlus->cd();
  hist_pMuPlus->SetTitle("p #mu^{+}");
  hist_pMuPlus->GetXaxis()->SetTitle("p #mu^{+}");
  hist_pMuPlus->Draw("hist");
  c_pMuPlus->SaveAs((plotPath + "/" + hist_pMuPlus->GetName() + ".png"));
  c_pMuPlus->SaveAs((plotPath + "/" + hist_pMuPlus->GetName() + ".pdf"));

  TCanvas* c_pMuMinus = new TCanvas("c_pMuMinus","c_pMuMinus");
  c_pMuMinus->cd();
  hist_pMuMinus->SetTitle("p #mu^{-}");
  hist_pMuMinus->GetXaxis()->SetTitle("p #mu^{-}");
  hist_pMuMinus->Draw("hist");
  c_pMuMinus->SaveAs((plotPath + "/" + hist_pMuMinus->GetName() + ".png"));
  c_pMuMinus->SaveAs((plotPath + "/" + hist_pMuMinus->GetName() + ".pdf"));
  
  TCanvas* c_pTot = new TCanvas("c_pTot","c_pTot");
  c_pTot->cd();
  hist_pTot->SetTitle("p #mu^{+} + p #mu^{-}");
  hist_pTot->GetXaxis()->SetTitle("p #mu^{+} + p #mu^{-}");
  hist_pTot->Draw("hist");
  c_pTot->SaveAs((plotPath + "/" + hist_pTot->GetName() + ".png"));
  c_pTot->SaveAs((plotPath + "/" + hist_pTot->GetName() + ".pdf"));

  TCanvas* c_chi2MuPlus = new TCanvas("c_chi2MuPlus","c_chi2MuPlus");
  c_chi2MuPlus->cd();
  hist_chi2MuPlus->SetTitle("#Chi^{2} #mu^{+}");
  hist_chi2MuPlus->GetXaxis()->SetTitle("#Chi^{2} #mu^{+}");
  hist_chi2MuPlus->Draw("hist");
  c_chi2MuPlus->SaveAs((plotPath + "/" + hist_chi2MuPlus->GetName() + ".png"));
  c_chi2MuPlus->SaveAs((plotPath + "/" + hist_chi2MuPlus->GetName() + ".pdf"));
  
  TCanvas* c_chi2MuMinus = new TCanvas("c_chi2MuMinus","c_chi2MuMinus");
  c_chi2MuMinus->cd();
  hist_chi2MuMinus->SetTitle("#Chi^{2} #mu^{-}");
  hist_chi2MuMinus->GetXaxis()->SetTitle("#Chi^{2} #mu^{-}");
  hist_chi2MuMinus->Draw("hist");
  c_chi2MuMinus->SaveAs((plotPath + "/" + hist_chi2MuMinus->GetName() + ".png"));
  c_chi2MuMinus->SaveAs((plotPath + "/" + hist_chi2MuMinus->GetName() + ".pdf"));
   
  // xh in det10
  TCanvas* c_det10 = new TCanvas("c_det10","c_det10");
  c_det10->cd();
  hist_xh_det10_MuPlus->SetTitle("xh in det10");
  hist_xh_det10_MuPlus->GetXaxis()->SetTitle("mm");
  hist_xh_det10_MuPlus->SetLineColor(kRed);
  hist_xh_det10_MuPlus->SetMaximum(1.05 * max(hist_xh_det10_MuMinus->GetMaximum(),hist_xh_det10_MuPlus->GetMaximum()));
  hist_xh_det10_MuPlus->Draw("hist");
  hist_xh_det10_MuMinus->SetTitle("");
  hist_xh_det10_MuMinus->SetLineColor(kBlue);
  hist_xh_det10_MuMinus->Draw("histsame");
  gStyle->SetOptStat(0);

  TLegend* l_det10 = new TLegend(0.84,0.76,0.94,0.87);
  l_det10->AddEntry(hist_xh_det10_MuPlus, "#mu^{+}","l");
  l_det10->AddEntry(hist_xh_det10_MuMinus,"#mu^{-}","l");
  l_det10->SetFillColor(kWhite);
  l_det10->SetLineColor(kBlack);
  l_det10->SetTextFont(43);
  l_det10->SetTextSize(20);
  l_det10->Draw();
  c_det10->Update();

  c_det10->SaveAs((plotPath + "/" + c_det10->GetName() + ".png"));    
  c_det10->SaveAs((plotPath + "/" + c_det10->GetName() + ".pdf"));     

  // xh in det20
  TCanvas* c_det20 = new TCanvas("c_det20","c_det20");
  c_det20->cd();
  hist_xh_det20_MuPlus->SetTitle("xh in det20");
  hist_xh_det20_MuPlus->GetXaxis()->SetTitle("mm");
  hist_xh_det20_MuPlus->SetLineColor(kRed);
  hist_xh_det20_MuPlus->SetMaximum(1.05 * max(hist_xh_det20_MuMinus->GetMaximum(),hist_xh_det20_MuPlus->GetMaximum()));
  hist_xh_det20_MuPlus->Draw("hist");
  hist_xh_det20_MuMinus->SetTitle("");
  hist_xh_det20_MuMinus->SetLineColor(kBlue);
  hist_xh_det20_MuMinus->Draw("histsame");
  gStyle->SetOptStat(0);

  TLegend* l_det20 = new TLegend(0.84,0.76,0.94,0.87);
  l_det20->AddEntry(hist_xh_det20_MuPlus, "#mu^{+}","l");
  l_det20->AddEntry(hist_xh_det20_MuMinus,"#mu^{-}","l");
  l_det20->SetFillColor(kWhite);
  l_det20->SetLineColor(kBlack);
  l_det20->SetTextFont(43);
  l_det20->SetTextSize(20);
  l_det20->Draw();
  c_det20->Update();

  c_det20->SaveAs((plotPath + "/" + c_det20->GetName() + ".png"));    
  c_det20->SaveAs((plotPath + "/" + c_det20->GetName() + ".pdf"));     

  // xh in det30
  TCanvas* c_det30 = new TCanvas("c_det30","c_det30");
  c_det30->cd();
  hist_xh_det30_MuPlus->SetTitle("xh in det30");
  hist_xh_det30_MuPlus->GetXaxis()->SetTitle("mm");
  hist_xh_det30_MuPlus->SetLineColor(kRed);
  hist_xh_det30_MuPlus->SetMaximum(1.05 * max(hist_xh_det30_MuMinus->GetMaximum(),hist_xh_det30_MuPlus->GetMaximum()));
  hist_xh_det30_MuPlus->Draw("hist");
  hist_xh_det30_MuMinus->SetTitle("");
  hist_xh_det30_MuMinus->SetLineColor(kBlue);
  hist_xh_det30_MuMinus->Draw("histsame");
  gStyle->SetOptStat(0);

  TLegend* l_det30 = new TLegend(0.84,0.76,0.94,0.87);
  l_det30->AddEntry(hist_xh_det30_MuPlus, "#mu^{+}","l");
  l_det30->AddEntry(hist_xh_det30_MuMinus,"#mu^{-}","l");
  l_det30->SetFillColor(kWhite);
  l_det30->SetLineColor(kBlack);
  l_det30->SetTextFont(43);
  l_det30->SetTextSize(20);
  l_det30->Draw();
  c_det30->Update();

  c_det30->SaveAs((plotPath + "/" + c_det30->GetName() + ".png"));    
  c_det30->SaveAs((plotPath + "/" + c_det30->GetName() + ".pdf")); 

  // xh in det31
  TCanvas* c_det31 = new TCanvas("c_det31","c_det31");
  c_det31->cd();
  hist_xh_det31_MuPlus->SetTitle("xh in det31");
  hist_xh_det31_MuPlus->GetXaxis()->SetTitle("mm");
  hist_xh_det31_MuPlus->SetLineColor(kRed);
  hist_xh_det31_MuPlus->SetMaximum(1.05 * max(hist_xh_det31_MuMinus->GetMaximum(),hist_xh_det31_MuPlus->GetMaximum()));
  hist_xh_det31_MuPlus->Draw("hist");
  hist_xh_det31_MuMinus->SetTitle("");
  hist_xh_det31_MuMinus->SetLineColor(kBlue);
  hist_xh_det31_MuMinus->Draw("histsame");
  gStyle->SetOptStat(0);

  TLegend* l_det31 = new TLegend(0.84,0.76,0.94,0.87);
  l_det31->AddEntry(hist_xh_det31_MuPlus, "#mu^{+}","l");
  l_det31->AddEntry(hist_xh_det31_MuMinus,"#mu^{-}","l");
  l_det31->SetFillColor(kWhite);
  l_det31->SetLineColor(kBlack);
  l_det31->SetTextFont(43);
  l_det31->SetTextSize(20);
  l_det31->Draw();
  c_det31->Update();

  c_det31->SaveAs((plotPath + "/" + c_det31->GetName() + ".png"));    
  c_det31->SaveAs((plotPath + "/" + c_det31->GetName() + ".pdf")); 

  // xh in det32
  TCanvas* c_det32 = new TCanvas("c_det32","c_det32");
  c_det32->cd();
  hist_xh_det32_MuPlus->SetTitle("xh in det32");
  hist_xh_det32_MuPlus->GetXaxis()->SetTitle("mm");
  hist_xh_det32_MuPlus->SetLineColor(kRed);
  hist_xh_det32_MuPlus->SetMaximum(1.05 * max(hist_xh_det32_MuMinus->GetMaximum(),hist_xh_det32_MuPlus->GetMaximum()));
  hist_xh_det32_MuPlus->Draw("hist");
  hist_xh_det32_MuMinus->SetTitle("");
  hist_xh_det32_MuMinus->SetLineColor(kBlue);
  hist_xh_det32_MuMinus->Draw("histsame");
  gStyle->SetOptStat(0);

  TLegend* l_det32 = new TLegend(0.84,0.76,0.94,0.87);
  l_det32->AddEntry(hist_xh_det32_MuPlus, "#mu^{+}","l");
  l_det32->AddEntry(hist_xh_det32_MuMinus,"#mu^{-}","l");
  l_det32->SetFillColor(kWhite);
  l_det32->SetLineColor(kBlack);
  l_det32->SetTextFont(43);
  l_det32->SetTextSize(20);
  l_det32->Draw();
  c_det32->Update();

  c_det32->SaveAs((plotPath + "/" + c_det32->GetName() + ".png"));    
  c_det32->SaveAs((plotPath + "/" + c_det32->GetName() + ".pdf")); 

  // xh in det33
  TCanvas* c_det33 = new TCanvas("c_det33","c_det33");
  c_det33->cd();
  hist_xh_det33_MuPlus->SetTitle("xh in det33");
  hist_xh_det33_MuPlus->GetXaxis()->SetTitle("mm");
  hist_xh_det33_MuPlus->SetLineColor(kRed);
  hist_xh_det33_MuPlus->SetMaximum(1.05 * max(hist_xh_det33_MuMinus->GetMaximum(),hist_xh_det33_MuPlus->GetMaximum()));
  hist_xh_det33_MuPlus->Draw("hist");
  hist_xh_det33_MuMinus->SetTitle("");
  hist_xh_det33_MuMinus->SetLineColor(kBlue);
  hist_xh_det33_MuMinus->Draw("histsame");
  gStyle->SetOptStat(0);

  TLegend* l_det33 = new TLegend(0.84,0.76,0.94,0.87);
  l_det33->AddEntry(hist_xh_det33_MuPlus, "#mu^{+}","l");
  l_det33->AddEntry(hist_xh_det33_MuMinus,"#mu^{-}","l");
  l_det33->SetFillColor(kWhite);
  l_det33->SetLineColor(kBlack);
  l_det33->SetTextFont(43);
  l_det33->SetTextSize(20);
  l_det33->Draw();
  c_det33->Update();

  c_det33->SaveAs((plotPath + "/" + c_det33->GetName() + ".png"));    
  c_det33->SaveAs((plotPath + "/" + c_det33->GetName() + ".pdf")); 

  // xh in det34
  TCanvas* c_det34 = new TCanvas("c_det34","c_det34");
  c_det34->cd();
  hist_xh_det34_MuPlus->SetTitle("xh in det34");
  hist_xh_det34_MuPlus->GetXaxis()->SetTitle("mm");
  hist_xh_det34_MuPlus->SetLineColor(kRed);
  hist_xh_det34_MuPlus->SetMaximum(1.05 * max(hist_xh_det34_MuMinus->GetMaximum(),hist_xh_det34_MuPlus->GetMaximum()));
  hist_xh_det34_MuPlus->Draw("hist");
  hist_xh_det34_MuMinus->SetTitle("");
  hist_xh_det34_MuMinus->SetLineColor(kBlue);
  hist_xh_det34_MuMinus->Draw("histsame");
  gStyle->SetOptStat(0);

  TLegend* l_det34 = new TLegend(0.84,0.76,0.94,0.87);
  l_det34->AddEntry(hist_xh_det34_MuPlus, "#mu^{+}","l");
  l_det34->AddEntry(hist_xh_det34_MuMinus,"#mu^{-}","l");
  l_det34->SetFillColor(kWhite);
  l_det34->SetLineColor(kBlack);
  l_det34->SetTextFont(43);
  l_det34->SetTextSize(20);
  l_det34->Draw();
  c_det34->Update();

  c_det34->SaveAs((plotPath + "/" + c_det34->GetName() + ".png"));    
  c_det34->SaveAs((plotPath + "/" + c_det34->GetName() + ".pdf")); 

  // xh in det35
  TCanvas* c_det35 = new TCanvas("c_det35","c_det35");
  c_det35->cd();
  hist_xh_det35_MuPlus->SetTitle("xh in det35");
  hist_xh_det35_MuPlus->GetXaxis()->SetTitle("mm");
  hist_xh_det35_MuPlus->SetLineColor(kRed);
  hist_xh_det35_MuPlus->SetMaximum(1.05 * max(hist_xh_det35_MuMinus->GetMaximum(),hist_xh_det35_MuPlus->GetMaximum()));
  hist_xh_det35_MuPlus->Draw("hist");
  hist_xh_det35_MuMinus->SetTitle("");
  hist_xh_det35_MuMinus->SetLineColor(kBlue);
  hist_xh_det35_MuMinus->Draw("histsame");
  gStyle->SetOptStat(0);

  TLegend* l_det35 = new TLegend(0.84,0.76,0.94,0.87);
  l_det35->AddEntry(hist_xh_det35_MuPlus, "#mu^{+}","l");
  l_det35->AddEntry(hist_xh_det35_MuMinus,"#mu^{-}","l");
  l_det35->SetFillColor(kWhite);
  l_det35->SetLineColor(kBlack);
  l_det35->SetTextFont(43);
  l_det35->SetTextSize(20);
  l_det35->Draw();
  c_det35->Update();

  c_det35->SaveAs((plotPath + "/" + c_det35->GetName() + ".png"));    
  c_det35->SaveAs((plotPath + "/" + c_det35->GetName() + ".pdf")); 

  // xh in det36
  TCanvas* c_det36 = new TCanvas("c_det36","c_det36");
  c_det36->cd();
  hist_xh_det36_MuPlus->SetTitle("xh in det36");
  hist_xh_det36_MuPlus->GetXaxis()->SetTitle("mm");
  hist_xh_det36_MuPlus->SetLineColor(kRed);
  hist_xh_det36_MuPlus->SetMaximum(1.05 * max(hist_xh_det36_MuMinus->GetMaximum(),hist_xh_det36_MuPlus->GetMaximum()));
  hist_xh_det36_MuPlus->Draw("hist");
  hist_xh_det36_MuMinus->SetTitle("");
  hist_xh_det36_MuMinus->SetLineColor(kBlue);
  hist_xh_det36_MuMinus->Draw("histsame");
  gStyle->SetOptStat(0);

  TLegend* l_det36 = new TLegend(0.84,0.76,0.94,0.87);
  l_det36->AddEntry(hist_xh_det36_MuPlus, "#mu^{+}","l");
  l_det36->AddEntry(hist_xh_det36_MuMinus,"#mu^{-}","l");
  l_det36->SetFillColor(kWhite);
  l_det36->SetLineColor(kBlack);
  l_det36->SetTextFont(43);
  l_det36->SetTextSize(20);
  l_det36->Draw();
  c_det36->Update();

  c_det36->SaveAs((plotPath + "/" + c_det36->GetName() + ".png"));    
  c_det36->SaveAs((plotPath + "/" + c_det36->GetName() + ".pdf")); 

  // xh in det37
  TCanvas* c_det37 = new TCanvas("c_det37","c_det37");
  c_det37->cd();
  hist_xh_det37_MuPlus->SetTitle("xh in det37");
  hist_xh_det37_MuPlus->GetXaxis()->SetTitle("mm");
  hist_xh_det37_MuPlus->SetLineColor(kRed);
  hist_xh_det37_MuPlus->SetMaximum(1.05 * max(hist_xh_det37_MuMinus->GetMaximum(),hist_xh_det37_MuPlus->GetMaximum()));
  hist_xh_det37_MuPlus->Draw("hist");
  hist_xh_det37_MuMinus->SetTitle("");
  hist_xh_det37_MuMinus->SetLineColor(kBlue);
  hist_xh_det37_MuMinus->Draw("histsame");
  gStyle->SetOptStat(0);

  TLegend* l_det37 = new TLegend(0.84,0.76,0.94,0.87);
  l_det37->AddEntry(hist_xh_det37_MuPlus, "#mu^{+}","l");
  l_det37->AddEntry(hist_xh_det37_MuMinus,"#mu^{-}","l");
  l_det37->SetFillColor(kWhite);
  l_det37->SetLineColor(kBlack);
  l_det37->SetTextFont(43);
  l_det37->SetTextSize(20);
  l_det37->Draw();
  c_det37->Update();

  c_det37->SaveAs((plotPath + "/" + c_det37->GetName() + ".png"));    
  c_det37->SaveAs((plotPath + "/" + c_det37->GetName() + ".pdf"));     

  cout<<" Plots done! =) "<<endl; 

}
