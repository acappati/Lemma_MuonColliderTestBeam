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


// fill histos of det36 and det37 function
// excluding part of the bkg 
void makeHistos(int runNumber)
{
  // define input file name
  TString inputFileName;
  inputFileName.Form("/afs/cern.ch/user/a/abertoli/public/lemma/si-%d.root",runNumber);

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


  // histos for x occupancy of det 36 and 37 
  TH1F* hist_xh_det36 = new TH1F("hist_xh_det36","hist_xh_det36",60,0.,22.);
  TH1F* hist_xh_det37 = new TH1F("hist_xh_det37","hist_xh_det37",60,0.,22.);

  // loop over tree entries 
  Long64_t entries = inputTree->GetEntries();
  for(Long64_t z=0; z<entries; ++z){

    inputTree->GetEntry(z);


    bool hitDet32 = false;
    bool hitDet33 = false;
    bool hitDet34 = false;
    bool hitDet35 = false;
    bool hitDet36 = false;
    bool hitDet37 = false;


    for(int i=0; i<nhits;i++){

      if(xh[i]>0.){
        if(subdet[i]==32) {hitDet32 = true; }
        if(subdet[i]==33) {hitDet33 = true; }
        if(subdet[i]==34) {hitDet34 = true; }
        if(subdet[i]==35) {hitDet35 = true; }
        if(subdet[i]==36) {hitDet36 = true; }
        if(subdet[i]==37) {hitDet37 = true; }
        
        //cout<< "***"<<subdet[i] <<" "<<hitDet32<<" "<<hitDet32<<" "<<hitDet33<<" "<<hitDet34<<" "<<hitDet35<<" "<<hitDet36<<" "<<hitDet37<<endl;

      }
    }// end loop over number of hits 


    // fill histo for det36 only if there are hits also in the dets 32, 34 and 36 
    if(hitDet32 && hitDet34 && hitDet36){
      for(int i=0; i<nhits;i++){
        if(xh[i]>0.){
          if(subdet[i]==36) {hist_xh_det36->Fill(xh[i]);}
        }
      }
    }
    // fill histo for det37 only if there are hits also in the dets 33, 35 and 37
    if(hitDet33 && hitDet35 && hitDet37){
      for(int i=0; i<nhits;i++){
        if(xh[i]>0.){
          if(subdet[i]==37) {hist_xh_det37->Fill(xh[i]);}
        }
      }
    }

  }//end loop over entries 


  // *** save histos in a root file
  TString outputFileName;
  outputFileName.Form("plotsDet36Det37_Run%d.root",runNumber);
  TFile* fOutHistos = new TFile(outputFileName,"recreate");
  fOutHistos->cd();

  hist_xh_det36->Write(hist_xh_det36->GetName());
  hist_xh_det37->Write(hist_xh_det37->GetName());

  fOutHistos->Close();
  delete fOutHistos;

  cout<<" root file filled and created!"<<endl;


}//end makeHistos function 



// do the fit function
void doTheFit()
{
  // number of runs to be analyzed
  // current of the magnet +580A
  int runNumberPlus580A[6]  = {500332,500327,500325,500320,500319,500314};
  TString srunNumberPlus580A[6]  = {"500332","500327","500325","500320","500319","500314"};
  // current of the magnet -580A
  int runNumberMinus580A[6] = {500330,500329,500323,500322,500317,500316};
  TString srunNumberMinus580A[6] = {"500330","500329","500323","500322","500317","500316"};

 
  // vectors for the TGraph 
  float beamEnergy[6] = {16.,18.,20.,22.,24.,26.};
  float beamEnergy_err[6] = {0.1,0.1,0.1,0.1,0.1,0.1};
  string sbeamEnergy[6] = {"16GeV","18GeV","20GeV","22GeV","24GeV","26GeV"};

  float meanDet37Plus580A[6];
  float meanDet37Plus580A_err[6];
  float sigmaDet37Plus580A[6];
  float sigmaDet37Plus580A_err[6];

  float meanDet36Minus580A[6];
  float meanDet36Minus580A_err[6];
  float sigmaDet36Minus580A[6];
  float sigmaDet36Minus580A_err[6];


  // define output path and make output directory
  TString plotPath = "plotFitsProfile";
  gSystem->Exec(("mkdir -p "+plotPath));  
 
  

  // current of the magnet +580A
  for(int i=0; i<6; i++){
  
    TString inputFileName;
    inputFileName.Form("plotsDet36Det37_Run%d.root",runNumberPlus580A[i]);

    TFile* inputFile = TFile::Open(inputFileName);
    TH1F* hist_det37 = (TH1F*)inputFile->Get("hist_xh_det37");

    float minFit1 = max(hist_det37->GetXaxis()->GetBinCenter(hist_det37->GetMaximumBin())-3.,0. );
    float maxFit1 = min(hist_det37->GetXaxis()->GetBinCenter(hist_det37->GetMaximumBin())+3.,22.);
    TF1* gfit1 = new TF1("gfit1","gaus",minFit1,maxFit1);
    gfit1->SetParameter(1,hist_det37->GetMean());
    gfit1->SetParameter(2,hist_det37->GetRMS());
    gfit1->SetLineColor(kBlue);
  
    hist_det37->Fit("gfit1","R");

    float minFit2 = max(gfit1->GetParameter(1)-1.5*gfit1->GetParameter(2),0. );
    float maxFit2 = min(gfit1->GetParameter(1)+1.5*gfit1->GetParameter(2),22.);
    TF1* gfit2 = new TF1("gfit2","gaus",minFit2,maxFit2);
    gfit2->SetParameter(1,hist_det37->GetMean());
    gfit2->SetParameter(2,hist_det37->GetRMS());
    gfit2->SetLineColor(kRed);
  
    hist_det37->Fit("gfit2","R");

    // save fit param
    meanDet37Plus580A[i]      = gfit2->GetParameter(1);
    meanDet37Plus580A_err[i]  = gfit2->GetParError(1);
    sigmaDet37Plus580A[i]     = gfit2->GetParameter(2);
    sigmaDet37Plus580A_err[i] = gfit2->GetParError(2);
    

    TCanvas* c = new TCanvas(); 
    c->cd();
    hist_det37->GetXaxis()->SetTitle("[cm]");
    hist_det37->SetTitle("x view - det37 run " + srunNumberPlus580A[i]);
    hist_det37->Draw("hist");
    //gfit1->Draw("samel");
    gfit2->Draw("samel");

    gStyle->SetOptStat(0);
    //gStyle->SetOptFit(111);

    TPaveText* pvtext1 = new TPaveText(0.74,0.75,0.98,0.95,"brNDC");
    pvtext1->AddText(Form("#chi^{2}/ndf    %.2f / %.0i",gfit2->GetChisquare(),gfit2->GetNDF()));
    pvtext1->AddText(Form("Const  %.2f #pm %.2f",gfit2->GetParameter(0),gfit2->GetParError(0)));
    pvtext1->AddText(Form("Mean     %.2f #pm %.2f",gfit2->GetParameter(1),gfit2->GetParError(1)));
    pvtext1->AddText(Form("Sigma    %.2f #pm %.2f",gfit2->GetParameter(2),gfit2->GetParError(2)));
    pvtext1->SetTextSize(0.03);
    pvtext1->SetFillColor(kWhite);
    pvtext1->SetBorderSize(1);
    pvtext1->SetTextFont(40);
    pvtext1->SetTextSize(0.037);
    pvtext1->SetTextFont(42);
    pvtext1->SetTextAlign(12);
    pvtext1->Draw(); 

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

    c->SaveAs(Form(plotPath + "/" + hist_det37->GetName() + "_run%d.png",runNumberPlus580A[i]));

    delete c;
    delete gfit1;
    delete gfit2;
    delete hist_det37;
    delete inputFile;
    delete pvtext;

  }



  // current of the magnet -580A
  for(int i=0; i<6; i++){
  
    TString inputFileName;
    inputFileName.Form("plotsDet36Det37_Run%d.root",runNumberMinus580A[i]);

    TFile* inputFile = TFile::Open(inputFileName);
    TH1F* hist_det36 = (TH1F*)inputFile->Get("hist_xh_det36");

    float minFit1 = max(hist_det36->GetXaxis()->GetBinCenter(hist_det36->GetMaximumBin())-3.,0. );
    float maxFit1 = min(hist_det36->GetXaxis()->GetBinCenter(hist_det36->GetMaximumBin())+3.,22.);
    TF1* gfit1 = new TF1("gfit1","gaus",minFit1,maxFit1);
    gfit1->SetParameter(1,hist_det36->GetMean());
    gfit1->SetParameter(2,hist_det36->GetRMS());
    gfit1->SetLineColor(kBlue);
  
    hist_det36->Fit("gfit1","R");

    float minFit2 = max(gfit1->GetParameter(1)-1.5*gfit1->GetParameter(2),0. );
    float maxFit2 = min(gfit1->GetParameter(1)+1.5*gfit1->GetParameter(2),22.);
    TF1* gfit2 = new TF1("gfit2","gaus",minFit2,maxFit2);
    gfit2->SetParameter(1,hist_det36->GetMean());
    gfit2->SetParameter(2,hist_det36->GetRMS());
    gfit2->SetLineColor(kRed);
  
    hist_det36->Fit("gfit2","R");

    // save fit param
    meanDet36Minus580A[i]      = gfit2->GetParameter(1);
    meanDet36Minus580A_err[i]  = gfit2->GetParError(1);
    sigmaDet36Minus580A[i]     = gfit2->GetParameter(2);
    sigmaDet36Minus580A_err[i] = gfit2->GetParError(2);
    

    TCanvas* c = new TCanvas(); 
    c->cd();
    hist_det36->GetXaxis()->SetTitle("[cm]");
    hist_det36->SetTitle("x view - det36 run " + srunNumberMinus580A[i]);
    hist_det36->Draw("hist");
    //gfit1->Draw("samel");
    gfit2->Draw("samel");

    gStyle->SetOptStat(0);
    //gStyle->SetOptFit(111);

    TPaveText* pvtext1 = new TPaveText(0.74,0.75,0.98,0.95,"brNDC");
    pvtext1->AddText(Form("#chi^{2}/ndf    %.2f / %.0i",gfit2->GetChisquare(),gfit2->GetNDF()));
    pvtext1->AddText(Form("Const  %.2f #pm %.2f",gfit2->GetParameter(0),gfit2->GetParError(0)));
    pvtext1->AddText(Form("Mean     %.2f #pm %.2f",gfit2->GetParameter(1),gfit2->GetParError(1)));
    pvtext1->AddText(Form("Sigma    %.2f #pm %.2f",gfit2->GetParameter(2),gfit2->GetParError(2)));
    pvtext1->SetTextSize(0.03);
    pvtext1->SetFillColor(kWhite);
    pvtext1->SetBorderSize(1);
    pvtext1->SetTextFont(40);
    pvtext1->SetTextSize(0.037);
    pvtext1->SetTextFont(42);
    pvtext1->SetTextAlign(12);
    pvtext1->Draw(); 

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

    c->SaveAs(Form(plotPath + "/" + hist_det36->GetName() + "_run%d.png",runNumberMinus580A[i]));

    delete c;
    delete gfit1;
    delete gfit2;
    delete hist_det36;
    delete inputFile;
    delete pvtext;

  }


  // graphs 
  //mean vs beam en - det37
  TCanvas* c_meanDet37Plus580A = new TCanvas("c_meanDet37Plus580A","c_meanDet37Plus580A");
  TGraphErrors* g_meanDet37Plus580A = new TGraphErrors(6,beamEnergy,meanDet37Plus580A,beamEnergy_err,meanDet37Plus580A_err);
  g_meanDet37Plus580A->SetTitle("hits distrib mean on det37 - x view");
  g_meanDet37Plus580A->GetXaxis()->SetTitle("e^{+} beam energy [GeV]");
  g_meanDet37Plus580A->GetYaxis()->SetTitle("[cm]");
  g_meanDet37Plus580A->SetMarkerStyle(20);
  g_meanDet37Plus580A->SetMarkerColor(kBlue);
  g_meanDet37Plus580A->SetMarkerSize(0.9);
  c_meanDet37Plus580A->cd();
  g_meanDet37Plus580A->Draw("AP");
  c_meanDet37Plus580A->SaveAs(plotPath + "/" + c_meanDet37Plus580A->GetName() + ".png");

  //sigma vs beam en - det37
  TCanvas* c_sigmaDet37Plus580A = new TCanvas("c_sigmaDet37Plus580A","c_sigmaDet37Plus580A");
  TGraphErrors* g_sigmaDet37Plus580A = new TGraphErrors(6,beamEnergy,sigmaDet37Plus580A,beamEnergy_err,sigmaDet37Plus580A_err);
  g_sigmaDet37Plus580A->SetTitle("hits distrib sigma on det37 - x view");
  g_sigmaDet37Plus580A->GetXaxis()->SetTitle("e^{+} beam energy [GeV]");
  g_sigmaDet37Plus580A->GetYaxis()->SetTitle("[cm]");
  g_sigmaDet37Plus580A->SetMarkerStyle(20);
  g_sigmaDet37Plus580A->SetMarkerColor(kBlue);
  g_sigmaDet37Plus580A->SetMarkerSize(0.9);
  c_sigmaDet37Plus580A->cd();
  g_sigmaDet37Plus580A->Draw("AP");
  c_sigmaDet37Plus580A->SaveAs(plotPath + "/" + c_sigmaDet37Plus580A->GetName() + ".png");

  //mean vs beam en - det36
  TCanvas* c_meanDet36Minus580A = new TCanvas("c_meanDet36Minus580A","c_meanDet36Minus580A");
  TGraphErrors* g_meanDet36Minus580A = new TGraphErrors(6,beamEnergy,meanDet36Minus580A,beamEnergy_err,meanDet36Minus580A_err);
  g_meanDet36Minus580A->SetTitle("hits distrib mean on det37 - x view");
  g_meanDet36Minus580A->GetXaxis()->SetTitle("e^{+} beam energy [GeV]");
  g_meanDet36Minus580A->GetYaxis()->SetTitle("[cm]");
  g_meanDet36Minus580A->SetMarkerStyle(20);
  g_meanDet36Minus580A->SetMarkerColor(kBlue);
  g_meanDet36Minus580A->SetMarkerSize(0.9);
  c_meanDet36Minus580A->cd();
  g_meanDet36Minus580A->Draw("AP");
  c_meanDet36Minus580A->SaveAs(plotPath + "/" + c_meanDet36Minus580A->GetName() + ".png");

  //sigma vs beam en - det36
  TCanvas* c_sigmaDet36Minus580A = new TCanvas("c_sigmaDet36Minus580A","c_sigmaDet36Minus580A");
  TGraphErrors* g_sigmaDet36Minus580A = new TGraphErrors(6,beamEnergy,sigmaDet36Minus580A,beamEnergy_err,sigmaDet36Minus580A_err);
  g_sigmaDet36Minus580A->SetTitle("hits distrib sigma on det37 - x view");
  g_sigmaDet36Minus580A->GetXaxis()->SetTitle("e^{+} beam energy [GeV]");
  g_sigmaDet36Minus580A->GetYaxis()->SetTitle("[cm]");
  g_sigmaDet36Minus580A->SetMarkerStyle(20);
  g_sigmaDet36Minus580A->SetMarkerColor(kBlue);
  g_sigmaDet36Minus580A->SetMarkerSize(0.9);
  c_sigmaDet36Minus580A->cd();
  g_sigmaDet36Minus580A->Draw("AP");
  c_sigmaDet36Minus580A->SaveAs(plotPath + "/" + c_sigmaDet36Minus580A->GetName() + ".png");

  
  cout<<"Plots done!"<<endl;

}//end do the fit function 



// main function
void plotProfileTrend_noTargetRuns()
{

  bool redoHistos = true;

  if(redoHistos){
     makeHistos(500332);
     makeHistos(500327);
     makeHistos(500325);
     makeHistos(500320);
     makeHistos(500319);
     makeHistos(500314);
     makeHistos(500330);
     makeHistos(500329);
     makeHistos(500323);
     makeHistos(500322);
     makeHistos(500317);
     makeHistos(500316);
  }
 

  doTheFit();

}


