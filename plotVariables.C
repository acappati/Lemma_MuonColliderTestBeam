// ************************************
// 
// usage: root -l -b -q plotVariables.C
//
// ************************************

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
#include "TF1.h"
#include "TGraph.h"
#include "TPaveText.h"


using namespace std;


//Extrapolate track function: track points are refitted
Double_t extrapolate_track_x(Double_t z0 ,Double_t x_pos_mum[12], Double_t x_pos_mum_err[12], Double_t z_x_pos_mum[12], Double_t& x_ext, Double_t& x_ext_err, Double_t& dx_on_dz_ext, Double_t& dx_on_dz_ext_err){

   //Magnetic field (box)
   Double_t zM= 17193-846;
   Double_t z1= zM-1000.;
   Double_t z2= zM+1000.;
   Double_t B = 1.7476;

   Double_t chi2 = 99999;

   //3 cases: z > z_magnet_out=z2, z < z_magnet_entry=z1, z1<z<z2 

   if (z0>z2) {

      Int_t npoints=0;

      for (Int_t k=0; k<12; k++) {
 
         if (z_x_pos_mum[k]>z2) npoints++;

         }

      if (npoints < 2) return chi2;
   
      Double_t z_x_pos_mum_err[12];
      for (Int_t k=0; k<12; k++) {z_x_pos_mum_err[k]=0;}

      //Change TGraph with TGraphErrors to use errors
      //TGraphErrors* graph = new TGraphErrors(npoints, z_x_pos_mum, x_pos_mum, z_x_pos_mum_err, x_pos_mum_err);

      TGraph* graph = new TGraph(npoints, z_x_pos_mum, x_pos_mum);

      //Linear fit using points after magnet

      TF1* line = new TF1("line", "[0]+[1]*(x-[2])");

      Double_t theta_min = 0.03;
      Double_t theta_max = 0.07;

      line->FixParameter(2,z0);


      graph->Fit(line,"Q");


      chi2 =  line->GetChisquare()/line->GetNDF();

      x_ext = line->GetParameter(0);
      x_ext_err = line->GetParError(0);

      dx_on_dz_ext = line->GetParameter(1);
      dx_on_dz_ext_err = line->GetParError(1);

      delete graph;
      delete line;

   } else if (z0 < z1)
   {

      Int_t npoints=0;

      for (Int_t k=0; k<12; k++) {
 
         if ((z_x_pos_mum[k]<-0.001 || z_x_pos_mum[k]>0.001)  && z_x_pos_mum[k] < 22700) npoints++;

         }    

      if (npoints < 2) return chi2;

      Double_t z_x_pos_mum_err[12];
      for (Int_t k=0; k<12; k++) {z_x_pos_mum_err[k]=0;}

      //Change TGraph with TGraphErrors to use errors
      //TGraphErrors* graph = new TGraphErrors(npoints, z_x_pos_mum, x_pos_mum, z_x_pos_mum_err, x_pos_mum_err);

      TGraph* graph = new TGraph(npoints, z_x_pos_mum, x_pos_mum);

      //Fit to the complete track: line + parabula + line, 3 free parameters

      TF1* trajectory = new TF1("trajectory", "(x<[3])*([0]+[1]*(x-[3]))+(x>[4])*([0]-[2]*([4]-[3])*([4]-[3])+([1]+2*[2]*([4]-[3]))*(x-[3]))+(x>[3])*(x<[4])*([0]+[1]*(x-[3])+[2]*(x-[3])*(x-[3]))");

      trajectory->FixParameter(3,z1);
      trajectory->FixParameter(4,z2);
  
      graph->Fit(trajectory,"Q");

      Double_t R = 1./(2.*trajectory->GetParameter(2));
      Double_t p = -B/(1e+9/TMath::C()) * R;
  
      chi2 = trajectory->GetChisquare()/trajectory->GetNDF();

      x_ext = trajectory->Eval(z0);
      x_ext_err = trajectory->GetParError(0);

      dx_on_dz_ext = trajectory->GetParameter(1);
      dx_on_dz_ext_err = trajectory->GetParError(1);


      delete graph;
      delete trajectory;

   } else 
   {
   cout<< "Extrapolation not defined" << endl;
   }

   return chi2;

}


// doTheHistos function: read root file and do histos 
void doTheHistos(TString inputFileName, TString label){

  
  Double_t chi2m;
  Double_t x_pos_mum[12];
  Double_t x_pos_mum_err[12];
  Double_t z_x_pos_mum[12];
  Double_t x_pos_DT_mum[8];
  Double_t z_pos_DT_mum[8];
  Double_t p_mum;
  Double_t p_mup;
  Double_t chi2p;
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
  Double_t Calo_EnDep[25];
  Int_t    event_type;

  TFile* inputFile = TFile::Open(inputFileName);
  TTree* inputTree = (TTree*)inputFile->Get("lemma");

  inputTree->SetBranchAddress("chi2m",	        &chi2m);	     
  inputTree->SetBranchAddress("x_pos_mum",      &x_pos_mum[0]); 
  inputTree->SetBranchAddress("x_pos_mum_err",  &x_pos_mum_err[0]);
  inputTree->SetBranchAddress("z_x_pos_mum",    &z_x_pos_mum[0]);
  inputTree->SetBranchAddress("x_pos_DT_mum",   &x_pos_DT_mum[0]);
  inputTree->SetBranchAddress("z_pos_DT_mum",   &z_pos_DT_mum[0]);
  inputTree->SetBranchAddress("p_mum",          &p_mum);	     
  inputTree->SetBranchAddress("p_mup",          &p_mup);	     
  inputTree->SetBranchAddress("chi2p",          &chi2p);	       
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

  TH1F* hist_xh_det62_MuPlus  = new TH1F("hist_xh_det62_MuPlus", "hist_xh_det62_MuPlus", 50,-500.,500.);
  TH1F* hist_xh_det61_MuMinus = new TH1F("hist_xh_det61_MuMinus","hist_xh_det61_MuMinus",50,-500.,500.);

  TH1F* hist_xext_MuMinus = new TH1F("hist_xext_MuMinus","hist_xext_MuMinus",20,-50.,50.);
  TH1F* hist_xext_MuPlus  = new TH1F("hist_xext_MuPlus", "hist_xext_MuPlus", 20,-50.,50.);

  // loop over tree entries 
  Long64_t entries = inputTree->GetEntries();
  for(Long64_t z=0; z<entries; ++z){

    inputTree->GetEntry(z);

    if( p_mup > 0. && p_mum > 0. ) {

      hist_pMuPlus->Fill(p_mup);      //momentum for mu plus
      hist_chi2MuPlus->Fill(chi2p);   //chi2 for mu plus tracks

      hist_pMuMinus->Fill(p_mum);     //momentum for mu minus
      hist_chi2MuMinus->Fill(chi2m);  //chi2 for mu minus tracks

      hist_pTot->Fill(p_mum + p_mup); //total momentum 

      // histos for DTs
      hist_xh_det62_MuPlus->Fill(x_pos_DT_mup[0]);
      hist_xh_det61_MuMinus->Fill(x_pos_DT_mum[0]);


      // use extrapolate_track_x function 
      Double_t z0 = 10.*(569.5-84.6); // subdet 30
      Double_t x_ext = -9999;
      Double_t x_ext_err = -9999;
      Double_t dx_on_dz_ext = -9999;
      Double_t dx_on_dz_ext_err = -9999;
      Double_t chi2 = 99999;
      //inputs are: z0, x_pos_mum, x_pos_mum_err, z_pos_mum
      //outputs are saved in:
      //x_ext (extrapolated x position at z=z0), 
      //x_ext_err, dx_on_dz_ext (extrapolated dx/dz at z=z0), 
      //dx_on_dz_ext_err
      chi2 = extrapolate_track_x(z0, x_pos_mum, x_pos_mum_err, z_x_pos_mum, x_ext, x_ext_err, dx_on_dz_ext, dx_on_dz_ext_err);
      if( chi2 < 9999. ){
        // cout << x_ext << " " << x_ext_err << endl;
        hist_xext_MuMinus->Fill(x_ext);
      }
      chi2 = extrapolate_track_x(z0, x_pos_mup, x_pos_mup_err, z_x_pos_mup, x_ext, x_ext_err, dx_on_dz_ext, dx_on_dz_ext_err);
      if( chi2 < 9999. ){
        // cout << x_ext << " " << x_ext_err << endl;
        hist_xext_MuPlus->Fill(x_ext);
      }


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

      } // end loop over i

    } // end if (p_mup > 0. && p_mum > 0.)

  }//end over tree entries 



  // save histos into a root file 
  TString outFileName = "plotVariables_" + label + ".root";
  TFile* fOutHistos = new TFile(outFileName,"recreate");
  fOutHistos->cd();

  hist_pMuPlus->Write(hist_pMuPlus->GetName());
  hist_pMuMinus->Write(hist_pMuMinus->GetName());
  hist_pTot->Write(hist_pTot->GetName());      
  hist_chi2MuPlus->Write(hist_chi2MuPlus->GetName());
  hist_chi2MuMinus->Write(hist_chi2MuMinus->GetName());
                                
  hist_xh_det10_MuPlus->Write(hist_xh_det10_MuPlus->GetName());
  hist_xh_det20_MuPlus->Write(hist_xh_det20_MuPlus->GetName());
  hist_xh_det30_MuPlus->Write(hist_xh_det30_MuPlus->GetName());
  hist_xh_det31_MuPlus->Write(hist_xh_det31_MuPlus->GetName());
  hist_xh_det32_MuPlus->Write(hist_xh_det32_MuPlus->GetName());
  hist_xh_det33_MuPlus->Write(hist_xh_det33_MuPlus->GetName());
  hist_xh_det34_MuPlus->Write(hist_xh_det34_MuPlus->GetName());
  hist_xh_det35_MuPlus->Write(hist_xh_det35_MuPlus->GetName());
  hist_xh_det36_MuPlus->Write(hist_xh_det36_MuPlus->GetName());
  hist_xh_det37_MuPlus->Write(hist_xh_det37_MuPlus->GetName());
                              
  hist_xh_det10_MuMinus->Write(hist_xh_det10_MuMinus->GetName()); 
  hist_xh_det20_MuMinus->Write(hist_xh_det20_MuMinus->GetName()); 
  hist_xh_det30_MuMinus->Write(hist_xh_det30_MuMinus->GetName()); 
  hist_xh_det31_MuMinus->Write(hist_xh_det31_MuMinus->GetName()); 
  hist_xh_det32_MuMinus->Write(hist_xh_det32_MuMinus->GetName()); 
  hist_xh_det33_MuMinus->Write(hist_xh_det33_MuMinus->GetName()); 
  hist_xh_det34_MuMinus->Write(hist_xh_det34_MuMinus->GetName()); 
  hist_xh_det35_MuMinus->Write(hist_xh_det35_MuMinus->GetName()); 
  hist_xh_det36_MuMinus->Write(hist_xh_det36_MuMinus->GetName()); 
  hist_xh_det37_MuMinus->Write(hist_xh_det37_MuMinus->GetName()); 
                              
  hist_xh_det62_MuPlus->Write(hist_xh_det62_MuPlus->GetName());  
  hist_xh_det61_MuMinus->Write(hist_xh_det61_MuMinus->GetName()); 
                              
  hist_xext_MuMinus->Write(hist_xext_MuMinus->GetName());
  hist_xext_MuPlus->Write(hist_xext_MuPlus->GetName()); 

  fOutHistos->Close();
  delete fOutHistos;

  cout<<" root file filled and created!"<<endl;

}



// data MC comparison function
void dataMCComparison(TString plotDataMCOutputPath){
  
  // read data file 
  TFile *inFile_Data = TFile::Open("plotVariables_DATA.root");

  TH1F* hist_pMuPlus_Data     = (TH1F*)inFile_Data->Get("hist_pMuPlus");
  TH1F* hist_pMuMinus_Data    = (TH1F*)inFile_Data->Get("hist_pMuMinus");
  TH1F* hist_pTot_Data        = (TH1F*)inFile_Data->Get("hist_pTot");      
  TH1F* hist_chi2MuPlus_Data  = (TH1F*)inFile_Data->Get("hist_chi2MuPlus");
  TH1F* hist_chi2MuMinus_Data = (TH1F*)inFile_Data->Get("hist_chi2MuMinus");
                                
  TH1F* hist_xh_det10_MuPlus_Data = (TH1F*)inFile_Data->Get("hist_xh_det10_MuPlus");
  TH1F* hist_xh_det20_MuPlus_Data = (TH1F*)inFile_Data->Get("hist_xh_det20_MuPlus");
  TH1F* hist_xh_det30_MuPlus_Data = (TH1F*)inFile_Data->Get("hist_xh_det30_MuPlus");
  TH1F* hist_xh_det31_MuPlus_Data = (TH1F*)inFile_Data->Get("hist_xh_det31_MuPlus");
  TH1F* hist_xh_det32_MuPlus_Data = (TH1F*)inFile_Data->Get("hist_xh_det32_MuPlus");
  TH1F* hist_xh_det33_MuPlus_Data = (TH1F*)inFile_Data->Get("hist_xh_det33_MuPlus");
  TH1F* hist_xh_det34_MuPlus_Data = (TH1F*)inFile_Data->Get("hist_xh_det34_MuPlus");
  TH1F* hist_xh_det35_MuPlus_Data = (TH1F*)inFile_Data->Get("hist_xh_det35_MuPlus");
  TH1F* hist_xh_det36_MuPlus_Data = (TH1F*)inFile_Data->Get("hist_xh_det36_MuPlus");
  TH1F* hist_xh_det37_MuPlus_Data = (TH1F*)inFile_Data->Get("hist_xh_det37_MuPlus");
                              
  TH1F* hist_xh_det10_MuMinus_Data = (TH1F*)inFile_Data->Get("hist_xh_det10_MuMinus"); 
  TH1F* hist_xh_det20_MuMinus_Data = (TH1F*)inFile_Data->Get("hist_xh_det20_MuMinus"); 
  TH1F* hist_xh_det30_MuMinus_Data = (TH1F*)inFile_Data->Get("hist_xh_det30_MuMinus"); 
  TH1F* hist_xh_det31_MuMinus_Data = (TH1F*)inFile_Data->Get("hist_xh_det31_MuMinus"); 
  TH1F* hist_xh_det32_MuMinus_Data = (TH1F*)inFile_Data->Get("hist_xh_det32_MuMinus"); 
  TH1F* hist_xh_det33_MuMinus_Data = (TH1F*)inFile_Data->Get("hist_xh_det33_MuMinus"); 
  TH1F* hist_xh_det34_MuMinus_Data = (TH1F*)inFile_Data->Get("hist_xh_det34_MuMinus"); 
  TH1F* hist_xh_det35_MuMinus_Data = (TH1F*)inFile_Data->Get("hist_xh_det35_MuMinus"); 
  TH1F* hist_xh_det36_MuMinus_Data = (TH1F*)inFile_Data->Get("hist_xh_det36_MuMinus"); 
  TH1F* hist_xh_det37_MuMinus_Data = (TH1F*)inFile_Data->Get("hist_xh_det37_MuMinus"); 
                              
  TH1F* hist_xh_det62_MuPlus_Data  = (TH1F*)inFile_Data->Get("hist_xh_det62_MuPlus");  
  TH1F* hist_xh_det61_MuMinus_Data = (TH1F*)inFile_Data->Get("hist_xh_det61_MuMinus"); 
                              
  TH1F* hist_xext_MuMinus_Data = (TH1F*)inFile_Data->Get("hist_xext_MuMinus");
  TH1F* hist_xext_MuPlus_Data  = (TH1F*)inFile_Data->Get("hist_xext_MuPlus"); 


  // read MC file 
  TFile *inFile_MC = TFile::Open("plotVariables_MC.root");

  TH1F* hist_pMuPlus_MC     = (TH1F*)inFile_MC->Get("hist_pMuPlus");
  TH1F* hist_pMuMinus_MC    = (TH1F*)inFile_MC->Get("hist_pMuMinus");
  TH1F* hist_pTot_MC        = (TH1F*)inFile_MC->Get("hist_pTot");      
  TH1F* hist_chi2MuPlus_MC  = (TH1F*)inFile_MC->Get("hist_chi2MuPlus");
  TH1F* hist_chi2MuMinus_MC = (TH1F*)inFile_MC->Get("hist_chi2MuMinus");
                                
  TH1F* hist_xh_det10_MuPlus_MC = (TH1F*)inFile_MC->Get("hist_xh_det10_MuPlus");
  TH1F* hist_xh_det20_MuPlus_MC = (TH1F*)inFile_MC->Get("hist_xh_det20_MuPlus");
  TH1F* hist_xh_det30_MuPlus_MC = (TH1F*)inFile_MC->Get("hist_xh_det30_MuPlus");
  TH1F* hist_xh_det31_MuPlus_MC = (TH1F*)inFile_MC->Get("hist_xh_det31_MuPlus");
  TH1F* hist_xh_det32_MuPlus_MC = (TH1F*)inFile_MC->Get("hist_xh_det32_MuPlus");
  TH1F* hist_xh_det33_MuPlus_MC = (TH1F*)inFile_MC->Get("hist_xh_det33_MuPlus");
  TH1F* hist_xh_det34_MuPlus_MC = (TH1F*)inFile_MC->Get("hist_xh_det34_MuPlus");
  TH1F* hist_xh_det35_MuPlus_MC = (TH1F*)inFile_MC->Get("hist_xh_det35_MuPlus");
  TH1F* hist_xh_det36_MuPlus_MC = (TH1F*)inFile_MC->Get("hist_xh_det36_MuPlus");
  TH1F* hist_xh_det37_MuPlus_MC = (TH1F*)inFile_MC->Get("hist_xh_det37_MuPlus");
                              
  TH1F* hist_xh_det10_MuMinus_MC = (TH1F*)inFile_MC->Get("hist_xh_det10_MuMinus"); 
  TH1F* hist_xh_det20_MuMinus_MC = (TH1F*)inFile_MC->Get("hist_xh_det20_MuMinus"); 
  TH1F* hist_xh_det30_MuMinus_MC = (TH1F*)inFile_MC->Get("hist_xh_det30_MuMinus"); 
  TH1F* hist_xh_det31_MuMinus_MC = (TH1F*)inFile_MC->Get("hist_xh_det31_MuMinus"); 
  TH1F* hist_xh_det32_MuMinus_MC = (TH1F*)inFile_MC->Get("hist_xh_det32_MuMinus"); 
  TH1F* hist_xh_det33_MuMinus_MC = (TH1F*)inFile_MC->Get("hist_xh_det33_MuMinus"); 
  TH1F* hist_xh_det34_MuMinus_MC = (TH1F*)inFile_MC->Get("hist_xh_det34_MuMinus"); 
  TH1F* hist_xh_det35_MuMinus_MC = (TH1F*)inFile_MC->Get("hist_xh_det35_MuMinus"); 
  TH1F* hist_xh_det36_MuMinus_MC = (TH1F*)inFile_MC->Get("hist_xh_det36_MuMinus"); 
  TH1F* hist_xh_det37_MuMinus_MC = (TH1F*)inFile_MC->Get("hist_xh_det37_MuMinus"); 
                              
  TH1F* hist_xh_det62_MuPlus_MC  = (TH1F*)inFile_MC->Get("hist_xh_det62_MuPlus");  
  TH1F* hist_xh_det61_MuMinus_MC = (TH1F*)inFile_MC->Get("hist_xh_det61_MuMinus"); 
                              
  TH1F* hist_xext_MuMinus_MC = (TH1F*)inFile_MC->Get("hist_xext_MuMinus");
  TH1F* hist_xext_MuPlus_MC  = (TH1F*)inFile_MC->Get("hist_xext_MuPlus"); 


  gStyle->SetOptStat(0);
 
  //plot histos

  // pMuPlus plot
  TCanvas* c_pMuPlus = new TCanvas("c_pMuPlus","c_pMuPlus");
  c_pMuPlus->cd();
  hist_pMuPlus_MC->SetTitle("p #mu^{+}");
  hist_pMuPlus_MC->GetXaxis()->SetTitle("p #mu^{+}");
  hist_pMuPlus_MC->SetLineColor(kRed);
  hist_pMuPlus_MC->Draw("hist");
  hist_pMuPlus_Data->SetMarkerStyle(20);
  hist_pMuPlus_Data->SetMarkerColor(kBlack);
  hist_pMuPlus_Data->SetLineColor(kBlack);
  hist_pMuPlus_Data->Draw("samepe");
  TPaveText* pv_pMuPlus = new TPaveText(0.74,0.73,0.98,0.95,"brNDC");
  pv_pMuPlus->AddText("MC: "); ((TText*)pv_pMuPlus->GetListOfLines()->Last())->SetTextColor(kRed);
  pv_pMuPlus->AddText(Form("entries: %.2f",hist_pMuPlus_MC->GetEntries())); ((TText*)pv_pMuPlus->GetListOfLines()->Last())->SetTextColor(kRed);
  pv_pMuPlus->AddText(Form("mean:    %.2f",hist_pMuPlus_MC->GetMean()));    ((TText*)pv_pMuPlus->GetListOfLines()->Last())->SetTextColor(kRed);
  pv_pMuPlus->AddText("Data: "); ((TText*)pv_pMuPlus->GetListOfLines()->Last())->SetTextColor(kBlack);
  pv_pMuPlus->AddText(Form("entries: %.2f",hist_pMuPlus_Data->GetEntries())); ((TText*)pv_pMuPlus->GetListOfLines()->Last())->SetTextColor(kBlack);
  pv_pMuPlus->AddText(Form("mean:    %.2f",hist_pMuPlus_Data->GetMean()));    ((TText*)pv_pMuPlus->GetListOfLines()->Last())->SetTextColor(kBlack);
  pv_pMuPlus->SetTextSize(0.03);
  pv_pMuPlus->SetFillColor(kWhite);
  pv_pMuPlus->SetBorderSize(1);
  pv_pMuPlus->SetTextFont(40);
  pv_pMuPlus->SetTextSize(0.037);
  pv_pMuPlus->SetTextFont(42);
  pv_pMuPlus->SetTextAlign(12);
  pv_pMuPlus->Draw(); 
  c_pMuPlus->Update();
  c_pMuPlus->SaveAs((plotDataMCOutputPath + "/" + c_pMuPlus->GetName() + ".png"));

  // pMuMinus plot
  TCanvas* c_pMuMinus = new TCanvas("c_pMuMinus","c_pMuMinus");
  c_pMuMinus->cd();
  hist_pMuMinus_MC->SetTitle("p #mu^{+}");
  hist_pMuMinus_MC->GetXaxis()->SetTitle("p #mu^{+}");
  hist_pMuMinus_MC->SetLineColor(kBlue);
  hist_pMuMinus_MC->Draw("hist");
  hist_pMuMinus_Data->SetMarkerStyle(20);
  hist_pMuMinus_Data->SetMarkerColor(kBlack);
  hist_pMuMinus_Data->SetLineColor(kBlack);
  hist_pMuMinus_Data->Draw("samepe");
  TPaveText* pv_pMuMinus = new TPaveText(0.74,0.73,0.98,0.95,"brNDC");
  pv_pMuMinus->AddText("MC: "); ((TText*)pv_pMuMinus->GetListOfLines()->Last())->SetTextColor(kBlue);
  pv_pMuMinus->AddText(Form("entries: %.2f",hist_pMuMinus_MC->GetEntries())); ((TText*)pv_pMuMinus->GetListOfLines()->Last())->SetTextColor(kBlue);
  pv_pMuMinus->AddText(Form("mean:    %.2f",hist_pMuMinus_MC->GetMean()));    ((TText*)pv_pMuMinus->GetListOfLines()->Last())->SetTextColor(kBlue);
  pv_pMuMinus->AddText("Data: "); ((TText*)pv_pMuMinus->GetListOfLines()->Last())->SetTextColor(kBlack);
  pv_pMuMinus->AddText(Form("entries: %.2f",hist_pMuMinus_Data->GetEntries())); ((TText*)pv_pMuMinus->GetListOfLines()->Last())->SetTextColor(kBlack);
  pv_pMuMinus->AddText(Form("mean:    %.2f",hist_pMuMinus_Data->GetMean()));    ((TText*)pv_pMuMinus->GetListOfLines()->Last())->SetTextColor(kBlack);
  pv_pMuMinus->SetTextSize(0.03);
  pv_pMuMinus->SetFillColor(kWhite);
  pv_pMuMinus->SetBorderSize(1);
  pv_pMuMinus->SetTextFont(40);
  pv_pMuMinus->SetTextSize(0.037);
  pv_pMuMinus->SetTextFont(42);
  pv_pMuMinus->SetTextAlign(12);
  pv_pMuMinus->Draw(); 
  c_pMuMinus->Update();
  c_pMuMinus->SaveAs((plotDataMCOutputPath + "/" + c_pMuMinus->GetName() + ".png"));
 
  
  

  // TCanvas* c_pTot = new TCanvas("c_pTot","c_pTot");
  // c_pTot->cd();
  // hist_pTot->SetTitle("p #mu^{+} + p #mu^{-}");
  // hist_pTot->GetXaxis()->SetTitle("p #mu^{+} + p #mu^{-}");
  // hist_pTot->Draw("hist");
  // c_pTot->SaveAs((plotDataMCOutputPath + "/" + hist_pTot->GetName() + ".png"));
 

  // TCanvas* c_chi2MuPlus = new TCanvas("c_chi2MuPlus","c_chi2MuPlus");
  // c_chi2MuPlus->cd();
  // hist_chi2MuPlus->SetTitle("#Chi^{2} #mu^{+}");
  // hist_chi2MuPlus->GetXaxis()->SetTitle("#Chi^{2} #mu^{+}");
  // hist_chi2MuPlus->Draw("hist");
  // c_chi2MuPlus->SaveAs((plotDataMCOutputPath + "/" + hist_chi2MuPlus->GetName() + ".png"));
 

  // TCanvas* c_chi2MuMinus = new TCanvas("c_chi2MuMinus","c_chi2MuMinus");
  // c_chi2MuMinus->cd();
  // hist_chi2MuMinus->SetTitle("#Chi^{2} #mu^{-}");
  // hist_chi2MuMinus->GetXaxis()->SetTitle("#Chi^{2} #mu^{-}");
  // hist_chi2MuMinus->Draw("hist");
  // c_chi2MuMinus->SaveAs((plotDataMCOutputPath + "/" + hist_chi2MuMinus->GetName() + ".png"));
 

  


  // // xh in det30
  // TCanvas* c_det30 = new TCanvas("c_det30","c_det30");
  // c_det30->cd();
  // hist_xh_det30_MuPlus->SetTitle("xh in det30");
  // hist_xh_det30_MuPlus->GetXaxis()->SetTitle("mm");
  // hist_xh_det30_MuPlus->SetLineColor(kRed);
  // hist_xh_det30_MuPlus->SetMaximum(1.05 * max(hist_xh_det30_MuMinus->GetMaximum(),hist_xh_det30_MuPlus->GetMaximum()));
  // hist_xh_det30_MuPlus->Draw("hist");
  // hist_xh_det30_MuMinus->SetTitle("");
  // hist_xh_det30_MuMinus->SetLineColor(kBlue);
  // hist_xh_det30_MuMinus->Draw("histsame");

  // TLegend* l_det30 = new TLegend(0.84,0.76,0.94,0.87);
  // l_det30->AddEntry(hist_xh_det30_MuPlus, "#mu^{+}","l");
  // l_det30->AddEntry(hist_xh_det30_MuMinus,"#mu^{-}","l");
  // l_det30->SetFillColor(kWhite);
  // l_det30->SetLineColor(kBlack);
  // l_det30->SetTextFont(43);
  // l_det30->SetTextSize(20);
  // l_det30->Draw();
  // c_det30->Update();
  // c_det30->SaveAs((plotDataMCOutputPath + "/" + c_det30->GetName() + ".png"));    
  

  // // xh in det31
  // TCanvas* c_det31 = new TCanvas("c_det31","c_det31");
  // c_det31->cd();
  // hist_xh_det31_MuPlus->SetTitle("xh in det31");
  // hist_xh_det31_MuPlus->GetXaxis()->SetTitle("mm");
  // hist_xh_det31_MuPlus->SetLineColor(kRed);
  // hist_xh_det31_MuPlus->SetMaximum(1.05 * max(hist_xh_det31_MuMinus->GetMaximum(),hist_xh_det31_MuPlus->GetMaximum()));
  // hist_xh_det31_MuPlus->Draw("hist");
  // hist_xh_det31_MuMinus->SetTitle("");
  // hist_xh_det31_MuMinus->SetLineColor(kBlue);
  // hist_xh_det31_MuMinus->Draw("histsame");

  // TLegend* l_det31 = new TLegend(0.84,0.76,0.94,0.87);
  // l_det31->AddEntry(hist_xh_det31_MuPlus, "#mu^{+}","l");
  // l_det31->AddEntry(hist_xh_det31_MuMinus,"#mu^{-}","l");
  // l_det31->SetFillColor(kWhite);
  // l_det31->SetLineColor(kBlack);
  // l_det31->SetTextFont(43);
  // l_det31->SetTextSize(20);
  // l_det31->Draw();
  // c_det31->Update();
  // c_det31->SaveAs((plotDataMCOutputPath + "/" + c_det31->GetName() + ".png"));    
 

  // // xh in det32
  // TCanvas* c_det32 = new TCanvas("c_det32","c_det32");
  // c_det32->cd();
  // hist_xh_det32_MuPlus->SetTitle("xh in det32");
  // hist_xh_det32_MuPlus->GetXaxis()->SetTitle("mm");
  // hist_xh_det32_MuPlus->SetLineColor(kRed);
  // hist_xh_det32_MuPlus->SetMaximum(1.05 * max(hist_xh_det32_MuMinus->GetMaximum(),hist_xh_det32_MuPlus->GetMaximum()));
  // hist_xh_det32_MuPlus->Draw("hist");
  // hist_xh_det32_MuMinus->SetTitle("");
  // hist_xh_det32_MuMinus->SetLineColor(kBlue);
  // hist_xh_det32_MuMinus->Draw("histsame");

  // TLegend* l_det32 = new TLegend(0.84,0.76,0.94,0.87);
  // l_det32->AddEntry(hist_xh_det32_MuPlus, "#mu^{+}","l");
  // l_det32->AddEntry(hist_xh_det32_MuMinus,"#mu^{-}","l");
  // l_det32->SetFillColor(kWhite);
  // l_det32->SetLineColor(kBlack);
  // l_det32->SetTextFont(43);
  // l_det32->SetTextSize(20);
  // l_det32->Draw();
  // c_det32->Update();
  // c_det32->SaveAs((plotDataMCOutputPath + "/" + c_det32->GetName() + ".png"));    
 

  // // xh in det33
  // TCanvas* c_det33 = new TCanvas("c_det33","c_det33");
  // c_det33->cd();
  // hist_xh_det33_MuPlus->SetTitle("xh in det33");
  // hist_xh_det33_MuPlus->GetXaxis()->SetTitle("mm");
  // hist_xh_det33_MuPlus->SetLineColor(kRed);
  // hist_xh_det33_MuPlus->SetMaximum(1.05 * max(hist_xh_det33_MuMinus->GetMaximum(),hist_xh_det33_MuPlus->GetMaximum()));
  // hist_xh_det33_MuPlus->Draw("hist");
  // hist_xh_det33_MuMinus->SetTitle("");
  // hist_xh_det33_MuMinus->SetLineColor(kBlue);
  // hist_xh_det33_MuMinus->Draw("histsame");

  // TLegend* l_det33 = new TLegend(0.84,0.76,0.94,0.87);
  // l_det33->AddEntry(hist_xh_det33_MuPlus, "#mu^{+}","l");
  // l_det33->AddEntry(hist_xh_det33_MuMinus,"#mu^{-}","l");
  // l_det33->SetFillColor(kWhite);
  // l_det33->SetLineColor(kBlack);
  // l_det33->SetTextFont(43);
  // l_det33->SetTextSize(20);
  // l_det33->Draw();
  // c_det33->Update();
  // c_det33->SaveAs((plotDataMCOutputPath + "/" + c_det33->GetName() + ".png"));    
  

  // // xh in det34
  // TCanvas* c_det34 = new TCanvas("c_det34","c_det34");
  // c_det34->cd();
  // hist_xh_det34_MuPlus->SetTitle("xh in det34");
  // hist_xh_det34_MuPlus->GetXaxis()->SetTitle("mm");
  // hist_xh_det34_MuPlus->SetLineColor(kRed);
  // hist_xh_det34_MuPlus->SetMaximum(1.05 * max(hist_xh_det34_MuMinus->GetMaximum(),hist_xh_det34_MuPlus->GetMaximum()));
  // hist_xh_det34_MuPlus->Draw("hist");
  // hist_xh_det34_MuMinus->SetTitle("");
  // hist_xh_det34_MuMinus->SetLineColor(kBlue);
  // hist_xh_det34_MuMinus->Draw("histsame");

  // TLegend* l_det34 = new TLegend(0.84,0.76,0.94,0.87);
  // l_det34->AddEntry(hist_xh_det34_MuPlus, "#mu^{+}","l");
  // l_det34->AddEntry(hist_xh_det34_MuMinus,"#mu^{-}","l");
  // l_det34->SetFillColor(kWhite);
  // l_det34->SetLineColor(kBlack);
  // l_det34->SetTextFont(43);
  // l_det34->SetTextSize(20);
  // l_det34->Draw();
  // c_det34->Update();
  // c_det34->SaveAs((plotDataMCOutputPath + "/" + c_det34->GetName() + ".png"));    
  

  // // xh in det35
  // TCanvas* c_det35 = new TCanvas("c_det35","c_det35");
  // c_det35->cd();
  // hist_xh_det35_MuPlus->SetTitle("xh in det35");
  // hist_xh_det35_MuPlus->GetXaxis()->SetTitle("mm");
  // hist_xh_det35_MuPlus->SetLineColor(kRed);
  // hist_xh_det35_MuPlus->SetMaximum(1.05 * max(hist_xh_det35_MuMinus->GetMaximum(),hist_xh_det35_MuPlus->GetMaximum()));
  // hist_xh_det35_MuPlus->Draw("hist");
  // hist_xh_det35_MuMinus->SetTitle("");
  // hist_xh_det35_MuMinus->SetLineColor(kBlue);
  // hist_xh_det35_MuMinus->Draw("histsame");

  // TLegend* l_det35 = new TLegend(0.84,0.76,0.94,0.87);
  // l_det35->AddEntry(hist_xh_det35_MuPlus, "#mu^{+}","l");
  // l_det35->AddEntry(hist_xh_det35_MuMinus,"#mu^{-}","l");
  // l_det35->SetFillColor(kWhite);
  // l_det35->SetLineColor(kBlack);
  // l_det35->SetTextFont(43);
  // l_det35->SetTextSize(20);
  // l_det35->Draw();
  // c_det35->Update();
  // c_det35->SaveAs((plotDataMCOutputPath + "/" + c_det35->GetName() + ".png"));    
   

  // // xh in det36
  // TCanvas* c_det36 = new TCanvas("c_det36","c_det36");
  // c_det36->cd();
  // hist_xh_det36_MuPlus->SetTitle("xh in det36");
  // hist_xh_det36_MuPlus->GetXaxis()->SetTitle("mm");
  // hist_xh_det36_MuPlus->SetLineColor(kRed);
  // hist_xh_det36_MuPlus->SetMaximum(1.05 * max(hist_xh_det36_MuMinus->GetMaximum(),hist_xh_det36_MuPlus->GetMaximum()));
  // hist_xh_det36_MuPlus->Draw("hist");
  // hist_xh_det36_MuMinus->SetTitle("");
  // hist_xh_det36_MuMinus->SetLineColor(kBlue);
  // hist_xh_det36_MuMinus->Draw("histsame");

  // TLegend* l_det36 = new TLegend(0.84,0.76,0.94,0.87);
  // l_det36->AddEntry(hist_xh_det36_MuPlus, "#mu^{+}","l");
  // l_det36->AddEntry(hist_xh_det36_MuMinus,"#mu^{-}","l");
  // l_det36->SetFillColor(kWhite);
  // l_det36->SetLineColor(kBlack);
  // l_det36->SetTextFont(43);
  // l_det36->SetTextSize(20);
  // l_det36->Draw();
  // c_det36->Update();
  // c_det36->SaveAs((plotDataMCOutputPath + "/" + c_det36->GetName() + ".png"));    
  

  // // xh in det37
  // TCanvas* c_det37 = new TCanvas("c_det37","c_det37");
  // c_det37->cd();
  // hist_xh_det37_MuPlus->SetTitle("xh in det37");
  // hist_xh_det37_MuPlus->GetXaxis()->SetTitle("mm");
  // hist_xh_det37_MuPlus->SetLineColor(kRed);
  // hist_xh_det37_MuPlus->SetMaximum(1.05 * max(hist_xh_det37_MuMinus->GetMaximum(),hist_xh_det37_MuPlus->GetMaximum()));
  // hist_xh_det37_MuPlus->Draw("hist");
  // hist_xh_det37_MuMinus->SetTitle("");
  // hist_xh_det37_MuMinus->SetLineColor(kBlue);
  // hist_xh_det37_MuMinus->Draw("histsame");

  // TLegend* l_det37 = new TLegend(0.84,0.76,0.94,0.87);
  // l_det37->AddEntry(hist_xh_det37_MuPlus, "#mu^{+}","l");
  // l_det37->AddEntry(hist_xh_det37_MuMinus,"#mu^{-}","l");
  // l_det37->SetFillColor(kWhite);
  // l_det37->SetLineColor(kBlack);
  // l_det37->SetTextFont(43);
  // l_det37->SetTextSize(20);
  // l_det37->Draw();
  // c_det37->Update();
  // c_det37->SaveAs((plotDataMCOutputPath + "/" + c_det37->GetName() + ".png"));    
       

  // TCanvas* c_det6x = new TCanvas("c_det6x","c_det6x");
  // c_det6x->cd();
  // hist_xh_det62_MuPlus->SetLineColor(kRed); hist_xh_det62_MuPlus->SetMarkerColor(kRed);
  // hist_xh_det61_MuMinus->SetLineColor(kBlue); hist_xh_det61_MuMinus->SetMarkerColor(kBlue);
  // hist_xh_det62_MuPlus->Draw("pe");
  // hist_xh_det61_MuMinus->Draw("samepe");
  // c_det6x->SaveAs((plotDataMCOutputPath + "/" + c_det6x->GetName() + ".png"));   


  // TCanvas* c_ExtTest = new TCanvas("c_ExtTest","c_ExtTest");
  // c_ExtTest->cd();
  // hist_xext_MuPlus->SetLineColor(kRed); hist_xext_MuPlus->SetMarkerColor(kRed);
  // hist_xext_MuMinus->SetLineColor(kBlue); hist_xext_MuMinus->SetMarkerColor(kBlue);
  // hist_xext_MuPlus->Draw("pe");
  // hist_xext_MuMinus->Draw("samepe");


  cout<<" Plots done! =) "<<endl; 


}




// main function 
void plotVariables(){

  // define input files 
  TString inputFile_Data = "/afs/cern.ch/user/a/abertoli/public/lemma/reco/reco-333to337.root";
  TString inputFile_MC   = "/afs/cern.ch/user/a/abertoli/public/lemma/reco/reco-mupmum.root"; 

  // call do the histos function
  doTheHistos(inputFile_Data, "DATA");
  doTheHistos(inputFile_MC, "MC");

  // define output path and make output directory for data/MC comparison
  TString plotDataMCOutputPath = "LemmaVariables_DataMCComparison_reco-333to337";
  gSystem->Exec(("mkdir -p "+plotDataMCOutputPath));

  // call data/MC comparison function
  dataMCComparison(plotDataMCOutputPath);


}
