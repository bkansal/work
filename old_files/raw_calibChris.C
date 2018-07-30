
#include <vector>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TF1.h>
#include <TF2.h>
#include "TH2F.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TGraph.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "Math/SMatrix.h"
#include "Math/SVector.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLine.h"
#include <string>
#include <iostream>
#include <math.h>
#include "TDirectory.h"

//#include "CrystalBall.C"

using namespace std;

unsigned sampleRangeHigh = 500;
bool old_logic = false;

vector<float> assignvalues(vector<float> *pfcID_, vector<float> *Ecalenergy_, 
			   vector<float> *Hcalenergy_, vector<float> *dr) {

  vector<float> energies;
  float e = 0.0 , h = 0.0;
  cout<<"pfcID size ----> "<<pfcID_->size()<<endl;//<<"pfcID_->at(00) ----->"<<pfcID_->at(00)<<endl;
  for(unsigned ii = 0; ii < pfcID_->size(); ii++) {
        cout<<" pfcID_:" << pfcID_->at(ii) << endl;
    if (old_logic) {
      if (pfcID_->at(ii) == 4 && dr->at(ii) < 0.2) e += Ecalenergy_->at(ii);
      if (pfcID_->at(ii) == 5 && dr->at(ii) < 0.4) h += Hcalenergy_->at(ii);
      
    }
    else {
      if (pfcID_->at(ii) == 5 && dr->at(ii) < 0.4) {
	e += Ecalenergy_->at(ii);
	h += Hcalenergy_->at(ii);
      }
    }
  }
  energies.push_back(e);
  energies.push_back(h);
  return energies;

}


void getValuesFromTree(TTree* tree, vector<double>& ETrueEnergies, 
                       vector<double>& ecalEnergies, 
                       vector<double>& hcalEnergies, vector<double>& etas, 
                       vector<double>& phis)
{
   Float_t         true_;
   Float_t         p_;
   Float_t         ecal_;
   Float_t         hcal_;
   Float_t         eta_;
   Float_t         phi_;
   TBranch        *b_true; 
   TBranch        *b_p;   
   TBranch        *b_ecal;   
   TBranch        *b_hcal;   
   TBranch        *b_eta;    
   TBranch        *b_phi;    

   vector<float>        *pfcID_;
   vector<float>        *E_ecal_;
   vector<float>        *E_hcal_;
   vector<float>        *dr_;
   TBranch        *b_pfcID;   
   TBranch        *b_E_ecal;
   TBranch        *b_E_hcal;
   TBranch        *b_dr;

   pfcID_ = 0;
   E_ecal_ = 0;
   E_hcal_ = 0;
   dr_ = 0;

   tree->SetMakeClass(1);
   
   if(tree->GetBranchStatus("true"))
      tree->SetBranchAddress("true", &true_, &b_true);
   tree->SetBranchAddress("p", &p_, &b_p);
   tree->SetBranchAddress("ecal", &ecal_, &b_ecal);
   tree->SetBranchAddress("hcal", &hcal_, &b_hcal);
   tree->SetBranchAddress("eta", &eta_, &b_eta);
   tree->SetBranchAddress("phi", &phi_, &b_phi);
   tree->SetBranchAddress("pfcID", &pfcID_, &b_pfcID);
   tree->SetBranchAddress("Eecal", &E_ecal_, &b_E_ecal);
   tree->SetBranchAddress("Ehcal", &E_hcal_, &b_E_hcal);
   tree->SetBranchAddress("dr", &dr_, &b_dr);

   double sigmaEcalHcal=1;
   long veto = 0 ;
   //int count = 0;
   bool flag[10] = {0,0,0,0,0,0,0,0,0,0};
   for( unsigned entry = 0; entry < std::min((unsigned)50000000,(unsigned)(tree->GetEntriesFast()) ); entry++)
     {
       tree->GetEntry(entry);
       
       if (fabs(eta_) < 2.4 && p_ == 0) continue;
       //if (true_ > 2.5 ) continue;
       if(tree->GetBranchStatus("true"))
	 ETrueEnergies.push_back(true_);
       else
	 ETrueEnergies.push_back(p_);
       if(pfcID_->size() != 0) {
	 vector<float> tmp = assignvalues(pfcID_, E_ecal_, E_hcal_, dr_);
	 ecalEnergies.push_back(tmp.at(0));
	 hcalEnergies.push_back(tmp.at(1));
       }
       else {
	 ecalEnergies.push_back(ecal_);
	 hcalEnergies.push_back(hcal_);
       }
       etas.push_back(eta_);
       phis.push_back(phi_);
       
       unsigned N = tree->GetEntriesFast();
       int frac = ((double)entry/N)*100;
       switch(frac) {
       case 10 : if (!flag[0]) { cout<<"10%"<<endl; flag[0] = 1; } break;
       case 20 : if (!flag[1]) { cout<<"20%"<<endl; flag[1] = 1; } break;
       case 30 : if (!flag[2]) { cout<<"30%"<<endl; flag[2] = 1; } break;
       case 40 : if (!flag[3]) { cout<<"40%"<<endl; flag[3] = 1; } break;
       case 50 : if (!flag[4]) { cout<<"50%"<<endl; flag[4] = 1; } break;
       case 60 : if (!flag[5]) { cout<<"60%"<<endl; flag[5] = 1; } break;
       case 70 : if (!flag[6]) { cout<<"70%"<<endl; flag[6] = 1; } break;
       case 80 : if (!flag[7]) { cout<<"80%"<<endl; flag[7] = 1; } break;
       case 90 : if (!flag[8]) { cout<<"90%"<<endl; flag[8] = 1; } break;
       case 99 : if (!flag[9]) { cout<<"100%"<<endl; flag[9] = 1; } break;
       default : break;
       
       }
       
     }

   cout<<" Entries "<<ecalEnergies.size()<<endl;
   cout<<" Vetoed events "<<veto<<endl;
   //exit(0);

}
TTree* sTree;

vector<double> ETrueEnergies;  //The values that are taken from the root file
vector<double> ecalEnergies;
vector<double> hcalEnergies;
vector<double> etas;
vector<double> phis;


void calibChris() 
{
   gROOT->Reset();
   gStyle->SetCanvasColor(0);

   // InitBarrelAlpha();
   // LoadOldThresholds();
   // //LoadNewThresholds();

   gStyle->SetOptFit(0);

   
   TChain* chain= new TChain("s");
   
   //chain->Add("./../PGun_2_500_10_0_3_upgrade2018_ECAL_pfB.root");
   chain->Add("./step3.root");
   sTree = (TTree*)chain;
   cout<<"Reading input tree..."<<endl;
   getValuesFromTree(sTree, ETrueEnergies, ecalEnergies, 
                     hcalEnergies, etas, phis);


   int n=30,m=50;

TH2F* rawecalvsrawhcal = new TH2F("rawecalvsrawhcal","uncorrected Ecal energy vs uncorrected Hcal energy",500,-2,sampleRangeHigh,500,-2,sampleRangeHigh);
TH1F* uncorr_ecal = new TH1F("uncorr_ecal","Uncorrected Ecal energy ",500,-2,sampleRangeHigh);
TH1F* uncorr_hcal = new TH1F("uncorr_hcal","Uncorrected Hcal energy ",500,-2,sampleRangeHigh);
TH1F* response_EH = new TH1F("response_EH", "1D Response distribution",125,-1.5,2);
TH1F* response_H = new TH1F("response_H", "1D Response distribution",125,-1.5,2);
TH1F* response_EH_0_1 = new TH1F("response_EH_0_1", "1D Response distribution",125,-1.5,2);
 TH2F* rawhcal_response = new TH2F("rawhcal_response","Response vs Etrue",500,0,sampleRangeHigh,150,-1.5,1.5); 
 TH2F* rawecalhcal_response = new TH2F("rawecalhcal_response","Response vs Etrue",500,0,sampleRangeHigh,150,-1.5,1.5); 


TH2F* rawetrue_vs_etas_EH = new TH2F("rawetrue_vs_etas_EH","Raw etrue energy vs etas with response (avg hcalenergy+ecalenergy/etrue) in z axis",n,0,3.0,m,0,500);
 TH2F* rawetrue_vs_etas_H = new TH2F("rawetrue_vs_etas_H","Raw etrue energy vs etas with response (avg hcalenergy/etrue) in z axis",n,0,3.0,m,0,500);

 TH2F* rawecal_vs_etas = new TH2F("rawecal_vs_etas","Raw ecal energy vs etas",50,0,3.0,500,0,500);
 TH2F* rawhcal_vs_etas = new TH2F("rawhcal_vs_etas","Raw hcal energy vs etas",50,0,3.0,500,0,500);
 TH2F* rawecal_vs_etas_EH = new TH2F("rawecal_vs_etas_EH","Raw ecal energy vs etas",50,0,3.0,500,0,500);
 TH2F* rawhcal_vs_etas_EH = new TH2F("rawhcal_vs_etas_EH","Raw hcal energy vs etas",50,0,3.0,500,0,500);
 TH2F* rawhcal_vs_etas_H = new TH2F("rawhcal_vs_etas_H","Raw hcal energy  vs etas",50,0,3.0,500,0,500);
 TH2F* rawetrue_vs_etas_w_1_EH = new TH2F("rawetrue_vs_etas_w_1_EH","Raw etrue energy vs etas",n,0,3.0,m,0,500);
 TH2F* rawetrue_vs_etas_w_1_H = new TH2F("rawetrue_vs_etas_w_1_H","Raw etrue energy vs etas",n,0,3.0,m,0,500);

//plots summary
   for(double i=2; i <= ETrueEnergies.size() ; i++)
     {       // if (ETrueEnergies[i] >= 10  && ETrueEnergies[i] < 20)
	 {
	    // if(fabs(etas[i]) > 0.5) continue;// delete me                                                                                 
         
	   rawecal_vs_etas->Fill(fabs(etas[i]),ecalEnergies[i]);
	   rawhcal_vs_etas->Fill(fabs(etas[i]),hcalEnergies[i]);
	      
	    // if (fabs(etas[i]) < 0 ) continue;
	  if (ecalEnergies[i] != 0 && hcalEnergies[i] != 0 )
	    {if(ecalEnergies[i] > 0)
		 {rawecalvsrawhcal->Fill(ecalEnergies[i] ,hcalEnergies[i]);
	       uncorr_ecal->Fill(ecalEnergies[i]);
	       uncorr_hcal->Fill(hcalEnergies[i]);
	       response_EH->Fill((ecalEnergies[i] + hcalEnergies[i]-ETrueEnergies[i])/ETrueEnergies[i]);
	       rawecal_vs_etas_EH->Fill(fabs(etas[i]),ecalEnergies[i]);
	       rawhcal_vs_etas_EH->Fill(fabs(etas[i]),hcalEnergies[i]);
	       rawetrue_vs_etas_w_1_EH->Fill(fabs(etas[i]),ETrueEnergies[i]);
	       rawetrue_vs_etas_EH->Fill(fabs(etas[i]),ETrueEnergies[i],(hcalEnergies[i]+ecalEnergies[i])/ETrueEnergies[i]);   
	       rawecalhcal_response->Fill(ETrueEnergies[i],(ecalEnergies[i] + hcalEnergies[i]-ETrueEnergies[i])/ETrueEnergies[i]);	
	 }
	      if (ecalEnergies[i] < 1)
		response_EH_0_1->Fill((ecalEnergies[i] + hcalEnergies[i]-ETrueEnergies[i])/ETrueEnergies[i]);
	      
	      
	    }
	   if (ecalEnergies[i] == 0 && hcalEnergies[i] !=0)
	     { response_H->Fill((ecalEnergies[i] + hcalEnergies[i]-ETrueEnergies[i])/ETrueEnergies[i]);
	       rawhcal_response->Fill(ETrueEnergies[i],(ecalEnergies[i] + hcalEnergies[i]-ETrueEnergies[i])/ETrueEnergies[i]);
	       rawhcal_vs_etas_H->Fill(fabs(etas[i]),hcalEnergies[i]);
	      rawetrue_vs_etas_w_1_H->Fill(fabs(etas[i]),ETrueEnergies[i]);
	       rawetrue_vs_etas_H->Fill(fabs(etas[i]),ETrueEnergies[i],(hcalEnergies[i]+ecalEnergies[i])/ETrueEnergies[i]);  
	  
	     }	
	 }

   }
     response_EH->Scale(1.0/response_EH->Integral());
     response_H->Scale(1.0/response_H->Integral());
     uncorr_ecal->Scale(1.0/uncorr_ecal->Integral());
     uncorr_hcal->Scale(1.0/uncorr_hcal->Integral());
     response_EH_0_1->Scale(1.0/response_EH_0_1->Integral());

     rawecalvsrawhcal->SetXTitle("uncorrected ecal energy GeV");
     rawecalvsrawhcal->SetYTitle("uncorrected hcal energy GeV");
     response_EH->SetLineColor(4);
     response_EH_0_1->SetLineColor(2);
     response_H->SetLineColor(1);
     TH2F* h_frac_EH = (TH2F*)rawetrue_vs_etas_EH->Clone("fraction_EH");
     h_frac_EH->GetZaxis()->SetTitle("z frac");
     h_frac_EH->Divide(h_frac_EH, rawetrue_vs_etas_w_1_EH, 1., 1., "B");

     TH2F* h_frac_H = (TH2F*)rawetrue_vs_etas_H->Clone("fraction_H");
     h_frac_H->GetZaxis()->SetTitle("z frac");
     h_frac_H->Divide(h_frac_H, rawetrue_vs_etas_w_1_H, 1., 1., "B");
     // TPaveStats *s = (TPaveStats*)response_H->GetListOfFunctions()->FindObject("stats");
     // s->SetY1NDC (0.8);
     // s->SetY2NDC (1);
     // s->SetX1NDC (0.7);
     // s->SetX1NDC (1);
     // s->SetTextColor(response_H->GetLineColor());

     // TPaveStats *t = (TPaveStats*)response_EH->GetListOfFunctions()->FindObject("stats");
     // t->SetY1NDC (0.6);
     // t->SetY2NDC (0.8);
     // t->SetX1NDC (0.7);
     // t->SetX1NDC (1);
     // t->SetTextColor(response_EH->GetLineColor());

     // TPaveStats *u = (TPaveStats*)response_EH_0_1->GetListOfFunctions()->FindObject("stats");
     // u->SetY1NDC (0.4);
     // u->SetY2NDC (0.6);
     // u->SetX1NDC (0.7);
     // u->SetX1NDC (1);
     // u->SetTextColor(response_EH_0_1->GetLineColor());

     // response_EH->Draw("hist");
     // response_EH_0_1->Draw("hist sames");
     // response_H->Draw("hist sames");
      rawecalhcal_response->Draw();
     //  uncorr_hcal->Draw("sames");
     // rawecalvsrawhcal->Draw("colz");
     
      TFile* file=new TFile("output.root","recreate");
      
      rawecalvsrawhcal->Write();
      uncorr_ecal->Write();
      uncorr_hcal->Write();
      response_EH->Write();
      response_H->Write();
      response_EH_0_1->Write();
      rawhcal_response->Write();
      rawetrue_vs_etas_EH->Write();
      rawecal_vs_etas->Write();
      rawhcal_vs_etas->Write();
      rawecal_vs_etas_EH->Write();
      rawhcal_vs_etas_EH->Write();
      rawhcal_vs_etas_H->Write();
      rawetrue_vs_etas_w_1_EH->Write();
      rawetrue_vs_etas_w_1_H->Write();
      rawetrue_vs_etas_H->Write();
      
      h_frac_EH->Write();
      h_frac_H->Write();
      file->Close();
}
