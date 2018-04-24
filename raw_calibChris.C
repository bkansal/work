
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

//#include "CrystalBall.C"

using namespace std;

unsigned sampleRangeHigh = 20;
bool old_logic = false;

vector<float> assignvalues(vector<float> *pfcID_, vector<float> *Ecalenergy_, 
			   vector<float> *Hcalenergy_, vector<float> *dr) {

  vector<float> energies;
  float e = 0.0 , h = 0.0;
  for(unsigned ii = 0; ii < pfcID_->size(); ii++) {
    //cout<<" pfcID_:" << pfcID_->at(ii) << endl;
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
       
       // if(ecal_<0.4) continue; //FIXME MM
       if (fabs(eta_) < 2.4 && p_ == 0) continue;
       //if (fabs(eta_) > 2.5 && (true_/cosh(eta_) < 5)) { continue;}
       //if (true_ < 48 || true_ > 52 ) continue;
       //if (true_ < 178 || true_ > 182 ) continue;
       //if (true_/cosh(eta_) < 5 ) continue;
       //////  HEP17 Veto
       //if (phi_ < -0.4 && phi_ > -1.0 && eta_ < 3.0 && eta_ > 1.5) { veto++; continue; } 
       //if (fabs(eta_) < 1.5)
       // if(true_!=60    ) continue;
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
   
   chain->Add("./../PGun_2_500_10_0_3_upgrade2018_ECAL_pfB.root");

   sTree = (TTree*)chain;
   cout<<"Reading input tree..."<<endl;
   getValuesFromTree(sTree, ETrueEnergies, ecalEnergies, 
                     hcalEnergies, etas, phis);



TH2F* rawecalvsrawhcal = new TH2F("rawecalvsrawhcal","uncorrected Ecal energy vs uncorrected Hcal energy",500,-2,sampleRangeHigh,500,-2,sampleRangeHigh);
TH1F* uncorr_ecal = new TH1F("uncorr_ecal","Uncorrected Ecal energy ",500,-2,sampleRangeHigh);
TH1F* uncorr_hcal = new TH1F("uncorr_hcal","Uncorrected Hcal energy ",500,-2,sampleRangeHigh);
TH1F* response_EH = new TH1F("response_EH", "1D Response distribution",500,-1.5,2);
TH1F* response_H = new TH1F("response_H", "1D Response distribution",500,-1.5,2);
TH1F* response_EH_0_1 = new TH1F("response_EH_0_1", "1D Response distribution",500,-1.5,2);





//plots summary
   for(double i=2; i <= ETrueEnergies.size() ; i++)
     {        if (ETrueEnergies[i] >= 10  && ETrueEnergies[i] < 20)
	 {
	   // if(fabs(etas[i]) > 1.5) continue;// delete me                                                                                               
	   // if (fabs(etas[i]) < 0 ) continue;
	  if (ecalEnergies[i] != 0 && hcalEnergies[i] != 0 )
	    {if(ecalEnergies[i] > 0)
		 {rawecalvsrawhcal->Fill(ecalEnergies[i] ,hcalEnergies[i]);
	       uncorr_ecal->Fill(ecalEnergies[i]);
	       uncorr_hcal->Fill(hcalEnergies[i]);
	       response_EH->Fill((ecalEnergies[i] + hcalEnergies[i]-ETrueEnergies[i])/ETrueEnergies[i]);
		 }
	      if (ecalEnergies[i] < 1)
		response_EH_0_1->Fill((ecalEnergies[i] + hcalEnergies[i]-ETrueEnergies[i])/ETrueEnergies[i]);

	      
	    }
	   if (ecalEnergies[i] == 0 && hcalEnergies[i] !=0)
	     response_H->Fill((ecalEnergies[i] + hcalEnergies[i]-ETrueEnergies[i])/ETrueEnergies[i]);

	 }
   }
     response_EH->Scale(1.0/response_EH->Integral());
     response_H->Scale(1.0/response_H->Integral());
     uncorr_ecal->Scale(1.0/uncorr_ecal->Integral());
     uncorr_hcal->Scale(1.0/uncorr_hcal->Integral());
     response_EH_0_1->Scale(1.0/response_EH_0_1->Integral());

     rawecalvsrawhcal->SetXTitle("uncorrected ecal energy GeV");
     rawecalvsrawhcal->SetYTitle("uncorrected hcal energy GeV");
     // response_EH->Draw("hist");
     // response_H->Draw("hist sames");
     //uncorr_ecal->Draw();
     //  uncorr_hcal->Draw("sames");
     // rawecalvsrawhcal->Draw("colz");
     
      TFile* file=new TFile("output.root","recreate");
      rawecalvsrawhcal->Write();
      uncorr_ecal->Write();
      uncorr_hcal->Write();
      response_EH->Write();
      response_H->Write();
      response_EH_0_1->Write();

      file->Close();
}
