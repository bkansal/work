
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

double sigC_ = 5.;

// unsigned sampleRangeHigh = 200;
unsigned sampleRangeHigh = 500;
bool freezeparameters = true;
bool useMean = false;
bool useMedian = false;
bool changeRange =false;
bool old_logic = false;
bool drawpT = false;
bool drawResoFit = false;
bool saveCanvas = true;
char* _region_ = (char*)"EC_outside_tracker";
//char* _region_ = (char*)"EC_within_tracker";
//char* _region_ = (char*)"barrel";
//char* _region_ = (char*)"Full";

float _etaMin_ = 0.0;
float _etaMax_ = 0.0;


//threshold
/////////////MM
double aEH = 3.5;
double aE = 3.5;

// double aEH = 4.5;
// double aE = 4.5;

double aH = 2.5;//3.0;

double aEHe = 3.5;
double aEe = 3.5;

//spandey Feb 23 2018
// double aEHe = 4.5;
// double aEe = 4.5;

double aHe = 2.5;
/////////////MM


/*
//threshold
/////////////SP
double aEH = 3.675;
double aE = 3.675;
double aH = 2.625;//3.0;

double aEHe = 3.675;
double aEe = 3.675;
double aHe = 2.625;
/////////////SP
*/


double lBs = 0.25; //0.25 B //0.5 E
double mBs=1.; //1 B
double hBs=5.0; //10 40 B 50 E

double RBE = 4; //rebinning for alpha-beta curves

//Eta factor : ecal*factor + hcal
double factorB = 1.3;//Optimal //A factor put in by hand to make the eta dependence agree
double factorE = 1.3; //better. Should fit for the value later. 
//double factorE = 1.7; //better. Should fit for the value later.  //shubham

vector<double> sigmas;

void LoadOldThresholds() {
  aEH = 3.5; aE = 3.5; aH = 2.5; //3.5 2.5
  aEHe = 3.5; aEe = 3.5; aHe = 2.5;
}

//spandey
void LoadNewThresholds() {
  aEH = 3.8; aE = 3.8; aH = 3.0; //3.5 2.5
  aEHe = 3.8; aEe = 3.8; aHe = 3.0;
}


///////////////////////////////////////////////////////////////////////////////
//Container class that holds all the coefficients for a particular Etrue bin
///////////////////////////////////////////////////////////////////////////////

class ABC
{
   private:

      vector<double> ETrueEnergies_; //Variables that fall in this ETrue bin, 
      vector<double> ecalEnergies_;  //which is defined by binLowEdge and 
      vector<double> hcalEnergies_;  //binHighEdge. sigmaEcalHcal is simply the
      vector<double> etas_;          //uncertainty from the ecal and hcal(not
      vector<double> sigmaEcalHcal_; //uncertainty in the calibration consts)
      double binLowEdge_;
      double binHighEdge_;
      double etaMinFit_;     //From what etas to fit the B and C constants.
      double etaMaxFit_;     
      double etaMinEtaFit_;  //From what etas to fit the Alpha and Beta 
      double etaMaxEtaFit_;  //constants.
      bool isBarrel_;
      

      double ETrueAverage_; //Average and RMS of the ETrue values in the bin.
      double ETrueRMS_;
      double a_;  //Constant values.
      double b_;
      double c_;

      double sigmaB_;  //Uncertainty in the constants
      double sigmaC_;

   public:
      ABC(double binLowEdge, double binHighEdge, bool isBarrel); 

      bool addEntry(double ETrueEnergy, double ecalEnergy, double hcalEnergy,
                    double eta);  //Adds an event to the ETrue bin
      double getBinLowEdge();
      double getBinHighEdge();
      bool isBarrel();  //Checks if it is a barrel-type constant storage
      bool isEmpty();   //Checks if its empty
      bool isEmptyInFitRange();   //Checks if its empty in eta fit range
      unsigned getSize();        //Returns the various stored variables in the
      double getETrueAverage();  //ABC object
      double getETrueRMS();
      double getA();
      double getB();
      double getC();
      double getSigmaB();
      double getSigmaC();
       
      double getETrue(unsigned i);
      double getEcal(unsigned i); //Returns b*ecal for entry i
      double getHcal(unsigned i); //Returns c*hcal for entry i
      double getEta(unsigned i);
      double getNEntries();

      void computeETrueAverage();  //Computes the various calibration constants
      void computeETrueRMS();      //and other stored elements in the object.
      void computeA(double a);     //Right now the constant "a" is not computed
      void computeB();             //but rather just set.
      void computeC();
      bool computeBC();
      void clear();
};

ABC::ABC(double binLowEdge, double binHighEdge, 
	 bool _isBarrel) 

{
  binLowEdge_ = binLowEdge;
  binHighEdge_ = binHighEdge;
  isBarrel_ = _isBarrel;

   if(isBarrel_)
   {
      etaMinFit_ = 0.0;
      etaMaxFit_ = 1.0;
      etaMinEtaFit_ = 0.0;
      etaMaxEtaFit_ = 1.3;
   }
   else
   {
      etaMinFit_ = 1.6;
      etaMaxFit_ = 2.2; //FIXME 2.2
      etaMinEtaFit_ = 1.6;
      etaMaxEtaFit_ = 2.8; //FIXME 2.8
   }
   
   a_ = 0;
   b_ = 0;
   c_ = 0;
   ETrueAverage_ = 0;
   ETrueRMS_ = 0;
   sigmaB_ = 0;
   sigmaC_ = 0;
}
bool ABC::addEntry(double ETrue, double ecalEnergy, double hcalEnergy, double eta)
{

   
   double sigmaEcalHcal;

   if(isBarrel_) 
      sigmaEcalHcal = sqrt(0.08*0.08 + 1.04*1.04*(std::max(ecalEnergy + 
                                                           hcalEnergy, 1.0)));
   else
      sigmaEcalHcal = sqrt(0.04*0.04 + 1.80*1.80*(std::max(ecalEnergy + 
                                                           hcalEnergy, 1.0)));

   if( fabs(ecalEnergy + hcalEnergy + a_ - ETrue) > sigC_*sigmaEcalHcal || 
      (ecalEnergy + hcalEnergy) < 0.5 || ETrue < 1.0 || ETrue < binLowEdge_ || 
      ETrue> binHighEdge_ || fabs(eta) > etaMaxFit_ || fabs(eta) < etaMinFit_ )
      return false;


   ETrueEnergies_.push_back(ETrue);
   ecalEnergies_.push_back(ecalEnergy);
   hcalEnergies_.push_back(hcalEnergy);
   etas_.push_back(eta);
   sigmaEcalHcal_.push_back(sigmaEcalHcal);

   // if(ecalEnergy==0 && isBarrel_)
   //   cout<<hcalEnergy<<"  "<<fabs(ecalEnergy + hcalEnergy + a_ - ETrue)<<"   "<<sigmaEcalHcal<<endl;

   return true;
}

double ABC::getBinLowEdge() {return binLowEdge_;}
double ABC::getBinHighEdge() {return binHighEdge_;}
double ABC::getETrueAverage() {return ETrueAverage_;}
double ABC::getETrueRMS() {return ETrueRMS_;}
double ABC::getA() {return a_;}
double ABC::getB() {return b_;}
double ABC::getC() {return c_;}
double ABC::getSigmaB() {return sigmaB_;}
double ABC::getSigmaC() {return sigmaC_;}

double ABC::getETrue(unsigned i) {return ETrueEnergies_[i];}
double ABC::getEcal(unsigned i) {return ecalEnergies_[i];}
double ABC::getHcal(unsigned i) {return hcalEnergies_[i];}
double ABC::getEta(unsigned i) {return etas_[i];}
double ABC::getNEntries(){return ETrueEnergies_.size();}
bool ABC::isBarrel() {return isBarrel_;}
bool ABC::isEmpty() 
{
   if(ETrueEnergies_.size() == 0) return true;
   else return false;
}
bool ABC::isEmptyInFitRange() 
{
  //  cout<<"etas_.size():"<<etas_.size()<<endl;
   for(unsigned i = 0; i < etas_.size(); i++)
   {
      if(fabs(etas_[i]) < etaMaxFit_ && fabs(etas_[i]) > etaMinFit_)
         return false;
   }
   return true;
}
unsigned ABC::getSize()
{
   return ETrueEnergies_.size();
}
void ABC::computeETrueAverage()
{
   double totalETrue = 0;
   int numberSkipped = 0;

   for(unsigned i = 0; i < ETrueEnergies_.size(); i++)
   {
     if(fabs(etas_[i]) > etaMaxFit_ || fabs(etas_[i]) < etaMinFit_ ) 
      {
         numberSkipped++;
         continue;
      }   

     totalETrue += ETrueEnergies_[i];
   }
   
   ETrueAverage_ = totalETrue/(ETrueEnergies_.size() - numberSkipped);
} 
void ABC::computeETrueRMS()
{
   double totalETrueSquared = 0;
   int numberSkipped = 0;

   for(unsigned i = 0; i<ETrueEnergies_.size(); i++)
   {
      if(fabs(etas_[i]) > etaMaxFit_ || fabs(etas_[i]) < etaMinFit_ ) 
      {
         numberSkipped++;
         continue;
      }
       
      totalETrueSquared += ETrueEnergies_[i]*ETrueEnergies_[i];
   }
   
   ETrueRMS_ = sqrt(
      (totalETrueSquared/(ETrueEnergies_.size() - numberSkipped) -
       ETrueAverage_*ETrueAverage_)/(ETrueEnergies_.size() - numberSkipped)); 
} 
void ABC::computeA(double a) 
{
   a_ = a; 
}
void ABC::computeB() 
{
   double totalEcalSquared = 0;
   double totalEMinusATimesEcal = 0;
   
   for (unsigned i = 0; i < ETrueEnergies_.size(); i++)
   {
      if(fabs(etas_[i]) > etaMaxFit_ || fabs(etas_[i]) < etaMinFit_ ) continue;

      totalEcalSquared += 2* ecalEnergies_[i]*ecalEnergies_[i]/
      (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);
      totalEMinusATimesEcal += 2*(ETrueEnergies_[i] - a_)*ecalEnergies_[i]/
      (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);
   }
   
   b_ = totalEMinusATimesEcal/totalEcalSquared;
   sigmaB_ = sqrt(1/totalEcalSquared);
}
void ABC::computeC() 
{
   double totalHcalSquared = 0;
   double totalEMinusATimesHcal = 0;
   
   //cout<<" size C : "<<ETrueEnergies_.size()<<endl;

   for (unsigned i = 0; i < ETrueEnergies_.size(); i++)
   {
     totalHcalSquared += 2*hcalEnergies_[i]*hcalEnergies_[i]/
       (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);
     totalEMinusATimesHcal += 2*(ETrueEnergies_[i] - a_)*hcalEnergies_[i]/
       (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);
   }

   c_ = totalEMinusATimesHcal/totalHcalSquared;
   //cout<<totalEMinusATimesHcal<<"   "<<totalHcalSquared<<"   "<<c_<<endl;
   //cout<<"totalEMinusATimesHcal: "<<totalEMinusATimesHcal<<",   totalHcalSquared: "<<totalHcalSquared<<",   c_: "<<c_<<endl;

   // //FIXME MM
   // TH1F* cplot=new TH1F("cplot","",3000,-1, 5); 
   // float median=c_;
   // for (unsigned i = 0; i < ETrueEnergies_.size(); i++)
   //   {
   //     cplot->Fill( (ETrueEnergies_[i] - a_)/hcalEnergies_[i] );
   //   }

   // //  cout<<ETrueEnergies_[0]<<"  "<<" c_ "<<c_<<"   "<<(c_ + median1(cplot))/2.<<"   "<<median1(cplot)<<endl;
   // // if(ETrueEnergies_[0]>20)
   // //   c_ = (c_ + median1(cplot))/2.;
   // if(ETrueEnergies_[0]<=20)
   //   c_ = median1(cplot);

   // delete cplot;
   
   //================================


   sigmaC_ = sqrt(1/totalHcalSquared);
}
bool ABC::computeBC()
{
   ROOT::Math::SMatrix<double,2, 2, ROOT::Math::MatRepStd<double,2> > coeffs;
   ROOT::Math::SVector<double, 2> consts;
   ROOT::Math::SVector<double, 2> values;
   bool isInverted;
 
   coeffs(0, 0) = 0;
   coeffs(0, 1) = 0;
   coeffs(1, 0) = 0;
   coeffs(1, 1) = 0;
   consts(0) = 0;
   consts(1) = 0;
   //Create the matrix that will be inverted and the vector which will multiply
   //that matrix to find the b and c calibration constants.

   //cout<<" size AB : "<<ETrueEnergies_.size()<<endl;

   for(unsigned i = 0; i < ETrueEnergies_.size(); ++i)
   {
      if(fabs(etas_[i]) > etaMaxFit_ || fabs(etas_[i]) < etaMinFit_ ) continue;
      
      // if(i<10)
      // 	cout<<ETrueEnergies_[i]<<"   "<<ecalEnergies_[i]<<"   "<<hcalEnergies_[i]<<"   "<<etas_[i]<<"  "<<sigmaEcalHcal_[i]<<"   "<<a_;


      coeffs(0, 0) += 2*ecalEnergies_[i]*ecalEnergies_[i]/
         (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);
      
      coeffs(0, 1) += 2*ecalEnergies_[i]*hcalEnergies_[i]/
         (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);
      
      coeffs(1, 0) += 2*ecalEnergies_[i]*hcalEnergies_[i]/
         (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);
      
      coeffs(1, 1) += 2*hcalEnergies_[i]*hcalEnergies_[i]/
         (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);
      
      consts(0) += 2*(ETrueEnergies_[i] - a_)*ecalEnergies_[i] /
         (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);
      
      consts(1) += 2*(ETrueEnergies_[i] - a_)*hcalEnergies_[i]/
         (sigmaEcalHcal_[i]*sigmaEcalHcal_[i]);

      // if(i<10)
      // 	cout<<" ---> "<<coeffs(0, 0)<<"  "<<coeffs(0, 0)
      // 	    <<"  "<<coeffs(0, 0)<<"  "<<coeffs(0, 0)
      // 	    <<"  "<<consts(0)<<"  "<<consts(1)<<endl;

   }
   isInverted = coeffs.Invert();
   // coeffs.Print(cout); cout<<endl;
   // consts.Print(cout);


   if(isInverted && sqrt(coeffs(0,0)) <  100000) //Make sure it inverted successfully (i.e. det != 0)
   {

      values = coeffs*consts;
      
      b_ = values(0);
      c_ = values(1);
      sigmaB_ = sqrt(coeffs(0,0));
      sigmaC_ = sqrt(coeffs(1,1));

      //cout<<" b: "<<b_<<"  c: "<<c_<<endl;

      // //FIXME MM
      // TH1F* cplot=new TH1F("cplot","",3000,-1, 5); 
      //   for(unsigned i = 0; i < ETrueEnergies_.size(); ++i)
      // 	  {
      // 	    if(fabs(etas_[i]) > etaMaxFit_ || fabs(etas_[i]) < etaMinFit_ ) continue;
	    
      // 	    cplot->Fill( (ETrueEnergies_[i] - a_ - b_*ecalEnergies_[i])/hcalEnergies_[i] ); 
      // 	  }

      // 	//	cout<<ETrueEnergies_[0]<<"  "<<" c_ "<<c_<<"   "<<(c_ + median1(cplot))/2.<<"   "<<median1(cplot)<<endl;
	
      // 	c_ = (c_ + median1(cplot))/2.;
      // 	//	c_ = median1(cplot);

      // 	delete cplot;
	//================================
	
      return true;


   }
   else return false;
}

void ABC::clear()
{
   ETrueEnergies_.clear();
   ecalEnergies_.clear();
   hcalEnergies_.clear();
   etas_.clear();
   
   delete &binLowEdge_;
   delete &binHighEdge_;
   delete &isBarrel_;
   delete &ETrueAverage_;
   delete &ETrueRMS_;
   delete &a_;
   delete &b_;
   delete &c_;
   delete &sigmaB_;
   delete &sigmaC_;
}



///////////////////////////////////////////////////////////////////////////////
//The class that holds the calibration information over all values of ETrue. 
//It takes in a collection of ABC objects and fits each calibration 
//constant to a function. These functions are then used to find the calibrated
//energy.
///////////////////////////////////////////////////////////////////////////////
class Calibration
{
   private:
      double ETrueMax_;
      bool isBarrel_;
      
      vector<double> ETrueMeansABC_;
      vector<double> ETrueRMSsABC_;

      vector<double> as_;
      vector<double> bs_;
      vector<double> cs_;
 
      vector<double> sigmaBs_;
      vector<double> sigmaCs_;
      
      TGraph *graphA_ ;
      TGraphErrors *graphB_;
      TGraphErrors *graphC_;

      TGraph *graphBError_;
      TGraph *graphCError_;
      TF1* functionA_;
      TF1* functionB_;
      TF1* functionC_;
      
   public:
      Calibration();
      Calibration(double ETrueMax, bool isBarrel);
      //Adds a graph point by hand, i.e. putting in the individual values.
      void addGraphPoints(double ETrueAverage, double ETrueRMS,
                          double a, double b, double sigmaB, double C,
                          double sigmaC ); 
      //Adds a graph point by taking apart an ABC object.
      void addGraphPoints(ABC* abc);
  //Creates the graphs after the points have been added.
  void initializeGraphs(string option);  
  double getETrueMax();      //Returns the various objects  that the 
  TGraphErrors* getGraph();  //calibration class holds.
  TF1* getFunctionA();
  TF1* getFunctionB();
  TF1* getFunctionC();
  //Returns calibrated energy without any eta dependence.
  double getCalibratedEnergy(double ETrue, double ecalEnergy, 
			     double hcalEnergy); 
 
  void setETrueMax(double ETrueMax);
  bool fitAsToFunction(TF1 *functionA);  //Fits the functions to their 
  bool fitAsToFunction();                //graph points. One that takes in
  bool fitBsToFunction(TF1 *functionB);  //a function which then fits it to
  bool fitBsToFunction();                //that function, and one that is
      bool fitCsToFunction(TF1 *functionC);  //used to simply improve the fit.
      bool fitCsToFunction();
  
      void drawCoeffGraph(string graph, string tag);  //Makes and draws a graph of a 
      void drawSigmaGraph(string graph);  //coefficient. Takes in as an 
                                          //argument a, b, c, alpha, or beta.
      void printBs();
      void printCs();

};
Calibration::Calibration()
{
   ETrueMax_ = 1000;
   isBarrel_ = true;
}
Calibration::Calibration(double ETrueMax, bool isBarrel) 
          
{
   ETrueMax_ = ETrueMax;
   isBarrel_ = isBarrel;
}

void Calibration::addGraphPoints(double ETrueAverage, double ETrueRMS,
                                 double a, double b, double sigmaB,   
				 double c,double sigmaC)
{
   ETrueMeansABC_.push_back(ETrueAverage);
   ETrueRMSsABC_.push_back(ETrueRMS);
   as_.push_back(a);
   bs_.push_back(b);
   cs_.push_back(c);
   
   sigmaBs_.push_back(sigmaB);
   sigmaCs_.push_back(sigmaC);
 
}

void Calibration::addGraphPoints(ABC* abc)
{
   if(abc->isEmpty() || (abc->getSigmaC() == 0 && abc->getC() == 0)) return;

   ETrueMeansABC_.push_back(abc->getETrueAverage());
   ETrueRMSsABC_.push_back(abc->getETrueRMS());
   as_.push_back(abc->getA());
   bs_.push_back(abc->getB());
   cs_.push_back(abc->getC());


   sigmaBs_.push_back(abc->getSigmaB());
   sigmaCs_.push_back(abc->getSigmaC());
}

void Calibration::initializeGraphs(string option)
{
   vector<double> x;
   vector<double> sigmaX;
   vector<double> y;
   vector<double> sigmaY;
   
   if(option == "abc" || option == "ABC" || option == "all")
   {
      for(unsigned i = 0; i < ETrueMeansABC_.size(); i++)
      {
         if( bs_[i] == 0 && sigmaBs_[i] == 0.0) continue;
         
         x.push_back(ETrueMeansABC_[i]);
         sigmaX.push_back(ETrueRMSsABC_[i]);
         y.push_back(bs_[i]);
         sigmaY.push_back(sigmaBs_[i]);
         
      }
      
      graphB_ = new TGraphErrors(x.size(), &x[0], &y[0], &sigmaX[0], 
                                 &sigmaY[0]);
      
      x.clear();
      sigmaX.clear();
      y.clear();
      sigmaY.clear();
      
      for(unsigned i = 0; i < ETrueMeansABC_.size(); i++)
      {
         if( cs_[i] == 0 && sigmaCs_[i] == 0.0) continue;
         
         x.push_back(ETrueMeansABC_[i]);
         sigmaX.push_back(ETrueRMSsABC_[i]);
         y.push_back(cs_[i]);
         sigmaY.push_back(sigmaCs_[i]);
         
      }
      
      graphC_ = new TGraphErrors(x.size(), &x[0], &y[0], &sigmaX[0], 
                                 &sigmaY[0]);
      
      x.clear();
      sigmaX.clear();
      y.clear();
      sigmaY.clear();
   }
 
}

double Calibration::getETrueMax() {return ETrueMax_;}
TGraphErrors* Calibration::getGraph() {return graphB_;}
TF1* Calibration::getFunctionA() {return functionA_;}
TF1* Calibration::getFunctionB() {return functionB_;}
TF1* Calibration::getFunctionC() {return functionC_;}
double Calibration::getCalibratedEnergy(double ETrue, double ecalEnergy, 
                                        double hcalEnergy)
{
  double a = functionA_->Eval(ETrue);
  double b = functionB_->Eval(ETrue);
  double c = functionC_->Eval(ETrue);
  
   return a+ b*ecalEnergy + c*hcalEnergy;
}

void Calibration::setETrueMax(double ETrueMax){ETrueMax_ = ETrueMax;}

bool Calibration::fitAsToFunction(TF1 *functionA)
{
   functionA_ = functionA;
//   graphA_->Fit(functionA_->GetName(), "Q", "", 0, ETrueMax_);
   return true;
}
bool Calibration::fitAsToFunction()
{
   graphA_->Fit(functionA_->GetName(), "Q", "", 0, ETrueMax_);
   return true;
}

bool Calibration::fitBsToFunction(TF1 *functionB)
{
   functionB_ = functionB;
   graphB_->Fit(functionB_->GetName(), "Q", "", 1.5, ETrueMax_);
   return true;
}
bool Calibration::fitBsToFunction()
{
   graphB_->Fit(functionB_->GetName(), "Q", "", 1.5, ETrueMax_);
   return true;
}
bool Calibration::fitCsToFunction(TF1 *functionC)
{
   functionC_ = functionC;
   graphC_->Fit(functionC_->GetName(), "Q", "", 1.5, ETrueMax_);

   return true;
}
bool Calibration::fitCsToFunction()
{
   graphC_->Fit(functionC_->GetName(), "Q", "", 1.5, ETrueMax_);
   return true;
}

void Calibration::drawCoeffGraph(string graph, string tag)
{
  

  //  cout<<" Tag = "<<tag<<endl;

   string saveString;
   //TCanvas* canvas = new TCanvas( (graph+tag).c_str() , graph.c_str(), 1200,600 );
   TCanvas* canvas = new TCanvas( (graph+tag).c_str() , graph.c_str(), 500,300 );
   TH2F* histo = new TH2F("histoCG", "", sampleRangeHigh, 0, sampleRangeHigh, sampleRangeHigh,  -2.0, 2.0); 

   canvas->cd();
   histo->SetStats(0);
   histo->GetXaxis()->SetTitle("True Energy(GeV)");
   histo->GetYaxis()->SetTitle("Value");
   histo->Draw();
   
   gPad->SetGridx();
   gPad->SetGridy();
   
   // TLegend* leg=new TLegend(0.53,0.23,0.91,0.44);
   // leg->SetFillColor(0); leg->SetShadowColor(0);
   // TLine* line=new TLine(0,0,0,0);
   if(graph == "a" || graph == "A") 
     {
       histo->SetTitle("A Parameter");
       graphB_->SetTitle("A vs True Energy");
       histo->GetYaxis()->SetTitle( ("A coefficient ("+tag+")").c_str() );
       graphB_->SetMarkerStyle(22);
       graphB_->SetMarkerSize(1);
       graphB_->SetFillColor(0);

       graphB_->Draw("P");
       // if(tag=="EH") {
       // 	faBarrel->Draw("Lsame+");
       //   	faBarrel52x->Draw("Lsame+");
       // }

       //  leg->AddEntry(

       saveString = "ACoefficient" + tag + ".gif";
       canvas->SaveAs(saveString.c_str());
     }
   if(graph == "b" || graph == "B") 
     {
       histo->SetTitle("B Parameter");
       graphC_->SetTitle("B vs True Energy");
       histo->GetYaxis()->SetTitle( ("B coefficient ("+tag+")").c_str() );
       graphC_->SetMarkerStyle(22);
       graphC_->SetMarkerSize(1);
       graphC_->SetFillColor(0);

      graphC_->Draw("P");
      // if(tag=="EH") {
      // 	faBarrel->Draw("Lsame+");
	//   	faBarrel52x->Draw("Lsame+");
      //}

      //  leg->AddEntry(

      saveString = "BCoefficient" + tag + ".gif";
      canvas->SaveAs(saveString.c_str());
   }
   else if(graph == "c" || graph == "C") 
   {
      histo->SetTitle("C parameter");
      graphC_->SetTitle("C vs True Energy");
      histo->GetYaxis()->SetTitle( ("C coefficient ("+tag+")").c_str() );
      graphC_->SetMarkerStyle(22);
      graphC_->SetMarkerSize(1);
      //graphC_->SetMarkerColor(1);
      graphC_->SetFillColor(0);
      
      graphC_->Draw("P");
      // if(tag=="EH") {
      // 	fbBarrel->Draw("Lsame+");
      // 	//   	fbBarrel52x->Draw("Lsame+");
      // }
      // else {
      // 	fcBarrel->Draw("Lsame+");
      // 	//    	fcBarrel52x->Draw("Lsame+");
      // }

      saveString = "CCoefficient" + tag + ".gif";
      canvas->SaveAs(saveString.c_str());
   }
   else cout << "No graph with that name" <<endl;
 
  TLine* line=new TLine(1,2,1,2);
  line->SetLineColor(kRed+1);
  line->SetLineWidth(2);
}


void Calibration::drawSigmaGraph(string graph)
{
  
   TCanvas* canvas2 = new TCanvas();
   TH2F* histo2 = new TH2F("histoS", "", 100, 0, ETrueMax_, 100,  0, .1); 
   
   canvas2->cd();
   histo2->SetStats(0);
   histo2->Draw();
   
   gPad->SetGridx();
   gPad->SetGridy();

   if(graph == "b" || graph =="B")
   {
      graphBError_->SetMarkerStyle(22);
      graphBError_->SetMarkerColor(2);
      graphBError_->SetMarkerSize(.5);
      graphBError_->Draw("P");

   }
   else if(graph == "c" || graph =="C")
   {
      graphCError_->SetMarkerStyle(22);
      graphCError_->SetMarkerColor(2);
      graphCError_->SetMarkerSize(.5);
      graphCError_->Draw("P");
   }

   else cout<<"No graph with that name"<<endl;
   
   delete histo2;

}

void Calibration::printBs()
{
   
   for(unsigned i = 0; i < as_.size(); i++)
   {
     cout<<bs_[i];
   }
}

void Calibration::printCs()
{
   
   for(unsigned i = 0; i < as_.size(); i++)
   {
      cout<<cs_[i]<<endl;
   }
}


double CalculateMedian(TH1F* h)
{
  // TFile f("histos.root");
  // TH1F *h = (TH1F*)f.Get("hgaus");

  int numBins = h->GetXaxis()->GetNbins();
  Double_t* x = new Double_t[numBins];
  Double_t* y = new Double_t[numBins];
  for (int i = 0; i < numBins; i++) {
    x[i] = h->GetBinCenter(i);
    y[i] = h->GetBinContent(i);
  }
  double MedianOfHisto = TMath::Median(numBins, &x[0], &y[0]);

  return MedianOfHisto;

  // cout<<"Median -----------> \t"<<MedianOfHisto;
  // return 0;
}


///////////////////////////////////////////////////////////////////////////////
//Global functions used in the main of the code.
///////////////////////////////////////////////////////////////////////////////


void drawGausFit(TH2F* inHisto, TGraph& response, TGraph& resolution)
{
  
   if(inHisto->GetEntries() == 0) return;
   
   vector<TH1F*> ETrueBin;
   TF1* gaus; 
   string name;
   char num[4];
   float rebin = 1.;
   TCanvas* canvas;
   TCanvas* temp = new TCanvas();
   //TLine *line = new TLine(0.0,0.0,sampleRangeHigh,0.0);
   int rangelow_ = 2, rangehigh_ = sampleRangeHigh, bins_ = sampleRangeHigh;

   if (drawpT) {
     rangehigh_ = 100;
     bins_ = 100;
   }
   TLine *line = new TLine(0.0,0.0,rangehigh_,0.0);   
   TH2F* respHisto = new TH2F("respHisto", "", bins_, rangelow_, rangehigh_, 100, -0.5, 0.5);
   TH2F* resoHisto = new TH2F("resoHisto", "", bins_, rangelow_, rangehigh_, 200, 0.0, 1.0);



   TGraph averages;
   TGraph rmss;

   vector<double> ETrue;
   vector<double> gausMean; 
   vector<double> gausSigma;
   vector<double> average;
   vector<double> rms;
   
   TFile* file1=new TFile("projections.root","recreate");

   
   temp->cd();//This TCanvas is only used since when we do the fit down below 
              //it creates an unwanted TCanvas. We will get rid of it later on 
              //in the function.

   // TCanvas* cccc=new TCanvas("balda","bacla");
   //cout<<"ETrue.back(), gausMean[0].back()"<<endl;
   //cout<<"**********Draw Gaus**********"<<endl;
   for(unsigned bin = 2; bin < sampleRangeHigh; )
   {
     double x_min = -1.0, x_max = 1.0;
     if(strcmp(inHisto->GetName(),"corrEtaEndcapEcalHcal") == 0 && strcmp(_region_,"EC_outside_tracker") == 0 && false)
       x_min = -0.70;
      name = "histcorhybrid";
      sprintf(num, "%i", bin);
      name += num;
      //Split up the TH2 into many TH1's for each ETrue bin.
     
      ETrueBin.push_back((TH1F*)inHisto->ProjectionY(name.c_str(),bin, 
                                                     bin + 4*rebin));
      //cout <<"bin to  bin + 4*rebin: "<<bin<<" to "<<(bin + 4*rebin)<<endl;//"   "<<ETrueBin.back()->GetEntries()<<endl;
      if(ETrueBin.back()->GetEntries() > 5)
	{
	  //Fit each ETrue bin to a gaus (iteratively done to get better fit)
	  //cout<<"ETrueBin.back()->GetEntries():"<<ETrueBin.back()->GetEntries()<<endl;
	  if(bin > 2) {

	    gaus =new TF1("gaus","gaus(0)",-3,3);
	    gaus->SetParameters(500.,0.,0.2);
	    //ETrueBin.back()->Fit("gaus", "Q", "", -1.0, 1.0);
	    ETrueBin.back()->Fit("gaus", "Q", "", x_min, x_max);
	    //ETrueBin.back()->Fit("gaus", "Q", "", -0.7, 0.7);
	    
	    
	    gaus = ETrueBin.back()->GetFunction("gaus");
	    //cout<<name<<" "<<gaus->GetParameter(1)<<"   "<<gaus->GetParameter(2)<<endl;
	    	    
	    if(gaus->GetParameter(2) < 0)
	      goto here1;
	    x_min = gaus->GetParameter(1)- gaus->GetParameter(2);

	    if(strcmp(inHisto->GetName(),"corrEtaEndcapEcalHcal") == 0 && strcmp(_region_,"EC_outside_tracker") == 0 && false)
	      x_min = (gaus->GetParameter(1) - gaus->GetParameter(2)) < -0.7 ? -0.7 : (gaus->GetParameter(1) - gaus->GetParameter(2));

		 
	    //x_max = (gaus->GetParameter(1)+ 2* gaus->GetParameter(2))>1.0 ? 1.0 : (gaus->GetParameter(1)+ 2* gaus->GetParameter(2));
	    x_max = (gaus->GetParameter(1)+ 2* gaus->GetParameter(2));
	    
	    ETrueBin.back()->Fit("gaus", "Q", "", x_min, x_max);

	    // ETrueBin.back()->Fit("gaus", "Q", "",
	    // 			 gaus->GetParameter(1) - 2*gaus->
	    // 			 GetParameter(2), 1.0);  

	    gaus = ETrueBin.back()->GetFunction("gaus");

	    x_min = gaus->GetParameter(1)- gaus->GetParameter(2);
	    if(strcmp(inHisto->GetName(),"corrEtaEndcapEcalHcal") == 0 && strcmp(_region_,"EC_outside_tracker") == 0 && false)
	      x_min = (gaus->GetParameter(1) - gaus->GetParameter(2)) < -0.7 ? -0.7 : (gaus->GetParameter(1) - gaus->GetParameter(2));

	    //x_max = (gaus->GetParameter(1)+ 2* gaus->GetParameter(2))>1.0 ? 1.0 : (gaus->GetParameter(1)+ 2* gaus->GetParameter(2));
	    x_max = (gaus->GetParameter(1)+ 2* gaus->GetParameter(2));
	    
	    
	    ETrueBin.back()->Fit("gaus", "Q", "", x_min, x_max);


            // ETrueBin.back()->Fit("gaus", "Q", "",
            //                      gaus->GetParameter(1) - 2*gaus->
            //                      GetParameter(2), 1.0);
            gaus = ETrueBin.back()->GetFunction("gaus");

	  here1:
	    // if(bin<=16)
	    // gausMean.push_back(ETrueBin.back()->GetMean());
	    // else

	    
            //gausSigma.push_back(gaus->GetParameter(2)/(1.0 + min(0.0, gaus->GetParameter(1))));

	    if (useMedian) {
	      gausMean.push_back(CalculateMedian(ETrueBin.back()));
	      gausSigma.push_back(CalculateMedian(ETrueBin.back())/(1.0 + min(0.0, ETrueBin.back()->GetMean())));

	    }
	    else if (useMean) {
	      gausMean.push_back(ETrueBin.back()->GetMean());
	      gausSigma.push_back(ETrueBin.back()->GetRMS()/(1.0 + min(0.0, ETrueBin.back()->GetMean())));
	    }

	    else {
	      gausMean.push_back(gaus->GetParameter(1));
	      gausSigma.push_back(gaus->GetParameter(2)/(1.0 + min(0.0, gaus->GetParameter(1))));
	    }
	    //gausMean.push_back(gaus->GetParameter(1));
	    // if (bin > 24 || bin == 2) {
	    //   gausMean.push_back(gaus->GetParameter(1));
	    // }
	    // else {
	    //   gausMean.push_back(ETrueBin.back()->GetMean());
	    //   gaus->Delete();
	    // }

	  }
	   else {
	   
	     gaus =new TF1("gaus","gaus",-3,3);
	     gaus->SetParameters( 500, 10, 5, 0, 0.20 );
	     gaus->FixParameter(2,5);
	     ETrueBin.back()->Fit("gaus", "QN0", "", -1.0, 1.0);
	     ETrueBin.back()->Fit("gaus", "QN0", "", -1.0, 1.0);
	     ETrueBin.back()->Fit("gaus", "Q", "", -1.0, 1.0);

	     gausMean.push_back(gaus->GetParameter(3));
	     gausSigma.push_back(gaus->GetParameter(4)/
				 (1.0 + min(0.0, gaus->GetParameter(3))));
	   }


	  // cout<<bin<<"   "<<median1(ETrueBin.back())<<endl;

	  // TFile oFile( ("tmp/"+name+".root").c_str() ,"RECREATE");
	  // ETrueBin.back()->Write();
	  // gaus->Write();
	  // oFile.Close();

	  // cccc->cd();
	  // ETrueBin.back()->Draw();
	  // // //   gaus->Draw("same");
	  // cccc->SaveAs( ("tmp/"+name+".png").c_str() );
	  // cccc->SaveAs( ("tmp/"+name+".C").c_str() );

            ETrue.push_back(bin + 2.0*rebin);
	    //cout<<"bin:"<<bin<<", rebin:"<<rebin<<", bin + 2*rebin:"<<(bin + 2*rebin)<<endl;
	    //cout<<ETrue.back()<<", "<<gausMean.back()<<endl;
	    //shubham
	    //cout<<ETrue.back()<<" ";
	    //  if(bin<=16)
	    //  gausMean.push_back(ETrueBin.back()->GetMean());
	    // else
	    //   gausMean.push_back(gaus->GetParameter(1));

	    
	   
            average.push_back(ETrueBin.back()->GetMean());
            rms.push_back(ETrueBin.back()->GetRMS());
	    

	    //cout<<bin<<"   "<<ETrue.back()<<"   "<<ETrueBin.back()->GetMean()<<" <> "<<gausMean.back()<<"   "<<ETrueBin.back()->GetMeanError()<<"   "<<gaus->GetParError(1)<<endl;
	    //cout<<bin<<"   "<<ETrue.back()<<"   "<<ETrueBin.back()->GetMean()<<" <> "<<gausMean.back()<<"   "<<ETrueBin.back()->GetMeanError()<<endl;

	    if (false)
	      gaus->Delete();

	    (ETrueBin.back())->Write();


	}

      
      bin += 2*rebin;
      
      //Increase bin size with increasing ETrue since there are fewer high 
      //energy events than low energy ones.
      if(bin > 10) rebin = 2.0;
      if(bin > 100) rebin = 5.0; //20
      if(bin > 1000) rebin = 20.0; //50
      //delete gaus;
      
   }

   file1->Close();
   // delete cccc;

   response = TGraph(ETrue.size(), &ETrue[0], &gausMean[0]); //Fill the graphs
   //response = TGraph(ETrue.size(), &ETrue[0], &average[0]); //Fill the graphs
   resolution = TGraph(ETrue.size(),&ETrue[0], &gausSigma[0]);
   averages =  TGraph(ETrue.size(), &ETrue[0], &average[0]);
   rmss = TGraph(ETrue.size(), &ETrue[0], &rms[0]);

   //Set up the graphs to look how you want them to.
   response.SetMarkerStyle(22);
   response.SetMarkerSize(0.8);
   response.SetMarkerColor(4);

   resolution.SetMarkerStyle(22);
   resolution.SetMarkerSize(0.8);
   resolution.SetMarkerColor(4);

   averages.SetMarkerStyle(22);
   averages.SetMarkerSize(0.8);
   averages.SetMarkerColor(4);

   rmss.SetMarkerStyle(22);
   rmss.SetMarkerSize(0.8);
   rmss.SetMarkerColor(4);

   line->SetLineStyle(1);
   line->SetLineWidth(2);
   line->SetLineColor(2);


   //  gStyle->SetOptStat(0); 
   //gStyle->SetOptFit(0);
   //canvas = new TCanvas(("canvas "+ (string)(inHisto->GetName()) ).c_str(), ("Response and Resolution "+ (string)(inHisto->GetName())).c_str(), 1000, 500);
   //spandey
   //canvas = new TCanvas(("canvas "+ (string)(inHisto->GetName()) ).c_str(), ("Response and Resolution "+ (string)(inHisto->GetName())).c_str(), 800, 400);
   canvas = new TCanvas(("canvas "+ (string)(inHisto->GetName()) ).c_str(), ("Response and Resolution "+ (string)(inHisto->GetName())).c_str(), 500, 300);


 
   canvas->Divide(2, 1);
   temp->~TCanvas();  //destroy the TCanvas 

   canvas->cd(1);
   gPad->SetGridx();
   gPad->SetGridy();
   respHisto->SetStats(0);
   respHisto->SetTitle("Response");
   respHisto->Draw();
   response.Draw("P");
   line->Draw();

   canvas->cd(2);
   gPad->SetGridx();
   gPad->SetGridy();
   resoHisto->SetStats(0);
   resoHisto->SetTitle("Resolution");
   resoHisto->Draw();
   resolution.Draw("P");


   respHisto->GetYaxis()->SetTitle("(E_{cor}-E_{true})/E_{true}");
   respHisto->GetXaxis()->SetTitle("E_{true} [GeV]");

   resoHisto->GetYaxis()->SetTitle("#sigma(E)/E_{true}");
   resoHisto->GetXaxis()->SetTitle("E_{true} [GeV]");

   if(drawResoFit) {
     //MM Fit Resolution
     TF1* f=new TF1( ("ResoFit"+ (string)(inHisto->GetName())).c_str(),"sqrt([0]*[0]+[1]*[1]/x+[2]*[2]/(x*x))",20,1000);// 3.5*4
     f->SetParameters(0.06,1.20,0.);
     f->SetParLimits(0,0,10);
     f->SetParLimits(1,0,10);
     f->SetParLimits(2,0,10);
     resolution.Fit(("ResoFit"+ (string)(inHisto->GetName())).c_str(),"QR");
     resolution.Fit(("ResoFit"+ (string)(inHisto->GetName())).c_str(),"QR");
     resolution.Fit(("ResoFit"+ (string)(inHisto->GetName())).c_str(),"R");
     
     
     
     string legend;
     int fres0 = (int)(f->GetParameter(0)*100.);
     int fres1 = (int)(10.*(f->GetParameter(0)*100.-fres0));
     int fres2 = (int)(f->GetParameter(1)*100.);
     // char text[100];
     // sprintf(text,"#sigma/E = %i%/#sqrt{E} + %i.%i%",fres2,fres0,fres1);
     TString text = "#sigma/E = ";
     text+=(int)fres2;
     text+="%/#sqrt{E} + ";
     text+=(int)fres0;
     text+=".";
     text+=(int)fres1;
     
     legend += text;
     TLegend *leg=new TLegend(0.30,0.75,0.85,0.85);
     leg->AddEntry((&resolution),legend.c_str(),"lp");
     leg->SetTextSize(0.04);
     leg->Draw();
   }
   if (saveCanvas) {
     //string cname = ((string)(inHisto->GetName()) ) + ".gif";
     char  cname[200];
     sprintf(cname,  "%s_updatedCode.gif",inHisto->GetName());
     canvas->Print(cname);
   }

}



vector<float> assignvalues(vector<float> *pfcID_, vector<float> *Ecalenergy_, 
			   vector<float> *Hcalenergy_, vector<float> *dr) {

  vector<float> energies;
  float e = 0.0 , h = 0.0;
  //  cout<<"pfcID size ----> "<<pfcID_->size()<<endl<<"pfcID_->at(00) ----->"<<pfcID_->at(00)<<endl;
  for(unsigned ii = 0; ii < pfcID_->size(); ii++) {
    //    cout<<" pfcID_:" << pfcID_->at(ii) << endl;
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



/// Bin Manager =========================================
vector<double> BinsETrue;
vector<double> BinsETrueEta;

unsigned int GetETrueBin(double etrue) {

  bool find=false;

  int bm=0;
  int bM=BinsETrue.size()-1;
  if(etrue< BinsETrue[ bm ]) return -1;
  if(etrue> BinsETrue[ bM ]) return -1;
  int n=0;
  while(!find) {
    
    if(etrue <= BinsETrue[ bm+(bM-bm)/2 ] ) {
      bM = bm+(bM-bm)/2;
    }
    else {
      bm = bm+(bM-bm)/2;
    }
    if( fabs(bm-bM)==1 )  {
      return bm;
    }
   
    if(n>(int)BinsETrue.size()) return -1;
    n++;
  }
  return -1;
}



///======================================================


///////////////////////////////////////////////////////////////////////////////
//All the needed variables for the main (calibChris) function.
///////////////////////////////////////////////////////////////////////////////

TFile* inputFile;
TTree* sTree;

vector<double> ETrueEnergies;  //The values that are taken from the root file
vector<double> ecalEnergies;
vector<double> hcalEnergies;
vector<double> etas;
vector<double> phis;

vector<ABC*> barrelABCEcalHcal; //Vectors of the ABC objects
vector<ABC*> barrelABCEcal;     //which hold all the calibration 
vector<ABC*> barrelABCHcal;     //constants for an individual bin.
vector<ABC*> endcapABCEcalHcal;
vector<ABC*> endcapABCEcal;
vector<ABC*> endcapABCHcal;

TF1* functionBarrelEcalHcalA;     //Functions that the calibration equations
TF1* functionBarrelEcalHcalB;     //are fit to
TF1* functionBarrelEcalHcalC;  
TF1* functionEndcapEcalHcalA;
TF1* functionEndcapEcalHcalB;
TF1* functionEndcapEcalHcalC;  

TF1* functionBarrelHcalA;
TF1* functionBarrelHcalB;
TF1* functionBarrelHcalC;  
TF1* functionEndcapHcalA;
TF1* functionEndcapHcalB;
TF1* functionEndcapHcalC;  

//Calibration objects which hold the all the calibration costants as functions
//of ETrue. 
Calibration* barrelWithEcalHcalCalib = new Calibration(sampleRangeHigh, true);
Calibration* barrelWithEcalCalib = new Calibration(sampleRangeHigh, true);
Calibration* barrelWithHcalCalib = new Calibration(sampleRangeHigh, true);
Calibration* endcapWithEcalHcalCalib = new Calibration(sampleRangeHigh, false);
Calibration* endcapWithEcalCalib = new Calibration(sampleRangeHigh, false);
Calibration* endcapWithHcalCalib = new Calibration(sampleRangeHigh, false);

//Temporary varibles that will be used in for loops just to cut down on the 
//length of the lines of code.
double etrue;
double ecal;
double hcal;
double eta;
double bpar;
double cpar;
double etrueMax;
double barrelEcalHcalB;
double barrelEcalHcalC;
double barrelHcalC;
double endcapEcalHcalB;
double endcapEcalHcalC;
double endcapHcalC;
double correctedE;


const char* functionEndcapEcalHcalB_e;
const char* functionEndcapEcalHcalC_e;
const char* functionEndcapHcalC_e;
const char* functionBarrelEcalHcalB_e;
const char* functionBarrelEcalHcalC_e;
const char* functionBarrelHcalC_e;

//All the differenct types of TH2's that will be filled in order to make resolution and response plots

TH2F* raw = new TH2F("raw","", 1000, 0, 1000, 150, -1.5, 1.5);
TH2F* corrEta = new TH2F("corrEta", "", 1000, 0, 1000, 150, -1.5, 1.5);

TH2F* rawBarrel = new TH2F("rawBarrel","", 1000, 0, 1000, 150, -1.5, 1.5);
TH2F* corrBarrel = new TH2F("corrBarrel", "", 1000, 0, 1000, 150, -1.5, 1.5);
TH2F* corrEtaBarrel = new TH2F("corrEtaBarrel", "", 1000, 0, 1000, 150, -1.5, 
                               1.5);
TH2F* rawBarrelEcalHcal = new TH2F("rawBarrelEcalHcal","", 1000, 0, 1000, 150, 
                                   -1.5, 1.5);
TH2F* corrBarrelEcalHcal = new TH2F("corrBarrelEcalHcal", "", 1000, 0, 1000, 
                                    150, -1.5, 1.5);
TH2F* corrEtaBarrelEcalHcal = new TH2F("corrEtaBarrelEcalHcal", "", 1000, 0, 
                                       1000, 150, -1.5, 1.5);
TH2F* rawBarrelHcal = new TH2F("rawBarrelHcal","", 1000, 0, 1000, 150, -1.5, 1.5 );
			       //1000, 0.0, 5.0);
TH2F* corrBarrelHcal = new TH2F("corrBarrelHcal", "", 1000, 0, 1000, 150, -1.5,
                                1.5);
TH2F* corrEtaBarrelHcal = new TH2F("corrEtaBarrelHcal", "", 1000, 0, 1000, 150,
                                   -1.5, 1.5);

TH2F* rawEndcap = new TH2F("rawEndcap","", 1000, 0, 1000, 150, -1.5, 1.5);
TH2F* corrEndcap = new TH2F("corrEndcap", "", 1000, 0, 1000, 150, -1.5, 1.5);
TH2F* corrEtaEndcap = new TH2F("corrEtaEndcap", "", 1000, 0, 1000, 150, -1.5, 
                               1.5);
//TH2F* rawEndcapEcalHcal = new TH2F("rawEndcapEcalHcal","", 1000, 0, 1000, 150, -1.5, 1.5);
TH2F* rawEndcapEcalHcal = new TH2F("rawEndcapEcalHcal","", 1000, 0, 1000, 500, -1.5, 10.0);

TH2F* corrEndcapEcalHcal = new TH2F("corrEndcapEcalHcal", "", 1000, 0, 1000, 
                                    150, -1.5, 1.5);
//TH2F* corrEtaEndcapEcalHcal = new TH2F("corrEtaEndcapEcalHcal", "", 1000, 0, 
//                                     1000, 150, -1.5, 1.5);
TH2F* corrEtaEndcapEcalHcal = new TH2F("corrEtaEndcapEcalHcal", "", 1000, 0, 1000, 500, -1.5, 10.0);

TH2F* rawEndcapHcal = new TH2F("rawEndcapHcal","", 1000, 0, 1000, 150, -1.5, 
                               1.5);
TH2F* corrEndcapHcal = new TH2F("corrEndcapHcal", "", 1000, 0, 1000, 150, -1.5,
                                1.5);
TH2F* corrEtaEndcapHcal = new TH2F("corrEtaEndcapHcal", "", 1000, 0, 1000, 150,
                                   -1.5, 1.5);

//Temporary TGraphs to passed drawGausFit
TGraph response;
TGraph resolution;
TGraph responseRaw;
TGraph resolutionRaw;
TGraph responseCor;
TGraph resolutionCor;
TGraph responseEta;
TGraph resolutionEta;
TGraph responseEtaEtaEH;
TGraph responseEtaHCorrEtaEH;
TGraph responseEtaEtaH;
TGraph responseEtaHCorrEtaH;


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

   if (strcmp(_region_, "barrel") == 0) {
     _etaMin_ = 0.0;
     _etaMax_ = 1.5;
   }
   else if (strcmp(_region_, "EC_within_tracker") == 0 ) {
     _etaMin_ = 1.55;
     _etaMax_ = 2.5;
   }
   
   else if (strcmp(_region_, "EC_outside_tracker") == 0 ) {
     _etaMin_ = 2.5;
     //_etaMax_ = 3.0;
     _etaMax_ = 2.75;
   }
   
   else if (strcmp(_region_,  "Full") == 0 ) {
     _etaMin_ = 1.55;
     _etaMax_ = 3.0;
   }

   cout<<"Creating abc objects..."<<endl;
  
   BinsETrue.clear();
   BinsETrueEta.clear();

   for(double bin = 0.0; bin < 10.0; bin = bin + lBs)
     {
       barrelABCEcalHcal.push_back(new ABC(bin, bin + lBs, true));
       barrelABCHcal.push_back(new ABC(bin, bin + lBs, true));    
       endcapABCEcalHcal.push_back(new ABC(bin, bin + lBs, false));
       endcapABCHcal.push_back(new ABC(bin, bin + lBs, false));
       BinsETrue.push_back(bin);
     }
   
   
   
   for(double bin = 10.0; bin < 100.0 ; bin = bin + mBs) //2
     {
       barrelABCEcalHcal.push_back(new ABC(bin, bin + mBs, true));
       barrelABCHcal.push_back(new ABC(bin, bin + mBs, true));
       endcapABCEcalHcal.push_back(new ABC(bin, bin + mBs, false));
       endcapABCHcal.push_back(new ABC(bin, bin + mBs, false));
       BinsETrue.push_back(bin);
   }
   
   
   for(double bin = 100.0; bin < 1000.0 ; bin = bin + hBs) //10
   {
     barrelABCEcalHcal.push_back(new ABC(bin, bin + hBs, true));
     barrelABCHcal.push_back(new ABC(bin, bin + hBs, true));
     endcapABCEcalHcal.push_back(new ABC(bin, bin + hBs, false));
     endcapABCHcal.push_back(new ABC(bin, bin + hBs, false));  
     BinsETrue.push_back(bin);
   }
   BinsETrue.push_back( BinsETrue.back() + hBs );

   //Fill all the ABC Objects with their respective events. They are all 
   //divided up into the six possible case ( (endcap or barrel)x(ecalhcal or 
   //ecal or hcal))
   

   TH1F* EcalSpectrum=new TH1F("EcalSpectrum","EcalSpectrum",1000,0,100);

   cout<<"Filling abc objects..."<<endl;
   for( unsigned bin = 0; bin < barrelABCEcal.size(); ++bin)
     {
       barrelABCEcalHcal[bin]->computeA(aEH);
       barrelABCHcal[bin]->computeA(aH);

       endcapABCEcalHcal[bin]->computeA(aEHe);
       endcapABCHcal[bin]->computeA(aHe);
     }
   

   {

     unsigned bin = 0;
     for(unsigned entry = 0; entry < ETrueEnergies.size(); entry++)
       {
	 etrue = ETrueEnergies[entry];
	 ecal = ecalEnergies[entry];
	 hcal = hcalEnergies[entry];
	 eta = etas[entry];

	 if(hcal == 0.0) continue;
	 if( etrue <1 ) continue;
	 if( etrue >sampleRangeHigh ) continue;
	 // if( etrue <10 ) continue;
	 // if( etrue >12 ) continue;

	 bin = GetETrueBin( etrue );

	 if( ecal > 0.0 && hcal > 0.0)
	   {
	     barrelABCEcalHcal[bin]->addEntry(etrue, ecal, hcal, eta);
	     endcapABCEcalHcal[bin]->addEntry(etrue, ecal, hcal, eta);
            
	     if(eta<1.3)
	       EcalSpectrum->Fill(ecal);

	   }
	 else if(ecal > 0.0)
	   {
	     barrelABCEcal[bin]->addEntry(etrue, ecal, hcal ,eta);
	     endcapABCEcal[bin]->addEntry(etrue, ecal, hcal ,eta);
	   }
	 else if(hcal > 0.0)
	   {
	     barrelABCHcal[bin]->addEntry(etrue, ecal, hcal, eta);
	     endcapABCHcal[bin]->addEntry(etrue, ecal, hcal, eta);
	   }
         
       }
   }


   //Compute the calibration constants along with their uncertainties for each
   //ETrue bin, then add their values to a Calibration object.
   cout<<"Computing a, b, c coefficients..."<<endl;

   for(unsigned bin = 2; bin < barrelABCEcalHcal.size() - 1; ++bin)
   {
      
      if(!barrelABCEcalHcal[bin]->isEmptyInFitRange())
      { 
         barrelABCEcalHcal[bin]->computeETrueAverage();
         barrelABCEcalHcal[bin]->computeETrueRMS();
         barrelABCEcalHcal[bin]->computeA(aEH);
         barrelABCEcalHcal[bin]->computeBC();
	 //exit(0);
      }
      
      if(!barrelABCEcal[bin]->isEmptyInFitRange())
      { 
         barrelABCEcal[bin]->computeETrueAverage();
         barrelABCEcal[bin]->computeETrueRMS();
         barrelABCEcal[bin]->computeA(aEH);
         barrelABCEcal[bin]->computeB();
      }
      if(!barrelABCHcal[bin]->isEmptyInFitRange())
      { 
         barrelABCHcal[bin]->computeETrueAverage();
         barrelABCHcal[bin]->computeETrueRMS();
         barrelABCHcal[bin]->computeA(aH);
         barrelABCHcal[bin]->computeC();
      }
      if(!endcapABCEcalHcal[bin]->isEmptyInFitRange())
      {
         endcapABCEcalHcal[bin]->computeETrueAverage();
         endcapABCEcalHcal[bin]->computeETrueRMS();
         endcapABCEcalHcal[bin]->computeA(aEHe);
         endcapABCEcalHcal[bin]->computeBC();
      }
      if(!endcapABCEcal[bin]->isEmptyInFitRange())
      {
         endcapABCEcal[bin]->computeETrueAverage();
         endcapABCEcal[bin]->computeETrueRMS();
         endcapABCEcal[bin]->computeA(aEe);
         endcapABCEcal[bin]->computeB();
      }
      if(!endcapABCHcal[bin]->isEmptyInFitRange())
      {
         endcapABCHcal[bin]->computeETrueAverage();
         endcapABCHcal[bin]->computeETrueRMS();
         endcapABCHcal[bin]->computeA(aHe);
         endcapABCHcal[bin]->computeC();
      }
      

      if(!barrelABCEcalHcal[bin]->isEmpty() && 
         barrelABCEcalHcal[bin]->getBinHighEdge() >
         barrelWithEcalHcalCalib->getETrueMax())
      {
         barrelWithEcalHcalCalib->setETrueMax(
            barrelABCEcalHcal[bin]->getBinHighEdge());
      }
      if(!barrelABCEcal[bin]->isEmpty() && 
         barrelABCEcal[bin]->getBinHighEdge() >
         barrelWithEcalCalib->getETrueMax())
      {
         barrelWithEcalCalib->setETrueMax(
            barrelABCEcal[bin]->getBinHighEdge());
      }
      if(!barrelABCHcal[bin]->isEmpty() && 
         barrelABCHcal[bin]->getBinHighEdge() >
         barrelWithHcalCalib->getETrueMax())
      {
         barrelWithHcalCalib->setETrueMax(
            barrelABCHcal[bin]->getBinHighEdge());
      }
      if(!endcapABCEcalHcal[bin]->isEmpty() && 
         endcapABCEcalHcal[bin]->getBinHighEdge() >
         endcapWithEcalHcalCalib->getETrueMax())
      {
         endcapWithEcalHcalCalib->setETrueMax(
            endcapABCEcalHcal[bin]->getBinHighEdge());
      }
      if(!endcapABCEcal[bin]->isEmpty() && 
         endcapABCEcal[bin]->getBinHighEdge() >
         endcapWithEcalCalib->getETrueMax())
      {
         endcapWithEcalCalib->setETrueMax(
            endcapABCEcal[bin]->getBinHighEdge());
      }
      if(!endcapABCHcal[bin]->isEmpty() && 
         endcapABCHcal[bin]->getBinHighEdge() >
         endcapWithHcalCalib->getETrueMax())
      {
         endcapWithHcalCalib->setETrueMax(
            endcapABCHcal[bin]->getBinHighEdge());
      }
 

      barrelWithEcalHcalCalib->addGraphPoints(barrelABCEcalHcal[bin]); 
      barrelWithEcalCalib->addGraphPoints(barrelABCEcal[bin]); 
      barrelWithHcalCalib->addGraphPoints(barrelABCHcal[bin]); 
      endcapWithEcalHcalCalib->addGraphPoints(endcapABCEcalHcal[bin]); 
      endcapWithEcalCalib->addGraphPoints(endcapABCEcal[bin]); 
      endcapWithHcalCalib->addGraphPoints(endcapABCHcal[bin]);                 


   }
   
   cout<<"Fitting a, b, c coefficients..."<<endl;
   //Initialize all the ABC graphs in the calibration objects.
   barrelWithEcalHcalCalib->initializeGraphs("abc");

   barrelWithEcalCalib->initializeGraphs("abc");

   barrelWithHcalCalib->initializeGraphs("abc");   

   endcapWithEcalHcalCalib->initializeGraphs("abc");

   endcapWithEcalCalib->initializeGraphs("abc");

   endcapWithHcalCalib->initializeGraphs("abc");   


   //Define the functions that you will fit your ABC calibration constants to.
   functionBarrelEcalHcalA = new TF1("functionBarrelEcalHcalA","[0]", 0, 1000);
   // functionBarrelEcalHcalB = new TF1("functionBarrelEcalHcalB","[0]+([1]+[2]/sqrt(x))*exp(-x/[3])-[4]*exp(-x*x/[5])", 0, 1000);
   functionBarrelEcalHcalB = new TF1("functionBarrelEcalHcalB","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))", 0, 1000);
   //functionBarrelEcalHcalC = new TF1("functionBarrelEcalHcalC","[0]+(([1]+([2]/sqrt(x)))*exp(-(x^[4]/[3])))",0,1000); //[0]+([1]+[2]/sqrt(x))*exp(-x/[3])-[4]*exp(-x*x/[5])", 0, 1000);
   functionBarrelEcalHcalC = new TF1("functionBarrelEcalHcalC","([0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5])))))",0,1000);
  
   functionEndcapEcalHcalA = new TF1("functionEndcapEcalHcalA","[0]", 0, 1000);
   //functionEndcapEcalHcalB = new TF1("functionEndcapEcalHcalB","[0]+([1]+[2]/sqrt(x))*exp(-x/[3])-[4]*exp(-x*x/[5])", 0, 1000);
   functionEndcapEcalHcalB = new TF1("functionEndcapEcalHcalB","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))", 0, 1000);
   functionEndcapEcalHcalC = new TF1("functionEndcapEcalHcalC","([0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5])))))", 0, 1000);

   functionBarrelHcalA = new TF1("functionBarrelHcalA","[0]", 0, 1000);
   functionBarrelHcalB = new TF1("functionBarrelHcalB","[0]", 0, 1000);
   // functionBarrelHcalC = new TF1("functionBarrelHcalC","[0]+([1]+[2]/sqrt(x))*exp(-x/[3])-[4]*exp(-x*x/[5])", 0, 1000);
   functionBarrelHcalC = new TF1("functionBarrelHcalC","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))", 0, 1000);
   //spandey
   //functionBarrelHcalC = new TF1("functionBarrelHcalC","1.03*([0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5])))))", 0, 1000);
  
   functionEndcapHcalA = new TF1("functionEndcapHcalA","[0]", 0, 1000);
   functionEndcapHcalB = new TF1("functionEndcapHcalB","[0]", 0, 1000);
   functionEndcapHcalC = new TF1("functionEndcapHcalC","([0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5])))))", 0, 1000); //[0]+([1]+[2]/sqrt(x))*exp(-x/[3])-[4]*exp(-x*x/[5])


   if(freezeparameters) {
   

     //Set the parameters of the function you just defined.
     functionBarrelEcalHcalA->FixParameter(0, aEH);
     // functionBarrelEcalHcalB->SetParameters(-1.38681e+01,1.49678e+01,3.45153e+00,1.04212e+00,-2.00910e-02, 9.41444e-16,-1.31641e+00,-7.07963e+00);

     //faBarrel
     // functionBarrelEcalHcalB->SetParameter(0, -13.9219);
     // functionBarrelEcalHcalB->SetParameter(1, 14.9124);
     // functionBarrelEcalHcalB->SetParameter(2, 5.38578);
     // functionBarrelEcalHcalB->SetParameter(3, 0.861981);
     // functionBarrelEcalHcalB->SetParameter(4, -0.00759275);
     // functionBarrelEcalHcalB->SetParameter(5, 0.00373563);
     // functionBarrelEcalHcalB->SetParameter(6, -1.17946);
     // functionBarrelEcalHcalB->SetParameter(7, -1.69561);



     functionBarrelEcalHcalB->FixParameter(0, -30.7141);
     functionBarrelEcalHcalB->FixParameter(1, 31.7583);
     functionBarrelEcalHcalB->FixParameter(2, 4.40594);
     functionBarrelEcalHcalB->FixParameter(3, 1.70914);
     functionBarrelEcalHcalB->FixParameter(4, 0.0613696);
     functionBarrelEcalHcalB->FixParameter(5, 0.000104857);
     functionBarrelEcalHcalB->FixParameter(6, -1.38927);
     functionBarrelEcalHcalB->FixParameter(7, -0.743082);


     // functionBarrelEcalHcalC->SetParameters(1.70114,0.404676,-3.88962,1.2109e+06,
     // 					  0.970741,0.0527482,2.60552,-0.8956);

     //fbBarrel
     functionBarrelEcalHcalC->FixParameter(0,2.253661);
     functionBarrelEcalHcalC->FixParameter(1,0.537715);
     functionBarrelEcalHcalC->FixParameter(2,-4.813746);
     functionBarrelEcalHcalC->FixParameter(3,12.109);
     functionBarrelEcalHcalC->FixParameter(4,1.805775);
     functionBarrelEcalHcalC->FixParameter(5,0.187919);
     functionBarrelEcalHcalC->FixParameter(6,-6.26234);
     functionBarrelEcalHcalC->FixParameter(7,-0.607392);


     //fcBarrel
     functionBarrelHcalC->FixParameter(0,1.5125962);
     functionBarrelHcalC->FixParameter(1,0.855057);
     functionBarrelHcalC->FixParameter(2,-6.041990);
     functionBarrelHcalC->FixParameter(3,2.08229);
     functionBarrelHcalC->FixParameter(4,0.592266);
     functionBarrelHcalC->FixParameter(5,0.0291232);
     functionBarrelHcalC->FixParameter(6,0.364802);
     functionBarrelHcalC->FixParameter(7,-1.50142);




     
     functionEndcapEcalHcalA->FixParameter(0, aEHe);
  

     //     functionEndcapEcalHcalB->SetParameters(0.930193,11.9536,-30.0337,0.76133,
     //					    0.077633,7.3809e-10,0.158734,-6.92163);

     //faEndcap
     //spandey
     functionEndcapEcalHcalB->FixParameter(0,1.17227);	
     functionEndcapEcalHcalB->FixParameter(1,13.1489);
     functionEndcapEcalHcalB->FixParameter(2,-29.1672);
     functionEndcapEcalHcalB->FixParameter(3,0.604223);
     functionEndcapEcalHcalB->FixParameter(4,0.0426363);
     functionEndcapEcalHcalB->FixParameter(5,3.30898e-15);
     functionEndcapEcalHcalB->FixParameter(6,0.165293);
     functionEndcapEcalHcalB->FixParameter(7,-7.56786);
 

     // functionEndcapEcalHcalC->SetParameters(-0.436687,2.73698,-3.1509,1.20536,
     //					  -1.39685,0.0180331,0.270058,-2.30372);


     // fbEndcap
     functionEndcapEcalHcalC->SetParameter(0,-0.974251);
     functionEndcapEcalHcalC->SetParameter(1,1.61733);
     functionEndcapEcalHcalC->SetParameter(2,0.0629183);
     functionEndcapEcalHcalC->FixParameter(3,7.78495);
     functionEndcapEcalHcalC->SetParameter(4,-0.774289);
     functionEndcapEcalHcalC->FixParameter(5,7.81399e-05);
     functionEndcapEcalHcalC->FixParameter(6,0.139116);
     functionEndcapEcalHcalC->FixParameter(7,-4.25551);
     

     functionBarrelHcalA->FixParameter(0, aH);
     functionBarrelHcalB->FixParameter(0, 0.0);

     functionEndcapHcalA->FixParameter(0, aHe);
     functionEndcapHcalB->FixParameter(0, 0.0);

     //fcEndcap 
     functionEndcapHcalC->FixParameter(0,1.01863);
     functionEndcapHcalC->FixParameter(1,1.29787);
     functionEndcapHcalC->FixParameter(2,-3.97293);
     functionEndcapHcalC->FixParameter(3,21.7805);
     functionEndcapHcalC->FixParameter(4,0.810195);
     functionEndcapHcalC->FixParameter(5,0.234134);
     functionEndcapHcalC->FixParameter(6,1.42226);
     functionEndcapHcalC->FixParameter(7,-0.0997326);




   
   }
   else {

        //Set the parameters of the function you just defined.
   functionBarrelEcalHcalA->FixParameter(0, aEH);
   // functionBarrelEcalHcalB->SetParameters(1.32991e+00,1.52538e+02,-7.89866e+02,2.56441e-01,
   // 					  -1.13519e+00,2.71646e+00,1.86857e-01,4.68538e-01);
   functionBarrelEcalHcalB->SetParameters(-1.38681e+01,1.49678e+01,3.45153e+00,1.04212e+00,-2.00910e-02,
					  9.41444e-16,-1.31641e+00,-7.07963e+00);
   //** begin modified by seema
   functionBarrelEcalHcalB->FixParameter(0, -13.9219);
   functionBarrelEcalHcalB->FixParameter(1, 14.9124);
   functionBarrelEcalHcalB->FixParameter(2, 5.38578);
   functionBarrelEcalHcalB->FixParameter(3, 0.861981);
   functionBarrelEcalHcalB->FixParameter(4, -0.00759275);
   functionBarrelEcalHcalB->FixParameter(5, 3.73563e-23);
   functionBarrelEcalHcalB->FixParameter(6, -1.17946);
   functionBarrelEcalHcalB->FixParameter(7, -13.3644);
   //** end modified by seema
   
   functionBarrelEcalHcalC->SetParameters(1.70114,0.404676,-3.88962,1.2109e+06,
					  0.970741,0.0527482,2.60552,-0.8956);

   //** begin modified by seema
   //functionBarrelEcalHcalC->FixParameter(0, 1.70114 );
   //functionBarrelEcalHcalC->FixParameter(1, 0.404676 );
   //functionBarrelEcalHcalC->FixParameter(2, -3.88962);
   functionBarrelEcalHcalC->FixParameter(3, 1.2109e+06); 
   //functionBarrelEcalHcalC->FixParameter(4, 0.970741 );
   //functionBarrelEcalHcalC->FixParameter(5, 0.0527482);
   //functionBarrelEcalHcalC->FixParameter(6, 2.60552);
   //functionBarrelEcalHcalC->FixParameter(7, -0.8956);	
   //** end modified by seema  
   functionEndcapEcalHcalA->FixParameter(0, aEHe);
  
   functionEndcapEcalHcalB->SetParameters(0.930193,11.9536,-30.0337,0.76133,
					  0.0776373,7.3809e-10,0.158734,-6.92163);
   
   //** begin modified by seema  
   //functionEndcapEcalHcalB->FixParameter(0, 0.930193);
   ///functionEndcapEcalHcalB->FixParameter(0, 0.962468);
   ///functionEndcapEcalHcalB->FixParameter(1, 11.9536);
   //functionEndcapEcalHcalB->FixParameter(2, -30.0337);
   //functionEndcapEcalHcalB->FixParameter(2, -28.6722);
   ///functionEndcapEcalHcalB->FixParameter(2, -27.7088);
   //functionEndcapEcalHcalB->FixParameter(3, 0.76133);
   //functionEndcapEcalHcalB->FixParameter(3, 0.757575);   
   //functionEndcapEcalHcalB->FixParameter(3, 0.758274);   
   ///functionEndcapEcalHcalB->FixParameter(3, 0.755474);
   //functionEndcapEcalHcalB->FixParameter(4, 0.0776373);
   ///functionEndcapEcalHcalB->FixParameter(4, 0.0791012);
   //functionEndcapEcalHcalB->FixParameter(5, 7.3809e-10);
   //functionEndcapEcalHcalB->FixParameter(5, 2.6901e-11);
   ///functionEndcapEcalHcalB->FixParameter(5, 2.6901e-3);
   ///functionEndcapEcalHcalB->FixParameter(6, 0.158734);
   //functionEndcapEcalHcalB->FixParameter(7, -6.92163);
   ///functionEndcapEcalHcalB->FixParameter(7, -0.92163);
   //** end modified by seema  

   functionEndcapEcalHcalC->SetParameters(-0.436687,2.73698,-3.1509,1.20536,
					  -1.39685,0.0180331,0.270058,-2.30372);
   //** begin modified by seema  
   //functionEndcapEcalHcalC->FixParameter(0, -0.436687 );
   ///functionEndcapEcalHcalC->FixParameter(0, -0.43671);
   //functionEndcapEcalHcalC->FixParameter(1, 2.73698);
   ///functionEndcapEcalHcalC->FixParameter(1, 2.90096);
   //functionEndcapEcalHcalC->FixParameter(2, -3.1509);
   ///functionEndcapEcalHcalC->FixParameter(2, -5.10099);
   //functionEndcapEcalHcalC->FixParameter(3, 1.20536);
   ///functionEndcapEcalHcalC->FixParameter(3, 1.20771);
   //functionEndcapEcalHcalC->FixParameter(4, -1.39685);
   ///functionEndcapEcalHcalC->FixParameter(4, -1.30656);
   //functionEndcapEcalHcalC->FixParameter(5, 0.0180331);
   ///functionEndcapEcalHcalC->FixParameter(5, 0.0189607);
   //functionEndcapEcalHcalC->FixParameter(6, 0.270058);
   ///functionEndcapEcalHcalC->FixParameter(6, 0.270027);
   ///functionEndcapEcalHcalC->FixParameter(7, -2.30372);			  
   //** end modified by seema  
  

   functionBarrelHcalA->FixParameter(0, aH);
   functionBarrelHcalB->FixParameter(0, 0.0);
   functionBarrelHcalC->SetParameters(1.58827e+00,4.06865e-01,-3.69939e+00,1.28926e+03,
					7.13400e-01,2.21925e-02,1.47842e+00,-1.22041e+00);
      
   functionEndcapHcalA->FixParameter(0, aHe);
   functionEndcapHcalB->FixParameter(0, 0.0);
   
   functionEndcapHcalC->SetParameters(1.13795,1.21698,-3.81192,115.409,
				      0.673456,0.217077,1.95596,-0.252215);

   //CHANGED 30 Apr

   // functionEndcapHcalC->FixParameter(0, 1.13795);
   // functionEndcapHcalC->FixParameter(1, 1.21698);
   // functionEndcapHcalC->FixParameter(2, -3.81192);
   // //functionEndcapHcalC->FixParameter(3, 115.409);
   // functionEndcapHcalC->FixParameter(3, 60.0406);
   // functionEndcapHcalC->FixParameter(4, 0.673456);
   // functionEndcapHcalC->FixParameter(5, 0.217077);
   // functionEndcapHcalC->FixParameter(6, 1.95596);
   // functionEndcapHcalC->FixParameter(7, -0.252215);
   }

   barrelWithEcalHcalCalib->fitAsToFunction(functionBarrelEcalHcalA);
   //Printing parameters:
   barrelWithEcalHcalCalib->fitBsToFunction(functionBarrelEcalHcalB);
   barrelWithEcalHcalCalib->fitBsToFunction();



   barrelWithEcalHcalCalib->fitCsToFunction(functionBarrelEcalHcalC);


   barrelWithEcalHcalCalib->fitCsToFunction();


   barrelWithEcalHcalCalib->fitCsToFunction();

   endcapWithEcalHcalCalib->fitAsToFunction(functionEndcapEcalHcalA);

   // cout<<"********************************************"<<endl;
   // cout<<"Fit Parameters, functionEndcapEcalHcalB, First"<<endl;
   // for ( unsigned i = 0; i < 10; ++i ) {
   //   double barrelEcalHcalC_spandey = functionEndcapEcalHcalB->GetParameter(i);
   //   if ( barrelEcalHcalC_spandey != 0. )
   //     cout<<"  functionEndcapEcalHcalB,Parameter("<<i<<","<<barrelEcalHcalC_spandey<<");"<<endl;
   // }

   endcapWithEcalHcalCalib->fitBsToFunction(functionEndcapEcalHcalB);

   // cout<<"********************************************"<<endl;
   // cout<<"Fit Parameters, functionEndcapEcalHcalB, First"<<endl;
   // for ( unsigned i = 0; i < 10; ++i ) {
   //   double barrelEcalHcalC_spandey = functionEndcapEcalHcalB->GetParameter(i);
   //   if ( barrelEcalHcalC_spandey != 0. )
   //     cout<<"  functionEndcapEcalHcalB,Parameter("<<i<<","<<barrelEcalHcalC_spandey<<");"<<endl;
   // }

   endcapWithEcalHcalCalib->fitBsToFunction();

   // cout<<"********************************************"<<endl;
   // cout<<"Fit Parameters, functionEndcapEcalHcalB, First"<<endl;
   // for ( unsigned i = 0; i < 10; ++i ) {
   //   double barrelEcalHcalC_spandey = functionEndcapEcalHcalB->GetParameter(i);
   //   if ( barrelEcalHcalC_spandey != 0. )
   //     cout<<"  functionEndcapEcalHcalB,Parameter("<<i<<","<<barrelEcalHcalC_spandey<<");"<<endl;
   // }

   //   cout<<"********************************************"<<endl;
   endcapWithEcalHcalCalib->fitCsToFunction(functionEndcapEcalHcalC);
   endcapWithEcalHcalCalib->fitCsToFunction();
   endcapWithEcalHcalCalib->fitCsToFunction();
   //   cout<<"Fit check 1\n";
   barrelWithHcalCalib->fitAsToFunction(functionBarrelHcalA);
   barrelWithHcalCalib->fitBsToFunction(functionBarrelHcalB);
   //   cout<<"Fit check 1.a\n";
   barrelWithHcalCalib->fitBsToFunction();
   //   cout<<"Fit check 1.b\n";
   barrelWithHcalCalib->fitCsToFunction(functionBarrelHcalC);
   barrelWithHcalCalib->fitCsToFunction();
   
   barrelWithHcalCalib->fitCsToFunction();
   //   cout<<"Fit check 2\n";
   endcapWithHcalCalib->fitAsToFunction(functionEndcapHcalA);
   endcapWithHcalCalib->fitBsToFunction(functionEndcapHcalB);
   endcapWithHcalCalib->fitBsToFunction();

   //   cout<<"1 FITTING ENDCAP HCAL C#######"<<endl;
   endcapWithHcalCalib->fitCsToFunction(functionEndcapHcalC);
   //   cout<<"2 FITTING ENDCAP HCAL C#######"<<endl;
   endcapWithHcalCalib->fitCsToFunction();
   //   cout<<"3 FITTING ENDCAP HCAL C#######"<<endl;
   endcapWithHcalCalib->fitCsToFunction();
   //   cout<<"Fit check 3\n";



   //Fill all the TH2's that can be put into drawGausFit in order to produce 
   //response and resolution plots.
   cout<<"Making response and resolution plots..."<<endl;
   for(unsigned entry = 0; entry < ETrueEnergies.size(); ++entry)
   {

      etrue = ETrueEnergies[entry];
      ecal = ecalEnergies[entry];
      hcal = hcalEnergies[entry];
      eta = abs(etas[entry]);
      double phi = phis[entry];
      if((ecal + hcal) < 0.5) continue;
      if( etrue < 1.0) continue;
      if( hcal == 0) continue;
      //if( etrue/cosh(eta) > 10.0) continue;
      // if( ecal > 0) continue;
      //if(fabs(eta) < 2.5) continue; // delete me
 if(fabs(eta) < 1.5) //b, c fit range //shubham Mar 27
         {
            rawBarrel->Fill(etrue, (ecal + hcal - etrue)/etrue);
            
            if(ecal > 0 && hcal > 0)
            {
               correctedE = barrelWithEcalHcalCalib->
                  getCalibratedEnergy(etrue, ecal, hcal);


	   if(drawpT) {
	     etrue = etrue/cosh(eta);
	     // correctedEta = correctedEta/cosh(eta);
	     ecal = ecal/cosh(eta);
	     hcal = hcal/cosh(eta);
	     correctedE = correctedE/cosh(eta);
	   }

	   //rawBarrelEcalHcal->Fill(etrue, (ecal + hcal  - etrue)/etrue );
  	       rawBarrelEcalHcal->Fill(etrue, (ecal + hcal  - etrue)/etrue );
               corrBarrel->Fill(etrue, (correctedE - etrue)/etrue);
               corrBarrelEcalHcal->Fill(etrue, (correctedE - etrue)/etrue);
               
	       // hcorrEtaDependence->Fill(eta, (correctedE - etrue)/etrue);

               // rawEtaDependence->Fill(eta, (ecal + hcal - etrue)/etrue);
               // corrEtaDependence->Fill(eta, (correctedEta - etrue)/etrue);

	       // if(entry<5000) 
	       // 	 cout<<entry<<"   "<<eta<<"   "<<etrue<<"   "<<ecal+hcal<<"   "<<correctedE<<"   "<<correctedEta<<endl;

	       
	       // h_response_vs_phi_barrel_EH->Fill(phi, (correctedE - etrue)/etrue); //shuham
	     

            }
            else if(ecal == 0 && hcal > 0)
	      {
		correctedE = barrelWithHcalCalib->
                  getCalibratedEnergy(etrue, ecal, hcal);


	   if(drawpT) {
	     etrue = etrue/cosh(eta);
	     // correctedEta = correctedEta/cosh(eta);
	     ecal = ecal/cosh(eta);
	     hcal = hcal/cosh(eta);
	     correctedE = correctedE/cosh(eta);
	   }

		rawBarrelHcal->Fill(etrue, ( ecal + hcal - etrue)/etrue );// (etrue-3.0)/(ecal+hcal) );//, 11936/(3917*sigmas[entry]*sigmas[entry]) );// ( ecal + hcal - etrue)/etrue );

		corrBarrel->Fill(etrue, (correctedE- etrue )/etrue);// 
		corrBarrelHcal->Fill(etrue, (correctedE- etrue )/etrue);

		// h_response_vs_phi_barrel_H->Fill(phi, (correctedE - etrue)/etrue); //shuham

	    
	      }
         }
   
      //if(fabs(eta) < 2.5 && fabs(eta) > 1.55) //WITHIN TRACKER alpha beta fit range for endcap 
      //if(fabs(eta) < 3.0 && fabs(eta) > 1.55) //FULL EndCap alpha beta fit range for endcap   //shubham
	  //if(fabs(eta) < 3.0 && fabs(eta) > 2.5) //OUTSIDE TRACKER alpha beta fit range for endcap   //shubham
      if(fabs(eta) < _etaMax_ && fabs(eta) > _etaMin_) 
      {
	//if (fabs(eta) > 2.7) cout<<"yolo "<<fabs(eta)<<endl;
         raw->Fill(etrue, (ecal + hcal - etrue)/etrue);

	 ////////////////////////
	 // RAW Proxy
	 double etrue_proxy;
	 if (fabs(eta) > 2.5)
	   etrue_proxy = etrue;//ecal + hcal;
	 else
	   etrue_proxy = etrue;

         if(ecal > 0)
         {
           
	    correctedE = endcapWithEcalHcalCalib->
	      getCalibratedEnergy(etrue_proxy, ecal, hcal);

	   if(drawpT) {
	     etrue = etrue/cosh(eta);
	     etrue_proxy = etrue_proxy/cosh(eta);
	     ecal = ecal/cosh(eta);
	     hcal = hcal/cosh(eta);
	     correctedE = correctedE/cosh(eta);
	   }

          
         }
         else
         {
          
	    correctedE = endcapWithHcalCalib->
	      getCalibratedEnergy(etrue_proxy, ecal, hcal);


	   if(drawpT) {
	     etrue = etrue/cosh(eta);
	     etrue_proxy = etrue_proxy/cosh(eta);
	     ecal = ecal/cosh(eta);
	     hcal = hcal/cosh(eta);
	     correctedE = correctedE/cosh(eta);
	   }

              
         }
         //if(fabs(eta) < 2.2) //b, c fi trange
	 if(fabs(eta) < 3.0) //b, c fi trange   //shubham
         {
            rawEndcap->Fill(etrue, (ecal + hcal - etrue)/etrue);
            
            if(ecal > 0)
            {

	   
               correctedE = endcapWithEcalHcalCalib->
                  getCalibratedEnergy(etrue_proxy, ecal, hcal);


	   if(drawpT) {
	     etrue = etrue/cosh(eta);
	     etrue_proxy = etrue_proxy/cosh(eta);
	     ecal = ecal/cosh(eta);
	     hcal = hcal/cosh(eta);
	     correctedE = correctedE/cosh(eta);
	   }

	       rawEndcapEcalHcal->Fill(etrue, (ecal + hcal - etrue)/etrue);
	       corrEndcap->Fill(etrue, (correctedE - etrue)/etrue);
	       corrEndcapEcalHcal->Fill(etrue, (correctedE - etrue)/etrue);

	   

	       // rawEtaDependence->Fill(eta, (ecal + hcal - etrue)/etrue);
	       // corrEtaDependence->Fill(eta, (correctedEta - etrue)/etrue);
	       // hcorrEtaDependence->Fill(eta, (correctedE - etrue)/etrue);

	       //cout<<"yolo, eta:"<<eta<<endl;
	    
            }
            else 
            {
               correctedE = endcapWithHcalCalib->
                  getCalibratedEnergy(etrue, ecal, hcal);


	   if(drawpT) {
	     etrue = etrue/cosh(eta);
	     etrue_proxy = etrue_proxy/cosh(eta);
	     ecal = ecal/cosh(eta);
	     hcal = hcal/cosh(eta);
	     correctedE = correctedE/cosh(eta);
	   }
               
               rawEndcapHcal->Fill(etrue, (ecal + hcal - etrue)/etrue);
               corrEndcap->Fill(etrue, (correctedE - etrue)/etrue);
               corrEndcapHcal->Fill(etrue, (correctedE - etrue)/etrue);

	    

             
            }
         }
	 else  {   //shubham
	   if(ecal > 0) 
	     corrEndcapEcalHcal->Fill(etrue, (correctedE - etrue)/etrue);
	 }
      }
   }
   	 
   cout<<" Now Summary "<<endl;
   drawGausFit(rawBarrelEcalHcal,responseRaw,resolutionRaw);

TH2F* rawecalvsrawhcal = new TH2F("rawecalvsrawhcal","uncorrected Ecal energy vs uncorrected Hcal energy",500,-2,sampleRangeHigh,500,-2,sampleRangeHigh);
TH1F* uncorr_ecal = new TH1F("uncorr_ecal","Uncorrected Ecal energy ",500,-2,sampleRangeHigh);
TH1F* uncorr_hcal = new TH1F("uncorr_hcal","Uncorrected Hcal energy ",500,-2,sampleRangeHigh);
TH1F* response_EH = new TH1F("response_EH", "1D Response distribution",125,-1.5,2);
TH1F* response_H = new TH1F("response_H", "1D Response distribution",125,-1.5,2);
TH1F* response_EH_0_1 = new TH1F("response_EH_0_1", "1D Response distribution",125,-1.5,2);





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
     response_EH->SetLineColor(4);
     response_EH_0_1->SetLineColor(2);
     response_H->SetLineColor(1);
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

     response_EH->Draw("hist");
     response_EH_0_1->Draw("hist sames");
     response_H->Draw("hist sames");
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
