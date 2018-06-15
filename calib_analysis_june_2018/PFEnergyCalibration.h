#ifndef RecoParticleFlow_PFClusterTools_PFEnergyCalibration_h
#define RecoParticleFlow_PFClusterTools_PFEnergyCalibration_h 

#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "CondFormats/PhysicsToolsObjects/interface/PerformancePayloadFromTFormula.h"

#include "CondFormats/ESObjects/interface/ESEEIntercalibConstants.h"

class TF1;

// -*- C++ -*-
//
// Package:    PFClusterTools
// Class:      PFEnergyCalibration
// 
/**\class

 Description: An auxiliary class of the Particle-Flow algorithm,
              for the Calibration of Energy Deposits in ECAL and HCAL

 Implementation:
     Original Implementation of Calibration functions in PFAlgo/PFBlock 
     by Colin Bernet;
     Code moved into separate Class by Christian Veelken 

 Modification 
     To include energy-dependent and angular-dependent calibration
     By Patrick Janot
 
*/
//
// Original Author:  Christian Veelken
//         Created:  Tue Aug  8 16:26:18 CDT 2006
//
//

#include <iostream>

//#include "FWCore/ParameterSet/interface/ParameterSet.h"

class PFEnergyCalibration 
{
 public:
  PFEnergyCalibration(); 

  ~PFEnergyCalibration();

  // ecal calibration for photons
  double energyEm(const reco::PFCluster& clusterEcal,
		  std::vector<double> &EclustersPS1,
		  std::vector<double> &EclustersPS2,
		  bool crackCorrection = true) const;
  double energyEm(const reco::PFCluster& clusterEcal,
		  double ePS1,  double ePS2,
		  bool crackCorrection = true) const;

  double energyEm(const reco::PFCluster& clusterEcal,
		  std::vector<double> &EclustersPS1,
		  std::vector<double> &EclustersPS2,
		  double &ps1,double&ps2,
		  bool crackCorrection=true) const;
  double energyEm(const reco::PFCluster& clusterEcal,
		  double ePS1, double ePS2,
		  double &ps1,double&ps2,
		  bool crackCorrection=true) const;

  // ECAL+HCAL (abc) calibration, with E and eta dependent coefficients, for hadrons
  void energyEmHad(double t, double& e, double&h, double eta, double phi) const;
  
  // Initialize default E- and eta-dependent coefficient functional form
  void initializeCalibrationFunctions();

  // Set the run-dependent calibration functions from the global tag
  void setCalibrationFunctions(const PerformancePayloadFromTFormula *thePFCal) {
    pfCalibrations = thePFCal;
  }

  void initAlphaGamma_ESplanes_fromDB(const ESEEIntercalibConstants* esEEInterCalib){
    esEEInterCalib_ = esEEInterCalib;
  }


  friend std::ostream& operator<<(std::ostream& out, 
				  const PFEnergyCalibration& calib);

 protected:

  // Calibration functions from global tag
  const PerformancePayloadFromTFormula *pfCalibrations;
  const ESEEIntercalibConstants* esEEInterCalib_;
  
  // eta 0 -> 1
  std::unique_ptr<TF1> fa_0;
  std::unique_ptr<TF1> fb_0; 
  std::unique_ptr<TF1> fc_0; 
 
  // eta 1 -> 1.3
  std::unique_ptr<TF1> fa_1;
  std::unique_ptr<TF1> fb_1; 
  std::unique_ptr<TF1> fc_1;

  // eta 1.3 -> 1.6
  std::unique_ptr<TF1> fa_2;
  std::unique_ptr<TF1> fb_2; 
  std::unique_ptr<TF1> fc_2;
 
  // eta 1.6 -> 2
  std::unique_ptr<TF1> fa_3;
  std::unique_ptr<TF1> fb_3; 
  std::unique_ptr<TF1> fc_3;

  // eta 2 -> 2.5
  std::unique_ptr<TF1> fa_4;
  std::unique_ptr<TF1> fb_4; 
  std::unique_ptr<TF1> fc_4;

  // eta 2.5 -> 3
  std::unique_ptr<TF1> fa_5;
  std::unique_ptr<TF1> fb_5; 
  std::unique_ptr<TF1> fc_5;



 private:
  
  double minimum(double a,double b) const;
  double dCrackPhi(double phi, double eta) const;
  double CorrPhi(double phi, double eta) const;
  double CorrEta(double eta) const;
  double CorrBarrel(double E, double eta) const;
  double Alpha(double eta) const;
  double Beta(double E, double eta) const;
  double Gamma(double etaEcal) const;
  double EcorrBarrel(double E, double eta, double phi, bool crackCorrection=true) const;
  double EcorrZoneBeforePS(double E, double eta) const;
  double EcorrPS(double eEcal,double ePS1,double ePS2,double etaEcal) const;
  double EcorrPS(double eEcal,double ePS1,double ePS2,double etaEcal,double&, double&) const;
  double EcorrPS_ePSNil(double eEcal,double eta) const;
  double EcorrZoneAfterPS(double E, double eta) const;
  double Ecorr(double eEcal,double ePS1,double ePS2,double eta,double phi,bool crackCorrection=true) const;
  double Ecorr(double eEcal,double ePS1,double ePS2,double eta,double phi,double&,double&,bool crackCorrection=true) const;

  // The calibration functions
  double a_0(double x) const;
  double b_0(double x) const;
  double c_0(double x) const;
  double a_1(double x) const;
  double b_1(double x) const;
  double c_1(double x) const;
  double a_2(double x) const;
  double b_2(double x) const;
  double c_2(double x) const;
  double a_3(double x) const;
  double b_3(double x) const;
  double c_3(double x) const;
  double a_4(double x) const;
  double b_4(double x) const;
  double c_4(double x) const;
  double a_5(double x) const;
  double b_5(double x) const;
  double c_5(double x) const;
 
  // Threshold correction (offset)
  double threshE, threshH;

};

#endif


