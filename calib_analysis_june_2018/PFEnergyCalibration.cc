#include "RecoParticleFlow/PFClusterTools/interface/PFEnergyCalibration.h"
#include "CondFormats/PhysicsToolsObjects/interface/PerformancePayloadFromTFormula.h"

#include "CondFormats/ESObjects/interface/ESEEIntercalibConstants.h"

#include <TMath.h>
#include <cmath>
#include <vector>
#include <TF1.h>
#include <map>
#include <algorithm>
#include <numeric>

using namespace std;

//float eta [6] = {0,1,1.3,1.6,2.0,2.5,3.0};
PFEnergyCalibration::PFEnergyCalibration() : pfCalibrations(nullptr), esEEInterCalib_(nullptr)
{
  initializeCalibrationFunctions();
}

PFEnergyCalibration::~PFEnergyCalibration() 
{

}

void
PFEnergyCalibration::initializeCalibrationFunctions() {

  threshE = 3.5;
  threshH = 2.5;

  //calibChrisClean.C calibration parameters bhumika march 2018
  //change in code to make it eta independent --->  fa_Barrel, fa_Endcap etc to fa_0, fa_1, fa_2 etc where 0,1,2,3,4,5 are finer eta regions. 
  fa_0 = std::make_unique<TF1>("fa_0","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",1.,1000.);
  fa_0->SetParameter(0,-13.9219);
  fa_0->SetParameter(1,14.9124);
  fa_0->SetParameter(2,5.38578);
  fa_0->SetParameter(3,0.861981);
  fa_0->SetParameter(4,-0.00759275);
  fa_0->SetParameter(5,0.00373563);
  fa_0->SetParameter(6,-1.17946);
  fa_0->SetParameter(7,-1.69561);
  fb_0 = std::make_unique<TF1>("fb_0","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",1.,1000.);
  fb_0->SetParameter(0,2.25366);
  fb_0->SetParameter(1,0.537715);
  fb_0->SetParameter(2,-4.81375);
  fb_0->SetParameter(3,12.109);
  fb_0->SetParameter(4,1.80577);
  fb_0->SetParameter(5,0.187919);
  fb_0->SetParameter(6,-6.26234);
  fb_0->SetParameter(7,-0.607392);
  fc_0 = std::make_unique<TF1>("fc_0","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",1.,1000.);
  fc_0->SetParameter(1,0.855057);
  fc_0->SetParameter(2,-6.04199);
  fc_0->SetParameter(3,2.08229);
  fc_0->SetParameter(4,0.592266);
  fc_0->SetParameter(5,0.0291232);
  fc_0->SetParameter(6,0.364802);
  fc_0->SetParameter(7,-1.50142);
 
  fa_1 = std::make_unique<TF1>("fa_1","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",1.,1000.);
  fa_1->SetParameter(0,-13.9219);
  fa_1->SetParameter(1,14.9124);
  fa_1->SetParameter(2,5.38578);
  fa_1->SetParameter(3,0.861981);
  fa_1->SetParameter(4,-0.00759275);
  fa_1->SetParameter(5,0.00373563);
  fa_1->SetParameter(6,-1.17946);
  fa_1->SetParameter(7,-1.69561);
  fb_1 = std::make_unique<TF1>("fb_1","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",1.,1000.);
  fb_1->SetParameter(0,2.25366);
  fb_1->SetParameter(1,0.537715);
  fb_1->SetParameter(2,-4.81375);
  fb_1->SetParameter(3,12.109);
  fb_1->SetParameter(4,1.80577);
  fb_1->SetParameter(5,0.187919);
  fb_1->SetParameter(6,-6.26234);
  fb_1->SetParameter(7,-0.607392);
  fc_1 = std::make_unique<TF1>("fc_1","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",1.,1000.);
  fc_1->SetParameter(1,0.855057);
  fc_1->SetParameter(2,-6.04199);
  fc_1->SetParameter(3,2.08229);
  fc_1->SetParameter(4,0.592266);
  fc_1->SetParameter(5,0.0291232);
  fc_1->SetParameter(6,0.364802);
  fc_1->SetParameter(7,-1.50142);

  fa_2 = std::make_unique<TF1>("fa_1","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",1.,1000.);
  fa_2->SetParameter(0,-13.9219);
  fa_2->SetParameter(1,14.9124);
  fa_2->SetParameter(2,5.38578);
  fa_2->SetParameter(3,0.861981);
  fa_2->SetParameter(4,-0.00759275);
  fa_2->SetParameter(5,0.00373563);
  fa_2->SetParameter(6,-1.17946);
  fa_2->SetParameter(7,-1.69561);
  fb_2 = std::make_unique<TF1>("fb_1","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",1.,1000.);
  fb_2->SetParameter(0,2.25366);
  fb_2->SetParameter(1,0.537715);
  fb_2->SetParameter(2,-4.81375);
  fb_2->SetParameter(3,12.109);
  fb_2->SetParameter(4,1.80577);
  fb_2->SetParameter(5,0.187919);
  fb_2->SetParameter(6,-6.26234);
  fb_2->SetParameter(7,-0.607392);
  fc_2 = std::make_unique<TF1>("fc_1","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",1.,1000.);
  fc_2->SetParameter(1,0.855057);
  fc_2->SetParameter(2,-6.04199);
  fc_2->SetParameter(3,2.08229);
  fc_2->SetParameter(4,0.592266);
  fc_2->SetParameter(5,0.0291232);
  fc_2->SetParameter(6,0.364802);
  fc_2->SetParameter(7,-1.50142);
 
  fa_3 = std::make_unique<TF1>("fa_3","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",1.,1000.);
  fa_3->SetParameter(0,-13.9219);
  fa_3->SetParameter(1,14.9124);
  fa_3->SetParameter(2,5.38578);
  fa_3->SetParameter(3,0.861981);
  fa_3->SetParameter(4,-0.00759275);
  fa_3->SetParameter(5,0.00373563);
  fa_3->SetParameter(6,-1.17946);
  fa_3->SetParameter(7,-1.69561);
  fb_3 = std::make_unique<TF1>("fb_3","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",3.,1000.);
  fb_3->SetParameter(0,2.25366);
  fb_3->SetParameter(1,0.537715);
  fb_3->SetParameter(2,-4.81375);
  fb_3->SetParameter(3,12.109);
  fb_3->SetParameter(4,1.80577);
  fb_3->SetParameter(5,0.187919);
  fb_3->SetParameter(6,-6.26234);
  fb_3->SetParameter(7,-0.607392);
  fc_3 = std::make_unique<TF1>("fc_3","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",1.,1000.);
  fc_3->SetParameter(1,0.855057);
  fc_3->SetParameter(2,-6.04199);
  fc_3->SetParameter(3,2.08229);
  fc_3->SetParameter(4,0.592266);
  fc_3->SetParameter(5,0.0291232);
  fc_3->SetParameter(6,0.364802);
  fc_3->SetParameter(7,-1.50142);

  fa_4 = std::make_unique<TF1>("fa_4","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",1.,1000.);
  fa_4->SetParameter(0,-13.9219);
  fa_4->SetParameter(1,14.9124);
  fa_4->SetParameter(2,5.38578);
  fa_4->SetParameter(3,0.861981);
  fa_4->SetParameter(4,-0.00759275);
  fa_4->SetParameter(5,0.00373563);
  fa_4->SetParameter(6,-1.17946);
  fa_4->SetParameter(7,-1.69561);
  fb_4 = std::make_unique<TF1>("fb_4","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",1.,1000.);
  fb_4->SetParameter(0,2.25366);
  fb_4->SetParameter(1,0.537715);
  fb_4->SetParameter(2,-4.81375);
  fb_4->SetParameter(3,12.109);
  fb_4->SetParameter(4,1.80577);
  fb_4->SetParameter(5,0.187919);
  fb_4->SetParameter(6,-6.26234);
  fb_4->SetParameter(7,-0.607392);
  fc_4 = std::make_unique<TF1>("fc_4","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",1.,1000.);
  fc_4->SetParameter(1,0.855057);
  fc_4->SetParameter(2,-6.04199);
  fc_4->SetParameter(3,2.08229);
  fc_4->SetParameter(4,0.592266);
  fc_4->SetParameter(5,0.0291232);
  fc_4->SetParameter(6,0.364802);
  fc_4->SetParameter(7,-1.50142);

  fa_5 = std::make_unique<TF1>("fa_5","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",1.,1000.);
  fa_5->SetParameter(0,-13.9219);
  fa_5->SetParameter(1,14.9124);
  fa_5->SetParameter(2,5.38578);
  fa_5->SetParameter(3,0.861981);
  fa_5->SetParameter(4,-0.00759275);
  fa_5->SetParameter(5,0.00373563);
  fa_5->SetParameter(6,-1.17946);
  fa_5->SetParameter(7,-1.69561);
  fb_5 = std::make_unique<TF1>("fb_5","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",1.,1000.);
  fb_5->SetParameter(0,2.25366);
  fb_5->SetParameter(1,0.537715);
  fb_5->SetParameter(2,-4.81375);
  fb_5->SetParameter(3,12.109);
  fb_5->SetParameter(4,1.80577);
  fb_5->SetParameter(5,0.187919);
  fb_5->SetParameter(6,-6.26234);
  fb_5->SetParameter(7,-0.607392);
  fc_5 = std::make_unique<TF1>("fc_5","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",1.,1000.);
  fc_5->SetParameter(1,0.855057);
  fc_5->SetParameter(2,-6.04199);
  fc_5->SetParameter(3,2.08229);
  fc_5->SetParameter(4,0.592266);
  fc_5->SetParameter(5,0.0291232);
  fc_5->SetParameter(6,0.364802);
  fc_5->SetParameter(7,-1.50142);

}

void 
PFEnergyCalibration::energyEmHad(double t, double& e, double&h, double eta, double phi) const { 

  // Use calorimetric energy as true energy for neutral particles
  double tt = t;
  double ee = e;
  double hh = h;
  double a = 1.;
  double b = 1.;
  // double etaCorrE = 1.;
  // double etaCorrH = 1.;
  auto absEta = std::abs(eta);
  t = min(999.9,max(tt,e+h));
  if ( t < 1. ) return;

  // FOR eta region 0 to 1
  if ( absEta <= 1.0 ) { 
    // The energy correction
    a = e>0. ? a_0(t) : 1.;
    b = e>0. ? b_0(t) : c_0(t);
    double thresh = e > 0. ? threshE : threshH;

    // Protection against negative calibration
    if ( a < -0.25 || b < -0.25 ) { 
      a = 1.;
      b = 1.;
      thresh = 0.;
    }

    // The new estimate of the true energy
    t = min(999.9,max(tt, thresh+a*e+b*h));

    // The angular correction 
    // if ( e > 0. && thresh > 0. ) {
    //   etaCorrE = 1.0 + aEtaBarrelEH(t) + 1.3*bEtaBarrelEH(t)*absEta*absEta;
    //   etaCorrH = 1.0;
    // } else {
    //   etaCorrE = 1.0 + aEtaBarrelH(t) + 1.3*bEtaBarrelH(t)*absEta*absEta; 
    //   etaCorrH = 1.0 + aEtaBarrelH(t) + bEtaBarrelH(t)*absEta*absEta;
    // }
    if ( e > 0. && thresh > 0. ) 
      e = h > 0. ? threshE-threshH +  a * e : threshE +  a * e;
    if ( h > 0. && thresh > 0. ) {
      h = threshH + b * h;
    }

     
  }

// FOR eta region 1 to 1.3
  if ( absEta > 1 && absEta <= 1.3 ) { 
    // The energy correction
    a = e>0. ? a_1(t) : 1.;
    b = e>0. ? b_1(t) : c_1(t);
    double thresh = e > 0. ? threshE : threshH;

    // Protection against negative calibration
    if ( a < -0.25 || b < -0.25 ) { 
      a = 1.;
      b = 1.;
      thresh = 0.;
    }

    // The new estimate of the true energy
    t = min(999.9,max(tt, thresh+a*e+b*h));

    // The angular correction 
    // if ( e > 0. && thresh > 0. ) {
    //   etaCorrE = 1.0 + aEtaBarrelEH(t) + 1.3*bEtaBarrelEH(t)*absEta*absEta;
    //   etaCorrH = 1.0;
    // } else {
    //   etaCorrE = 1.0 + aEtaBarrelH(t) + 1.3*bEtaBarrelH(t)*absEta*absEta; 
    //   etaCorrH = 1.0 + aEtaBarrelH(t) + bEtaBarrelH(t)*absEta*absEta;
    // }
    if ( e > 0. && thresh > 0. ) 
      e = h > 0. ? threshE-threshH +  a * e : threshE +  a * e;
    if ( h > 0. && thresh > 0. ) {
      h = threshH + b * h;
    }

  }
// FOR eta region 1.3 to 1.6
  if ( absEta > 1.3 && absEta <= 1.6 ) { 
    // The energy correction
    a = e>0. ? a_2(t) : 1.;
    b = e>0. ? b_2(t) : c_2(t);
    double thresh = e > 0. ? threshE : threshH;

    // Protection against negative calibration
    if ( a < -0.25 || b < -0.25 ) { 
      a = 1.;
      b = 1.;
      thresh = 0.;
    }

    // The new estimate of the true energy
    t = min(999.9,max(tt, thresh+a*e+b*h));

    // The angular correction 
    // if ( e > 0. && thresh > 0. ) {
    //   etaCorrE = 1.0 + aEtaBarrelEH(t) + 1.3*bEtaBarrelEH(t)*absEta*absEta;
    //   etaCorrH = 1.0;
    // } else {
    //   etaCorrE = 1.0 + aEtaBarrelH(t) + 1.3*bEtaBarrelH(t)*absEta*absEta; 
    //   etaCorrH = 1.0 + aEtaBarrelH(t) + bEtaBarrelH(t)*absEta*absEta;
    // }
    if ( e > 0. && thresh > 0. ) 
      e = h > 0. ? threshE-threshH +  a * e : threshE +  a * e;
    if ( h > 0. && thresh > 0. ) {
      h = threshH + b * h;
    }

  }
 

// FOR eta region 1.6 to 2.0
  if ( absEta > 1.6 && absEta <= 2.0 ) { 
    // The energy correction
    a = e>0. ? a_3(t) : 1.;
    b = e>0. ? b_3(t) : c_3(t);
    double thresh = e > 0. ? threshE : threshH;

    // Protection against negative calibration
    if ( a < -0.25 || b < -0.25 ) { 
      a = 1.;
      b = 1.;
      thresh = 0.;
    }

    // The new estimate of the true energy
    t = min(999.9,max(tt, thresh+a*e+b*h));

    // The angular correction 
    // if ( e > 0. && thresh > 0. ) {
    //   etaCorrE = 1.0 + aEtaBarrelEH(t) + 1.3*bEtaBarrelEH(t)*absEta*absEta;
    //   etaCorrH = 1.0;
    // } else {
    //   etaCorrE = 1.0 + aEtaBarrelH(t) + 1.3*bEtaBarrelH(t)*absEta*absEta; 
    //   etaCorrH = 1.0 + aEtaBarrelH(t) + bEtaBarrelH(t)*absEta*absEta;
    // }
    if ( e > 0. && thresh > 0. ) 
      e = h > 0. ? threshE-threshH +  a * e : threshE +  a * e;
    if ( h > 0. && thresh > 0. ) {
      h = threshH + b * h;
    }

  }
// FOR eta region 2 to 2.5
  if ( absEta > 2 && absEta <= 2.5 ) { 
    // The energy correction
    a = e>0. ? a_4(t) : 1.;
    b = e>0. ? b_4(t) : c_4(t);
    double thresh = e > 0. ? threshE : threshH;

    // Protection against negative calibration
    if ( a < -0.25 || b < -0.25 ) { 
      a = 1.;
      b = 1.;
      thresh = 0.;
    }

    // The new estimate of the true energy
    t = min(999.9,max(tt, thresh+a*e+b*h));

    // The angular correction 
    // if ( e > 0. && thresh > 0. ) {
    //   etaCorrE = 1.0 + aEtaBarrelEH(t) + 1.3*bEtaBarrelEH(t)*absEta*absEta;
    //   etaCorrH = 1.0;
    // } else {
    //   etaCorrE = 1.0 + aEtaBarrelH(t) + 1.3*bEtaBarrelH(t)*absEta*absEta; 
    //   etaCorrH = 1.0 + aEtaBarrelH(t) + bEtaBarrelH(t)*absEta*absEta;
    // }
    if ( e > 0. && thresh > 0. ) 
      e = h > 0. ? threshE-threshH +  a * e : threshE +  a * e;
    if ( h > 0. && thresh > 0. ) {
      h = threshH + b * h;
    }
  }

// FOR eta region 2.5 to 3.0
  if ( absEta > 2.5 && absEta <= 3 ) { 
    // The energy correction
    a = e>0. ? a_5(t) : 1.;
    b = e>0. ? b_5(t) : c_5(t);
    double thresh = e > 0. ? threshE : threshH;

    // Protection against negative calibration
    if ( a < -0.25 || b < -0.25 ) { 
      a = 1.;
      b = 1.;
      thresh = 0.;
    }

    // The new estimate of the true energy
    t = min(999.9,max(tt, thresh+a*e+b*h));

    // The angular correction 
    // if ( e > 0. && thresh > 0. ) {
    //   etaCorrE = 1.0 + aEtaBarrelEH(t) + 1.3*bEtaBarrelEH(t)*absEta*absEta;
    //   etaCorrH = 1.0;
    // } else {
    //   etaCorrE = 1.0 + aEtaBarrelH(t) + 1.3*bEtaBarrelH(t)*absEta*absEta; 
    //   etaCorrH = 1.0 + aEtaBarrelH(t) + bEtaBarrelH(t)*absEta*absEta;
    // }
    if ( e > 0. && thresh > 0. ) 
      e = h > 0. ? threshE-threshH +  a * e : threshE +  a * e;
    if ( h > 0. && thresh > 0. ) {
      h = threshH + b * h;
    }

  }

  // Protection                                                                                                                                
  if ( e < 0. || h < 0. ) {

    // Some protection against crazy calibration                                                                                               
    if ( e < 0. ) e = ee;
    if ( h < 0. ) h = hh;
  }



  // And that's it !

  
}
  
// The calibration functions
  
 
  double 
    PFEnergyCalibration::a_0(double x) const { 

    if ( pfCalibrations ) { 
      BinningPointByMap point;
      point.insert(BinningVariables::JetEt, x);
      return pfCalibrations->getResult(PerformanceResult::PFfa_0,point); 

    } else { 

      return fa_0->Eval(x); 

    }
  }

double 
PFEnergyCalibration::b_0(double x) const { 

  if ( pfCalibrations ) { 

    BinningPointByMap point;
    point.insert(BinningVariables::JetEt, x);
    return pfCalibrations->getResult(PerformanceResult::PFfb_0,point); 

  } else { 

    return fb_0->Eval(x); 

  }
}

double 
PFEnergyCalibration::c_0(double x) const { 

  if ( pfCalibrations ) { 

    BinningPointByMap point;
    point.insert(BinningVariables::JetEt, x);
    return pfCalibrations->getResult(PerformanceResult::PFfc_0,point); 

  } else { 

    return fc_0->Eval(x); 

  }
}


double 
PFEnergyCalibration::a_1(double x) const { 

  if ( pfCalibrations ) { 

    BinningPointByMap point;
    point.insert(BinningVariables::JetEt, x);
    return pfCalibrations->getResult(PerformanceResult::PFfa_1,point); 

  } else { 

    return fa_1->Eval(x); 

  }
}

double 
PFEnergyCalibration::b_1(double x) const { 

  if ( pfCalibrations ) { 

    BinningPointByMap point;
    point.insert(BinningVariables::JetEt, x);
    return pfCalibrations->getResult(PerformanceResult::PFfb_1,point); 

  } else { 

    return fb_1->Eval(x); 

  }
}

double 
PFEnergyCalibration::c_1(double x) const { 

  if ( pfCalibrations ) { 

    BinningPointByMap point;
    point.insert(BinningVariables::JetEt, x);
    return pfCalibrations->getResult(PerformanceResult::PFfc_1,point); 

  } else { 

    return fc_1->Eval(x); 

  }
}

double 
PFEnergyCalibration::a_2(double x) const { 

  if ( pfCalibrations ) { 
    BinningPointByMap point;
    point.insert(BinningVariables::JetEt, x);
    return pfCalibrations->getResult(PerformanceResult::PFfa_2,point); 

  } else { 

    return fa_2->Eval(x); 

  }
}

double 
PFEnergyCalibration::b_2(double x) const { 

  if ( pfCalibrations ) { 

    BinningPointByMap point;
    point.insert(BinningVariables::JetEt, x);
    return pfCalibrations->getResult(PerformanceResult::PFfb_2,point); 

  } else { 

    return fb_2->Eval(x); 

  }
}

double 
PFEnergyCalibration::c_2(double x) const { 

  if ( pfCalibrations ) { 

    BinningPointByMap point;
    point.insert(BinningVariables::JetEt, x);
    return pfCalibrations->getResult(PerformanceResult::PFfc_2,point); 

  } else { 

    return fc_2->Eval(x); 

  }
}
double 
PFEnergyCalibration::a_3(double x) const { 

  if ( pfCalibrations ) { 
    BinningPointByMap point;
    point.insert(BinningVariables::JetEt, x);
    return pfCalibrations->getResult(PerformanceResult::PFfa_3,point); 

  } else { 

    return fa_3->Eval(x); 

  }
}

double 
PFEnergyCalibration::b_3(double x) const { 

  if ( pfCalibrations ) { 

    BinningPointByMap point;
    point.insert(BinningVariables::JetEt, x);
    return pfCalibrations->getResult(PerformanceResult::PFfb_3,point); 

  } else { 

    return fb_3->Eval(x); 

  }
}

double 
PFEnergyCalibration::c_3(double x) const { 

  if ( pfCalibrations ) { 

    BinningPointByMap point;
    point.insert(BinningVariables::JetEt, x);
    return pfCalibrations->getResult(PerformanceResult::PFfc_3,point); 

  } else { 

    return fc_3->Eval(x); 

  }
}

double 
PFEnergyCalibration::a_4(double x) const { 

  if ( pfCalibrations ) { 
    BinningPointByMap point;
    point.insert(BinningVariables::JetEt, x);
    return pfCalibrations->getResult(PerformanceResult::PFfa_4,point); 

  } else { 

    return fa_4->Eval(x); 

  }
}

double 
PFEnergyCalibration::b_4(double x) const { 

  if ( pfCalibrations ) { 

    BinningPointByMap point;
    point.insert(BinningVariables::JetEt, x);
    return pfCalibrations->getResult(PerformanceResult::PFfb_4,point); 

  } else { 

    return fb_4->Eval(x); 

  }
}

double 
PFEnergyCalibration::c_4(double x) const { 

  if ( pfCalibrations ) { 

    BinningPointByMap point;
    point.insert(BinningVariables::JetEt, x);
    return pfCalibrations->getResult(PerformanceResult::PFfc_4,point); 

  } else { 

    return fc_4->Eval(x); 

  }
}

double 
PFEnergyCalibration::a_5(double x) const { 

  if ( pfCalibrations ) { 
    BinningPointByMap point;
    point.insert(BinningVariables::JetEt, x);
    return pfCalibrations->getResult(PerformanceResult::PFfa_5,point); 

  } else { 

    return fa_5->Eval(x); 

  }
}

double 
PFEnergyCalibration::b_5(double x) const { 

  if ( pfCalibrations ) { 

    BinningPointByMap point;
    point.insert(BinningVariables::JetEt, x);
    return pfCalibrations->getResult(PerformanceResult::PFfb_5,point); 

  } else { 

    return fb_5->Eval(x); 

  }
}

double 
PFEnergyCalibration::c_5(double x) const { 

  if ( pfCalibrations ) { 

    BinningPointByMap point;
    point.insert(BinningVariables::JetEt, x);
    return pfCalibrations->getResult(PerformanceResult::PFfc_5,point); 

  } else { 

    return fc_5->Eval(x); 

  }
}

double
PFEnergyCalibration::energyEm(const reco::PFCluster& clusterEcal,
			      std::vector<double> &EclustersPS1,
			      std::vector<double> &EclustersPS2,
			      bool crackCorrection ) const {
  double ePS1(std::accumulate(EclustersPS1.begin(), EclustersPS1.end(), 0.0));
  double ePS2(std::accumulate(EclustersPS2.begin(), EclustersPS2.end(), 0.0));
  return energyEm(clusterEcal, ePS1, ePS2, crackCorrection);
}

double
PFEnergyCalibration::energyEm(const reco::PFCluster& clusterEcal,
			      double ePS1,
			      double ePS2,
			      bool crackCorrection ) const {
  double eEcal = clusterEcal.energy();
  //temporaty ugly fix
  reco::PFCluster myPFCluster=clusterEcal;
  myPFCluster.calculatePositionREP();
  double eta = myPFCluster.positionREP().eta();
  double phi = myPFCluster.positionREP().phi();

  double calibrated = Ecorr(eEcal,ePS1,ePS2,eta,phi, crackCorrection);
  // if(eEcal!=0 && calibrated==0) std::cout<<"Eecal = "<<eEcal<<"  eta = "<<eta<<"  phi = "<<phi<<std::endl; 
  return calibrated; 
}

double PFEnergyCalibration::energyEm(const reco::PFCluster& clusterEcal,
				     std::vector<double> &EclustersPS1,
				     std::vector<double> &EclustersPS2,
				     double& ps1,double& ps2,
				     bool crackCorrection) const {
  double ePS1(std::accumulate(EclustersPS1.begin(), EclustersPS1.end(), 0.0));
  double ePS2(std::accumulate(EclustersPS2.begin(), EclustersPS2.end(), 0.0));
  return energyEm(clusterEcal, ePS1, ePS2, ps1, ps2, crackCorrection);
}
double PFEnergyCalibration::energyEm(const reco::PFCluster& clusterEcal,
				     double ePS1, double ePS2,
				     double& ps1,double& ps2,
				     bool crackCorrection) const {
  double eEcal = clusterEcal.energy();
  //temporaty ugly fix
  reco::PFCluster myPFCluster=clusterEcal;
  myPFCluster.calculatePositionREP();
  double eta = myPFCluster.positionREP().eta();
  double phi = myPFCluster.positionREP().phi();

  double calibrated = Ecorr(eEcal,ePS1,ePS2,eta,phi,ps1,ps2,crackCorrection);
  // if(eEcal!=0 && calibrated==0) std::cout<<"Eecal = "<<eEcal<<"  eta = "<<eta<<"  phi = "<<phi<<std::endl; 
  return calibrated; 
}


std::ostream& operator<<(std::ostream& out,
			 const PFEnergyCalibration& calib) {

  if(!out ) return out;

  out<<"PFEnergyCalibration -- "<<endl;

  if ( calib.pfCalibrations ) { 

    static std::map<std::string, PerformanceResult::ResultType> functType;

   
    functType["PFfa_0"] = PerformanceResult::PFfa_0;
    functType["PFfb_0"] = PerformanceResult::PFfb_0;
    functType["PFfc_0"] = PerformanceResult::PFfc_0;
    functType["PFfa_1"] = PerformanceResult::PFfa_1;
    functType["PFfb_1"] = PerformanceResult::PFfb_1;
    functType["PFfc_1"] = PerformanceResult::PFfc_1;
    functType["PFfa_2"] = PerformanceResult::PFfa_2;
    functType["PFfb_2"] = PerformanceResult::PFfb_2;
    functType["PFfc_2"] = PerformanceResult::PFfc_2;
    functType["PFfa_3"] = PerformanceResult::PFfa_3;
    functType["PFfb_3"] = PerformanceResult::PFfb_3;
    functType["PFfc_3"] = PerformanceResult::PFfc_3;
    functType["PFfa_4"] = PerformanceResult::PFfa_4;
    functType["PFfb_4"] = PerformanceResult::PFfb_4;
    functType["PFfc_5"] = PerformanceResult::PFfc_4;
    functType["PFfa_5"] = PerformanceResult::PFfa_5;
    functType["PFfb_5"] = PerformanceResult::PFfb_5;
    functType["PFfc_5"] = PerformanceResult::PFfc_5;
    
    for(std::map<std::string,PerformanceResult::ResultType>::const_iterator 
	  func = functType.begin(); 
        func != functType.end(); 
        ++func) {    
      
      cout << "Function: " << func->first << endl;
      PerformanceResult::ResultType fType = func->second;
      calib.pfCalibrations->printFormula(fType);
    }

  } else { 
    
    std::cout << "Default calibration functions : " << std::endl;
    
    calib.fa_0->Print();
    calib.fb_0->Print();
    calib.fc_0->Print();
    calib.fa_1->Print();
    calib.fb_1->Print();
    calib.fc_1->Print();
    calib.fa_2->Print();
    calib.fb_2->Print();
    calib.fc_2->Print();
    calib.fa_3->Print();
    calib.fb_3->Print();
    calib.fc_3->Print();
    calib.fa_4->Print();
    calib.fb_4->Print();
    calib.fc_4->Print();
    calib.fa_5->Print();
    calib.fb_5->Print();
    calib.fc_5->Print();
   
  }
    
  return out;
}




///////////////////////////////////////////////////////////////
////                                                       ////  
////             CORRECTION OF PHOTONS' ENERGY             ////
////                                                       ////
////              Material effect: No tracker              ////
////       Tuned on CMSSW_2_1_0_pre4, Full Sim events      ////
////                                                       ////
///////////////////////////////////////////////////////////////
////                                                       ////
////            Jonathan Biteau - June 2008                ////
////                                                       ////
///////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////
////                                                       ////  
////  USEFUL FUNCTIONS FOR THE CORRECTION IN THE BARREL    ////
////                                                       ////
///////////////////////////////////////////////////////////////


//useful to compute the signed distance to the closest crack in the barrel
double
PFEnergyCalibration::minimum(double a,double b) const {
  if(TMath::Abs(b)<TMath::Abs(a)) a=b;
  return a;
}

namespace {
  constexpr double pi= M_PI;// 3.14159265358979323846;


  std::vector<double> fillcPhi() {
    std::vector<double> retValue;
    retValue.resize(18,0);
    retValue[0]=2.97025;
    for(unsigned i=1;i<=17;++i) retValue[i]=retValue[0]-2*i*pi/18;
    
    return retValue;
  }
  
  //Location of the 18 phi-cracks
  const std::vector<double> cPhi = fillcPhi();
}

//compute the unsigned distance to the closest phi-crack in the barrel
double
PFEnergyCalibration::dCrackPhi(double phi, double eta) const {

  
  //Shift of this location if eta<0
  constexpr double delta_cPhi=0.00638;

  double m; //the result

  //the location is shifted
  if(eta<0) phi +=delta_cPhi;

  if (phi>=-pi && phi<=pi){

    //the problem of the extrema
    if (phi<cPhi[17] || phi>=cPhi[0]){
      if (phi<0) phi+= 2*pi;
      m = minimum(phi -cPhi[0],phi-cPhi[17]-2*pi);        	
    }

    //between these extrema...
    else{
      bool OK = false;
      unsigned i=16;
      while(!OK){
	if (phi<cPhi[i]){
	  m=minimum(phi-cPhi[i+1],phi-cPhi[i]);
	  OK=true;
	}
	else i-=1;
      }
    }
  }
  else{
    m=0.;        //if there is a problem, we assum that we are in a crack
    std::cout<<"Problem in dminphi"<<std::endl;
  }
  if(eta<0) m=-m;   //because of the disymetry
  return m;
}

// corrects the effect of phi-cracks
double
PFEnergyCalibration::CorrPhi(double phi, double eta) const {

  // we use 3 gaussians to correct the phi-cracks effect
  constexpr double p1=   5.59379e-01;
  constexpr double p2=   -1.26607e-03;
  constexpr double p3=  9.61133e-04;

  constexpr double p4=   1.81691e-01;
  constexpr double p5=   -4.97535e-03;
  constexpr double p6=   1.31006e-03;

  constexpr double p7=   1.38498e-01;
  constexpr double p8=   1.18599e-04;
  constexpr double p9= 2.01858e-03;
  

  double dminphi = dCrackPhi(phi,eta);
  
  double result = (1+p1*TMath::Gaus(dminphi,p2,p3)+p4*TMath::Gaus(dminphi,p5,p6)+p7*TMath::Gaus(dminphi,p8,p9));

  return result;
}   


// corrects the effect of  |eta|-cracks
double
PFEnergyCalibration::CorrEta(double eta) const {
  
  // we use a gaussian with a screwness for each of the 5 |eta|-cracks
  constexpr double a[] = {6.13349e-01, 5.08146e-01, 4.44480e-01, 3.3487e-01, 7.65627e-01}; // amplitude
  constexpr double m[] = {-1.79514e-02, 4.44747e-01, 7.92824e-01, 1.14090e+00, 1.47464e+00}; // mean
  constexpr double s[] = {7.92382e-03, 3.06028e-03, 3.36139e-03, 3.94521e-03, 8.63950e-04}; // sigma
  constexpr double sa[] = {1.27228e+01, 3.81517e-02, 1.63507e-01, -6.56480e-02, 1.87160e-01}; // screwness amplitude
  constexpr double ss[] = {5.48753e-02, -1.00223e-02, 2.22866e-03, 4.26288e-04, 2.67937e-03}; // screwness sigma
  double result = 1;

  for(unsigned i=0;i<=4;i++) result+=a[i]*TMath::Gaus(eta,m[i],s[i])*(1+sa[i]*TMath::Sign(1.,eta-m[i])*TMath::Exp(-TMath::Abs(eta-m[i])/ss[i]));

  return result;
}


//corrects the global behaviour in the barrel
double
PFEnergyCalibration::CorrBarrel(double E, double eta) const {

  //Energy dependency
  /*
  //YM Parameters 52XX:
  constexpr double p0=1.00000e+00;
  constexpr double p1=3.27753e+01;
  constexpr double p2=2.28552e-02;
  constexpr double p3=3.06139e+00;
  constexpr double p4=2.25135e-01;
  constexpr double p5=1.47824e+00;
  constexpr double p6=1.09e-02;
  constexpr double p7=4.19343e+01;
  */
  constexpr double p0 = 0.9944;
  constexpr double p1 = 9.827;
  constexpr double p2 = 1.503;
  constexpr double p3 = 1.196;
  constexpr double p4 = 0.3349;
  constexpr double p5 = 0.89;
  constexpr double p6 = 0.004361;
  constexpr double p7 = 51.51;
  //Eta dependency
  constexpr double p8=2.705593e-03;
  
  double result = (p0+1/(p1+p2*TMath::Power(E,p3))+p4*TMath::Exp(-E/p5)+p6*TMath::Exp(-E*E/(p7*p7)))*(1+p8*eta*eta);

  return result;
}



///////////////////////////////////////////////////////////////
////                                                       ////  
////  USEFUL FUNCTIONS FOR THE CORRECTION IN THE ENDCAPS   ////
////  Parameters tuned for:                                ////
////          dR(ClustersPS1,ClusterEcal) < 0.08           ////
////          dR(ClustersPS2,ClusterEcal) < 0.13           ////
////                                                       ////
///////////////////////////////////////////////////////////////


//Alpha, Beta, Gamma give the weight of each sub-detector (PS layer1, PS layer2 and Ecal) in the areas of the endcaps where there is a PS
// Etot = Beta*eEcal + Gamma*(ePS1 + Alpha*ePS2) 

double
PFEnergyCalibration::Alpha(double eta) const {

  //Energy dependency
  constexpr double p0 = 5.97621e-01;

  //Eta dependency
  constexpr double p1 =-1.86407e-01;
  constexpr double p2 = 3.85197e-01; 

  //so that <feta()> = 1
  constexpr double norm = (p1+p2*(2.6+1.656)/2);

  double result = p0*(p1+p2*eta)/norm;

  return result;
}

double
PFEnergyCalibration::Beta(double E, double eta) const {

 //Energy dependency
  constexpr double p0 = 0.032;
  constexpr double p1 = 9.70394e-02;
  constexpr double p2 = 2.23072e+01;
  constexpr double p3 = 100;

  //Eta dependency
  constexpr double p4 = 1.02496e+00 ;
  constexpr double p5 = -4.40176e-03 ;

  //so that <feta()> = 1
  constexpr double norm = (p4+p5*(2.6+1.656)/2);

  double result = (1.0012+p0*TMath::Exp(-E/p3)+p1*TMath::Exp(-E/p2))*(p4+p5*eta)/norm;			  
  return result;
}


double
PFEnergyCalibration::Gamma(double etaEcal) const {

 //Energy dependency
  constexpr double p0 = 2.49752e-02;

  //Eta dependency
  constexpr double p1 = 6.48816e-02;
  constexpr double p2 = -1.59517e-02; 
 
  //so that <feta()> = 1
  constexpr double norm = (p1+p2*(2.6+1.656)/2);

  double result = p0*(p1+p2*etaEcal)/norm;					  

  return result;
}



///////////////////////////////////////////////////////////////
////                                                       ////  
////   THE CORRECTIONS IN THE BARREL AND IN THE ENDCAPS    ////
////                                                       ////
///////////////////////////////////////////////////////////////


// returns the corrected energy in the barrel (0,1.48)
// Global Behaviour, phi and eta cracks are taken into account
double
PFEnergyCalibration::EcorrBarrel(double E, double eta, double phi,
				 bool crackCorrection ) const {

  // double result = E*CorrBarrel(E,eta)*CorrEta(eta)*CorrPhi(phi,eta);
  double correction = crackCorrection ? std::max(CorrEta(eta),CorrPhi(phi,eta)) : 1.;
  double result = E * CorrBarrel(E,eta) * correction;

  return result;
}


// returns the corrected energy in the area between the barrel and the PS (1.48,1.65)
double
PFEnergyCalibration::EcorrZoneBeforePS(double E, double eta) const {

 //Energy dependency
  constexpr double p0 =1; 
  constexpr double p1 =0.18;
  constexpr double p2 =8.;

  //Eta dependency
  constexpr double p3 =0.3;
  constexpr double p4 =1.11;
  constexpr double p5 =0.025;
  constexpr double p6 =1.49;
  constexpr double p7 =0.6;

  //so that <feta()> = 1
  constexpr double norm = 1.21;

  double result = E*(p0+p1*TMath::Exp(-E/p2))*(p3+p4*TMath::Gaus(eta,p6,p5)+p7*eta)/norm;

  return result;
}


// returns the corrected energy in the PS (1.65,2.6)
// only when (ePS1>0)||(ePS2>0)
double
PFEnergyCalibration::EcorrPS(double eEcal,double ePS1,double ePS2,double etaEcal) const {

  // gives the good weights to each subdetector
  double E = Beta(1.0155*eEcal+0.025*(ePS1+0.5976*ePS2)/9e-5,etaEcal)*eEcal+Gamma(etaEcal)*(ePS1+Alpha(etaEcal)*ePS2)/9e-5 ;

  //Correction of the residual energy dependency
  constexpr double p0 = 1.00;
  constexpr double p1 = 2.18;
  constexpr double p2 =1.94;
  constexpr double p3 =4.13;
  constexpr double p4 =1.127;

  double result = E*(p0+p1*TMath::Exp(-E/p2)-p3*TMath::Exp(-E/p4));

  return result;
} 

// returns the corrected energy in the PS (1.65,2.6)
// only when (ePS1>0)||(ePS2>0)
double
PFEnergyCalibration::EcorrPS(double eEcal,double ePS1,double ePS2,double etaEcal,double & outputPS1, double & outputPS2) const {

  // gives the good weights to each subdetector
  double gammaprime=Gamma(etaEcal)/9e-5;

  if(outputPS1 == 0 && outputPS2 == 0 && esEEInterCalib_ != nullptr){
    // both ES planes working
    // scaling factor accounting for data-mc                                                                                 
    outputPS1=gammaprime*ePS1 * esEEInterCalib_->getGammaLow0();
    outputPS2=gammaprime*Alpha(etaEcal)*ePS2 * esEEInterCalib_->getGammaLow3();
  }
  else if(outputPS1 == 0 && outputPS2 == -1 && esEEInterCalib_ != nullptr){
    // ESP1 only working
    double corrTotES = gammaprime*ePS1 * esEEInterCalib_->getGammaLow0() * esEEInterCalib_->getGammaLow1();
    outputPS1 = gammaprime*ePS1 * esEEInterCalib_->getGammaLow0();
    outputPS2 = corrTotES - outputPS1;
  }
  else if(outputPS1 == -1 && outputPS2 == 0 && esEEInterCalib_ != nullptr){
    // ESP2 only working
    double corrTotES = gammaprime*Alpha(etaEcal)*ePS2 * esEEInterCalib_->getGammaLow3() * esEEInterCalib_->getGammaLow2();
    outputPS2 = gammaprime*Alpha(etaEcal)*ePS2 * esEEInterCalib_->getGammaLow3();
    outputPS1 = corrTotES - outputPS2;
  }
  else{
    // none working
    outputPS1 = gammaprime*ePS1;
    outputPS2 = gammaprime*Alpha(etaEcal)*ePS2;
  }

  double E = Beta(1.0155*eEcal+0.025*(ePS1+0.5976*ePS2)/9e-5,etaEcal)*eEcal+outputPS1+outputPS2;

  //Correction of the residual energy dependency
  constexpr double p0 = 1.00;
  constexpr double p1 = 2.18;
  constexpr double p2 =1.94;
  constexpr double p3 =4.13;
  constexpr double p4 =1.127;
  
  double corrfac=(p0+p1*TMath::Exp(-E/p2)-p3*TMath::Exp(-E/p4));
  outputPS1*=corrfac;
  outputPS2*=corrfac;
  double result = E*corrfac;

  return result;
} 


// returns the corrected energy in the PS (1.65,2.6)
// only when (ePS1=0)&&(ePS2=0)
double 
PFEnergyCalibration::EcorrPS_ePSNil(double eEcal,double eta) const {

  //Energy dependency
  constexpr double p0= 1.02;
  constexpr double p1= 0.165;
  constexpr double p2= 6.5 ;
  constexpr double p3=  2.1 ;

  //Eta dependency
  constexpr double p4 = 1.02496e+00 ;
  constexpr double p5 = -4.40176e-03 ;

  //so that <feta()> = 1
  constexpr double norm = (p4+p5*(2.6+1.656)/2);

  double result = eEcal*(p0+p1*TMath::Exp(-TMath::Abs(eEcal-p3)/p2))*(p4+p5*eta)/norm;
		  
  return result;
}


// returns the corrected energy in the area between the end of the PS and the end of the endcap (2.6,2.98)
double
PFEnergyCalibration::EcorrZoneAfterPS(double E, double eta) const {

  //Energy dependency
  constexpr double p0 =1; 
  constexpr double p1 = 0.058;
  constexpr double p2 =12.5;
  constexpr double p3 =-1.05444e+00;
  constexpr double p4 =-5.39557e+00;
  constexpr double p5 =8.38444e+00;
  constexpr double p6 = 6.10998e-01  ;

  //Eta dependency
  constexpr double p7 =1.06161e+00;
  constexpr double p8 = 0.41;
  constexpr double p9 =2.918;
  constexpr double p10 =0.0181;
  constexpr double p11= 2.05;
  constexpr double p12 =2.99;
  constexpr double p13=0.0287;

  //so that <feta()> = 1
  constexpr double norm=1.045;

  double result = E*(p0+p1*TMath::Exp(-(E-p3)/p2)+1/(p4+p5*TMath::Power(E,p6)))*(p7+p8*TMath::Gaus(eta,p9,p10)+p11*TMath::Gaus(eta,p12,p13))/norm;
  return result;
}




// returns the corrected energy everywhere
// this work should be improved between 1.479 and 1.52 (junction barrel-endcap)
double
PFEnergyCalibration::Ecorr(double eEcal,double ePS1,double ePS2,
			   double eta,double phi,
			   bool crackCorrection ) const {

  constexpr double endBarrel=1.48;
  constexpr double beginingPS=1.65;
  constexpr double endPS=2.6;
  constexpr double endEndCap=2.98;
 
  double result=0;

  eta=TMath::Abs(eta);

  if(eEcal>0){
    if(eta <= endBarrel)                         result = EcorrBarrel(eEcal,eta,phi,crackCorrection);
    else if(eta <= beginingPS)                   result = EcorrZoneBeforePS(eEcal,eta);
    else if((eta < endPS) && ePS1==0 && ePS2==0) result = EcorrPS_ePSNil(eEcal,eta);
    else if(eta < endPS)                         result = EcorrPS(eEcal,ePS1,ePS2,eta);
    else if(eta < endEndCap)                     result = EcorrZoneAfterPS(eEcal,eta); 
    else result =eEcal;
  }
  else result = eEcal;// useful if eEcal=0 or eta>2.98
  //protection
  if(result<eEcal) result=eEcal;
  return result;
}

// returns the corrected energy everywhere
// this work should be improved between 1.479 and 1.52 (junction barrel-endcap)
double
PFEnergyCalibration::Ecorr(double eEcal,double ePS1,double ePS2,double eta,double phi,double& ps1,double&ps2,bool crackCorrection)  const {

  constexpr double endBarrel=1.48;
  constexpr double beginingPS=1.65;
  constexpr double endPS=2.6;
  constexpr double endEndCap=2.98;
 
  double result=0;

  eta=TMath::Abs(eta);

  if(eEcal>0){
    if(eta <= endBarrel)                         result = EcorrBarrel(eEcal,eta,phi,crackCorrection);
    else if(eta <= beginingPS)                   result = EcorrZoneBeforePS(eEcal,eta);
    else if((eta < endPS) && ePS1==0 && ePS2==0) result = EcorrPS_ePSNil(eEcal,eta);
    else if(eta < endPS)                         result = EcorrPS(eEcal,ePS1,ePS2,eta,ps1,ps2);
    else if(eta < endEndCap)                     result = EcorrZoneAfterPS(eEcal,eta); 
    else result =eEcal;
  }
  else result = eEcal;// useful if eEcal=0 or eta>2.98
  // protection
  if(result<eEcal) result=eEcal;
  return result;
}
