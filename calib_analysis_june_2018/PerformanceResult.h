#ifndef PerformanceResult_h
#define PerformanceResult_h

#include "CondFormats/Serialization/interface/Serializable.h"

class PerformanceResult {
 public:
  enum ResultType {
    //
    // BTAG
    //
    BTAGBEFF=1001, BTAGBERR=1002, BTAGCEFF=1003, 
    BTAGCERR=1004, BTAGLEFF=1005, BTAGLERR=1006, BTAGNBEFF=1007, BTAGNBERR=1008,
    //
    // add corrections in case the table is for weights and not efficiencies
    //
    BTAGBEFFCORR=1009, BTAGBERRCORR=1010, BTAGCEFFCORR=1011, 
    BTAGCERRCORR=1012, BTAGLEFFCORR=1013, BTAGLERRCORR=1014, BTAGNBEFFCORR=1015, BTAGNBERRCORR=1016,
    //
    // MUONS
    //
    MUEFF=2001, MUERR=2002, MUFAKE=2003, MUEFAKE=2004,
    //
    // PF - calibrations
    //
    /* PFfa_BARREL = 3001, PFfa_ENDCAP = 3002, */
    /* PFfb_BARREL = 3003, PFfb_ENDCAP = 3004, */
    /* PFfc_BARREL = 3005, PFfc_ENDCAP = 3006, */
    /* PFfaEta_BARREL = 3007, PFfaEta_ENDCAP = 3008,  */
    /* PFfbEta_BARREL = 3009, PFfbEta_ENDCAP = 3010,  */
    /* PFfaEta_BARRELH = 3011, PFfaEta_ENDCAPH = 3012,  */
    /* PFfbEta_BARRELH = 3013, PFfbEta_ENDCAPH = 3014,  */
    /* PFfaEta_BARRELEH = 3015, PFfaEta_ENDCAPEH = 3016,  */
    /* PFfbEta_BARRELEH = 3017, PFfbEta_ENDCAPEH = 3018  */
    PFfa_0 = 3001 , PFfb_0 = 3002 , PFfc_0 = 3003,
    PFfa_1 = 3004 , PFfb_1 = 3005 , PFfc_1 = 3006,
    PFfa_2 = 3007 , PFfb_2 = 3008 , PFfc_2 = 3009,
    PFfa_3 = 3010 , PFfb_3 = 3011 , PFfc_3 = 3012,
    PFfa_4 = 3013 , PFfb_4 = 3014 , PFfc_4 = 3015,
    PFfa_5 = 3016 , PFfb_5 = 3017 , PFfc_5 = 3018,
   
    
};
  

 COND_SERIALIZABLE;
};

#endif
