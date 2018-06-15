import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")
#process.load("CondCore.DBCommon.CondDBCommon_cfi")
process.load("CondCore.CondDB.CondDB_cfi")
#process.CondDBCommon.connect = 'sqlite_file:PhysicsPerformance.db'


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)
process.source = cms.Source("EmptySource",
                            firstRun = cms.untracked.uint32(10)
                            )

# process.PoolDBOutputService.DBParameters.messageLevel = 3


process.mywriter = cms.EDAnalyzer(
  "ProducePFCalibrationObject",
  write = cms.untracked.bool(False),
  toWrite = cms.VPSet(
            cms.PSet(fType      = cms.untracked.string("PFfa_0"),
                     formula    = cms.untracked.string("[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))"),
                     limits     = cms.untracked.vdouble(1., 1000.),
                     parameters = cms.untracked.vdouble(-13.9219, 14.9124, 5.38578, 0.861981, -0.00759275, 0.00373563, -1.17946, -1.69561, )
                    ),
            cms.PSet(fType      = cms.untracked.string("PFfb_0"),
                     formula    = cms.untracked.string("[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))"),
                     limits     = cms.untracked.vdouble(1., 1000.),
                     parameters = cms.untracked.vdouble(2.25366, 0.537715, -4.81375, 12.109, 1.80577, 0.187919, -6.26234, -0.607392, )
                    ),
            cms.PSet(fType      = cms.untracked.string("PFfc_0"),
                     formula    = cms.untracked.string("[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))"),
                     limits     = cms.untracked.vdouble(1., 1000.),
                     parameters = cms.untracked.vdouble(1.5126, 0.855057, -6.04199, 2.08229, 0.592266, 0.0291232, 0.364802, -1.50142, )
                    ),
            cms.PSet(fType      = cms.untracked.string("PFfa_1"),
                     formula    = cms.untracked.string("[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))"),
                     limits     = cms.untracked.vdouble(1., 1000.),
                     parameters = cms.untracked.vdouble(-13.9219, 14.9124, 5.38578, 0.861981, -0.00759275, 0.00373563, -1.17946, -1.69561, )
                    ),
            cms.PSet(fType      = cms.untracked.string("PFfb_1"),
                     formula    = cms.untracked.string("[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))"),
                     limits     = cms.untracked.vdouble(1., 1000.),
                     parameters = cms.untracked.vdouble(2.25366, 0.537715, -4.81375, 12.109, 1.80577, 0.187919, -6.26234, -0.607392, )
                    ),
            cms.PSet(fType      = cms.untracked.string("PFfc_1"),
                     formula    = cms.untracked.string("[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))"),
                     limits     = cms.untracked.vdouble(1., 1000.),
                     parameters = cms.untracked.vdouble(1.5126, 0.855057, -6.04199, 2.08229, 0.592266, 0.0291232, 0.364802, -1.50142, )
                    ),
            cms.PSet(fType      = cms.untracked.string("PFfa_2"),
                     formula    = cms.untracked.string("[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))"),
                     limits     = cms.untracked.vdouble(1., 1000.),
                     parameters = cms.untracked.vdouble(-13.9219, 14.9124, 5.38578, 0.861981, -0.00759275, 0.00373563, -1.17946, -1.69561, )
                    ),
            cms.PSet(fType      = cms.untracked.string("PFfb_2"),
                     formula    = cms.untracked.string("[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))"),
                     limits     = cms.untracked.vdouble(1., 1000.),
                     parameters = cms.untracked.vdouble(2.25366, 0.537715, -4.81375, 12.109, 1.80577, 0.187919, -6.26234, -0.607392, )
                    ),
            cms.PSet(fType      = cms.untracked.string("PFfc_2"),
                     formula    = cms.untracked.string("[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))"),
                     limits     = cms.untracked.vdouble(1., 1000.),
                     parameters = cms.untracked.vdouble(1.5126, 0.855057, -6.04199, 2.08229, 0.592266, 0.0291232, 0.364802, -1.50142, )
                    ),
            cms.PSet(fType      = cms.untracked.string("PFfa_3"),
                     formula    = cms.untracked.string("[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))"),
                     limits     = cms.untracked.vdouble(1., 1000.),
                     parameters = cms.untracked.vdouble(-13.9219, 14.9124, 5.38578, 0.861981, -0.00759275, 0.00373563, -1.17946, -1.69561, )
                    ),
            cms.PSet(fType      = cms.untracked.string("PFfb_3"),
                     formula    = cms.untracked.string("[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))"),
                     limits     = cms.untracked.vdouble(1., 1000.),
                     parameters = cms.untracked.vdouble(2.25366, 0.537715, -4.81375, 12.109, 1.80577, 0.187919, -6.26234, -0.607392, )
                    ),
            cms.PSet(fType      = cms.untracked.string("PFfc_3"),
                     formula    = cms.untracked.string("[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))"),
                     limits     = cms.untracked.vdouble(1., 1000.),
                     parameters = cms.untracked.vdouble(1.5126, 0.855057, -6.04199, 2.08229, 0.592266, 0.0291232, 0.364802, -1.50142, )
                    ),
            cms.PSet(fType      = cms.untracked.string("PFfa_4"),
                     formula    = cms.untracked.string("[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))"),
                     limits     = cms.untracked.vdouble(1., 1000.),
                     parameters = cms.untracked.vdouble(-13.9219, 14.9124, 5.38578, 0.861981, -0.00759275, 0.00373563, -1.17946, -1.69561, )
                    ),
            cms.PSet(fType      = cms.untracked.string("PFfb_4"),
                     formula    = cms.untracked.string("[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))"),
                     limits     = cms.untracked.vdouble(1., 1000.),
                     parameters = cms.untracked.vdouble(2.25366, 0.537715, -4.81375, 12.109, 1.80577, 0.187919, -6.26234, -0.607392, )
                    ),
            cms.PSet(fType      = cms.untracked.string("PFfc_4"),
                     formula    = cms.untracked.string("[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))"),
                     limits     = cms.untracked.vdouble(1., 1000.),
                     parameters = cms.untracked.vdouble(1.5126, 0.855057, -6.04199, 2.08229, 0.592266, 0.0291232, 0.364802, -1.50142, )
                    ),
            cms.PSet(fType      = cms.untracked.string("PFfa_5"),
                     formula    = cms.untracked.string("[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))"),
                     limits     = cms.untracked.vdouble(1., 1000.),
                     parameters = cms.untracked.vdouble(-13.9219, 14.9124, 5.38578, 0.861981, -0.00759275, 0.00373563, -1.17946, -1.69561, )
                    ),
            cms.PSet(fType      = cms.untracked.string("PFfb_5"),
                     formula    = cms.untracked.string("[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))"),
                     limits     = cms.untracked.vdouble(1., 1000.),
                     parameters = cms.untracked.vdouble(2.25366, 0.537715, -4.81375, 12.109, 1.80577, 0.187919, -6.26234, -0.607392, )
                    ),
            cms.PSet(fType      = cms.untracked.string("PFfc_5"),
                     formula    = cms.untracked.string("[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))"),
                     limits     = cms.untracked.vdouble(1., 1000.),
                     parameters = cms.untracked.vdouble(1.5126, 0.855057, -6.04199, 2.08229, 0.592266, 0.0291232, 0.364802, -1.50142, )
                    ),
           
            ),


  read = cms.untracked.bool(True),
  toRead = cms.untracked.vstring("PFfa_0",
                                 "PFfb_0",
                                 "PFfc_0",
                                 "PFfa_1",
                                 "PFfb_1",
                                 "PFfc_1",
                                 "PFfa_2",
                                 "PFfb_2",
                                 "PFfc_2",
                                 "PFfa_3",
                                 "PFfb_3",
                                 "PFfc_3",
                                 "PFfa_4",
                                 "PFfb_4",
                                 "PFfc_4",
                                 "PFfa_5",
                                 "PFfb_5",
                                 "PFfc_5",
                                


                                 #### New Functions
                                 # "PFfaEta_BARRELEH",
                                 # "PFfbEta_BARRELEH",
                                 # "PFfaEta_ENDCAPEH",
                                 # "PFfbEta_ENDCAPEH",
                                 # "PFfaEta_BARRELH",
                                 # "PFfbEta_BARRELH",
                                 # "PFfaEta_ENDCAPH",
                                 # "PFfbEta_ENDCAPH",
                                 # #### Left older functions untouched for backward compatibility
                                 # "PFfaEta_BARREL",
                                 # "PFfbEta_BARREL",
                                 # "PFfaEta_ENDCAP",
                                 # "PFfbEta_ENDCAP",

                                 ) # same strings as fType
)


process.p = cms.Path(process.mywriter)

from CondCore.DBCommon.CondDBCommon_cfi import CondDBCommon
CondDBCommon.connect = "sqlite_file:PFCalibration.db"

process.PoolDBOutputService = cms.Service("PoolDBOutputService",
                                  CondDBCommon,
                                  toPut = cms.VPSet(cms.PSet(record = cms.string('PFCalibrationRcd'),
                                                             tag = cms.string('PFCalibration_v7_mc'),
                                                             timetype   = cms.untracked.string('runnumber')
                                                             )
                                                             ),
                                  loadBlobStreamer = cms.untracked.bool(False),
                                  #    timetype   = cms.untracked.string('lumiid')
                                  #    timetype   = cms.untracked.string('runnumber')
                                  )

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = '90X_upgrade2017_realistic_v20'
#process.GlobalTag.connect   = 'sqlite_file:/afs/cern.ch/user/c/cerminar/public/Alca/GlobalTag/GR_R_311_V2.db'

process.GlobalTag.toGet = cms.VPSet(
    cms.PSet(record = cms.string("PFCalibrationRcd"),
             tag = cms.string("PFCalibration_v7_mc"),
             connect = cms.string("sqlite_file:PFCalibration.db")
             #connect = cms.untracked.string("sqlite_file:PFCalibration.db")
             )
    )
