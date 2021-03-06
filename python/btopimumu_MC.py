################
#  variables
################
runMC               = True
run2012not2011      = True  # True for 2012 and False for 2011

print "\n@@@ CMSSW run configuration flags @@@"

if (runMC==True):
    print " running MC : ", runMC
else:
    print "not running MC, please xCheck!"


if (run2012not2011 == True):
    print "==> 2012 MC running "
else:
    print "==> 2011 MC running "

print "========================================\n"


##################
# cmssw configs
################## 

import FWCore.ParameterSet.Config as cms
from btopimumu_cfi import process

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
#     '/store/mc/Summer12_DR53X/BuToMuMuPi_EtaPtFilter_8TeV-pythia6-evtgen/AODSIM/PU_S10_START53_V7C-v1/10000/022D253A-A85D-E211-A9AC-00266CFFCCC8.root',
     '/store/mc/Summer12_DR53X/BuToJpsiPi_EtaPtFilter_8TeV-pythia6-evtgen/AODSIM/PU_S10_START53_V7C-v1/20000/007929BA-EF65-E211-B216-00266CFFCCC8.root',
                            )
                        )

if (run2012not2011 == True):
    process.GlobalTag.globaltag = cms.string('START53_V7G::All')
    print "\nGlobalTag : START53_V7G::All\n"
else:
    process.GlobalTag.globaltag = cms.string('START53_LV6A1::All')
    print "\nGlobalTag : START53_LV6A1::All\n"


# do trigger matching for muons
triggerProcessName = 'HLT'


if (run2012not2011 == True):
    process.cleanMuonTriggerMatchHLT  = cms.EDProducer(
            # match by DeltaR only (best match by DeltaR)
            'PATTriggerMatcherDRLessByR',
                src                   = cms.InputTag('cleanPatMuons'),
                # default producer label as defined in
                # PhysicsTools/PatAlgos/python/triggerLayer1/triggerProducer_cfi.py
                matched               = cms.InputTag('patTrigger'),
                matchedCuts           = cms.string('path("HLT_DoubleMu3p5_LowMass_Displaced*",0,0)'),
                maxDeltaR             = cms.double(0.1),
                # only one match per trigger object
                resolveAmbiguities    = cms.bool(True),
                # take best match found per reco object (by DeltaR here, see above)
                resolveByMatchQuality = cms.bool(False))
else:
    process.cleanMuonTriggerMatchHLT0 = cms.EDProducer(
            # match by DeltaR only (best match by DeltaR)
            'PATTriggerMatcherDRLessByR',
                src = cms.InputTag('cleanPatMuons'),
                # default producer label as defined in
                # PhysicsTools/PatAlgos/python/triggerLayer1/triggerProducer_cfi.py
                matched = cms.InputTag('patTrigger'),
                matchedCuts = cms.string('path("HLT_Dimuon6p5_LowMass_Displaced_v*",0,0)'),
                maxDeltaR = cms.double(0.1),
                # only one match per trigger object
                resolveAmbiguities = cms.bool(True),
                # take best match found per reco object (by DeltaR here, see above)
                resolveByMatchQuality = cms.bool(False))

    process.cleanMuonTriggerMatchHLT1 = cms.EDProducer(
             'PATTriggerMatcherDRLessByR',
                 src = cms.InputTag('cleanPatMuons'),
                 matched = cms.InputTag('patTrigger'),
                 matchedCuts = cms.string('path("HLT_Dimuon7_LowMass_Displaced_v*")'),
                 maxDeltaR = cms.double(0.1),
                 resolveAmbiguities = cms.bool(True),
                 resolveByMatchQuality = cms.bool(False))

    process.cleanMuonTriggerMatchHLT2 = cms.EDProducer(
             'PATTriggerMatcherDRLessByR',
                 src = cms.InputTag('cleanPatMuons'),
                 matched = cms.InputTag('patTrigger'),
                 matchedCuts = cms.string('path("HLT_DoubleMu4_LowMass_Displaced_v*")'),
                 maxDeltaR = cms.double(0.1),
                 resolveAmbiguities = cms.bool(True),
                 resolveByMatchQuality = cms.bool(False))

    process.cleanMuonTriggerMatchHLT3 = cms.EDProducer(
              'PATTriggerMatcherDRLessByR',
                  src = cms.InputTag('cleanPatMuons'),
                  matched = cms.InputTag('patTrigger'),
                  matchedCuts = cms.string('path("HLT_DoubleMu4p5_LowMass_Displaced_v*")'),
                  maxDeltaR = cms.double(0.1),
                  resolveAmbiguities = cms.bool(True),
                  resolveByMatchQuality = cms.bool(False))

    process.cleanMuonTriggerMatchHLT4 = cms.EDProducer(
              'PATTriggerMatcherDRLessByR',
                  src = cms.InputTag('cleanPatMuons'),
                  matched = cms.InputTag('patTrigger'),
                  matchedCuts = cms.string('path("HLT_DoubleMu5_LowMass_Displaced_v*")'),
                  maxDeltaR = cms.double(0.1),
                  resolveAmbiguities = cms.bool(True),
                  resolveByMatchQuality = cms.bool(False))

      
      
 
 
             
from PhysicsTools.PatAlgos.tools.trigTools import *

if (run2012not2011 == True):
    switchOnTriggerMatchEmbedding(process, triggerMatchers = ['cleanMuonTriggerMatchHLT'], hltProcess = triggerProcessName, outputModule = '')

    g_TriggerNames_LastFilterNames = [
        ('HLT_DoubleMu3p5_LowMass_Displaced',  'hltDisplacedmumuFilterDoubleMu3p5LowMass')  #5E32, v8.1-v8.3
        ]

else:
    switchOnTriggerMatchEmbedding(process, triggerMatchers = ['cleanMuonTriggerMatchHLT0','cleanMuonTriggerMatchHLT1','cleanMuonTriggerMatchHLT2','cleanMuonTriggerMatchHLT3','cleanMuonTriggerMatchHLT4'], hltProcess = triggerProcessName, outputModule = '')

    g_TriggerNames_LastFilterNames = [
             ('HLT_Dimuon6p5_LowMass_Displaced',
              'hltDisplacedmumuFilterLowMass'), #5E32

             ('HLT_Dimuon7_LowMass_Displaced',
              'hltDisplacedmumuFilterLowMass'), #1E33, 1.4E33

             ('HLT_DoubleMu4_LowMass_Displaced',
              'hltDisplacedmumuFilterLowMass'), #2E33

             ('HLT_DoubleMu4p5_LowMass_Displaced',
              'hltDisplacedmumuFilterDoubleMu4p5LowMass'), #3E33, 5E33

             ('HLT_DoubleMu5_LowMass_Displaced',
              'hltDisplacedmumuFilterDoubleMu5LowMass'), #3E33, 5E33
             ]


    
g_TriggerNames = [i[0] for i in g_TriggerNames_LastFilterNames]
g_LastFilterNames = [i[1] for i in g_TriggerNames_LastFilterNames]


process.ntuple.TriggerNames = cms.vstring(g_TriggerNames)
process.ntuple.LastFilterNames = cms.vstring(g_LastFilterNames)
process.ntuple.IsMonteCarlo = cms.untracked.bool(True)
process.ntuple.KeepGENOnly = cms.untracked.bool(False)
