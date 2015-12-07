
##################
# cmssw configs
################## 

import FWCore.ParameterSet.Config as cms
from btopimumu_cfi import process

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(

#     '/store/user/nsahoo/PYTHIA6_Bu2MuMuPi_GenOnly_8TeV/crab_edm_BToPiMuMu_MC_GenOnly_8TeV_v2/151207_114850/0000/PYTHIA6_Bu2MuMuPi_TuneZ2star_8TeV_cff_py_GEN_NoFilter_1.root',
        'file:PYTHIA6_Bu2MuMuPi_TuneZ2star_8TeV_cff_py_GEN_NoFilter.root',
                            )
                        )


process.GlobalTag.globaltag = cms.string('START53_V7G::All')
print "\nGlobalTag : START53_V7G::All\n"


process.ntuple.IsMonteCarlo = cms.untracked.bool(True)
process.ntuple.KeepGENOnly = cms.untracked.bool(True)
process.p = cms.Path(process.ntuple)
