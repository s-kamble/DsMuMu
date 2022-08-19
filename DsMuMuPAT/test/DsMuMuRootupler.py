import FWCore.ParameterSet.Config as cms
process = cms.Process("Rootuple")

#-- GEOMETRY + B-FIELD --#
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')
process.load('Configuration.EventContent.EventContent_cff')

#-- GLOBAL TAG --#
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data')
#process.GlobalTag.globaltag = cms.string('102X_dataRun2_v12')
#process.GlobalTag.globaltag = cms.string('106X_upgrade2018_realistic_v4') #For  Private MC
#process.GlobalTag.globaltag = cms.string('102X_upgrade2018_realistic_v21')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)

#-- NUMBER OF EVENTS --#
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10000))
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

#-- SOURCE FILES --#	
process.source = cms.Source("PoolSource",
    			     fileNames = cms.untracked.vstring(
    								'/store/data/Run2018B/DoubleMuonLowMass/MINIAOD/17Sep2018-v1/60000/6EF2991A-8E4C-C94A-BE5E-3E5D0A17DE08.root',
#    								'/store/data/Run2016F/DoubleMuonLowMass/MINIAOD/17Jul2018-v1/40000/10AEEFA3-7D8C-E811-9F2E-0CC47A544E12.root',     
#    								'/store/data/Run2016C/Charmonium/MINIAOD/17Jul2018-v1/20000/9C03CBE2-4B8B-E811-9299-0CC47AC17678.root',
#     								'/store/data/Run2018C/Charmonium/MINIAOD/PromptReco-v2/000/319/756/00000/EEF6CEC1-698B-E811-8081-02163E00AF5F.root',
#'file:/eos/user/s/skamble/Bc2DsMuMu-MC-samples/BcToDsMuMu_MC_MINIAOD/BcToDsMuMu_MC_MINIAOD_5.root',
#'file:/BctoDsMuMu_MC_MINIAOD_5.root'


								)
)

process.load("slimmedMuonsTriggerMatcher_cfi")

#process.load("myAnalyzers.DsMuMuPAT.DsMuMuRootupler_cfi")

process.rootuple = cms.EDAnalyzer('DsMuMu',
                          dimuons = cms.InputTag("slimmedMuons"),
                          Track = cms.InputTag("packedPFCandidates"),        
                          GenParticles = cms.InputTag("prunedGenParticles"),
                          packedGenParticles = cms.InputTag("packedGenParticles"),
                          primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          secondaryVerticesPtr = cms.InputTag("slimmedKshortVertices"),
                         
                          bslabel = cms.InputTag("offlineBeamSpot"),
                          TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                          OnlyBest = cms.bool(False),                         
                          isMC = cms.bool(True), # True for MC, false for data 
                          OnlyGen = cms.bool(False),
                       
                          MuonMass = cms.untracked.double(0.10565837),   #pdg mass
                          MuonMassErr = cms.untracked.double(3.5e-9),
                          KaonMass = cms.untracked.double(0.493677),
                          KaonMassErr = cms.untracked.double(1.6e-5),
                          PionMass = cms.untracked.double(0.13957018),
                          PionMassErr = cms.untracked.double(3.5e-7),
                          MuonMinPt = cms.untracked.double(4.0),
                          MuonMaxEta = cms.untracked.double(2.4),



)

process.TFileService = cms.Service("TFileService",
				    fileName = cms.string('DsMuMu_Rootuple_test_DM18B10.root')
)

process.p = cms.Path(process.rootuple)

