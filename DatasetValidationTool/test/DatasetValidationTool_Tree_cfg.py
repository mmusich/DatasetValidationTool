import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")


process.load("FWCore.MessageService.MessageLogger_cfi")
#process.load('Configuration.StandardSequences.MagneticField_cff') # B-field map
process.load('Configuration.Geometry.GeometryRecoDB_cff') # Ideal geometry and interface
#process.load("Configuration.Geometry.GeometryDB_cff")
#process.load("RecoVertex.BeamSpotProducer.BeamSpot_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff") # Global tag
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag.globaltag = GlobalTag('102X_dataRun2_Prompt_v1')
process.GlobalTag = GlobalTag(process.GlobalTag,'auto:run2_data','')
#process.GlobalTag = GlobalTag(process.GlobalTag,'80X_dataRun2_2016LegacyRepro_v4 Pset','')

process.MessageLogger.cerr.FwkReport.reportEvery = 500
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'root://cms-xrd-global.cern.ch//store/data/Run2016A/Cosmics/ALCARECO/TkAlCosmics0T-07Dec2018-v1/20000/783C9E32-480B-E911-92B2-20040FE8ECAC.root'
    #     'root://cms-xrd-global.cern.ch//store/data/Run2016H/DoubleMuon/ALCARECO/TkAlZMuMu-07Aug17-v1/90000/FA413FBB-A599-E711-B7AC-0CC47A7C3472.root'
    #     'root://cms-xrd-global.cern.ch//store/data/Run2016H/SingleMuon/ALCARECO/TkAlMuonIsolated-07Aug17-v1/90000/FCB1CBE1-0790-E711-8359-3417EBE2F316.root'
    #     'root://cms-xrd-global.cern.ch//store/data/Run2016E/ZeroBias/ALCARECO/TkAlMinBias-07Aug17-v1/90000/FEA7EA10-338D-E711-8AB6-003048FFD75A.root'

    )
)

process.demo = cms.EDAnalyzer('DatasetValidationTool_Tree',
     tracks = cms.InputTag("ALCARECOTkAlCosmicsCTF0T")
#     tracks = cms.InputTag("ALCARECOTkAlZMuMu")
#     tracks = cms.InputTag("ALCARECOTkAlMuonIsolated")
#     tracks = cms.InputTag("ALCARECOTkAlMinBias")
)

process.TFileService = cms.Service("TFileService",
     fileName = cms.string('outputfileTree_Res1.root')
)


process.p = cms.Path(process.demo)
