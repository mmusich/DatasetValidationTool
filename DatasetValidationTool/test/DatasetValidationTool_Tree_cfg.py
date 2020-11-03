import FWCore.ParameterSet.Config as cms
import os
import sys

#arglist = sys.argv
#filename = arglist[2]
#filename='JPsimumuPromptRun2018Bv1.txt'
#basename = os.path.splitext(filename)[0]

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.MagneticField_cff') # B-field map
process.load('Configuration.Geometry.GeometryRecoDB_cff') # Ideal geometry and interface
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff") # Global tag
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag,'auto:run2_data','')

process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000))

import Alignment.CommonAlignment.tools.trackselectionRefitting as trackselRefit
process.seqTrackselRefit = trackselRefit.getSequence(process, 'ALCARECOTkAlCosmicsCTF0T',
                                                     isPVValidation=False, 
                                                     TTRHBuilder='WithAngleAndTemplate',
                                                     usePixelQualityFlag=True,
                                                     openMassWindow=False,
                                                     cosmicsDecoMode=True,
                                                     cosmicsZeroTesla=False,
                                                     momentumConstraint=None,
                                                     cosmicTrackSplitting=False,
                                                     use_d0cut=False)

#import FWCore.Utilities.FileUtils as FileUtils
#readFiles = cms.untracked.vstring()
#readFiles = cms.untracked.vstring( FileUtils.loadListFromFile (os.environ['CMSSW_BASE']+'/src/DatasetValidation/DatasetValidationTool/test/'+str(filename)) )
process.source = cms.Source("PoolSource",#fileNames = readFiles,
    fileNames = cms.untracked.vstring(
        'root://cms-xrd-global.cern.ch//store/data/Run2016A/Cosmics/ALCARECO/TkAlCosmics0T-07Dec2018-v1/20000/783C9E32-480B-E911-92B2-20040FE8ECAC.root'
#        'root://cms-xrd-global.cern.ch//store/data/Run2016H/DoubleMuon/ALCARECO/TkAlZMuMu-07Aug17-v1/90000/FA413FBB-A599-E711-B7AC-0CC47A7C3472.root'
#        'root://cms-xrd-global.cern.ch//store/data/Run2016H/SingleMuon/ALCARECO/TkAlMuonIsolated-07Aug17-v1/90000/FCB1CBE1-0790-E711-8359-3417EBE2F316.root'
#        'root://cms-xrd-global.cern.ch//store/data/Run2016E/ZeroBias/ALCARECO/TkAlMinBias-07Aug17-v1/90000/FEA7EA10-338D-E711-8AB6-003048FFD75A.root'
    )
)

process.demo = cms.EDAnalyzer('DatasetValidationTool_Tree',
#     tracks = cms.InputTag("ALCARECOTkAlCosmicsCTF0T"),
      BS = cms.InputTag("offlineBeamSpot"),
      tracks = cms.InputTag("FinalTrackRefitter"),
      vertices= cms.InputTag("offlinePrimaryVertices"),
      IsResonance=cms.bool(True)
#     tracks = cms.InputTag("ALCARECOTkAlMuonIsolated")
#     tracks = cms.InputTag("ALCARECOTkAlMinBias")
)

process.TFileService = cms.Service("TFileService",
     fileName = cms.string('outputfileTree_JPsimumu2.root')
)


process.p = cms.Path(process.seqTrackselRefit*process.demo)
