import FWCore.ParameterSet.Config as cms
import os
import sys

#arglist = sys.argv
#filename = arglist[2]
#filename='MC_MinBias1000to1400_2018.txt'
#basename = os.path.splitext(filename)[0]

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.MagneticField_cff') # B-field map
process.load('Configuration.Geometry.GeometryRecoDB_cff') # Ideal geometry and interface
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cff") # Beamspot

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff") # Global tag
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag,'105X_upgrade2018_realistic_v4','')
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')

process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100))

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

#For processing files in txt files in the test directory
#import FWCore.Utilities.FileUtils as FileUtils
#readFiles = cms.untracked.vstring()
#readFiles = cms.untracked.vstring( FileUtils.loadListFromFile (os.environ['CMSSW_BASE']+'/src/DatasetValidation/DatasetValidationTool/test/'+str(filename)) )

#For proceesing single or less number of root files
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#         'root://cms-xrd-global.cern.ch//store/mc/RunIIWinter19CosmicDR/TKCosmics_0T/ALCARECO/TkAlCosmics0T-0T_103X_upgrade2018cosmics_realistic_deco_v7-v3/110000/18E5DA07-134D-B247-BFE0-1E3BB01D48EA.root'
#'root://cms-xrd-global.cern.ch//store/mc/RunIIWinter19PFCalibDRPremix/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/ALCARECO/TkAlMuonIsolated-2016Conditions_newPixelConditions_105X_mcRun2_asymptotic_newPixCond_v2-v1/00000/05601FA6-D476-9049-8B3C-D26AEE830054.root'
#process.source = cms.Source("PoolSource", fileNames = readFiles
#    fileNames = cms.untracked.vstring(
        'root://cms-xrd-global.cern.ch//store/data/Run2016A/Cosmics/ALCARECO/TkAlCosmics0T-07Dec2018-v1/20000/783C9E32-480B-E911-92B2-20040FE8ECAC.root',
#        'root://cms-xrd-global.cern.ch//store/data/Run2016H/DoubleMuon/ALCARECO/TkAlZMuMu-07Aug17-v1/90000/FA413FBB-A599-E711-B7AC-0CC47A7C3472.root'
#        'root://cms-xrd-global.cern.ch//store/data/Run2016H/SingleMuon/ALCARECO/TkAlMuonIsolated-07Aug17-v1/90000/FCB1CBE1-0790-E711-8359-3417EBE2F316.root'
#        'root://cms-xrd-global.cern.ch//store/data/Run2016E/ZeroBias/ALCARECO/TkAlMinBias-07Aug17-v1/90000/FEA7EA10-338D-E711-8AB6-003048FFD75A.root'
         'root://cms-xrd-global.cern.ch//store/data/Run2016A/Cosmics/ALCARECO/TkAlCosmics0T-07Dec2018-v1/20000/783C9E32-480B-E911-92B2-20040FE8ECAC.root',
         'root://cms-xrd-global.cern.ch//store/data/Run2016D/Cosmics/ALCARECO/TkAlCosmics0T-07Dec2018-v1/20000/8A715B05-EF0A-E911-B731-782BCB20ED64.root',
#         'root://cms-xrd-global.cern.ch//store/data/Run2016C/Cosmics/ALCARECO/TkAlCosmics0T-07Dec2018-v1/20000/F443086A-3C0B-E911-AA8F-14187741013C.root',
#         'root://cms-xrd-global.cern.ch//store/data/Run2016C/Cosmics/ALCARECO/TkAlCosmics0T-07Dec2018-v1/20000/10E8FF4A-530B-E911-9150-549F3525BBCC.root'
    )
)

process.demo = cms.EDAnalyzer('DatasetValidationTool', # Use the name of the C file or class name here
#     tracks = cms.InputTag("ALCARECOTkAlCosmicsCTF0T"),
      BS = cms.InputTag("offlineBeamSpot"),
      tracks = cms.InputTag("FinalTrackRefitter"),
      vertices= cms.InputTag("offlinePrimaryVertices"),
      IsResonance=cms.bool(False)
#     tracks = cms.InputTag("ALCARECOTkAlMuonIsolated")
#     tracks = cms.InputTag("ALCARECOTkAlMinBias")
)

process.TFileService = cms.Service("TFileService",
     fileName = cms.string('Test_Reference.root')
)

#process.p = cms.Path(process.demo)
process.p = cms.Path(process.seqTrackselRefit*process.demo)
