import FWCore.ParameterSet.Config as cms

process = cms.Process("Fullsim")

process.load("Configuration.StandardSequences.Services_cff")

process.load("Configuration.StandardSequences.Geometry_cff")

process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("Configuration.StandardSequences.FakeConditions_cff")

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.destinations = cms.untracked.vstring("detailedInfo_fullstdgeommu50")

# this config frament brings you the generator information
process.load("Configuration.StandardSequences.Generator_cff")

# this config frament brings you 3 steps of the detector simulation:
# -- vertex smearing (IR modeling)
# -- G4-based hit level detector simulation
# -- digitization (electronics readout modeling)
# it returns 2 sequences : 
# -- psim (vtx smearing + G4 sim)
# -- pdigi (digitization in all subsystems, i.e. tracker=pix+sistrips,
#           cal=ecal+ecal-0-suppression+hcal), muon=csc+dt+rpc)
#
process.load("Configuration.StandardSequences.Simulation_cff")

# please note the IMPORTANT: 
# in order to operate Digis, one needs to include Mixing module 
# (pileup modeling), at least in the 0-pileup mode
#
# There're 3 possible configurations of the Mixing module :
# no-pileup, low luminosity pileup, and high luminosity pileup
#
# they come, respectively, through the 3 config fragments below
#
# *each* config returns label "mix"; thus you canNOT have them
# all together in the same configuration, but only one !!!
#
process.load("Configuration.StandardSequences.MixingNoPileUp_cff")

#include "Configuration/StandardSequences/data/MixingLowLumiPileUp.cff" 
#include "Configuration/StandardSequences/data/MixingHighLumiPileUp.cff" 
process.load("Configuration.StandardSequences.L1Emulator_cff")

process.load("Configuration.StandardSequences.DigiToRaw_cff")

process.load("Configuration.StandardSequences.RawToDigi_cff")

#process.load("Configuration.StandardSequences.VtxSmearedBetafuncEarlyCollision_cff")
process.load("Configuration.StandardSequences.VtxSmearedGauss_cff")

process.load("Configuration.StandardSequences.Reconstruction_cff")

process.load("SimTracker.Configuration.SimTracker_cff")
process.simSiPixelDigis.MissCalibrate = False
process.simSiPixelDigis.AddPixelInefficiency = -1
process.simSiPixelDigis.LorentzAngle_DB = False
process.simSiPixelDigis.killModules = False

process.siPixelClusters.src = 'simSiPixelDigis'
process.siPixelClusters.MissCalibrate = False
process.siStripZeroSuppression.RawDigiProducersList[0].RawDigiProducer = 'simSiStripDigis'
process.siStripZeroSuppression.RawDigiProducersList[1].RawDigiProducer = 'simSiStripDigis'
process.siStripZeroSuppression.RawDigiProducersList[2].RawDigiProducer = 'simSiStripDigis'
process.siStripClusters.DigiProducersList[0].DigiProducer= 'simSiStripDigis'

# Event output
process.load("Configuration.EventContent.EventContent_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

process.load("FastSimulation/Configuration/FlatPtMuonGun_cfi")
process.FlatRandomPtGunSource.PGunParameters.PartID[0] = 13
process.FlatRandomPtGunSource.PGunParameters.MinPt = 0.9
process.FlatRandomPtGunSource.PGunParameters.MaxPt = 50.0
process.FlatRandomPtGunSource.PGunParameters.MinEta = -2.4
process.FlatRandomPtGunSource.PGunParameters.MaxEta = 2.4
process.FlatRandomPtGunSource.AddAntiParticle = cms.untracked.bool(True)

process.FEVT = cms.OutputModule("PoolOutputModule",
    process.FEVTSIMEventContent,
    fileName = cms.untracked.string('/uscms_data/d2/cheung/slhc/testfullstdg_muon_50GeV.root')
)

process.load("Validation.RecoTrack.cutsTPEffic_cfi")
process.load("Validation.RecoTrack.cutsTPFake_cfi")

process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")
process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")

process.load("Validation.RecoTrack.MultiTrackValidator_cff")
#process.multiTrackValidator.label = ['generalTracks']
### if using simple (non-iterative) or old (as in 1_8_4) tracking
#process.multiTrackValidator.label = ['ctfWithMaterialTracks']
process.multiTrackValidator.label = ['cutsRecoTracks']
process.multiTrackValidator.label_tp_effic = cms.InputTag("cutsTPEffic")
process.multiTrackValidator.label_tp_fake = cms.InputTag("cutsTPFake")
process.multiTrackValidator.associators = ['TrackAssociatorByHits']
process.multiTrackValidator.UseAssociators = True
process.multiTrackValidator.outputFile = "validfullstdg_muon_50GeV.root"
process.multiTrackValidator.nint = cms.int32(20)
process.multiTrackValidator.nintpT = cms.int32(25)
process.multiTrackValidator.maxpT = cms.double(50.0)

##### with John's changes ##############################
process.load("SLHCUpgradeSimulations.Geometry.oldTracking_wtriplets")
# restrict vertex fining in trackingtruthprod to smaller volume (note: these numbers in mm)
process.mergedtruth.volumeRadius = cms.double(100.0)
process.mergedtruth.volumeZ = cms.double(900.0)
process.mergedtruth.discardOutVolume = cms.bool(True)

#process.cutsTPEffic.ptMin = cms.double(2.5)
#process.cutsTPFake.ptMin = cms.double(2.0)
process.cutsTPFake.tip = cms.double(10.0)
process.cutsTPFake.lip = cms.double(90.0)
#NB: tracks are already filtered by the generalTracks sequence
#for additional cuts use the cutsRecoTracks filter:
process.load("Validation.RecoTrack.cutsRecoTracks_cfi")
process.cutsRecoTracks.src = cms.InputTag("ctfWithMaterialTracks")
process.cutsRecoTracks.quality = cms.string('')
process.cutsRecoTracks.minHit = cms.int32(3)
#process.cutsRecoTracks.minHit = cms.int32(8)
#process.cutsRecoTracks.minHit = cms.int32(6)
############ end John's changes ###########################

### produce an ntuple with pixel hits for analysis
process.ReadLocalMeasurement = cms.EDAnalyzer("StdHitNtuplizer",
   src = cms.InputTag("siPixelRecHits"),
   stereoRecHits = cms.InputTag("siStripMatchedRecHits","stereoRecHit"),
   rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
   matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
   #trackProducer = cms.InputTag("generalTracks"),
   ### if using simple (non-iterative) or old (as in 1_8_4) tracking
   trackProducer = cms.InputTag("ctfWithMaterialTracks"),
   OutputFile = cms.string("stdgrechitfullstdg_ntuple.root"),
   ### for using track hit association
   associatePixel = cms.bool(True),
   associateStrip = cms.bool(False),
   associateRecoTracks = cms.bool(False),
   ROUList = cms.vstring('g4SimHitsTrackerHitsPixelBarrelLowTof',
                         'g4SimHitsTrackerHitsPixelBarrelHighTof',
                         'g4SimHitsTrackerHitsPixelEndcapLowTof',
                         'g4SimHitsTrackerHitsPixelEndcapHighTof')
)

process.Timing =  cms.Service("Timing")
process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
     oncePerEventMode = cms.untracked.bool(False),
     showMallocInfo = cms.untracked.bool(False),
     moduleMemorySummary = cms.untracked.bool(True),
     ignoreTotal = cms.untracked.int32(5)
)

# need validation packages

process.p0 = cms.Path(process.pgen)
process.p1 = cms.Path(process.psim)
process.p2 = cms.Path(process.pdigi)
process.p3 = cms.Path(process.L1Emulator)
#process.p4 = cms.Path(process.DigiToRaw)
#process.p5 = cms.Path(process.RawToDigi)
process.p6 = cms.Path(process.trackerlocalreco)
process.p7 = cms.Path(process.offlineBeamSpot+process.oldTracking_wtriplets)
#process.p7 = cms.Path(process.reconstruction)
process.p8 = cms.Path(process.cutsTPEffic*process.cutsTPFake*process.cutsRecoTracks*process.multiTrackValidator)
process.p9 = cms.Path(process.ReadLocalMeasurement)
process.outpath = cms.EndPath(process.FEVT)
#process.schedule = cms.Schedule(process.p0,process.p1,process.p2,process.p3,process.p4,process.p5,process.p6,process.p7,process.outpath)
process.schedule = cms.Schedule(process.p0,process.p1,process.p2,process.p3,process.p6,process.p7,process.p8)
