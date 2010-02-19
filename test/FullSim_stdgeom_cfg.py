import FWCore.ParameterSet.Config as cms

process = cms.Process("Fullsim")

process.load("Configuration.StandardSequences.Services_cff")

process.load("Configuration.StandardSequences.Geometry_cff")

process.load("Configuration.StandardSequences.MagneticField_cff")

#process.load("Configuration.StandardSequences.FakeConditions_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'MC_31X_V8::All'
#process.GlobalTag.globaltag = 'MC_31X_V9::All'
#process.GlobalTag.globaltag = 'STARTUP3X_V8L::All'

process.siPixelFakeGainOfflineESSource = cms.ESSource("SiPixelFakeGainOfflineESSource",
    file = cms.FileInPath('SLHCUpgradeSimulations/Geometry/data/stdgeom/PixelSkimmedGeometry_stdgeom.txt')
)
process.es_prefer_fake_gain = cms.ESPrefer("SiPixelFakeGainOfflineESSource","siPixelFakeGainOfflineESSource")

process.siPixelFakeLorentzAngleESSource = cms.ESSource("SiPixelFakeLorentzAngleESSource",
    file = cms.FileInPath('SLHCUpgradeSimulations/Geometry/data/stdgeom/PixelSkimmedGeometry_stdgeom.txt')
)
process.es_prefer_fake_lorentz = cms.ESPrefer("SiPixelFakeLorentzAngleESSource","siPixelFakeLorentzAngleESSource")

process.load("CalibTracker.SiStripESProducers.fake.SiStripNoisesFakeESSource_cfi")
process.SiStripNoisesGenerator.NoiseStripLengthSlope=51. #dec mode
process.SiStripNoisesGenerator.NoiseStripLengthQuote=630.

process.siStripNoisesFakeESSource  = cms.ESSource("SiStripNoisesFakeESSource")
process.es_prefer_fake_strip_noise = cms.ESPrefer("SiStripNoisesFakeESSource",
                                                  "siStripNoisesFakeESSource")

process.load("CalibTracker.SiStripESProducers.fake.SiStripQualityFakeESSource_cfi")

process.siStripQualityFakeESSource  = cms.ESSource("SiStripQualityFakeESSource")
process.es_prefer_fake_strip_quality = cms.ESPrefer("SiStripQualityFakeESSource",
                                                     "siStripQualityFakeESSource")

process.load("CalibTracker.SiStripESProducers.fake.SiStripPedestalsFakeESSource_cfi")

process.siStripPedestalsFakeESSource  = cms.ESSource("SiStripPedestalsFakeESSource")
process.es_prefer_fake_strip_pedestal = cms.ESPrefer("SiStripPedestalsFakeESSource",
                                                     "siStripPedestalsFakeESSource")

process.load("CalibTracker.SiStripESProducers.fake.SiStripLorentzAngleFakeESSource_cfi")

process.siStripLorentzAngleFakeESSource  = cms.ESSource("SiStripLorentzAngleFakeESSource")
process.es_prefer_fake_strip_LA = cms.ESPrefer("SiStripLorentzAngleFakeESSource",
                                               "siStripLorentzAngleFakeESSource")

process.siStripLorentzAngleSimFakeESSource  = cms.ESSource("SiStripLorentzAngleSimFakeESSource")
process.es_prefer_fake_strip_LA_sim = cms.ESPrefer("SiStripLorentzAngleSimFakeESSource",
                                                   "siStripLorentzAngleSimFakeESSource")

process.load("CalibTracker.SiStripESProducers.fake.SiStripApvGainFakeESSource_cfi")
process.SiStripApvGainGenerator.MeanGain=1.0
process.SiStripApvGainGenerator.SigmaGain=0.0
process.SiStripApvGainGenerator.genMode = cms.string("default")

process.myStripApvGainFakeESSource = cms.ESSource("SiStripApvGainFakeESSource")
process.es_prefer_myStripApvGainFakeESSource  = cms.ESPrefer("SiStripApvGainFakeESSource",
                                                  "myStripApvGainFakeESSource")

process.myStripApvGainSimFakeESSource  = cms.ESSource("SiStripApvGainSimFakeESSource")
process.es_prefer_myStripApvGainSimFakeESSource = cms.ESPrefer("SiStripApvGainSimFakeESSource",
                                                               "myStripApvGainSimFakeESSource")

process.load("CalibTracker.SiStripESProducers.fake.SiStripThresholdFakeESSource_cfi")

process.siStripThresholdFakeESSource  = cms.ESSource("SiStripThresholdFakeESSource")
process.es_prefer_fake_strip_threshold = cms.ESPrefer("SiStripThresholdFakeESSource",
                                                     "siStripThresholdFakeESSource")

process.TrackerDigiGeometryESModule.applyAlignment = False

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
process.simSiPixelDigis.LorentzAngle_DB = False
process.simSiPixelDigis.killModules = False
process.simSiPixelDigis.useDB = False
process.simSiPixelDigis.DeadModules_DB = False

process.simSiPixelDigis.AddPixelInefficiency = -1

process.siPixelClusters.src = 'simSiPixelDigis'
process.siPixelClusters.MissCalibrate = False
#process.siStripZeroSuppression.RawDigiProducersList[0].RawDigiProducer = 'simSiStripDigis'
#process.siStripZeroSuppression.RawDigiProducersList[1].RawDigiProducer = 'simSiStripDigis'
#process.siStripZeroSuppression.RawDigiProducersList[2].RawDigiProducer = 'simSiStripDigis'
#process.siStripClusters.DigiProducersList[0].DigiProducer= 'simSiStripDigis'
process.siStripZeroSuppression.RawDigiProducersList = cms.VInputTag( cms.InputTag('simSiStripDigis','VirginRaw'), 
                                                                     cms.InputTag('simSiStripDigis','ProcessedRaw'),
                                                                     cms.InputTag('simSiStripDigis','ScopeMode'))
process.siStripClusters.DigiProducersList = cms.VInputTag(cms.InputTag('simSiStripDigis','ZeroSuppressed'),
                                                          cms.InputTag('siStripZeroSuppression','VirginRaw'),
                                                          cms.InputTag('siStripZeroSuppression','ProcessedRaw'),
                                                          cms.InputTag('siStripZeroSuppression','ScopeMode'))
## more needed for strips from RecoTracker/MeasurementDet/python/MeasurementTrackerESProducer_cfi.py
#process.MeasurementTracker.inactiveStripDetectorLabels = cms.VInputTag('siStripDigis')
#process.MeasurementTracker.inactiveStripDetectorLabels = cms.VInputTag('')
#process.MeasurementTracker.stripClusterProducer=cms.string('')
## needed to avoid RAW2DIGI and DIGI2Raw
#process.MeasurementTracker.inactiveStripDetectorLabels = cms.VInputTag()
#
#process.MeasurementTracker.UseStripModuleQualityDB     = cms.bool(False)
#process.MeasurementTracker.UseStripAPVFiberQualityDB   = cms.bool(False)
#process.MeasurementTracker.MaskBadAPVFibers            = cms.bool(False)
#process.MeasurementTracker.UseStripStripQualityDB      = cms.bool(False)
#process.MeasurementTracker.UsePixelModuleQualityDB   = cms.bool(False)
#process.MeasurementTracker.UsePixelROCQualityDB      = cms.bool(False)
#process.MeasurementTracker.PixelCPE = cms.string('PixelCPEGeneric')
##Ivan's version
#process.MeasurementTracker.stripClusterProducer=cms.string('')
process.MeasurementTracker.inactiveStripDetectorLabels = cms.VInputTag()
process.MeasurementTracker.UseStripModuleQualityDB     = cms.bool(False)
process.MeasurementTracker.UseStripAPVFiberQualityDB   = cms.bool(False)


# Event output
process.load("Configuration.EventContent.EventContent_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

process.load("FastSimulation/Configuration/FlatPtMuonGun_cfi")
#process.FlatRandomPtGunSource.PGunParameters.PartID[0] = 13
process.generator.PGunParameters.MinPt = 0.9
process.generator.PGunParameters.MaxPt = 50.0
process.generator.PGunParameters.MinEta = -2.4
process.generator.PGunParameters.MaxEta = 2.4
process.generator.Verbosity = 1
process.generator.AddAntiParticle = True

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
process.multiTrackValidator.label = ['ctfWithMaterialTracks']
#process.multiTrackValidator.label = ['cutsRecoTracks']
#process.multiTrackValidator.label_tp_effic = cms.InputTag("cutsTPEffic")
#process.multiTrackValidator.label_tp_fake = cms.InputTag("cutsTPFake")
process.multiTrackValidator.associators = ['TrackAssociatorByHits']
process.multiTrackValidator.UseAssociators = True
process.multiTrackValidator.outputFile = "validfullstdg_muon_50GeV.root"
process.multiTrackValidator.nint = cms.int32(20)
process.multiTrackValidator.nintpT = cms.int32(25)
process.multiTrackValidator.maxpT = cms.double(50.0)
process.multiTrackValidator.skipHistoFit = False
##### with John's changes ##############################
process.load("SLHCUpgradeSimulations.Geometry.oldTracking_wtriplets")
# restrict vertex fining in trackingtruthprod to smaller volume (note: these numbers in mm)
process.mergedtruth.volumeRadius = cms.double(100.0)
process.mergedtruth.volumeZ = cms.double(900.0)
process.mergedtruth.discardOutVolume = cms.bool(True)

process.cutsTPFake.tip = cms.double(10.0)
process.cutsTPFake.lip = cms.double(90.0)
#NB: tracks are already filtered by the generalTracks sequence
#for additional cuts use the cutsRecoTracks filter:
#process.load("Validation.RecoTrack.cutsRecoTracks_cfi")
#process.cutsRecoTracks.src = cms.InputTag("ctfWithMaterialTracks")
#process.cutsRecoTracks.quality = cms.string('')
#process.cutsRecoTracks.minHit = cms.int32(3)
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

process.anal = cms.EDAnalyzer("EventContentAnalyzer")

# need validation packages

process.p0 = cms.Path(process.generator+process.pgen)
process.p1 = cms.Path(process.psim)
process.p2 = cms.Path(process.pdigi)
process.p3 = cms.Path(process.L1Emulator)
process.p4 = cms.Path(process.DigiToRaw)
process.p5 = cms.Path(process.RawToDigi)
process.p6 = cms.Path(process.trackerlocalreco)
process.p7 = cms.Path(process.offlineBeamSpot+process.oldTracking_wtriplets)
#process.p7 = cms.Path(process.offlineBeamSpot)
#process.p7a = cms.Path(process.pixelTriplets)
#process.p7b = cms.Path(process.ckfTrackCandidates)
#process.p7c = cms.Path(process.ctfWithMaterialTracks)
#process.p7 = cms.Path(process.reconstruction)
#process.p8 = cms.Path(process.cutsTPEffic*process.cutsTPFake*process.cutsRecoTracks*process.multiTrackValidator)
process.p8 = cms.Path(process.cutsTPEffic*process.cutsTPFake*process.multiTrackValidator)
#process.p9 = cms.Path(process.anal*process.ReadLocalMeasurement)
process.p9 = cms.Path(process.ReadLocalMeasurement)
process.outpath = cms.EndPath(process.FEVT)
#process.schedule = cms.Schedule(process.p0,process.p1,process.p2,process.p3,process.p4,process.p5,process.p6,process.p7,process.outpath)
process.schedule = cms.Schedule(process.p0,process.p1,process.p2,process.p3,process.p6,process.p7,process.p8,process.p9)
