import FWCore.ParameterSet.Config as cms

process = cms.Process("Fastsimwdigi")

# Number of events to be generated
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

# Include the RandomNumberGeneratorService definition
process.load("FastSimulation.Configuration.RandomServiceInitialization_cff")
process.RandomNumberGeneratorService.simSiStripDigis = cms.PSet(
      initialSeed = cms.untracked.uint32(1234567),
      engineName = cms.untracked.string('HepJamesRandom'))
process.RandomNumberGeneratorService.simSiPixelDigis = cms.PSet(
      initialSeed = cms.untracked.uint32(1234567),
      engineName = cms.untracked.string('HepJamesRandom'))

# Generate H -> ZZ -> l+l- l'+l'- (l,l'=e or mu), with mH=180GeV/c2
#  process.load("FastSimulation.Configuration.HZZllll_cfi")
# Generate ttbar events
#  process.load("FastSimulation/Configuration/ttbar_cfi")
# Generate multijet events with different ptHAT bins
#  process.load("FastSimulation/Configuration/QCDpt80-120_cfi")
#  process.load("FastSimulation/Configuration/QCDpt600-800_cfi")
# Generate Minimum Bias Events
#  process.load("FastSimulation/Configuration/MinBiasEvents_cfi")
# Generate muons with a flat pT particle gun, and with pT=10.
process.load("FastSimulation/Configuration/FlatPtMuonGun_cfi")
process.FlatRandomPtGunSource.PGunParameters.PartID[0] = 13
## for 4 muons to test with vertex
#process.FlatRandomPtGunSource.PGunParameters.PartID = cms.untracked.vint32(13,-13,13,-13)
process.FlatRandomPtGunSource.PGunParameters.MinPt = 50.0
process.FlatRandomPtGunSource.PGunParameters.MaxPt = 50.0
process.FlatRandomPtGunSource.PGunParameters.MinEta = -2.5
process.FlatRandomPtGunSource.PGunParameters.MaxEta = 2.5

# Generate di-electrons with pT=35 GeV
# process.load("FastSimulation/Configuration/DiElectrons_cfi")

# from std full sim (eventsetup for digis?)
process.load("Configuration.StandardSequences.FakeConditions_cff")

# Famos sequences (fake conditions)
process.load("FastSimulation.Configuration.CommonInputsFake_cff")
process.load("FastSimulation.Configuration.FamosSequences_cff")
# replace with strawmanB geometry
process.load("SLHCUpgradeSimulations.Geometry.hybrid_cmsIdealGeometryXML_cff")
# does using an empty PixelSkimmedGeometry.txt file speeds up job with lots more channels?

process.SiPixelFakeGainOfflineESSource.file = 'SLHCUpgradeSimulations/Geometry/data/hybrid/EmptyPixelSkimmedGeometry.txt'
#process.SiPixelFakeGainOfflineESSource.file = 'SLHCUpgradeSimulations/Geometry/data/hybrid/PixelSkimmedGeometry.txt'
process.SiPixelFakeLorentzAngleESSource.file = 'SLHCUpgradeSimulations/Geometry/data/hybrid/PixelSkimmedGeometry.txt'

# Parametrized magnetic field (new mapping, 4.0 and 3.8T)
#process.load("Configuration.StandardSequences.MagneticField_40T_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.VolumeBasedMagneticFieldESProducer.useParametrizedTrackerField = True

#process.load("Configuration.StandardSequences.VtxSmearedBetafuncEarlyCollision_cff")
process.load("Configuration.StandardSequences.VtxSmearedGauss_cff")

# If you want to turn on/off pile-up
process.famosPileUp.PileUpSimulator.averageNumber = 0.0
# You may not want to simulate everything for your study
process.famosSimHits.SimulateCalorimetry = True
process.famosSimHits.SimulateTracking = True

## make occupancies more similar to full simulation
#process.famosSimHits.ParticleFilter.etaMax = 3.0
#process.famosSimHits.ParticleFilter.pTMin = 0.05
#process.famosSimHits.TrackerSimHits.pTmin = 0.05
#process.famosSimHits.TrackerSimHits.firstLoop = False

process.load("SimTracker.Configuration.SimTracker_cff")
process.simSiPixelDigis.ROUList =  ['famosSimHitsTrackerHits']
process.simSiPixelDigis.MissCalibrate = False
process.simSiPixelDigis.AddPixelInefficiency = -1
process.simSiPixelDigis.LorentzAngle_DB = False
process.simSiPixelDigis.killModules = False
process.simSiPixelDigis.NumPixelBarrel = cms.int32(6)
process.simSiPixelDigis.NumPixelEndcap = cms.int32(3)

#process.simSiStripDigis.Noise = False
process.simSiStripDigis.ROUList =  ['famosSimHitsTrackerHits']

process.load("Configuration.StandardSequences.Reconstruction_cff")
process.siPixelClusters.src = 'simSiPixelDigis'
process.siPixelClusters.MissCalibrate = False
process.siStripZeroSuppression.RawDigiProducersList[0].RawDigiProducer = 'simSiStripDigis'
process.siStripZeroSuppression.RawDigiProducersList[1].RawDigiProducer = 'simSiStripDigis'
process.siStripZeroSuppression.RawDigiProducersList[2].RawDigiProducer = 'simSiStripDigis'
process.siStripClusters.DigiProducersList[0].DigiProducer= 'simSiStripDigis'

process.load("SimGeneral.TrackingAnalysis.trackingParticles_cfi")
process.mergedtruth.TrackerHitLabels = ['famosSimHitsTrackerHits']
process.mergedtruth.simHitLabel = 'famosSimHits'

process.load("Validation.RecoTrack.cutsTPEffic_cfi")
process.load("Validation.RecoTrack.cutsTPFake_cfi")

process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")
process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
process.TrackAssociatorByHits.ROUList = ['famosSimHitsTrackerHits']

process.load("Validation.RecoTrack.MultiTrackValidator_cff")
process.multiTrackValidator.label = ['generalTracks']
### if using simple (non-iterative) or old (as in 1_8_4) tracking
#process.multiTrackValidator.label = ['ctfWithMaterialTracks']
process.multiTrackValidator.sim = 'famosSimHits'
process.multiTrackValidator.associators = ['TrackAssociatorByHits']
process.multiTrackValidator.UseAssociators = True
process.multiTrackValidator.outputFile = "validhybrid_muon_50GeV.root"

### make sure the correct (modified) error routine is used
process.siPixelRecHits.CPE = 'PixelCPEfromTrackAngle'
process.MeasurementTracker.PixelCPE = 'PixelCPEfromTrackAngle'
process.ttrhbwr.PixelCPE = 'PixelCPEfromTrackAngle'
process.mixedlayerpairs.BPix.TTRHBuilder = cms.string('WithTrackAngle')
process.mixedlayerpairs.FPix.TTRHBuilder = cms.string('WithTrackAngle')
process.pixellayertriplets.BPix.TTRHBuilder = cms.string('WithTrackAngle')
process.pixellayertriplets.FPix.TTRHBuilder = cms.string('WithTrackAngle')
process.ctfWithMaterialTracks.TTRHBuilder = cms.string('WithTrackAngle')
#next may not be needed
process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
process.TrackRefitter.TTRHBuilder = cms.string('WithTrackAngle')

#next may not be needed
process.load("RecoTracker.SiTrackerMRHTools.SiTrackerMultiRecHitUpdator_cff")
process.siTrackerMultiRecHitUpdator.TTRHBuilder = cms.string('WithTrackAngle')

#replace with correct component in cloned version (replace with original TTRH producer)
process.preFilterFirstStepTracks.TTRHBuilder = cms.string('WithTrackAngle')
process.secPixelRecHits.CPE = cms.string('PixelCPEfromTrackAngle')
process.seclayertriplets.BPix.TTRHBuilder = cms.string('WithTrackAngle')
process.seclayertriplets.FPix.TTRHBuilder = cms.string('WithTrackAngle')
process.secMeasurementTracker.PixelCPE = cms.string('PixelCPEfromTrackAngle')
process.secWithMaterialTracks.TTRHBuilder = cms.string('WithTrackAngle')
process.thPixelRecHits.CPE = cms.string('PixelCPEfromTrackAngle')
process.thlayerpairs.BPix.TTRHBuilder = cms.string('WithTrackAngle')
process.thlayerpairs.FPix.TTRHBuilder = cms.string('WithTrackAngle')
process.thMeasurementTracker.PixelCPE = cms.string('PixelCPEfromTrackAngle')
process.thWithMaterialTracks.TTRHBuilder = cms.string('WithTrackAngle')

### produce an ntuple with hits for analysis
process.ReadLocalMeasurement = cms.EDAnalyzer("StdHitNtuplizer",
   src = cms.InputTag("siPixelRecHits"),
   stereoRecHits = cms.InputTag("siStripMatchedRecHits","stereoRecHit"),
   rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
   matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
   trackProducer = cms.InputTag("generalTracks"),
   ### if using simple (non-iterative) or old (as in 1_8_4) tracking
   #trackProducer = cms.InputTag("ctfWithMaterialTracks"),
   OutputFile = cms.string("stdgrechitHY_ntuple.root"),
   ### for using track hit association
   associatePixel = cms.bool(True),
   associateStrip = cms.bool(False),
   associateRecoTracks = cms.bool(False),
   ROUList = cms.vstring('famosSimHitsTrackerHits')
)

process.o1 = cms.OutputModule(
    "PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *',
                                           'drop *_mix_*_*'),
    fileName = cms.untracked.string('/uscms_data/d2/cheung/slhc/fastsimHY_50mu.root')
)

process.outpath = cms.EndPath(process.o1)

# Make the job crash in case of missing product
process.options = cms.untracked.PSet( Rethrow = cms.untracked.vstring('ProductNotFound') )

process.Timing =  cms.Service("Timing")
process.load("FWCore/MessageService/MessageLogger_cfi")
process.MessageLogger.destinations = cms.untracked.vstring("detailedInfo_hybrid_mu50")
### to output debug messages for particular modules
# process.MessageLogger.detailedInfo_strawb_mu50 = cms.untracked.PSet(threshold = cms.untracked.string('DEBUG'))
# process.MessageLogger.debugModules= cms.untracked.vstring("*")


# Famos with tracks
process.p1 = cms.Path(process.famosWithTrackerHits)
process.p2 = cms.Path(process.trDigi)
process.p3 = cms.Path(process.trackerlocalreco)
process.p6 = cms.Path(process.offlineBeamSpot+process.recopixelvertexing*process.ckftracks)
process.p8 = cms.Path(process.trackingParticles*process.cutsTPEffic*process.cutsTPFake*process.multiTrackValidator)
process.p9 = cms.Path(process.ReadLocalMeasurement)
#process.schedule = cms.Schedule(process.p1,process.p2,process.p3,process.p9,process.outpath)
process.schedule = cms.Schedule(process.p1,process.p2,process.p3,process.p6,process.p8,process.p9,process.outpath)
