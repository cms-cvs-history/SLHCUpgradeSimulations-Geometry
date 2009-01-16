import FWCore.ParameterSet.Config as cms

process = cms.Process("Fastsimwdigi")

# Number of events to be generated
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

# Include the RandomNumberGeneratorService definition
#process.load("FastSimulation.Configuration.RandomServiceInitialization_cff")
process.RandomNumberGeneratorService = cms.Service(

    "RandomNumberGeneratorService",

    # To save the status of the last event (useful for crashes)
    # Just give a name to the file you want the status to be saved
    # otherwise just put saveFileName = ""
    saveFileName = cms.untracked.string(''),

    # To restore the status of the last event, just un-comment the
    # following line (and comment the saveFileName line!)
    # restoreFileName = cms.string('RandomEngineState.log'),

    # To reproduce events using the RandomEngineStateProducer (source
    # excluded), comment the sourceSeed definition, and un-comment 
    # the restoreStateLabel
    # restoreStateLabel = cms.string('randomEngineStateProducer'),

    # This is to initialize the random engine of the source
    theSource = cms.PSet(
        initialSeed = cms.untracked.uint32(123456789),
        engineName = cms.untracked.string('TRandom3')
    ),

    # This is to initialize the random engines used for  Famos
    VtxSmeared = cms.PSet(
        initialSeed = cms.untracked.uint32(123456789),
        engineName = cms.untracked.string('TRandom3')
    ),


    famosPileUp = cms.PSet(
        initialSeed = cms.untracked.uint32(918273),
        engineName = cms.untracked.string('TRandom3')
    ),

    famosSimHits = cms.PSet(
        initialSeed = cms.untracked.uint32(13579),
        engineName = cms.untracked.string('TRandom3')
    ),

    siTrackerGaussianSmearingRecHits = cms.PSet(
        initialSeed = cms.untracked.uint32(24680),
        engineName = cms.untracked.string('TRandom3')
    ),

    caloRecHits = cms.PSet(
        initialSeed = cms.untracked.uint32(654321),
        engineName = cms.untracked.string('TRandom3')
    ),

    paramMuons = cms.PSet(
        initialSeed = cms.untracked.uint32(54525),
        engineName = cms.untracked.string('TRandom3')
    ),

    l1ParamMuons = cms.PSet(
        initialSeed = cms.untracked.uint32(6453209),
        engineName = cms.untracked.string('TRandom3')
    ),

    MuonSimHits = cms.PSet(
        initialSeed = cms.untracked.uint32(987346),
        engineName = cms.untracked.string('TRandom3')
    ),

    simMuonRPCDigis = cms.PSet(
        initialSeed = cms.untracked.uint32(524964),
        engineName = cms.untracked.string('TRandom3')
    ),

    simMuonCSCDigis = cms.PSet(
        initialSeed = cms.untracked.uint32(525432),
        engineName = cms.untracked.string('TRandom3')
    ),

    simMuonDTDigis = cms.PSet(
        initialSeed = cms.untracked.uint32(67673876),
        engineName = cms.untracked.string('TRandom3')
    ),

   simSiStripDigis = cms.PSet(
      initialSeed = cms.untracked.uint32(1234567),
      engineName = cms.untracked.string('HepJamesRandom')
   ),
   simSiPixelDigis = cms.PSet(
      initialSeed = cms.untracked.uint32(1234567),
      engineName = cms.untracked.string('HepJamesRandom')
   )
)

randomEngineStateProducer = cms.EDProducer("RandomEngineStateProducer")


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
# replace FlatRandomPtGunSource.PGunParameters.PartID={13}
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
process.load("SLHCUpgradeSimulations.Geometry.strawmanb_cmsIdealGeometryXML_cff")


# Parametrized magnetic field (new mapping, 4.0 and 3.8T)
#process.load("Configuration.StandardSequences.MagneticField_40T_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.VolumeBasedMagneticFieldESProducer.useParametrizedTrackerField = True

# If you want to turn on/off pile-up
process.famosPileUp.PileUpSimulator.averageNumber = 0.0
# You may not want to simulate everything for your study
process.famosSimHits.SimulateCalorimetry = True
process.famosSimHits.SimulateTracking = True

# Tracker Digis (Pixel + SiStrips)
# returns sequence "trDigi"
#
process.load("SimTracker.Configuration.SimTracker_cff")
process.simSiPixelDigis.ROUList =  ['famosSimHitsTrackerHits']
process.simSiPixelDigis.MissCalibrate = False
process.simSiStripDigis.ROUList =  ['famosSimHitsTrackerHits']

process.load("Configuration.StandardSequences.DigiToRaw_cff")

process.load("Configuration.StandardSequences.RawToDigi_cff")

process.load("Configuration.StandardSequences.Reconstruction_cff")

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
process.multiTrackValidator.associators = ['TrackAssociatorByHits']
process.multiTrackValidator.UseAssociators = True
process.multiTrackValidator.outputFile = "validstrawb_muon_50GeV.root"

# To write out events (not need: FastSimulation _is_ fast!)
process.o1 = cms.OutputModule(
    "PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *',
                                           'drop *_mix_*_*'),
    fileName = cms.untracked.string('/uscms_data/d1/cheung/slhc/strawb_muon_50GeV.root')
)
process.outpath = cms.EndPath(process.o1)

# Keep output to a nice level
process.Timing =  cms.Service("Timing")
process.load("FWCore/MessageService/MessageLogger_cfi")
process.MessageLogger.destinations = cms.untracked.vstring("detailedInfo_strawb_mu50.txt")

# Make the job crash in case of missing product
process.options = cms.untracked.PSet( Rethrow = cms.untracked.vstring('ProductNotFound') )

# Famos with tracks
process.p1 = cms.Path(process.famosWithTrackerHits)
process.p2 = cms.Path(process.trDigi)
process.p3 = cms.Path(process.siPixelRawData*process.SiStripDigiToRaw*process.rawDataCollector)
process.p4 = cms.Path(process.siPixelDigis*process.SiStripRawToDigis)
process.p5 = cms.Path(process.trackerlocalreco)
process.p6 = cms.Path(process.offlineBeamSpot+process.recopixelvertexing*process.ckftracks)
process.p7 = cms.Path(process.trackingParticles*process.cutsTPEffic*process.cutsTPFake*process.multiTrackValidator)
process.schedule = cms.Schedule(process.p1,process.p2,process.p3,process.p4,process.p5,process.p6,process.p7,process.outpath)
