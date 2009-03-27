import FWCore.ParameterSet.Config as cms

process = cms.Process("VALID")

# Number of events to be generated
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(20000)
)

#Timing = cms.Service( )

# Include the RandomNumberGeneratorService definition
process.load("FastSimulation.Configuration.RandomServiceInitialization_cff")

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
process.FlatRandomPtGunSource.PGunParameters.MinEta = -3.0
process.FlatRandomPtGunSource.PGunParameters.MaxEta = 3.0
# Generate di-electrons with pT=35 GeV
# process.load("FastSimulation/Configuration/DiElectrons_cfi")


# Famos sequences (fake conditions)
process.load("FastSimulation.Configuration.CommonInputsFake_cff")
process.load("FastSimulation.Configuration.FamosSequences_cff")
# replace with strawmanB geometry
process.load("SLHCUpgradeSimulations.Geometry.strawmana_cmsIdealGeometryXML_cff")

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

# taking these from FastSimulation/Validation/python/TrackValidation_HighPurity_cff.py
# and e.g. FastSimulation/Validation/test/valTK_muon_100GeV_cfg.py
process.load("SimGeneral.TrackingAnalysis.trackingParticles_cfi")
process.mergedtruth.TrackerHitLabels = ['famosSimHitsTrackerHits']
process.mergedtruth.simHitLabel = 'famosSimHits'

#process.load("SLHCUpgradeSimulations.Geometry.cutsTPEffic_cfi")
#process.load("SLHCUpgradeSimulations.Geometry.cutsTPFake_cfi")
process.load("Validation.RecoTrack.cutsTPEffic_cfi")
process.load("Validation.RecoTrack.cutsTPFake_cfi")

process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")
process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
process.TrackAssociatorByHits.associateStrip = False
process.TrackAssociatorByHits.associatePixel = False
#process.TrackAssociatorByHits.ROUList = "TrackerHits"
process.TrackAssociatorByHits.ROUList = ['famosSimHitsTrackerHits']

process.load("Validation.RecoTrack.MultiTrackValidator_cff")
process.multiTrackValidator.label = ['generalTracks']
process.multiTrackValidator.associators = ['TrackAssociatorByHits']
process.multiTrackValidator.UseAssociators = True
process.multiTrackValidator.outputFile = "valid_muon_50GeV.root"


# Famos with tracks
process.p1 = cms.Path(process.famosWithTracks*process.trackerGSRecHitTranslator*process.trackingParticles)
process.p1 *= process.cutsTPEffic*process.cutsTPFake*process.multiTrackValidator

# To write out events (not need: FastSimulation _is_ fast!)
process.o1 = cms.OutputModule(
    "PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *', 
                                           'drop *_mix_*_*'),
    fileName = cms.untracked.string('/uscms_data/d2/cheung/slhc/test_muon_50GeV.root')
)
process.outpath = cms.EndPath(process.o1)

# Keep output to a nice level
process.Timing =  cms.Service("Timing")
process.load("FWCore/MessageService/MessageLogger_cfi")
process.MessageLogger.destinations = cms.untracked.vstring("detailedInfo_mu50.txt")

# Make the job crash in case of missing product
process.options = cms.untracked.PSet( Rethrow = cms.untracked.vstring('ProductNotFound') )
