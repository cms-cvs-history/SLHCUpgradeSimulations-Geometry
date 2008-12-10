import FWCore.ParameterSet.Config as cms

process = cms.Process("VALID")

# Number of events to be generated
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
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
process.FlatRandomPtGunSource.PGunParameters.MinEta = -2.5
process.FlatRandomPtGunSource.PGunParameters.MaxEta = 2.5
# Generate di-electrons with pT=35 GeV
# process.load("FastSimulation/Configuration/DiElectrons_cfi")


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
process.famosPileUp.PileUpSimulator.averageNumber = 5.0
# You may not want to simulate everything for your study
process.famosSimHits.SimulateCalorimetry = False
process.famosSimHits.SimulateTracking = True

# Famos with everything !
process.p1 = cms.Path(process.famosWithTracks)

# To write out events (not need: FastSimulation _is_ fast!)
process.o1 = cms.OutputModule(
    "PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *', 
                                           'drop *_mix_*_*'),
    fileName = cms.untracked.string('/uscms_data/d1/cheung/slhc/test_muon_50GeV.root')
)
process.outpath = cms.EndPath(process.o1)

# Keep output to a nice level
process.Timing =  cms.Service("Timing")
process.load("FWCore/MessageService/MessageLogger_cfi")
process.MessageLogger.destinations = cms.untracked.vstring("detailedInfo_mu50.txt")

# Make the job crash in case of missing product
process.options = cms.untracked.PSet( Rethrow = cms.untracked.vstring('ProductNotFound') )
