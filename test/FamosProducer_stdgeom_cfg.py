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
#process.load("FastSimulation/Configuration/MinBiasEvents_cfi")
# Generate muons with a flat pT particle gun
process.load("FastSimulation/Configuration/FlatPtMuonGun_cfi")
process.FlatRandomPtGunSource.PGunParameters.PartID[0] = 13
#process.FlatRandomPtGunSource.PGunParameters.PartID[0] = 211
### for 4 muons to test with vertex
##process.FlatRandomPtGunSource.PGunParameters.PartID = cms.untracked.vint32(13,-13,13,-13)
### for opposite sign back-to-back dimuon pairs
process.FlatRandomPtGunSource.PGunParameters.MinPt = 0.9
process.FlatRandomPtGunSource.PGunParameters.MaxPt = 50.0
process.FlatRandomPtGunSource.PGunParameters.MinEta = -2.4
process.FlatRandomPtGunSource.PGunParameters.MaxEta = 2.4
process.FlatRandomPtGunSource.AddAntiParticle = cms.untracked.bool(True)
# Generate QCD dijet events
#from Configuration.Generator.PythiaUESettings_cfi import *
###process.load("Configuration.Generator.PythiaUESettings_cfi")
#process.source = cms.Source("PythiaSource",
#   pythiaHepMCVerbosity = cms.untracked.bool(False),
#   maxEventsToPrint = cms.untracked.int32(0),
#   pythiaPylistVerbosity = cms.untracked.int32(0),
#   filterEfficiency = cms.untracked.double(1.0),
#   comEnergy = cms.untracked.double(14000.0),
#   PythiaParameters = cms.PSet(
#      pythiaUESettingsBlock,
#      processParameters = cms.vstring('MSEL=1               ! QCD hight pT processes',
#              'MSEL=0          ! user control',
#              'MSUB(11)=1      ! qq to qq',
#              'MSUB(68)=1      ! gg to gg',
#              'MSUB(28)=1      ! qg to qg',
#              'MSUB(53)=1      ! gg to qq',
#              'CKIN(3)=100.    ! Pt low cut but also the Et jet required',
#              'CKIN(3)=120.    ! Pt high cut but also the Et jet required',
#              'CKIN(13)=0.     ! etamin',
#              'CKIN(14)=2.5    ! etamax',
#              'CKIN(15)=-2.5   ! -etamax',
#              'CKIN(16)=0.     ! -etamin'),
#      # This is a vector of ParameterSet names to be read, in this order
#      parameterSets = cms.vstring('pythiaUESettings','processParameters')
#   )
#)

# J/psi from Pythia
#from Configuration.Generator.PythiaUESettings_cfi import *
#process.source = cms.Source("PythiaSource",
#    Phimin = cms.untracked.double(0.0),
#    maxEventsToPrint = cms.untracked.int32(1),
#    pythiaPylistVerbosity = cms.untracked.int32(1),
#    #  possibility to run single or double back-to-back particles with PYTHIA
#    # if ParticleID = 0, run PYTHIA
#    ParticleID = cms.untracked.int32(443),
#    pythiaHepMCVerbosity = cms.untracked.bool(True),
#    Etamin = cms.untracked.double(0.0),
#    DoubleParticle = cms.untracked.bool(False),
#    Phimax = cms.untracked.double(360.0),
#    Ptmin = cms.untracked.double(20.0),
#    Ptmax = cms.untracked.double(40.0),
#    Etamax = cms.untracked.double(2.4),
#    PythiaParameters = cms.PSet(
#        pythiaUESettingsBlock,
#        pythiaJpsiDecays = cms.vstring('MDME(858,1)=1                 ! J/psi -> ee turned ON',
#            'MDME(859,1)=1                 ! J/psi -> mumu turned ON',
#            'MDME(860,1)=0                 ! J/psi -> random turned OFF'),
#        # This is a vector of ParameterSet names to be read, in this order
#        parameterSets = cms.vstring('pythiaUESettings',
#            'pythiaJpsiDecays')
#    )
#)



# Generate di-electrons with pT=35 GeV
# process.load("FastSimulation/Configuration/DiElectrons_cfi")

# from std full sim
process.load("Configuration.StandardSequences.FakeConditions_cff")

# Famos sequences (fake conditions)
process.load("FastSimulation.Configuration.CommonInputsFake_cff")
process.load("FastSimulation.Configuration.FamosSequences_cff")
# does using an empty PixelSkimmedGeometry.txt file speeds up job with lots more channels?
process.SiPixelFakeGainOfflineESSource.file = 'SLHCUpgradeSimulations/Geometry/data/hybrid/EmptyPixelSkimmedGeometry.txt'

# Parametrized magnetic field (new mapping, 4.0 and 3.8T)
#process.load("Configuration.StandardSequences.MagneticField_40T_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.VolumeBasedMagneticFieldESProducer.useParametrizedTrackerField = True

#process.load("Configuration.StandardSequences.VtxSmearedBetafuncEarlyCollision_cff")
process.load("Configuration.StandardSequences.VtxSmearedGauss_cff")

# Replace std 10 TeV with 14 TeV pileup files and set the vertex smearing like signal
import FastSimulation.Event.GaussianVertexGenerator_cfi as GaussSmearing
process.famosPileUp.VertexGenerator = cms.PSet( GaussSmearing.myVertexGenerator )
import FastSimulation.PileUpProducer.PileUpSimulator_cfi as Pileup14TeV
process.famosPileUp.PileUpSimulator = cms.PSet( Pileup14TeV.PileUpSimulatorBlock.PileUpSimulator )

# Make sure CoM energy is 14 TeV if we are using pythia for the signal source
#process.PythiaSource.comEnergy = cms.untracked.double(14000.0)
#process.PythiaSource.maxEventsToPrint = 1

# If you want to turn on/off pile-up
process.famosPileUp.PileUpSimulator.averageNumber = 5.0
# You may not want to simulate everything for your study
process.famosSimHits.SimulateCalorimetry = True
process.famosSimHits.SimulateTracking = True

## make occupancies more similar to full simulation
process.famosSimHits.ParticleFilter.etaMax = 3.0
process.famosSimHits.ParticleFilter.pTMin = 0.05
process.famosSimHits.TrackerSimHits.pTmin = 0.05
process.famosSimHits.TrackerSimHits.firstLoop = False

process.load("SimTracker.Configuration.SimTracker_cff")
process.simSiPixelDigis.ROUList =  ['famosSimHitsTrackerHits']
process.simSiPixelDigis.MissCalibrate = False
process.simSiPixelDigis.AddPixelInefficiency = -1
process.simSiPixelDigis.LorentzAngle_DB = False
process.simSiPixelDigis.killModules = False

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
#process.mergedtruth.mergedBremsstrahlung = cms.bool(False)
#process.mergedtruth.firstGeneratedPairOnly = cms.bool(False)

process.load("Validation.RecoTrack.cutsTPEffic_cfi")
process.load("Validation.RecoTrack.cutsTPFake_cfi")
## if mergedBremsstrahlung is False
#process.cutsTPEffic.src = cms.InputTag("mergedtruth")
#process.cutsTPFake.src = cms.InputTag("mergedtruth")

process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")
process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
process.TrackAssociatorByHits.ROUList = ['famosSimHitsTrackerHits']

process.load("Validation.RecoTrack.MultiTrackValidator_cff")
#process.multiTrackValidator.label = ['generalTracks']
### if using simple (non-iterative) or old (as in 1_8_4) tracking
#process.multiTrackValidator.label = ['ctfWithMaterialTracks']
process.multiTrackValidator.label = ['cutsRecoTracks']
process.multiTrackValidator.label_tp_effic = cms.InputTag("cutsTPEffic")
process.multiTrackValidator.label_tp_fake = cms.InputTag("cutsTPFake")
process.multiTrackValidator.sim = 'famosSimHits'
process.multiTrackValidator.associators = ['TrackAssociatorByHits']
process.multiTrackValidator.UseAssociators = True
process.multiTrackValidator.outputFile = "validstdgeom_muon_50GeV.root"
process.multiTrackValidator.nint = cms.int32(20)
process.multiTrackValidator.nintpT = cms.int32(25)
process.multiTrackValidator.maxpT = cms.double(50.0)
#process.multiTrackValidator.firstGeneratedPairOnly = cms.bool(True)
#process.multiTrackValidator.firstGeneratedPairOnly = cms.bool(False)
##### with John's changes ##############################
process.load("SLHCUpgradeSimulations.Geometry.oldTracking_wtriplets")
# restrict vertex fining in trackingtruthprod to smaller volume (note: these numbers in mm) 
process.mergedtruth.volumeRadius = cms.double(100.0)
process.mergedtruth.volumeZ = cms.double(900.0)
process.mergedtruth.discardOutVolume = cms.bool(True)

###process.cutsTPEffic.ptMin = cms.double(2.5)
###process.cutsTPFake.ptMin = cms.double(2.0)
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

### make sure the correct (modified) error routine is used
#process.siPixelRecHits.CPE = 'PixelCPEfromTrackAngle'
#process.MeasurementTracker.PixelCPE = 'PixelCPEfromTrackAngle'
#process.ttrhbwr.PixelCPE = 'PixelCPEfromTrackAngle'
#process.mixedlayerpairs.BPix.TTRHBuilder = cms.string('WithTrackAngle')
#process.mixedlayerpairs.FPix.TTRHBuilder = cms.string('WithTrackAngle')
#process.pixellayertriplets.BPix.TTRHBuilder = cms.string('WithTrackAngle')
#process.pixellayertriplets.FPix.TTRHBuilder = cms.string('WithTrackAngle')
#process.ctfWithMaterialTracks.TTRHBuilder = cms.string('WithTrackAngle')
#next may not be needed
#process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
#process.TrackRefitter.TTRHBuilder = cms.string('WithTrackAngle')

#next may not be needed
#process.load("RecoTracker.SiTrackerMRHTools.SiTrackerMultiRecHitUpdator_cff")
#process.siTrackerMultiRecHitUpdator.TTRHBuilder = cms.string('WithTrackAngle')

#replace with correct component in cloned version (replace with original TTRH producer)
#process.preFilterFirstStepTracks.TTRHBuilder = cms.string('WithTrackAngle')
#process.secPixelRecHits.CPE = cms.string('PixelCPEfromTrackAngle')
#process.seclayertriplets.BPix.TTRHBuilder = cms.string('WithTrackAngle')
#process.seclayertriplets.FPix.TTRHBuilder = cms.string('WithTrackAngle')
#process.secMeasurementTracker.PixelCPE = cms.string('PixelCPEfromTrackAngle')
#process.secWithMaterialTracks.TTRHBuilder = cms.string('WithTrackAngle')
#process.thPixelRecHits.CPE = cms.string('PixelCPEfromTrackAngle')
#process.thlayerpairs.BPix.TTRHBuilder = cms.string('WithTrackAngle')
#process.thlayerpairs.FPix.TTRHBuilder = cms.string('WithTrackAngle')
#process.thMeasurementTracker.PixelCPE = cms.string('PixelCPEfromTrackAngle')
#process.thWithMaterialTracks.TTRHBuilder = cms.string('WithTrackAngle')

### produce an ntuple with hits for analysis
process.ReadLocalMeasurement = cms.EDAnalyzer("StdHitNtuplizer",
   src = cms.InputTag("siPixelRecHits"),
   stereoRecHits = cms.InputTag("siStripMatchedRecHits","stereoRecHit"),
   rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
   matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
   simtrackHits = cms.InputTag("famosSimHits"),
   #trackProducer = cms.InputTag("generalTracks"),
   ### if using simple (non-iterative) or old (as in 1_8_4) tracking
   trackProducer = cms.InputTag("ctfWithMaterialTracks"),
   OutputFile = cms.string("stdgrechitstdgeom_ntuple.root"),
   ### for using track hit association
   associatePixel = cms.bool(True),
   associateStrip = cms.bool(False),
   associateRecoTracks = cms.bool(False),
   ROUList = cms.vstring('famosSimHitsTrackerHits')
)

TrackingParticleSelectionForTP = cms.PSet(
    lipTP = cms.double(30.0),
    chargedOnlyTP = cms.bool(True),
    stableOnlyTP = cms.bool(True),
    pdgIdTP = cms.vint32(),
    signalOnlyTP = cms.bool(True),
    minRapidityTP = cms.double(-2.4),
    minHitTP = cms.int32(0),
    ptMinTP = cms.double(0.9),
    maxRapidityTP = cms.double(2.4),
    tipTP = cms.double(3.5)
)
process.TPanal = cms.EDAnalyzer("TPNtuplizer",
   TrackingParticleSelectionForTP,
   label = cms.VInputTag(cms.InputTag("ctfWithMaterialTracks")),
   label_tp_effic = cms.InputTag("mergedtruth","MergedTrackTruth"),
   label_tp_fake = cms.InputTag("mergedtruth","MergedTrackTruth"),
   associators = cms.vstring('TrackAssociatorByHits'),
   UseAssociators = cms.bool(True),
   OutputFile = cms.string("tpanal_ntuple.root")
)
#process.TPanal.label_tp_effic = cms.InputTag("cutsTPEffic")
#process.TPanal.label_tp_fake = cms.InputTag("cutsTPFake")

process.o1 = cms.OutputModule(
    "PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *',
                                           'drop *_mix_*_*'),
    fileName = cms.untracked.string('/uscms_data/d2/cheung/slhc/fastsimStd_50mu.root')
)

process.outpath = cms.EndPath(process.o1)

# Make the job crash in case of missing product
process.options = cms.untracked.PSet( Rethrow = cms.untracked.vstring('ProductNotFound') )

process.Timing =  cms.Service("Timing")
process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
     oncePerEventMode = cms.untracked.bool(False),
     showMallocInfo = cms.untracked.bool(False),
     moduleMemorySummary = cms.untracked.bool(True),
     ignoreTotal = cms.untracked.int32(5)
)

process.load("FWCore/MessageService/MessageLogger_cfi")
#process.MessageLogger.destinations = cms.untracked.vstring("detailedInfo_stdgeom_mu50")
### to output debug messages for particular modules
# process.MessageLogger.detailedInfo_strawb_mu50 = cms.untracked.PSet(threshold = cms.untracked.string('DEBUG'))
# process.MessageLogger.debugModules= cms.untracked.vstring("*")

#process.anal = cms.EDAnalyzer("EventContentAnalyzer")

# Famos with tracks
process.p1 = cms.Path(process.famosWithTrackerHits)
process.p2 = cms.Path(process.trDigi)
process.p3 = cms.Path(process.trackerlocalreco)
process.p6 = cms.Path(process.oldTracking_wtriplets)
#process.p6 = cms.Path(process.offlineBeamSpot+process.recopixelvertexing*process.ckftracks)
process.p8 = cms.Path(process.trackingParticles*process.cutsTPEffic*process.cutsTPFake*process.cutsRecoTracks*process.multiTrackValidator)
#process.p9 = cms.Path(process.ReadLocalMeasurement)
process.p9 = cms.Path(process.TPanal)
#process.p9 = cms.Path(process.anal)
#process.schedule = cms.Schedule(process.p1,process.p2,process.p3,process.p6,process.p8,process.p9,process.outpath)
process.schedule = cms.Schedule(process.p1,process.p2,process.p3,process.p6,process.p8)
#process.schedule = cms.Schedule(process.p1,process.p2,process.p3,process.p6,process.p8,process.p9)

