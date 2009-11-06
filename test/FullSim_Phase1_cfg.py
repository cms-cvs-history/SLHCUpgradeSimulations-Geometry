import FWCore.ParameterSet.Config as cms

process = cms.Process("Fullsim")

process.load("Configuration.StandardSequences.Services_cff")

process.load("Configuration.StandardSequences.Geometry_cff")
#PhaseI Geometry
process.load("SLHCUpgradeSimulations.Geometry.PhaseI_cmsSimIdealGeometryXML_cff")

process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("Configuration.StandardSequences.FakeConditions_cff")
process.SiPixelFakeGainOfflineESSource.file = 'SLHCUpgradeSimulations/Geometry/data/PhaseI/PixelSkimmedGeometry_phase1.txt'
process.SiPixelFakeLorentzAngleESSource.file = 'SLHCUpgradeSimulations/Geometry/data/PhaseI/PixelSkimmedGeometry_phase1.txt'

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.destinations = cms.untracked.vstring("detailedInfo_fullph1geommu50")

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

#process.load("Configuration.StandardSequences.MixingLowLumiPileUp_cff")
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
process.simSiPixelDigis.NumPixelBarrel = cms.int32(4)
process.simSiPixelDigis.NumPixelEndcap = cms.int32(3)

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

## J/psi from pythia
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

process.FEVT = cms.OutputModule("PoolOutputModule",
    process.FEVTSIMEventContent,
    fileName = cms.untracked.string('/uscms_data/d2/cheung/slhc/testfullph1g_muon_50GeV.root')
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
process.multiTrackValidator.outputFile = "validfullph1g_muon_50GeV.root"
process.multiTrackValidator.nint = cms.int32(20)
process.multiTrackValidator.nintpT = cms.int32(25)
process.multiTrackValidator.maxpT = cms.double(50.0)

##### with John's changes ##############################
process.load("SLHCUpgradeSimulations.Geometry.oldTracking_wtriplets")
process.pixellayertriplets.layerList = cms.vstring('BPix1+BPix2+BPix3',
        'BPix1+BPix3+BPix4',
        'BPix2+BPix3+BPix4',
        'BPix1+BPix2+BPix4',
        'BPix1+BPix2+FPix1_pos',
        'BPix1+BPix2+FPix1_neg',
        'BPix1+FPix1_pos+FPix2_pos',
        'BPix1+FPix1_neg+FPix2_neg',
        'BPix1+FPix2_pos+FPix3_pos',
        'BPix1+FPix2_neg+FPix3_neg',
        'FPix1_pos+FPix2_pos+FPix3_pos',
        'FPix1_neg+FPix2_neg+FPix3_neg')
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
   OutputFile = cms.string("stdgrechitfullph1g_ntuple.root"),
   ### for using track hit association
   associatePixel = cms.bool(True),
   associateStrip = cms.bool(False),
   associateRecoTracks = cms.bool(False),
   ROUList = cms.vstring('g4SimHitsTrackerHitsPixelBarrelLowTof',
                         'g4SimHitsTrackerHitsPixelBarrelHighTof',
                         'g4SimHitsTrackerHitsPixelEndcapLowTof',
                         'g4SimHitsTrackerHitsPixelEndcapHighTof')
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
   label = cms.VInputTag(cms.InputTag("generalTracks")),
   #label = cms.VInputTag(cms.InputTag("ctfWithMaterialTracks")),
   label_tp_effic = cms.InputTag("mergedtruth","MergedTrackTruth"),
   label_tp_fake = cms.InputTag("mergedtruth","MergedTrackTruth"),
   associators = cms.vstring('TrackAssociatorByHits'),
   UseAssociators = cms.bool(True),
   OutputFile = cms.string("tpanalfull_ntuple.root")
)
#process.TPanal.label_tp_effic = cms.InputTag("cutsTPEffic")
#process.TPanal.label_tp_fake = cms.InputTag("cutsTPFake")

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
#process.p9 = cms.Path(process.TPanal)
process.outpath = cms.EndPath(process.FEVT)
#process.schedule = cms.Schedule(process.p0,process.p1,process.p2,process.p3,process.p4,process.p5,process.p6,process.p7,process.outpath)
#process.schedule = cms.Schedule(process.p0,process.p1,process.p2,process.p3,process.p4,process.p5,process.p6,process.p7,process.p8,process.p9)
#process.schedule = cms.Schedule(process.p0,process.p1,process.p2,process.p3,process.p4,process.p5,process.p6,process.p7,process.p8)
process.schedule = cms.Schedule(process.p0,process.p1,process.p2,process.p3,process.p6,process.p7,process.p8)