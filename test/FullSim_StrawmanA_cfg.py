import FWCore.ParameterSet.Config as cms

process = cms.Process("Fullsim")

process.load("Configuration.StandardSequences.Services_cff")

process.load("Configuration.StandardSequences.Geometry_cff")
# replace with strawmanA geometry
process.load("SLHCUpgradeSimulations.Geometry.strawmana_cmsIdealGeometryXML_cff")

process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("Configuration.StandardSequences.FakeConditions_cff")
process.SiPixelFakeGainOfflineESSource.file = 'SLHCUpgradeSimulations/Geometry/data/strawmana/PixelSkimmedGeometry.txt'
process.SiPixelFakeLorentzAngleESSource.file = 'SLHCUpgradeSimulations/Geometry/data/strawmana/PixelSkimmedGeometry.txt'

process.load("FWCore/MessageService/MessageLogger_cfi")
process.MessageLogger.destinations = cms.untracked.vstring("detailedInfo_fullstrawamu50")

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
process.load("SimTracker.Configuration.SimTracker_cff")
from SimTracker.Configuration.SimTracker_cff import *

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
from Configuration.StandardSequences.DigiToRaw_cff import *

process.load("Configuration.StandardSequences.RawToDigi_cff")
from Configuration.StandardSequences.RawToDigi_cff import *

#process.load("Configuration.StandardSequences.VtxSmearedBetafuncEarlyCollision_cff")
process.load("Configuration.StandardSequences.VtxSmearedGauss_cff")

process.load("Configuration.StandardSequences.Reconstruction_cff")

process.load("RecoLocalTracker.Configuration.RecoLocalTracker_cff")
process.load("RecoTracker.Configuration.RecoTracker_cff")
process.load("TrackingTools.Configuration.TrackingTools_cff")

process.load("TrackingTools.TrajectoryFiltering.TrajectoryFilter_cff")
from TrackingTools.TrajectoryFiltering.TrajectoryFilter_cff import *

from RecoLocalTracker.SiPixelClusterizer.SiPixelClusterizer_cfi import *
from RecoLocalTracker.SiStripClusterizer.SiStripClusterizer_SimData2_cfi import *
from RecoLocalTracker.SiStripZeroSuppression.SiStripZeroSuppression_SimData_cfi import *

from RecoLocalTracker.SiPixelRecHits.SiPixelRecHits_cfi import *
from RecoLocalTracker.SiStripRecHitConverter.SiStripRecHitConverter_cfi import *
from RecoLocalTracker.SiStripRecHitConverter.SiStripRecHitMatcher_cfi import *
from RecoLocalTracker.SiStripRecHitConverter.StripCPEfromTrackAngle_cfi import *

process.load("SimTracker.TrackAssociation.trackMCMatchSequence_cff")
from SimTracker.TrackAssociation.trackMCMatchSequence_cff import *

#process.load("SimTracker.Configuration.SimTracker_cff")
#from SimTracker.Configuration.SimTracker_cff import *
process.load("SimTracker.SiPixelDigitizer.PixelDigi_cfi")
process.simSiPixelDigis.MissCalibrate = False
process.simSiPixelDigis.AddPixelInefficiency = -1
process.simSiPixelDigis.LorentzAngle_DB = False
process.simSiPixelDigis.killModules = False
process.simSiPixelDigis.NumPixelBarrel = cms.int32(6)

process.load("RecoLocalTracker.SiPixelClusterizer.SiPixelClusterizer_cfi")
from RecoLocalTracker.SiPixelClusterizer.SiPixelClusterizer_cfi import *
process.siPixelClusters.src = 'simSiPixelDigis'
process.siPixelClusters.MissCalibrate = False

process.load("RecoLocalTracker.SiStripClusterizer.SiStripClusterizer_SimData2_cfi")
from RecoLocalTracker.SiStripClusterizer.SiStripClusterizer_SimData2_cfi import *

process.load("RecoLocalTracker.SiStripZeroSuppression.SiStripZeroSuppression_SimData_cfi")
from RecoLocalTracker.SiStripZeroSuppression.SiStripZeroSuppression_SimData_cfi import *

process.siStripZeroSuppression.RawDigiProducersList[0].RawDigiProducer = 'simSiStripDigis'
process.siStripZeroSuppression.RawDigiProducersList[1].RawDigiProducer = 'simSiStripDigis'
process.siStripZeroSuppression.RawDigiProducersList[2].RawDigiProducer = 'simSiStripDigis'
process.siStripClusters.DigiProducersList[0].DigiProducer= 'simSiStripDigis'

# Event output
process.load("Configuration.EventContent.EventContent_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(20000)
)

process.load("FastSimulation/Configuration/FlatPtMuonGun_cfi")
# replace FlatRandomPtGunSource.PGunParameters.PartID={13}
process.FlatRandomPtGunSource.PGunParameters.MinPt = 50.0
process.FlatRandomPtGunSource.PGunParameters.MaxPt = 50.0
process.FlatRandomPtGunSource.PGunParameters.MinEta = -2.5
process.FlatRandomPtGunSource.PGunParameters.MaxEta = 2.5

process.FEVT = cms.OutputModule("PoolOutputModule",
    process.FEVTSIMEventContent,
    fileName = cms.untracked.string('/uscms_data/d1/cheung/slhc/testfullstrawa_muon_50GeV.root')
)

from SimGeneral.Configuration.SimGeneral_cff import *
from Validation.RecoTrack.MultiTrackValidator_cff import *
from SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cff import *
process.load("SimGeneral.TrackingAnalysis.trackingParticles_cfi")
process.load("Validation.RecoTrack.cuts_cff")

process.load("Validation.RecoTrack.cutsTPEffic_cfi")
process.load("Validation.RecoTrack.cutsTPFake_cfi")

process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")
process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")

process.load("Validation.RecoTrack.MultiTrackValidator_cff")
process.multiTrackValidator.label = ['generalTracks']
process.multiTrackValidator.associators = ['TrackAssociatorByHits']
#process.multiTrackValidator.UseAssociators = True
process.multiTrackValidator.outputFile = "validfullstrawa_muon_50GeV.root"
if (process.multiTrackValidator.label[0] == 'generalTracks'):
    process.multiTrackValidator.UseAssociators = cms.bool(False)
else:
    process.multiTrackValidator.UseAssociators = cms.bool(True)
######
process.multiTrackValidator.maxpT = cms.double(150)
process.multiTrackValidator.nintpT = cms.int32(300)

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

## these are tighter than in iterative tracking (3 and 0.3)
#process.newTrajectoryFilter.filterPset.minimumNumberOfHits = 5
#process.newTrajectoryFilter.filterPset.minPt = 0.9
## keep all tracks from first step
#process.withLooseQuality.keepAllTracks = True

### produce an ntuple with pixel hits for analysis
process.ReadLocalMeasurement = cms.EDAnalyzer("StdHitNtuplizer",
   src = cms.InputTag("siPixelRecHits"),
   stereoRecHits = cms.InputTag("siStripMatchedRecHits","stereoRecHit"),
   rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
   matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
   trackProducer = cms.InputTag("generalTracks"),
   ### if using simple (non-iterative) or old (as in 1_8_4) tracking
   #trackProducer = cms.InputTag("ctfWithMaterialTracks"),
   OutputFile = cms.string("stdgrechitfulla_ntuple.root"),
   ### for using track hit association
   associatePixel = cms.bool(True),
   associateStrip = cms.bool(False),
   associateRecoTracks = cms.bool(False),
   ROUList = cms.vstring('g4SimHitsTrackerHitsPixelBarrelLowTof', 
                         'g4SimHitsTrackerHitsPixelBarrelHighTof', 
                         'g4SimHitsTrackerHitsPixelEndcapLowTof', 
                         'g4SimHitsTrackerHitsPixelEndcapHighTof')
)

### modules to write out PixelSkimmedGeometry.txt file
#process.writedet = cms.EDProducer("SiPixelDetInfoFileWriter",
#   FilePath = cms.untracked.string("PixelSkimmedGeometry_strawb.txt")
#)

### modules to write output navigational information for tracking
#process.Tracer = cms.Service("Tracer",
#    indentation = cms.untracked.string('$$')
#)
#process.navigationSchoolAnalyzer = cms.EDAnalyzer("NavigationSchoolAnalyzer",
#    navigationSchoolName = cms.string('SimpleNavigationSchool')
#)

process.Timing =  cms.Service("Timing")

# need validation packages

process.p0 = cms.Path(process.pgen)
process.p1 = cms.Path(process.psim)
process.mydigi = cms.Sequence(cms.SequencePlaceholder("randomEngineStateProducer")*cms.SequencePlaceholder("mix")+trDigi+process.trackingParticles)
process.p2 = cms.Path(process.pdigi)
## do not do cal and muon digis
#process.p2 = cms.Path(process.mydigi)
process.p3 = cms.Path(process.L1Emulator)
process.myDigiToRaw = cms.Sequence(siPixelRawData*SiStripDigiToRaw*rawDataCollector)
#process.p4 = cms.Path(process.DigiToRaw)
process.p4 = cms.Path(process.myDigiToRaw)
process.myRawToDigi = cms.Sequence(siPixelDigis+SiStripRawToDigis)
#process.p5 = cms.Path(process.RawToDigi)
process.p5= cms.Path(process.myRawToDigi)
pixeltrackerlocalreco = cms.Sequence(siPixelClusters*siPixelRecHits)
mystriptrackerlocalreco = cms.Sequence(siStripZeroSuppression*siStripClusters*siStripMatchedRecHits)
mytrackerlocalreco = cms.Sequence(pixeltrackerlocalreco*mystriptrackerlocalreco)
process.myreco = cms.Sequence(mytrackerlocalreco+process.offlineBeamSpot+process.recopixelvertexing+process.ckftracks)
#process.p6 = cms.Path(process.reconstruction)
process.p6 = cms.Path(process.myreco)
process.mytrackMCMatchSequence = cms.Sequence(trackMCMatch*trackingParticleRecoTrackAsssociation)
process.p7 = cms.Path(process.mytrackMCMatchSequence+process.cutsTPEffic*process.cutsTPFake*process.multiTrackValidator)
#process.p7 = cms.Path(process.cutsTPEffic*process.cutsTPFake*process.multiTrackValidator)
#process.p8 = cms.Path(process.writedet)
#process.p8 = cms.Path(process.navigationSchoolAnalyzer)
process.p8 = cms.Path(process.ReadLocalMeasurement)
process.outpath = cms.EndPath(process.FEVT)
#process.schedule = cms.Schedule(process.p0,process.p1,process.p2,process.p3,process.p4,process.p5,process.p6,process.p7,process.outpath)
#process.schedule = cms.Schedule(process.p0,process.p1,process.p2,process.p4,process.p5,process.p6,process.p7,process.p8,process.outpath)
process.schedule = cms.Schedule(process.p0,process.p1,process.p2,process.p6,process.p7,process.p8,process.outpath)
