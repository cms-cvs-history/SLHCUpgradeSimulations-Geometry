import FWCore.ParameterSet.Config as cms

process = cms.Process("Rec")
# this example configuration offers some minimum 
# annotation, to help users get through; please
# don't hesitate to read through the comments
# use MessageLogger to redirect/suppress multiple
# service messages coming from the system
#
# in this config below, we use the replace option to make
# the logger let out messages of severity ERROR (INFO level
# will be suppressed), and we want to limit the number to 10
#
process.load("Configuration.StandardSequences.Services_cff")

process.load("Configuration.StandardSequences.Geometry_cff")

#PhaseI Geometry
process.load("SLHCUpgradeSimulations.Geometry.PhaseI_cmsIdealGeometryXML_cff")

process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("Configuration.StandardSequences.FakeConditions_cff")

process.SiPixelFakeGainOfflineESSource.file = 'SLHCUpgradeSimulations/Geometry/data/PhaseI/PixelSkimmedGeometry.txt'
process.SiPixelFakeLorentzAngleESSource.file = 'SLHCUpgradeSimulations/Geometry/data/PhaseI/PixelSkimmedGeometry.txt'

process.load("FWCore.MessageService.MessageLogger_cfi")

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

    
# Event output
process.load("Configuration.EventContent.EventContent_cff")


process.Timing = cms.Service("Timing")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(20000)
)
process.source = cms.Source("FlatRandomPtGunSource",
    PGunParameters = cms.untracked.PSet(
        # you can request more than 1 particle
        # since PartID is a vector, you can place in as many 
        # PDG id's as you wish, comma seaparated
        #
        PartID = cms.untracked.vint32(13),
        MaxEta = cms.untracked.double(2.5),
        MaxPhi = cms.untracked.double(3.14159265359),
        MaxPt  = cms.untracked.double(100.0),
        MinEta = cms.untracked.double(-2.5),
        MinPhi = cms.untracked.double(-3.14159265359), ## in radians
        MinPt  = cms.untracked.double(100.0)
    ),
    Verbosity = cms.untracked.int32(1), ## set to 1 (or greater)  for printouts

    AddAntiParticle = cms.untracked.bool(True) ## back-to-back particles

)

# Full Reco
process.load("RecoLocalTracker.Configuration.RecoLocalTracker_cff")
process.load("RecoTracker.Configuration.RecoTracker_cff")
process.load("TrackingTools.Configuration.TrackingTools_cff")

# Stuff needed to produce efficiency plots
from SimGeneral.Configuration.SimGeneral_cff import *  ####
from Validation.RecoTrack.MultiTrackValidator_cff import *   ####
from SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cff import *   ####

process.load("SimGeneral.TrackingAnalysis.trackingParticles_cfi")  ####
process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")

# Uncomment only if digis are available
process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")

    # Filters
process.load("Validation.RecoTrack.cuts_cff")

    # Track Validator    
process.load("Validation.RecoTrack.MultiTrackValidator_cff")
#process.multiTrackValidator.label = ['cutsRecoTracks']
process.multiTrackValidator.label = ['generalTracks']      ####

if (process.multiTrackValidator.label[0] == 'generalTracks'):
    process.multiTrackValidator.UseAssociators = cms.bool(False)
else:
    process.multiTrackValidator.UseAssociators = cms.bool(True)
######
process.multiTrackValidator.outputFile = cms.string('TrackVal.root')
process.multiTrackValidator.maxpT = cms.double(150)
process.multiTrackValidator.nintpT = cms.int32(300)
#process.multiTrackValidator.associatormap = cms.InputTag("assoc2GsfTracks")
process.load("SimTracker.SiPixelDigitizer.PixelDigi_cfi")
process.simSiPixelDigis.MissCalibrate = cms.bool(False)
process.simSiPixelDigis.AddPixelInefficiency = cms.int32(-1)
process.simSiPixelDigis.LorentzAngle_DB = cms.bool(False)
process.simSiPixelDigis.killModules = cms.bool(False)
process.simSiPixelDigis.NumPixelBarrel = cms.int32(4)


process.load("RecoLocalTracker.SiPixelClusterizer.SiPixelClusterizer_cfi")
from RecoLocalTracker.SiPixelClusterizer.SiPixelClusterizer_cfi import *
process.siPixelClusters.MissCalibrate = cms.untracked.bool(False)
process.siPixelClusters.src = cms.InputTag("simSiPixelDigis")

process.load("RecoLocalTracker.SiStripClusterizer.SiStripClusterizer_SimData2_cfi")
from RecoLocalTracker.SiStripClusterizer.SiStripClusterizer_SimData2_cfi import *
#process.siStripClusters.DigiProducer = cms.string('simSiStripDigis')

process.load("RecoLocalTracker.SiStripZeroSuppression.SiStripZeroSuppression_SimData_cfi")
from RecoLocalTracker.SiStripZeroSuppression.SiStripZeroSuppression_SimData_cfi import *

from RecoLocalTracker.SiPixelRecHits.SiPixelRecHits_cfi import *
from RecoLocalTracker.SiStripRecHitConverter.SiStripRecHitConverter_cfi import *
from RecoLocalTracker.SiStripRecHitConverter.SiStripRecHitMatcher_cfi import *
from RecoLocalTracker.SiStripRecHitConverter.StripCPEfromTrackAngle_cfi import *

process.load("SimTracker.TrackAssociation.trackMCMatchSequence_cff")
from SimTracker.TrackAssociation.trackMCMatchSequence_cff import *

 # standard "prescription of what to keep in edm::Event upon output
   #

process.FEVT = cms.OutputModule("PoolOutputModule",
    process.FEVTSIMEventContent,
    fileName = cms.untracked.string('PhysVal.root')
)

process.p0 = cms.Path(process.pgen)
process.p1 = cms.Path(process.psim)
process.mydigi = cms.Sequence(cms.SequencePlaceholder("randomEngineStateProducer")*cms.SequencePlaceholder("mix")+trDigi+process.trackingParticles)
process.p2 = cms.Path(process.mydigi)
process.p3 = cms.Path(process.L1Emulator)
process.myDigiToRaw = cms.Sequence(siPixelRawData*SiStripDigiToRaw*rawDataCollector)
process.p4 = cms.Path(process.myDigiToRaw)
process.myRawToDigi = cms.Sequence(siPixelDigis+SiStripRawToDigis)
process.p5= cms.Path(process.myRawToDigi)
pixeltrackerlocalreco = cms.Sequence(siPixelClusters*siPixelRecHits)
mystriptrackerlocalreco = cms.Sequence(siStripZeroSuppression*siStripClusters*siStripMatchedRecHits)
mytrackerlocalreco = cms.Sequence(pixeltrackerlocalreco*mystriptrackerlocalreco)
process.myreco = cms.Sequence(mytrackerlocalreco+process.offlineBeamSpot+process.recopixelvertexing+process.ckftracks)
process.p6 = cms.Path(process.myreco)
process.mytrackMCMatchSequence = cms.Sequence(trackMCMatch*trackingParticleRecoTrackAsssociation)
process.p7 = cms.Path(process.mytrackMCMatchSequence+process.cutsTPEffic*process.cutsTPFake*process.multiTrackValidator)
process.outpath = cms.EndPath(process.FEVT)

process.schedule = cms.Schedule(process.p0,process.p1,process.p2,process.p4,process.p5,process.p6,process.p7)
