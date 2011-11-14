# Auto generated configuration file
# using: 
# Revision: 1.303.2.3 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: SLHCUpgradeSimulations/Configuration/python/FourMuPt_1_50_cfi.py -s GEN,FASTSIM,HLT:GRun --pileup=NoPileUp --geometry DB -n 10 --conditions auto:mc --eventcontent FEVTDEBUG --datatier GEN-SIM-DIGI-RECO --beamspot Gauss --no_exec --python_filename FASTSIM_4muons_cfg.py

import FWCore.ParameterSet.Config as cms

process = cms.Process('FASTSIMWDIGI')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('FastSimulation.Configuration.EventContent_cff')
process.load('SLHCUpgradeSimulations.Geometry.mixLowLumPU_FastSim14TeV_cff')
#process.load('FastSimulation.PileUpProducer.PileUpSimulator_NoPileUp_cff')
process.load('FastSimulation.Configuration.Geometries_MC_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('FastSimulation.Configuration.FamosSequences_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedParameters_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# Include the RandomNumberGeneratorService definition
process.load("FastSimulation.Configuration.RandomServiceInitialization_cff")

process.RandomNumberGeneratorService.simSiStripDigis = cms.PSet(
      initialSeed = cms.untracked.uint32(1234567),
      engineName = cms.untracked.string('HepJamesRandom'))
process.RandomNumberGeneratorService.simSiPixelDigis = cms.PSet(
      initialSeed = cms.untracked.uint32(1234567),
      engineName = cms.untracked.string('HepJamesRandom'))

process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Configuration.StandardSequences.EndOfProcess_cff')

process.load("SLHCUpgradeSimulations.Geometry.fakeConditions_stdgeom_cff")
process.load("SLHCUpgradeSimulations.Geometry.recoFromSimDigis_cff")
process.load("SLHCUpgradeSimulations.Geometry.upgradeTracking_stdgeom_cff")

process.ctfWithMaterialTracks.TTRHBuilder = 'WithTrackAngle'
process.PixelCPEGenericESProducer.UseErrorsFromTemplates = cms.bool(True)  #FG set True to use errors from templates
process.PixelCPEGenericESProducer.TruncatePixelCharge = cms.bool(False)
process.PixelCPEGenericESProducer.LoadTemplatesFromDB = cms.bool(True)  #FG set True to load the last version of the templates
process.PixelCPEGenericESProducer.IrradiationBiasCorrection = False
process.PixelCPEGenericESProducer.DoCosmics = False

## CPE for other steps
process.siPixelRecHits.CPE = cms.string('PixelCPEGeneric')
process.initialStepTracks.TTRHBuilder = cms.string('WithTrackAngle')
process.lowPtTripletStepTracks.TTRHBuilder = cms.string('WithTrackAngle')
process.pixelPairStepTracks.TTRHBuilder = cms.string('WithTrackAngle')
process.detachedTripletStepTracks.TTRHBuilder = cms.string('WithTrackAngle')
process.mixedTripletStepTracks.TTRHBuilder = cms.string('WithTrackAngle')
process.pixelLessStepTracks.TTRHBuilder = cms.string('WithTrackAngle')
process.tobTecStepTracks.TTRHBuilder = cms.string('WithTrackAngle')

# Need these lines to stop some errors about missing siStripDigis collections.
# should add them to fakeConditions_Phase1_cff
process.MeasurementTracker.inactiveStripDetectorLabels = cms.VInputTag()
process.MeasurementTracker.UseStripModuleQualityDB     = cms.bool(False)
process.MeasurementTracker.UseStripAPVFiberQualityDB   = cms.bool(False)
process.MeasurementTracker.UseStripStripQualityDB      = cms.bool(False)
process.MeasurementTracker.UsePixelModuleQualityDB     = cms.bool(False)
process.MeasurementTracker.UsePixelROCQualityDB        = cms.bool(False)
process.lowPtTripletStepMeasurementTracker.inactiveStripDetectorLabels = cms.VInputTag()
process.lowPtTripletStepMeasurementTracker.UseStripModuleQualityDB     = cms.bool(False)
process.lowPtTripletStepMeasurementTracker.UseStripAPVFiberQualityDB   = cms.bool(False)
process.lowPtTripletStepMeasurementTracker.UseStripStripQualityDB      = cms.bool(False)
process.lowPtTripletStepMeasurementTracker.UsePixelModuleQualityDB     = cms.bool(False)
process.lowPtTripletStepMeasurementTracker.UsePixelROCQualityDB        = cms.bool(False)
process.pixelPairStepMeasurementTracker.inactiveStripDetectorLabels = cms.VInputTag()
process.pixelPairStepMeasurementTracker.UseStripModuleQualityDB     = cms.bool(False)
process.pixelPairStepMeasurementTracker.UseStripAPVFiberQualityDB   = cms.bool(False)
process.pixelPairStepMeasurementTracker.UseStripStripQualityDB      = cms.bool(False)
process.pixelPairStepMeasurementTracker.UsePixelModuleQualityDB     = cms.bool(False)
process.pixelPairStepMeasurementTracker.UsePixelROCQualityDB        = cms.bool(False)
process.detachedTripletStepMeasurementTracker.inactiveStripDetectorLabels = cms.VInputTag()
process.detachedTripletStepMeasurementTracker.UseStripModuleQualityDB     = cms.bool(False)
process.detachedTripletStepMeasurementTracker.UseStripAPVFiberQualityDB   = cms.bool(False)
process.detachedTripletStepMeasurementTracker.UseStripStripQualityDB      = cms.bool(False)
process.detachedTripletStepMeasurementTracker.UsePixelModuleQualityDB     = cms.bool(False)
process.detachedTripletStepMeasurementTracker.UsePixelROCQualityDB        = cms.bool(False)
process.mixedTripletStepMeasurementTracker.inactiveStripDetectorLabels = cms.VInputTag()
process.mixedTripletStepMeasurementTracker.UseStripModuleQualityDB     = cms.bool(False)
process.mixedTripletStepMeasurementTracker.UseStripAPVFiberQualityDB   = cms.bool(False)
process.mixedTripletStepMeasurementTracker.UseStripStripQualityDB      = cms.bool(False)
process.mixedTripletStepMeasurementTracker.UsePixelModuleQualityDB     = cms.bool(False)
process.mixedTripletStepMeasurementTracker.UsePixelROCQualityDB        = cms.bool(False)
process.pixelLessStepMeasurementTracker.inactiveStripDetectorLabels = cms.VInputTag()
process.tobTecStepMeasurementTracker.inactiveStripDetectorLabels = cms.VInputTag()

## for fastsim we need these ################################
process.TrackerGeometricDetESModule.fromDDD=cms.bool(True)
process.TrackerDigiGeometryESModule.fromDDD=cms.bool(True)
process.simSiPixelDigis.ROUList =  ['famosSimHitsTrackerHits']
process.simSiStripDigis.ROUList =  ['famosSimHitsTrackerHits']
process.load("SimGeneral.TrackingAnalysis.trackingParticles_cfi")
process.mergedtruth.simHitCollections.tracker = ['famosSimHitsTrackerHits']
process.mergedtruth.simHitCollections.pixel = []
process.mergedtruth.simHitCollections.muon = []
process.mergedtruth.simHitLabel = 'famosSimHits'
## make occupancies more similar to full simulation
process.famosSimHits.ParticleFilter.etaMax = 3.0
process.famosSimHits.ParticleFilter.pTMin = 0.05
process.famosSimHits.TrackerSimHits.pTmin = 0.05
process.famosSimHits.TrackerSimHits.firstLoop = False
#############################################################
process.Timing =  cms.Service("Timing")

# If you want to turn on/off pile-up, default is no pileup
process.famosPileUp.PileUpSimulator.averageNumber = 50.00
### if doing inefficiency at <PU>=50
process.simSiPixelDigis.AddPixelInefficiency = 20

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.3 $'),
    annotation = cms.untracked.string('SLHCUpgradeSimulations/Configuration/python/FourMuPt_1_50_cfi.py nevts:10'),
    name = cms.untracked.string('PyReleaseValidation')
)

# Output definition
extendAOD = cms.untracked.vstring('keep *_MEtoEDMConverter_*_*')
process.AODSIMEventContent.outputCommands.extend(extendAOD)
process.output = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    #outputCommands = process.RECOSIMEventContent.outputCommands,
    outputCommands = process.AODSIMEventContent.outputCommands,
    #fileName = cms.untracked.string('file:reco.root'),
    fileName = cms.untracked.string('file:recoAODSIM.root'),
    dataset = cms.untracked.PSet(
        #dataTier = cms.untracked.string('GEN-SIM-RECO'),
        dataTier = cms.untracked.string('AODSIM'),
        filterName = cms.untracked.string('')
    )
)

# Other statements
process.GlobalTag.globaltag = 'DESIGN42_V11::All'

process.famosSimHits.SimulateCalorimetry = True
process.famosSimHits.SimulateTracking = True
process.GaussVtxSmearingParameters.type = cms.string("Gaussian")
process.famosSimHits.VertexGenerator = process.GaussVtxSmearingParameters
process.famosPileUp.VertexGenerator = process.GaussVtxSmearingParameters

####################
process.load("SLHCUpgradeSimulations.Configuration.TTbar_Tauola_14TeV_cfi")
#process.load("SLHCUpgradeSimulations.Configuration.HERWIGPP_POWHEG_H120_bbbar_Z_ll_14TeV_cff")
#process.load("SLHCUpgradeSimulations.Configuration.ZZ_MMorBB_TuneZ2_14TeV_pythia6_tauola_cff")
#process.load("SLHCUpgradeSimulations.Configuration.ZMM_14TeV_cfi")
#process.load("SLHCUpgradeSimulations.Configuration.FourMuPt_1_200_cfi")
##########################################################

# Make the job crash in case of missing product
process.options = cms.untracked.PSet( Rethrow = cms.untracked.vstring('ProductNotFound') )

process.anal = cms.EDAnalyzer("EventContentAnalyzer")


process.load('FastSimulation.CaloRecHitsProducer.CaloRecHits_cff')
from FastSimulation.CaloRecHitsProducer.CaloRecHits_cff import *
from RecoLocalCalo.HcalRecAlgos.hcalRecAlgoESProd_cfi import *
# Calo Towers
from RecoJets.Configuration.CaloTowersRec_cff import *

process.load('RecoTracker.Configuration.RecoTracker_cff')
# Calo RecHits producer (with no HCAL miscalibration by default)

# Muon RecHit sequence
from RecoLocalMuon.Configuration.RecoLocalMuon_cff import *
csc2DRecHits.stripDigiTag = cms.InputTag("simMuonCSCDigis","MuonCSCStripDigi")
csc2DRecHits.wireDigiTag = cms.InputTag("simMuonCSCDigis","MuonCSCWireDigi")
rpcRecHits.rpcDigiLabel = 'simMuonRPCDigis'
dt1DRecHits.dtDigiLabel = 'simMuonDTDigis'
dt1DCosmicRecHits.dtDigiLabel = 'simMuonDTDigis'

# Muon reconstruction sequence
from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
from RecoMuon.TrackingTools.MuonTrackLoader_cff import *
KFSmootherForMuonTrackLoader.Propagator = 'SmartPropagatorAny'
from RecoMuon.MuonSeedGenerator.standAloneMuonSeeds_cff import *
from RecoMuon.StandAloneMuonProducer.standAloneMuons_cff import *
from FastSimulation.Configuration.globalMuons_cff import *
globalMuons.GLBTrajBuilderParameters.TrackTransformer.TrackerRecHitBuilder = 'WithoutRefit'
globalMuons.GLBTrajBuilderParameters.TrackerRecHitBuilder = 'WithoutRefit'
globalMuons.GLBTrajBuilderParameters.TransformerOutPropagator = cms.string('SmartPropagatorAny')
globalMuons.GLBTrajBuilderParameters.MatcherOutPropagator = cms.string('SmartPropagator')

from RecoMuon.GlobalMuonProducer.tevMuons_cfi import *
GlobalMuonRefitter.TrackerRecHitBuilder = 'WithoutRefit'
GlobalMuonRefitter.Propagator = 'SmartPropagatorAny'
GlobalTrajectoryBuilderCommon.TrackerRecHitBuilder = 'WithoutRefit'
tevMuons.RefitterParameters.TrackerRecHitBuilder = 'WithoutRefit'
tevMuons.RefitterParameters.Propagator =  'SmartPropagatorAny'
KFSmootherForRefitInsideOut.Propagator = 'SmartPropagatorAny'
KFSmootherForRefitOutsideIn.Propagator = 'SmartPropagator'
KFFitterForRefitInsideOut.Propagator = 'SmartPropagatorAny'
KFFitterForRefitOutsideIn.Propagator = 'SmartPropagatorAny'


#from RecoEgamma.EgammaElectronProducers.electronSequence_cff import *
#from RecoEgamma.EgammaPhotonProducers.photonSequence_cff import *
#from RecoEgamma.EgammaPhotonProducers.conversionSequence_cff import *
#from RecoEgamma.EgammaPhotonProducers.conversionTrackSequence_cff import *
#from RecoEgamma.EgammaPhotonProducers.allConversionSequence_cff import *
#allConversions.src = 'gsfGeneralConversionTrackMerger'
########################################

# Famos with tracks
process.p0 = cms.Path(process.generator)
process.generation_step = cms.Path(process.pgen_genonly)
process.othergeneration_step = cms.Path(process.GeneInfo+process.genJetMET)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)

process.p1 = cms.Path(process.famosWithTrackerAndCaloHits)
process.p2 = cms.Path(process.trDigi*process.trackingParticles)
process.reconstruction_step     = cms.Path(process.trackerlocalreco*
                                           process.offlineBeamSpot+
                                           process.recopixelvertexing*
                                           process.ckftracks_wodEdX*process.trackExtrapolator*                                            
                                           process.particleFlowCluster*
                                           process.ecalClusters*
                                           process.caloTowersRec*
                                           process.vertexreco*
###                                           process.egammaGlobalReco*
                                           process.electronGsfTracking*process.conversionTrackSequence*process.conversionTrackSequenceNoEcalSeeded* 
                                           process.allConversionSequence*
                                           process.pfTrackingGlobalReco*
                                           process.jetGlobalReco*
                                           process.famosMuonSequence*
                                           process.famosMuonIdAndIsolationSequence*
###                                           process.highlevelreco
                                           process.egammaHighLevelRecoPrePF*
                                           process.particleFlowReco*
                                           process.egammaHighLevelRecoPostPF*
                                           process.jetHighLevelReco*
                                           process.tautagging*
###                                           process.metrecoPlusHCALNoise*
                                           process.btagging*
                                           process.recoPFMET*
                                           process.PFTau*
                                           process.regionalCosmicTracksSeq*
###                                           process.muoncosmichighlevelreco*
                                           process.reducedRecHits
)



process.p7 = cms.Path(process.anal)

process.endjob_step             = cms.Path(process.endOfProcess)
process.out_step                = cms.EndPath(process.output)

process.schedule = cms.Schedule(process.generation_step,process.othergeneration_step,process.genfiltersummary_step,process.p1,process.p2,process.reconstruction_step,process.endjob_step,process.out_step)
# filter all path with the production filter sequence
for path in process.paths:
        getattr(process,path)._seq = process.generator * getattr(process,path)._seq

