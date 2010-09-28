# Auto generated configuration file
# using: 
# Revision: 1.155 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: step2 -s RECO -n 100 --conditions MC_3XY_V9A::All --datatier GEN-SIM-RECO --eventcontent RECOSIM --beamspot Gauss --fileout file:reco.root --filein file:raw.root --python_filename RecoMuon_Fullsim_cfg.py --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('RECO')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.MixingNoPileUp_cff')
#process.load("SLHCUpgradeSimulations.Geometry.mixLowLumPU_Phase1_R39F16_cff")
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContent_cff')

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1.2.1 $'),
    annotation = cms.untracked.string('step2/phase1 nevts:100'),
    name = cms.untracked.string('PyReleaseValidation')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('OtherCMS', 
        'StdException', 
        'Unknown', 
        'BadAlloc', 
        'BadExceptionType', 
        'ProductNotFound', 
        'DictionaryNotFound', 
        'InsertFailure', 
        'Configuration', 
        'LogicError', 
        'UnimplementedFeature', 
        'InvalidReference', 
        'NullPointerError', 
        'NoProductSpecified', 
        'EventTimeout', 
        'EventCorruption', 
        'ScheduleExecutionFailure', 
        'EventProcessorFailure', 
        'FileInPathError', 
        'FileOpenError', 
        'FileReadError', 
        'FatalRootError', 
        'MismatchedInputFiles', 
        'ProductDoesNotSupportViews', 
        'ProductDoesNotSupportPtr', 
        'NotFound')
)
# Input source
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring('/store/user/cheung/phase1/r39v28/muon/TenMuon1-50GeV_RAWSIM_0pu.root'),
  duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    fileName = cms.untracked.string('file:valid_reco.root'),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    )
)
#I'm only interested in the validation stuff
process.output.outputCommands = cms.untracked.vstring('drop *','keep *_MEtoEDMConverter_*_*')


# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'MC_3XY_V9A::All'

### Phase 1 Geometry modifications ###############################################

process.load("SLHCUpgradeSimulations.Geometry.PhaseI_cmsSimIdealGeometryXML_cff")
process.Timing =  cms.Service("Timing")
process.mix.playback = True

#process.MessageLogger.destinations = cms.untracked.vstring("detailedInfo_fullph1geom")

process.load("SLHCUpgradeSimulations.Geometry.fakeConditions_Phase1_cff")
process.load("SLHCUpgradeSimulations.Geometry.recoFromSimDigis_cff")
#process.load("SLHCUpgradeSimulations.Geometry.upgradeTracking_phase1_cff")
process.load("SLHCUpgradeSimulations.Geometry.upgradeTracking_phase1_quad_wo45_cff")

process.ctfWithMaterialTracks.TTRHBuilder = 'WithTrackAngle'

process.load("RecoLocalTracker.SiPixelRecHits.PixelCPEGeneric_cfi")
process.PixelCPEGenericESProducer.Upgrade = True
process.PixelCPEGenericESProducer.SmallPitch = False
process.PixelCPEGenericESProducer.UseErrorsFromTemplates = False
process.PixelCPEGenericESProducer.TruncatePixelCharge = False
process.PixelCPEGenericESProducer.IrradiationBiasCorrection = False
process.PixelCPEGenericESProducer.DoCosmics = False
process.PixelCPEGenericESProducer.LoadTemplatesFromDB = False

## uncomment for changes for small pixel size
#process.load("SLHCUpgradeSimulations.Geometry.smallPixelSizeChanges_cff")

## CPE for other steps
process.siPixelRecHits.CPE = cms.string('PixelCPEGeneric')
process.newPixelRecHits.CPE = cms.string('PixelCPEGeneric')
process.secPixelRecHits.CPE = cms.string('PixelCPEGeneric')
process.thPixelRecHits.CPE = cms.string('PixelCPEGeneric')
process.preFilterZeroStepTracks.TTRHBuilder = cms.string('WithTrackAngle')
process.preFilterStepOneTracks.TTRHBuilder = cms.string('WithTrackAngle')
process.secWithMaterialTracks.TTRHBuilder = cms.string('WithTrackAngle')
process.thWithMaterialTracks.TTRHBuilder = cms.string('WithTrackAngle')
process.fourthWithMaterialTracks.TTRHBuilder = cms.string('WithTrackAngle')
process.fifthWithMaterialTracks.TTRHBuilder = cms.string('WithTrackAngle')

# Need these lines to stop some errors about missing siStripDigis collections. 
# should add them to fakeConditions_Phase1_cff
process.MeasurementTracker.UseStripStripQualityDB      = cms.bool(False)
process.MeasurementTracker.UsePixelModuleQualityDB     = cms.bool(False)
process.MeasurementTracker.UsePixelROCQualityDB        = cms.bool(False)
process.newMeasurementTracker.inactiveStripDetectorLabels = cms.VInputTag()
process.newMeasurementTracker.UseStripModuleQualityDB     = cms.bool(False)
process.newMeasurementTracker.UseStripAPVFiberQualityDB   = cms.bool(False)
process.newMeasurementTracker.UseStripStripQualityDB      = cms.bool(False)
process.newMeasurementTracker.UsePixelModuleQualityDB     = cms.bool(False)
process.newMeasurementTracker.UsePixelROCQualityDB        = cms.bool(False)
process.secMeasurementTracker.inactiveStripDetectorLabels = cms.VInputTag()
process.secMeasurementTracker.UseStripModuleQualityDB     = cms.bool(False)
process.secMeasurementTracker.UseStripAPVFiberQualityDB   = cms.bool(False)
process.secMeasurementTracker.UseStripStripQualityDB      = cms.bool(False)
process.secMeasurementTracker.UsePixelModuleQualityDB     = cms.bool(False)
process.secMeasurementTracker.UsePixelROCQualityDB        = cms.bool(False)
process.thMeasurementTracker.inactiveStripDetectorLabels = cms.VInputTag()
process.thMeasurementTracker.UseStripModuleQualityDB     = cms.bool(False)
process.thMeasurementTracker.UseStripAPVFiberQualityDB   = cms.bool(False)
process.thMeasurementTracker.UseStripStripQualityDB      = cms.bool(False)
process.thMeasurementTracker.UsePixelModuleQualityDB     = cms.bool(False)
process.thMeasurementTracker.UsePixelROCQualityDB        = cms.bool(False)
process.fourthMeasurementTracker.inactiveStripDetectorLabels = cms.VInputTag()
process.fifthMeasurementTracker.inactiveStripDetectorLabels = cms.VInputTag()

### Now Validation and other user functions #########################################
process.load("Validation.RecoTrack.cutsTPEffic_cfi")
process.load("Validation.RecoTrack.cutsTPFake_cfi")

process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")
process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")

process.load('Configuration.StandardSequences.Validation_cff')
### look look at OOTB generalTracks and high purity collections
### for high purity also look at 6 and 8 hit requirements
### some definitions in Validation/RecoTrack/python/TrackValidation_cff.py

import PhysicsTools.RecoAlgos.recoTrackSelector_cfi

process.cutsRecoTracksHpw6hits = PhysicsTools.RecoAlgos.recoTrackSelector_cfi.recoTrackSelector.clone()
process.cutsRecoTracksHpw6hits.quality=cms.vstring("highPurity")
process.cutsRecoTracksHpw6hits.minHit=cms.int32(6)

process.cutsRecoTracksHpw8hits = PhysicsTools.RecoAlgos.recoTrackSelector_cfi.recoTrackSelector.clone()
process.cutsRecoTracksHpw8hits.quality=cms.vstring("highPurity")
process.cutsRecoTracksHpw8hits.minHit=cms.int32(8)

process.trackValidator.label=cms.VInputTag(cms.InputTag("generalTracks"),
                                           cms.InputTag("cutsRecoTracksHp"),
                                           cms.InputTag("cutsRecoTracksHpw6hits"),
                                           cms.InputTag("cutsRecoTracksHpw8hits"),
                                           cms.InputTag("cutsRecoTracksZeroHp"),
                                           cms.InputTag("cutsRecoTracksFirstHp"),
                                           cms.InputTag("cutsRecoTracksSecondHp"),
                                           cms.InputTag("cutsRecoTracksThirdHp")
                                           )
process.trackValidator.associators = ['TrackAssociatorByHits']
process.trackValidator.UseAssociators = True
process.trackValidator.nint = cms.int32(20)
process.trackValidator.nintpT = cms.int32(25)
process.trackValidator.maxpT = cms.double(50.0)

process.slhcTracksValidation = cms.Sequence(process.cutsRecoTracksHp*
                                 process.cutsRecoTracksHpw6hits*
                                 process.cutsRecoTracksHpw8hits*
                                 process.cutsRecoTracksZeroHp*
                                 process.cutsRecoTracksFirstHp*
                                 process.cutsRecoTracksSecondHp*
                                 process.cutsRecoTracksThirdHp*
                                 process.trackValidator)
############ end John's changes ###########################
process.ReadLocalMeasurement = cms.EDAnalyzer("StdHitNtuplizer",
   src = cms.InputTag("siPixelRecHits"),
   stereoRecHits = cms.InputTag("siStripMatchedRecHits","stereoRecHit"),
   rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
   matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
   trackProducer = cms.InputTag("generalTracks"),
   ### if using simple (non-iterative) or old (as in 1_8_4) tracking
   #trackProducer = cms.InputTag("ctfWithMaterialTracks"),
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
process.anal = cms.EDAnalyzer("EventContentAnalyzer")
############# seed validator ##############################
process.load("Validation.RecoTrack.TrackerSeedValidator_cff")
process.trackerSeedValidator.min = -2.5
process.trackerSeedValidator.max = 2.5
process.trackerSeedValidator.nint = 25
process.trackerSeedValidator.useFabsEta = cms.bool(False)
process.trackerSeedValidator.maxpT = 50
process.trackerSeedValidator.nintpT = cms.int32(25)
#process.trackerSeedValidator.maxpT = cms.double(400)
#process.trackerSeedValidator.nintpT = cms.int32(800)
#process.trackerSeedValidator.label = cms.VInputTag(cms.InputTag("newSeedFromTriplets"))
#process.trackerSeedValidator.label = cms.VInputTag(cms.InputTag("pixelTriplets"))
#process.trackerSeedValidator.label = ['newSeedFromTriplets','newSeedFromPairs','newCombinedSeeds']
process.trackerSeedValidator.label = ['newCombinedSeeds','newSeedFromTriplets','stepOneSeedFromTriplets',
                                      'secTriplets','thSeedFromTriplets']
process.trackerSeedValidator.label_tp_effic = cms.InputTag("cutsTPEffic")
process.trackerSeedValidator.label_tp_fake = cms.InputTag("cutsTPFake")
process.trackerSeedValidator.useLogPt=cms.untracked.bool(False)
## missing (and not used apparently)? (should put this into TrackerSeedValidator_cff.py)
process.trackerSeedValidator.parametersDefiner = cms.string('LhcParametersDefinerForTP')
process.trackerSeedValidator.minVertpos = cms.double(0)
process.trackerSeedValidator.maxVertpos = cms.double(5)
process.trackerSeedValidator.nintVertpos = cms.int32(100)
process.trackerSeedValidator.minZpos = cms.double(-10)
process.trackerSeedValidator.maxZpos = cms.double(10)
process.trackerSeedValidator.nintZpos = cms.int32(100)

### back to standard commands
process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
process.Refitter = process.TrackRefitter.clone()
#process.Refitter.src = 'ctfWithMaterialTracks'
process.Refitter.src = 'generalTracks'
process.Refitter.TTRHBuilder= cms.string("WithTrackAngle")

process.load("SLHCUpgradeSimulations.Validation.SLHCPixelHitAnalyzer_cfi")
process.SLHCPixelHitAnalyzer.trajectoryInput = cms.string('Refitter')
process.SLHCPixelHitAnalyzer.useAllPixel = cms.bool(False)
process.SLHCPixelHitAnalyzer.isCosmic = cms.bool(False)
process.SLHCPixelHitAnalyzer.isSim = cms.bool(True)
process.SLHCPixelHitAnalyzer.OutputFile = cms.string('pixelhitdata.root')
process.SLHCPixelHitAnalyzer.trackProducer = cms.InputTag("Refitter")
#process.load("SLHCUpgradeSimulations.Validation.SLHCAllPixelHitAnalyzer_cfi")
#process.SLHCAllPixelHitAnalyzer.trajectoryInput = cms.string('Refitter')
#process.SLHCAllPixelHitAnalyzer.useAllPixel = cms.bool(False)
#process.SLHCAllPixelHitAnalyzer.isCosmic = cms.bool(False)
#process.SLHCAllPixelHitAnalyzer.isSim = cms.bool(True)
#process.SLHCAllPixelHitAnalyzer.trajectoryInput = cms.InputTag("ctfWithMaterialTracks")
#process.SLHCAllPixelHitAnalyzer.OutputFile = cms.string('allpixelhitdata.root')

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("pixelhitdata.root"),
                                   closeFileFast = cms.untracked.bool(True)
                                   )

process.load("RecoVertex.Configuration.RecoVertex_cff")
#process.offlinePrimaryVertices.TrackLabel = cms.InputTag("ctfWithMaterialTracks")
#process.offlinePrimaryVerticesWithBS.TrackLabel = cms.InputTag("ctfWithMaterialTracks")
#process.generalV0Candidates.trackRecoAlgorithm = cms.InputTag('ctfWithMaterialTracks')

process.load("SLHCUpgradeSimulations.Validation.SLHCPrimaryVertexAnalyzer_cfi")
#process.slhcSimpleVertexAnalysis.recoTrackProducer = cms.untracked.string('ctfWithMaterialTracks')
process.slhcSimpleVertexAnalysis.recoTrackProducer = cms.untracked.string('generalTracks')

process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
     oncePerEventMode = cms.untracked.bool(True),
     showMallocInfo = cms.untracked.bool(True),
     moduleMemorySummary = cms.untracked.bool(True),
     ignoreTotal = cms.untracked.int32(5)
)

# Path and EndPath definitions
#process.reconstruction_step = cms.Path(process.reconstruction)
process.mix_step = cms.Path(process.mix)
process.reconstruction_step = cms.Path(process.trackerlocalreco*process.offlineBeamSpot+process.recopixelvertexing*process.ckftracks_wodEdXandSteps4and5)
process.debug_step = cms.Path(process.anal)
#process.validation_step = cms.Path(process.cutsTPEffic*process.cutsTPFake*process.trackerSeedValidator*process.trackValidator)
process.validation_step = cms.Path(process.cutsTPEffic*process.cutsTPFake*process.trackerSeedValidator*process.slhcTracksValidation)
process.user_step0 = cms.Path(process.vertexreco*process.slhcSimpleVertexAnalysis*process.Refitter*process.SLHCPixelHitAnalyzer)
#process.user_step0 = cms.Path(process.Refitter*process.SLHCAllPixelHitAnalyzer)
process.user_step = cms.Path(process.ReadLocalMeasurement)
process.endjob_step = cms.Path(process.endOfProcess)
process.out_step = cms.EndPath(process.output)

# Schedule definition
#process.schedule = cms.Schedule(process.mix_step,process.reconstruction_step,process.user_step0,process.endjob_step)
#process.schedule = cms.Schedule(process.mix_step,process.reconstruction_step,process.user_step,process.endjob_step)
#process.schedule = cms.Schedule(process.mix_step,process.reconstruction_step,process.validation_step,process.endjob_step)
process.schedule = cms.Schedule(process.mix_step,process.reconstruction_step,process.validation_step,process.endjob_step,process.out_step)
