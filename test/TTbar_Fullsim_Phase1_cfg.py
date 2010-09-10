# Auto generated configuration file
# using: 
# Revision: 1.172.2.5 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: Configuration/Generator/python/SingleMuPt10_cfi.py -s GEN,SIM,DIGI,L1,DIGI2RAW,HLT,RAW2DIGI,L1Reco -n 100 --conditions DESIGN_36_V10::All --datatier GEN-SIM-RAW --eventcontent RAWSIM --beamspot Gauss --python_filename Muon_Fullsim_cfg.py --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('HLT')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.MixingNoPileUp_cff')
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('Configuration.StandardSequences.VtxSmearedGauss_cff')
process.load('Configuration.StandardSequences.Sim_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('HLTrigger.Configuration.HLT_1E31_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContent_cff')

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.172.2.5 $'),
    annotation = cms.untracked.string('Configuration/Generator/python/TTbar_cfi.py Altered nevts:100'),
    name = cms.untracked.string('PyReleaseValidation')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(25)
)
process.options = cms.untracked.PSet(

)
# Input source
process.source = cms.Source("EmptySource")

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('TTbar_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_RAW2DIGI_L1Reco.root'),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RAW'),
        filterName = cms.untracked.string('')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'DESIGN_36_V10::All'
process.generator = cms.EDFilter("Pythia6GeneratorFilter",
    pythiaPylistVerbosity = cms.untracked.int32(0),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(14000.0),
    maxEventsToPrint = cms.untracked.int32(0),
    PythiaParameters = cms.PSet(
        pythiaUESettings = cms.vstring('MSTJ(11)=3     ! Choice of the fragmentation function',
            'MSTJ(22)=2     ! Decay those unstable particles',
            'PARJ(71)=10 .  ! for which ctau  10 mm',
            'MSTP(2)=1      ! which order running alphaS',
            'MSTP(33)=0     ! no K factors in hard cross sections',
            'MSTP(51)=10042 ! structure function chosen (external PDF CTEQ6L1)',
            'MSTP(52)=2     ! work with LHAPDF',
            'MSTP(81)=1     ! multiple parton interactions 1 is Pythia default',
            'MSTP(82)=4     ! Defines the multi-parton model',
            'MSTU(21)=1     ! Check on possible errors during program execution',
            'PARP(82)=1.8387   ! pt cutoff for multiparton interactions',
            'PARP(89)=1960. ! sqrts for which PARP82 is set',
            'PARP(83)=0.5   ! Multiple interactions: matter distrbn parameter',
            'PARP(84)=0.4   ! Multiple interactions: matter distribution parameter',
            'PARP(90)=0.16  ! Multiple interactions: rescaling power',
            'PARP(67)=2.5    ! amount of initial-state radiation',
            'PARP(85)=1.0  ! gluon prod. mechanism in MI',
            'PARP(86)=1.0  ! gluon prod. mechanism in MI',
            'PARP(62)=1.25   ! ',
            'PARP(64)=0.2    ! ',
            'MSTP(91)=1      !',
            'PARP(91)=2.1   ! kt distribution',
            'PARP(93)=15.0  ! '),
        processParameters = cms.vstring('MSEL      = 0     ! User defined processes',
            'MSUB(81)  = 1     ! qqbar to QQbar',
            'MSUB(82)  = 1     ! gg to QQbar',
            'MSTP(7)   = 6     ! flavour = top',
            'PMAS(6,1) = 175.  ! top quark mass'),
        parameterSets = cms.vstring('pythiaUESettings',
            'processParameters')
    )
)

### PhaseI Geometry and modifications ###############################################
process.load("SLHCUpgradeSimulations.Geometry.PhaseI_cmsSimIdealGeometryXML_cff")
process.Timing =  cms.Service("Timing")
from SimTracker.Configuration.SimTracker_EventContent_cff import *
from SimGeneral.Configuration.SimGeneral_EventContent_cff import *

process.RAWSIMEventContent.outputCommands.extend(SimGeneralFEVTDEBUG.outputCommands)
process.RAWSIMEventContent.outputCommands.extend(SimTrackerFEVTDEBUG.outputCommands)

process.siPixelFakeGainOfflineESSource = cms.ESSource("SiPixelFakeGainOfflineESSource",
  file = cms.FileInPath('SLHCUpgradeSimulations/Geometry/data/PhaseI/EmptyPixelSkimmedGeometry_phase1.txt')
)
process.es_prefer_fake_gain = cms.ESPrefer("SiPixelFakeGainOfflineESSource","siPixelFakeGainOfflineESSource")

process.siPixelFakeLorentzAngleESSource = cms.ESSource("SiPixelFakeLorentzAngleESSource",
  file = cms.FileInPath('SLHCUpgradeSimulations/Geometry/data/PhaseI/PixelSkimmedGeometry_phase1.txt')
)
process.es_prefer_fake_lorentz = cms.ESPrefer("SiPixelFakeLorentzAngleESSource","siPixelFakeLorentzAngleESSource")

process.load("CalibTracker.SiStripESProducers.fake.SiStripNoisesFakeESSource_cfi")
# Next two lines are antiquated in form, and we use the proper values by default
#process.SiStripNoisesGenerator.NoiseStripLengthSlope=51. #dec mode
#process.SiStripNoisesGenerator.NoiseStripLengthQuote=630.

process.siStripNoisesFakeESSource  = cms.ESSource("SiStripNoisesFakeESSource")
process.es_prefer_fake_strip_noise = cms.ESPrefer("SiStripNoisesFakeESSource",
                                                  "siStripNoisesFakeESSource")

process.load("CalibTracker.SiStripESProducers.fake.SiStripQualityFakeESSource_cfi")

process.siStripQualityFakeESSource  = cms.ESSource("SiStripQualityFakeESSource")
process.es_prefer_fake_strip_quality = cms.ESPrefer("SiStripQualityFakeESSource",
                                                     "siStripQualityFakeESSource")

process.load("CalibTracker.SiStripESProducers.fake.SiStripPedestalsFakeESSource_cfi")

process.siStripPedestalsFakeESSource  = cms.ESSource("SiStripPedestalsFakeESSource")
process.es_prefer_fake_strip_pedestal = cms.ESPrefer("SiStripPedestalsFakeESSource",
                                                     "siStripPedestalsFakeESSource")

process.load("CalibTracker.SiStripESProducers.fake.SiStripLorentzAngleFakeESSource_cfi")

process.siStripLorentzAngleFakeESSource  = cms.ESSource("SiStripLorentzAngleFakeESSource")
process.es_prefer_fake_strip_LA = cms.ESPrefer("SiStripLorentzAngleFakeESSource",
                                               "siStripLorentzAngleFakeESSource")

process.siStripLorentzAngleSimFakeESSource  = cms.ESSource("SiStripLorentzAngleSimFakeESSource")
process.es_prefer_fake_strip_LA_sim = cms.ESPrefer("SiStripLorentzAngleSimFakeESSource",
                                                   "siStripLorentzAngleSimFakeESSource")

process.load("CalibTracker.SiStripESProducers.fake.SiStripApvGainFakeESSource_cfi")
process.SiStripApvGainGenerator.MeanGain=1.0
process.SiStripApvGainGenerator.SigmaGain=0.0
process.SiStripApvGainGenerator.genMode = cms.string("default")
process.myStripApvGainFakeESSource = cms.ESSource("SiStripApvGainFakeESSource")
process.es_prefer_myStripApvGainFakeESSource  = cms.ESPrefer("SiStripApvGainFakeESSource",
                                                  "myStripApvGainFakeESSource")

process.myStripApvGainSimFakeESSource  = cms.ESSource("SiStripApvGainSimFakeESSource")
process.es_prefer_myStripApvGainSimFakeESSource = cms.ESPrefer("SiStripApvGainSimFakeESSource",
                                                               "myStripApvGainSimFakeESSource")

process.load("CalibTracker.SiStripESProducers.fake.SiStripThresholdFakeESSource_cfi")

process.siStripThresholdFakeESSource  = cms.ESSource("SiStripThresholdFakeESSource")
process.es_prefer_fake_strip_threshold = cms.ESPrefer("SiStripThresholdFakeESSource",
                                                     "siStripThresholdFakeESSource")

process.TrackerDigiGeometryESModule.applyAlignment = False
process.simSiPixelDigis.MissCalibrate = False
process.simSiPixelDigis.LorentzAngle_DB = False
process.simSiPixelDigis.killModules = False
process.simSiPixelDigis.useDB = False
process.simSiPixelDigis.DeadModules_DB = False
process.simSiPixelDigis.NumPixelBarrel = cms.int32(4)
process.simSiPixelDigis.NumPixelEndcap = cms.int32(3)
process.simSiPixelDigis.AddPixelInefficiency = -1

process.MeasurementTracker.inactiveStripDetectorLabels = cms.VInputTag()
process.MeasurementTracker.UseStripModuleQualityDB     = cms.bool(False)
process.MeasurementTracker.UseStripAPVFiberQualityDB   = cms.bool(False)

process.mergedtruth.volumeRadius = cms.double(100.0)
process.mergedtruth.volumeZ = cms.double(900.0)
process.mergedtruth.discardOutVolume = cms.bool(True)

### back to standard job commands ##################################################

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.endjob_step = cms.Path(process.endOfProcess)
process.out_step = cms.EndPath(process.output)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.simulation_step,process.digitisation_step,process.L1simulation_step)
#process.schedule = cms.Schedule(process.generation_step,process.simulation_step,process.digitisation_step,process.L1simulation_step,process.digi2raw_step)
#process.schedule.extend(process.HLTSchedule)
#process.schedule.extend([process.raw2digi_step,process.L1Reco_step,process.endjob_step,process.out_step])
process.schedule.extend([process.endjob_step,process.out_step])
# special treatment in case of production filter sequence  
for path in process.paths: 
    getattr(process,path)._seq = process.generator*getattr(process,path)._seq
