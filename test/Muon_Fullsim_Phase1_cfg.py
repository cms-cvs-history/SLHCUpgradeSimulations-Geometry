# Auto generated configuration file
# using: 
# Revision: 1.172.2.12 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: Configuration/Generator/python/SingleMuPt10_cfi.py -s GEN,SIM,DIGI,L1 -n 25 --conditions DESIGN_36_V10::All --datatier GEN-SIM --eventcontent FEVTDEBUG --beamspot Gauss --python_filename Muon_Fullsim_Phase1_cfg.py --slhc R39F16 --no_exec
# Then altered for four muons from 0.9 to 50 GeV
import FWCore.ParameterSet.Config as cms

process = cms.Process('L1')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.MixingNoPileUp_cff')
process.load('SLHCUpgradeSimulations.Geometry.PhaseI_cmsSimIdealGeometryXML_R39F16_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('Configuration.StandardSequences.VtxSmearedGauss_cff')
process.load('Configuration.StandardSequences.Sim_cff')
process.load('SLHCUpgradeSimulations.Geometry.Digi_Phase1_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('SLHCUpgradeSimulations.Geometry.fakeConditions_Phase1_cff')
process.load('SLHCUpgradeSimulations.Geometry.fakeConditions_Phase1_R39F16_cff')
process.load('Configuration.EventContent.EventContent_cff')

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.172.2.12 $'),
    annotation = cms.untracked.string('Configuration/Generator/python/SingleMuPt10_cfi.py nevts:25'),
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
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    fileName = cms.untracked.string('FourMuPt1_50_cfi_py_GEN_SIM_DIGI_L1.root'),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM'),
        filterName = cms.untracked.string('')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'DESIGN_36_V10::All'
process.generator = cms.EDProducer("FlatRandomPtGunProducer",
    PGunParameters = cms.PSet(
        MaxPt = cms.double(50.0),
        MinPt = cms.double(0.9),
        PartID = cms.vint32(-13,-13),
        MaxEta = cms.double(2.5),
        MaxPhi = cms.double(3.14159265359),
        MinEta = cms.double(-2.5),
        MinPhi = cms.double(-3.14159265359)
    ),
    Verbosity = cms.untracked.int32(0),
    psethack = cms.string('Four mu pt 1 to 50'),
    AddAntiParticle = cms.bool(True),
    firstRun = cms.untracked.uint32(1)
)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.endjob_step = cms.Path(process.endOfProcess)
process.out_step = cms.EndPath(process.output)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.simulation_step,process.digitisation_step,process.L1simulation_step,process.endjob_step,process.out_step)

# special treatment in case of production filter sequence
for path in process.paths: 
    getattr(process,path)._seq = process.generator*getattr(process,path)._seq
