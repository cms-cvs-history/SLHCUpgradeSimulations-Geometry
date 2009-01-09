import FWCore.ParameterSet.Config as cms

process = cms.Process("GeometryTest")

# Number of events to be processed
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

process.load("FastSimulation.Configuration.CommonInputsFake_cff")
# replace with strawmanB geometry
process.load("SLHCUpgradeSimulations.Geometry.strawmanb_cmsIdealGeometryXML_cff")

process.source = cms.Source("EmptySource")

process.o1 = cms.OutputModule("AsciiOutputModule")
process.outpath = cms.EndPath(process.o1)

process.prod = cms.EDAnalyzer("ModuleInfo_StrawmanB",
   fromDDD = cms.bool(True)
)
process.p = cms.Path(process.prod)
