import FWCore.ParameterSet.Config as cms

process = cms.Process("ICALIB")
process.load("Configuration.StandardSequences.FakeConditions_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("SLHCUpgradeSimulations.Geometry.strawmana_cmsIdealGeometryXML_cff")
##process.load("Geometry.TrackerNumberingBuilder.trackerNumberingGeometry_cfi")

##process.load("Geometry.TrackerRecoData.trackerRecoGeometryXML_cfi")

##from Geometry.TrackerGeometryBuilder.trackerGeometry_cfi import *
##applyAlignment = cms.bool(False)


##process.load("Geometry.TrackerGeometryBuilder.trackerGeometry_cfi")
##process.load("Geometry.TrackerGeometryBuilder.idealForDigiTrackerGeometry_cff")


#from Alignment.CommonAlignmentProducer.FakeAlignmentSource_cfi import *

#process.fakeGlobalPositionSource = cms.ESSource("fakeGlobalPositionSource",
#    recordName = cms.string('GlobalPositionRcd'),
#    firstValid = cms.vuint32(1),
#    iovIsRunNotTime = cms.bool(True)
#)

process.source = cms.Source("EmptyIOVSource",
    firstValue = cms.uint64(1),
    lastValue = cms.uint64(1),
    timetype = cms.string('runnumber'),
    interval = cms.uint64(1)
)

process.MessageLogger = cms.Service("MessageLogger",
    insert_logfile = cms.untracked.PSet(
        threshold = cms.untracked.string('INFO')
    ),
    destinations = cms.untracked.vstring('./logfile.txt')
)

process.Timing = cms.Service("Timing")

process.prod = cms.EDFilter("SiStripDetInfoFileWriter",
    FilePathStrip = cms.untracked.string('myfile.txt'),
    FilePathPixel = cms.untracked.string('PixelSkimmedGeometry_strawmanA.txt')
)

process.asciiPrint = cms.OutputModule("AsciiOutputModule")

process.p = cms.Path(process.prod)
process.ep = cms.EndPath(process.asciiPrint)


