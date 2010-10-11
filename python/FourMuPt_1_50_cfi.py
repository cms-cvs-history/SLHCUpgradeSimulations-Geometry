import FWCore.ParameterSet.Config as cms

generator = cms.EDProducer("FlatRandomPtGunProducer",
    PGunParameters = cms.PSet(
        MaxPt = cms.double(50.0),
        MinPt = cms.double(0.9),
        PartID = cms.vint32(-13,-13),
        MaxEta = cms.double(2.5),
        MaxPhi = cms.double(3.14159265359),
        MinEta = cms.double(-2.5),
        MinPhi = cms.double(-3.14159265359) ## in radians

    ),
    Verbosity = cms.untracked.int32(0), ## set to 1 (or greater)  for printouts

    psethack = cms.string('Four mu pt 1 to 50'),
    AddAntiParticle = cms.bool(True),
    firstRun = cms.untracked.uint32(1)
)
