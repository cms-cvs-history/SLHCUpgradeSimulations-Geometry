
import FWCore.ParameterSet.Config as cms

def customise(process):
     process.load("SLHCUpgradeSimulations.Geometry.fakeConditions_stdgeom_cff")
     return (process)

