import FWCore.ParameterSet.Config as cms

# define all the changes for unganging ME1a
def me1a(process):
    #from Configuration.StandardSequences.GeometryDB_cff import *
    process.CSCGeometryESModule.useGangedStripsInME1a = False
    process.idealForDigiCSCGeometry.useGangedStripsInME1a = False
    #done
    return process
