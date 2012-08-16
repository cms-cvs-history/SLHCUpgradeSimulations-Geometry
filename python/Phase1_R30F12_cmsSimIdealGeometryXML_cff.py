import FWCore.ParameterSet.Config as cms

#  Tracking Geometry
from SLHCUpgradeSimulations.Geometry.Phase1_cmsSimIdealGeometryXML_cfi import *
from Geometry.CommonDetUnit.globalTrackingGeometry_cfi import *
#from Geometry.CommonDetUnit.globalTrackingGeometryDB_cfi import *

#hardwire these here
TrackerGeometricDetESModule.pixelGeometryConstants.layerNumberPXB=cms.uint32(18)
TrackerGeometricDetESModule.pixelGeometryConstants.totalBlade=cms.uint32(56)


#Tracker
from RecoTracker.GeometryESProducer.TrackerRecoGeometryESProducer_cfi import *
from Geometry.TrackerNumberingBuilder.trackerNumberingGeometry_cfi import *

#Muon
from RecoMuon.DetLayers.muonDetLayerGeometry_cfi import *
from Geometry.MuonNumbering.muonNumberingInitialization_cfi import *

#  Calorimeters
from Geometry.CaloEventSetup.CaloTopology_cfi import *
from Geometry.CaloEventSetup.AlignedCaloGeometryDBReader_cfi import *
from Geometry.CaloEventSetup.EcalTrigTowerConstituents_cfi import *
from Geometry.EcalMapping.EcalMapping_cfi import *
from Geometry.EcalMapping.EcalMappingRecord_cfi import *

#  Alignment
from Geometry.TrackerGeometryBuilder.idealForDigiTrackerGeometry_cff import *
#from Geometry.TrackerGeometryBuilder.idealForDigiTrackerGeometryDB_cff import *
from Geometry.CSCGeometryBuilder.idealForDigiCscGeometryDB_cff import *
from Geometry.DTGeometryBuilder.idealForDigiDtGeometryDB_cff import *

