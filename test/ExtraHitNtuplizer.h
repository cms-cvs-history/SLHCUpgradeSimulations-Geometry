#ifndef ExtraHitNtuplizer_h
#define ExtraHitNtuplizer_h

/** \class ExtraHitNtuplizer
 * 
 *
 ************************************************************/

#include "DataFormats/TrackerRecHit2D/interface/SiTrackerGSRecHit2DCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/EDProduct.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/TrackReco/interface/Track.h"

//#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h" 
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"

#include "CondFormats/SiPixelObjects/interface/SiPixelLorentzAngle.h"
#include "SimDataFormats/SLHC/interface/StackedTrackerTypes.h"
#include <TH1F.h>
#include <TH2F.h>

class TTree;
class TFile;
class RectangularPixelTopology;

class TransientInitialStateEstimator;
class MagneticField;
class TrackerGeometry;
class TrajectoryStateOnSurface;
class PTrajectoryStateOnDet;

class ExtraHitNtuplizer : public edm::EDAnalyzer
//   Add to Buildfile and compile
/*   To run include the following into your cfg.py file, then add ReadLocalMeasurement2 to your final process
process.ReadLocalMeasurement2 = cms.EDAnalyzer("ExtraHitNtuplizer",
   src = cms.InputTag("siPixelRecHits"),
   #stereoRecHits = cms.InputTag("siStripMatchedRecHits","stereoRecHit"),
   #rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
   trackProducer = cms.InputTag("ctfWithMaterialTracks"),
   ### if using simple (non-iterative) or old (as in 1_8_4) tracking
   #trackProducer = cms.InputTag("ctfWithMaterialTracks"),
   OutputFile = cms.string("extrarechit_ntuple.root"),
   ### for using track hit association
   associatePixel = cms.bool(True),
   associateStrip = cms.bool(False),
   associateRecoTracks = cms.bool(False),
   ROUList = cms.vstring('famosSimHitsTrackerHits')
)

*/
{
 public:
  
  explicit ExtraHitNtuplizer(const edm::ParameterSet& conf);
  virtual ~ExtraHitNtuplizer();
  virtual void beginJob(const edm::EventSetup& es);
  virtual void endJob();
  virtual void analyze(const edm::Event& e, const edm::EventSetup& es);

 protected:

  void fillEvt(const edm::Event& );
  void fillSRecHit(const int subid, SiTrackerGSRecHit2DCollection::const_iterator pixeliter,
                   const GeomDet* theGeom);
  //void fillPRecHit(const int subid, SiPixelRecHitCollection::const_iterator pixeliter,
  //                 const GeomDet* PixGeom);
  void fillPRecHit(const int subid,
                   SiPixelRecHitCollection::const_iterator pixeliter,
                   const int num_simhit,
                   std::vector<PSimHit>::const_iterator closest_simhit,
                   const GeomDet* PixGeom,
                   edmNew::DetSet<SiPixelCluster>::const_iterator cluster,
                   float cluster_adc,
                   const int layers_struck,
                   const int nlayers_struck,
  		   int nstrk_lay,int nstrk_lad, int nstrk_mod);
  void fillPRecHit(const int subid, trackingRecHit_iterator pixeliter,
                   const GeomDet* PixGeom);
  void getthetaetaphi(
                   const double X, const double Y, const double Z,
                   double& theta, double& eta, double& phi);
  void getlorentz( const PixelGeomDetUnit* pixDet,
                   float& xdrift, float& ydrift);
  void fillSimHit(std::vector<PSimHit>::const_iterator simhit_i,
                                   int nstrk_lay,int nstrk_lad, int nstrk_mod);

  void fillStubSimHit(  cmsUpgrades::GlobalStub_PSimHit_Collection::const_iterator stubiter  );
  void fillStubDigiHit( cmsUpgrades::GlobalStub_PixelDigi_Collection::const_iterator stubiter);

 private:
  TH2F* globalStubPositions_XY;
  TH2F* globalStubPositions_RZ;
  edm::ParameterSet conf_;
  const TrackerGeometry*  theGeometry;
  const MagneticField* themag;
  mutable const SiPixelLorentzAngle * lorentzAngle_;
  edm::InputTag src_;
  edm::InputTag src_cluster;
  edm::InputTag src_digi;
  void init();
  
  typedef GloballyPositioned<double> Frame;
  //--- Structures for ntupling:
  struct evt
  {
    int run;
    int evtnum;
    
    void init();
  } evt_,stripevt_;

  struct StubSimHit
  {
    float gx;           // X position of Sim Hit Stub in Global Coord. (cm)
    float gy;           // Y position of Sim Hit Stub in Global Coord. (cm)
    float gz;		// Z position of Sim Hit Stub in Global Coord. (cm)
    float gr;		// R position of Sim Hit Stub in Global Coord. (cm)
 


    void init();
  } stubSimHit_;

  struct StubDigiHit
  {
    float gx;           // X position of Digi Hit Stub in Global Coord. (cm)
    float gy;           // Y position of Digi Hit Stub in Global Coord. (cm)
    float gz;           // Z position of Digi Hit Stub in Global Coord. (cm)
    float gr;           // R position of Digi Hit Stub in Global Coord. (cm)



    void init();
  } stubDigiHit_;


  struct SimHit
  { // No guarantee that it will work with disks and strips ...
    float x;      // X position of SimHit in Local Coord. (cm)
    float y;      // Y position of SimHit in Local Coord. (cm)
    float gx;     // X position of SimHit in Global Coord. (cm)
    float gy;     // Y position of SimHit in Global Coord. (cm)
    float gz;     // Z position of SimHit in Global Coord. (cm)
    float gr;	  // R position of SimHit in Global Coord. (cm)
    float theta;  // Theta of SimHit
    float eta;    // Eta of SimHit
    float phi;    // Phi of SimHit
    float spreadx;// Size of SimHit in X (units of cm)
    float spready;// Size of SimHit in Y (units of cm)

    int subid;    // 1 = barrel, 2 = disk
    int layer;    // barrel layer, starting from interaction point
    int ladder;   // ladder number
    int module;   // module number
    int nstrk_layer;    // Number of SimHit strikes in the same layer
    int nstrk_ladder;   // Number of SimHit strikes in the same ladder
    int nstrk_module;	// Number of SimHit strikes in the same module

    void init();
  } simHit_;
  
  struct RecHit 
  { // No guarantee that it will work with disks and strips ...
    float x;      // X position of RecHit in Local Coord. (cm)
    float y;      // Y position of RecHit in Local Coord. (cm)
    float xx;     // Err**2 in X position of RecHit (cm**2)
    float xy;     // Err corelation, not used for most layers
    float yy;     // Err**2 in Y position of RecHit (cm**2)
    float row;    // X position of RecHit in Local Coord. (units of PitchX)
    float col;    // Y position of RecHit in Local Coord. (units of PitchY)
    float gx;     // X position of RecHit in Global Coord. (cm)
    float gy;     // Y position of RecHit in Global Coord. (cm)
    float gz;     // Z position of RecHit in Global Coord. (cm)
    float gr;     // R position of RecHit in Global Coord. (cm)
    int subid;    // 1 = barrel, 2 = disk
    int layer;    // barrel layer, starting from interaction point
    int ladder;   // ladder
    int module;   // module
    int nstrk_layer;    // Number of RecHit strikes in the same layer
    int nstrk_ladder;   // Number of RecHit strikes in the same ladder
    int nstrk_module;   // Number of RecHit strikes in the same module
    int nsimhit;  // Not really a valid variable...
    float hx, hy; // X,Y position of SimHit in Local Coord. (cm)
    float tx, ty; // Not really a valid variable...
    float theta, phi, eta;  // theta, phi, eta of RecHit
    int num_pix;            // Number of pixels that went into constructing the RecHit
    int spreadx;            // Size of RecHit in X (Number of pixels)
    int spready;            // Size of RecHit in Y (Number of pixels)
    float hspreadx;         // Size of SimHit in X (units of cm)
    float hspready;         // Size of SimHit in Y (units of cm)
    float pitchx;           // Pixel Pitch in X (cm)
    float pitchy;           // Pixel Pitch in Y (cm)
    float thickness;        // thickness of the pixel
    float cluster_adc;      // Total adc count for cluster (units of electrons)
    int hnlayersstruck;     // number of simhits in the same layer 
    float hphi, heta;       // phi, eta of SimHit
    float hgx, hgy, hgz;    // x,y,z position of SimHit in Global Coord. (cm)

    void init();
  } recHit_, striprecHit_;

  TFile * tfile_;
  TTree * pixeltree_;
  TTree * striptree_;
  TTree * pixeltree2_;
  TTree * pixelsimhittree_;
  TTree * pixelstubsimtree_;
  TTree * pixelstubdigitree_;
};

#endif
