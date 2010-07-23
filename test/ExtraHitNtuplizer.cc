// File: ExtraHitNtuplizer.cc
// Description: see ExtraHitNtuplizer.h
// Authors: Orig. H. Cheung as StdHitNtuplizer.cc in SLHC_2_2_3_005 --> Then heavily modified by E. Brownson
//--------------------------------------------------------------


#include "SLHCUpgradeSimulations/Geometry/test/ExtraHitNtuplizer.h"

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// DataFormats
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/OwnVector.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h" 
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h" 
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/Common/interface/Ref.h"

#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
// Geometry
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerTopology/interface/RectangularPixelTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelTopologyBuilder.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "SimDataFormats/SLHC/interface/StackedTrackerTypes.h"


// For ROOT
#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>
#include <TH2F.h>
#include <TH1F.h>

// STD
#include <memory>
#include <string>
#include <iostream>

using namespace std;
using namespace edm;
using namespace reco;

ExtraHitNtuplizer::ExtraHitNtuplizer(edm::ParameterSet const& conf) : 
  conf_(conf), 
  src_( conf.getParameter<edm::InputTag>( "src" ) ),
  src_cluster("siPixelClusters"),
  src_digi("simSiPixelDigis"),
  tfile_(0), 
  pixeltree_(0), 
  striptree_(0),
  pixeltree2_(0),
  pixelsimhittree_(0)
{
}


ExtraHitNtuplizer::~ExtraHitNtuplizer() { }  

void ExtraHitNtuplizer::endJob() 
{
  //globalStubPositions_XY->Print("test.eps");
  std::cout << " ExtraHitNtuplizer::endJob" << std::endl;
  tfile_->Write();
  tfile_->Close();
}



void ExtraHitNtuplizer::beginJob(const edm::EventSetup& es)
{
  std::cout << " ExtraHitNtuplizer::beginJob" << std::endl;
  std::string outputFile = conf_.getParameter<std::string>("OutputFile");
 
  tfile_ = new TFile ( outputFile.c_str() , "RECREATE" );
  pixeltree_ = new TTree("PixelNtuple","Pixel hit analyzer ntuple");
  striptree_ = new TTree("StripNtuple","Strip hit analyzer ntuple");
  pixeltree2_ = new TTree("Pixel2Ntuple","Track Pixel hit analyzer ntuple");
  pixelsimhittree_ = new TTree("PixelSimNtuple","Pixel simhit analyzer ntuple");
  pixelstubsimtree_= new TTree("PixelStubSimNtuple","Pixel stub simhit analyzer ntuple");
  pixelstubdigitree_= new TTree("PixelStubDigiNtuple","Pixel stub digihit analyzer ntuple");
  int bufsize = 64000;

  //Common Branch
  pixelstubsimtree_->Branch("evt",    &evt_,      "run/I:evtnum/I", bufsize);
  pixelstubsimtree_->Branch("pixel_stubSimHit", &stubSimHit_,
    "gx/F:gy:gz:gr", bufsize);

  pixelstubdigitree_->Branch("evt",    &evt_,      "run/I:evtnum/I", bufsize);
  pixelstubdigitree_->Branch("pixel_stubDigiHit", &stubDigiHit_,
    "gx/F:gy:gz:gr", bufsize);

  pixeltree_->Branch("evt",    &evt_,      "run/I:evtnum/I", bufsize);
  pixeltree_->Branch("pixel_recHit", &recHit_, 
    "x/F:y:xx:xy:yy:row:col:gx:gy:gz:gr:subid/I:layer:ladder:module:nstrk_layer:nstrk_ladder:nstrk_module:nsimhit:hx/F:hy:tx:ty:theta:phi:eta:num_pix/I:spreadx:spready:hspreadx/F:hspready:pitchx:pitchy:thickness:cluster_adc:hnlayersstruck/I:hphi/F:heta:hgx:hgy:hgz", bufsize);
  pixeltree2_->Branch("evt",    &evt_,      "run/I:evtnum/I", bufsize);
  pixeltree2_->Branch("pixel_recHit", &recHit_, 
    "x/F:y:xx:xy:yy:row:col:gx:gy:gz:gr:subid/I:layer:nsimhit:hx/F:hy:tx:ty:theta:phi", bufsize);
  
  // Strip Branches 
  striptree_->Branch("evt",    &evt_,      "run/I:evtnum/I", bufsize);
  striptree_->Branch("strip_recHit", &striprecHit_,
    "x/F:y:xx:xy:yy:row:col:gx:gy:gz:gr:subid/I:layer:nsimhit:hx/F:hy:tx:ty:theta:phi", bufsize);

  pixelsimhittree_->Branch("evt",    &evt_,      "run/I:evtnum/I", bufsize);
  pixelsimhittree_->Branch("pixel_simHit", &simHit_,  
    "x/F:y:gx:gy:gz:gr:theta:eta:phi:spreadx:spready:subid/I:layer:ladder:module:nstrk_layer:nstrk_ladder:nstrk_module", bufsize);

  // geometry setup
  edm::ESHandle<TrackerGeometry>        geometry;

  es.get<TrackerDigiGeometryRecord>().get(geometry);

  theGeometry = &(*geometry);

  edm::ESHandle<MagneticField> mag;
  es.get<IdealMagneticFieldRecord>().get(mag);
  themag = &(*mag);

}

// Functions that gets called by framework every event
void ExtraHitNtuplizer::analyze(const edm::Event& e, const edm::EventSetup& es)
{
 // Look at stubs stubs ... Probably should be moved later, but is fine here for now.
 // SimHit Stubs
 Handle< cmsUpgrades::GlobalStub_PSimHit_Collection > SimGlobalStubHandle;
 e.getByLabel( "GlobalStubsFromSimHits" , SimGlobalStubHandle);
 const cmsUpgrades::GlobalStub_PSimHit_Collection *tempSimGlobalStubs = SimGlobalStubHandle.product();
 cmsUpgrades::GlobalStub_PSimHit_Collection::const_iterator tempSimGlobalStubIter;
 for (  tempSimGlobalStubIter = tempSimGlobalStubs->begin(); tempSimGlobalStubIter != tempSimGlobalStubs->end() ; ++tempSimGlobalStubIter ) {
      //id = tempSimGlobalStubIter->Id();
      //globalStubPositions_XY->Fill(tempGlobalStubIter->position().x(), tempGlobalStubIter->position().y());
      //globalStubPositions_RZ->Fill(tempGlobalStubIter->position().z(), tempGlobalStubIter->position().perp());
      //std::cout <<"\nextrahit with global simstub at r,z = "<< tempSimGlobalStubIter->position().perp()<<" , "<<tempSimGlobalStubIter->position().z();
      fillStubSimHit(tempSimGlobalStubIter);
      fillEvt(e);
      pixelstubsimtree_->Fill();
      init();
     }
  // Digi Hit stubs
 Handle< cmsUpgrades::GlobalStub_PixelDigi_Collection > DigiGlobalStubHandle;
 e.getByLabel( "GlobalStubsFromPixelDigis" , DigiGlobalStubHandle);
 const cmsUpgrades::GlobalStub_PixelDigi_Collection *tempDigiGlobalStubs = DigiGlobalStubHandle.product();
 cmsUpgrades::GlobalStub_PixelDigi_Collection::const_iterator tempDigiGlobalStubIter;
 for (  tempDigiGlobalStubIter = tempDigiGlobalStubs->begin(); tempDigiGlobalStubIter != tempDigiGlobalStubs->end() ; ++tempDigiGlobalStubIter ) {
      //std::cout <<"\nextrahit with global digistub at r,z = "<< tempDigiGlobalStubIter->position().perp()<<" , "<<tempDigiGlobalStubIter->position().z();
      fillStubDigiHit(tempDigiGlobalStubIter);
      fillEvt(e);
      pixelstubdigitree_->Fill();
      init();
     }

// End of Stubs
  // fastsim rechits
  edm::Handle<SiPixelRecHitCollection> recHitColl;
  e.getByLabel( src_, recHitColl);

  // siPixelClusters before RecHits are made
  edm::Handle< edmNew::DetSetVector<SiPixelCluster> > clusterHitColl;
  e.getByLabel( src_cluster, clusterHitColl);

  // pixel digis
  edm::Handle< edm::DetSetVector<PixelDigi> >  digiHitColl;
  e.getByLabel( src_digi, digiHitColl);



  // for finding matched simhit
  TrackerHitAssociator associate( e, conf_ );

  //std::cout << " Step A: Standard RecHits found " << recHitColl->size() << std::endl;
  if(recHitColl->size() > 0) {
    //Loop over all rechits in SiPixelRecHitCollection (can also loop only over DetId)
    SiPixelRecHitCollection::const_iterator theRecHitRangeIteratorBegin = recHitColl->begin();
    SiPixelRecHitCollection::const_iterator theRecHitRangeIteratorEnd   = recHitColl->end();
    SiPixelRecHitCollection::const_iterator iterRecHit;

    std::string detname ;
    std::vector<PSimHit> matched;
    std::vector<PSimHit>::const_iterator closest_simhit;

    for ( iterRecHit = theRecHitRangeIteratorBegin; 
          iterRecHit != theRecHitRangeIteratorEnd; ++iterRecHit) {

      const DetId& detId =  iterRecHit->geographicalId();
      const GeomDet* geomDet( theGeometry->idToDet(detId) );

      unsigned int subdetId = detId.subdetId();
      int layerNumber=0;
      int ringNumber = 0;
      int stereo = 0;
      if ( subdetId == StripSubdetector::TIB) {
        detname = "TIB";
	TIBDetId tibid(detId.rawId());
	layerNumber = tibid.layer();
	stereo = tibid.stereo();
      } else if ( subdetId ==  StripSubdetector::TOB ) {
        detname = "TOB";
	TOBDetId tobid(detId.rawId());
	layerNumber = tobid.layer();
	stereo = tobid.stereo();
      } else if ( subdetId ==  StripSubdetector::TID) {
        detname = "TID";
	TIDDetId tidid(detId.rawId());
	layerNumber = tidid.wheel();
	ringNumber = tidid.ring();
	stereo = tidid.stereo();
      } else if ( subdetId ==  StripSubdetector::TEC ) {
        detname = "TEC";
	TECDetId tecid(detId.rawId());
	layerNumber = tecid.wheel();
	ringNumber = tecid.ring();
	stereo = tecid.stereo();
      } else if ( subdetId ==  PixelSubdetector::PixelBarrel ) {
        detname = "PXB";
	PXBDetId pxbid(detId.rawId());
	layerNumber = pxbid.layer();
	stereo = 1;
      } else if ( subdetId ==  PixelSubdetector::PixelEndcap ) {
        detname = "PXF";
	PXFDetId pxfid(detId.rawId());
	layerNumber = pxfid.disk();
	stereo = 1;
      }
      // *************************************************
      // get the digis and look at them...
      //   edm::Handle< edm::DetSetVector<PixelDigi> >  digiHitColl;
      //edm::DetSetVector<PixelDigi>::const_iterator digi_iter = digiHitColl.begin();
      //for( ; digi_iter != digiHitColl.end(); digi_iter++) {
      // DSViter->detId(); 
      // DSViter.detId();  

      //}

      // End of digis
      // *************************************************
      // get num_pix_hit that went into making the RecHit, use same style found in SiPixelRecHitConverter.cc
      // Note DetSetVector is a vector intrinsically sorted by detid...
      //int num_pix_hit =0, Spread_Y=0, Spread_X=0,
      int parity_checkA=0, looped_SiPixelCluster=0;
      float cluster_adc=0.0;
      edmNew::DetSet<SiPixelCluster>::const_iterator Cluster;
      const edmNew::DetSetVector<SiPixelCluster>& clusterHitColl_ = *clusterHitColl;
      edmNew::DetSetVector<SiPixelCluster>::const_iterator clusterHitColl_Iter=clusterHitColl_.begin();
      for ( ; clusterHitColl_Iter != clusterHitColl_.end() ; clusterHitColl_Iter++)
          { edmNew::DetSet<SiPixelCluster>::const_iterator clustIt = clusterHitColl_Iter->begin(), clustEnd = clusterHitColl_Iter->end();
            for ( ; clustIt != clustEnd; clustIt++)
            {
            // Will need to augment next line (and others also) for high pileup
            if (detId==clusterHitColl_Iter->detId())
               {
               Cluster = clustIt;parity_checkA=1;looped_SiPixelCluster++;
               const vector<SiPixelCluster::Pixel>& pixelsVec = clustIt->pixels();
               for (unsigned int i = 0;  i < pixelsVec.size(); ++i) 
                   {
 	           //float pixx = pixelsVec[i].x; // index as float=i+0.5
	           //float pixy = pixelsVec[i].y; // same
	           float adc = float(pixelsVec[i].adc); // In units of electrons
                   cluster_adc += adc ;
                   }
               }
            }
          }
      if (looped_SiPixelCluster>=2)  {std::cout <<"\nERROR in ExtraHitNtuplizer, SiPixelCluster had multiple clusters on a detId."
                                                <<"  This gives an internal logic problem\n";}
      if (parity_checkA==0) {std::cout <<"\nERROR in ExtraHitNtuplizer, SiPixelCluster not found will give error in final ntuple.\n";}
      // get matched simhit
      // *************************************************
        Handle< std::vector<PSimHit> > pixelSimHitHandle;  
        e.getByLabel("famosSimHits", "TrackerHits", pixelSimHitHandle);
        const std::vector<PSimHit>* pixelSimHits=pixelSimHitHandle.product();
      int sim_stacks_struck=0,num_stacks_struck=0;
      for (std::vector<PSimHit>::const_iterator simHitIter = pixelSimHits->begin(); simHitIter<pixelSimHits->end(); simHitIter++)
      {
          if (subdetId==PXBDetId(simHitIter->detUnitId()).subdetId()&&layerNumber==PXBDetId(simHitIter->detUnitId()).layer()) {
             num_stacks_struck++;
          }
      }
      // End of check for simhit in both stack layers
      matched.clear();
      matched = associate.associateHit(*iterRecHit);
      if ( !matched.empty() ) {
        float closest = 9999.9;
        std::vector<PSimHit>::const_iterator closestit = matched.begin();
        LocalPoint lp = iterRecHit->localPosition();
        float rechit_x = lp.x();
        float rechit_y = lp.y();
        //loop over simhits and find closest
        for (std::vector<PSimHit>::const_iterator m = matched.begin(); m<matched.end(); m++) 
        {
          // End stack hit confirmation.
          float sim_x1 = (*m).entryPoint().x();
          float sim_x2 = (*m).exitPoint().x();
          float sim_xpos = 0.5*(sim_x1+sim_x2);
          float sim_y1 = (*m).entryPoint().y();
          float sim_y2 = (*m).exitPoint().y();
          float sim_ypos = 0.5*(sim_y1+sim_y2);
            
          float x_res = fabs(sim_xpos - rechit_x);
          float y_res = fabs(sim_ypos - rechit_y);
          float dist = sqrt(x_res*x_res + y_res*y_res);
          if ( dist < closest ) {
                closest = dist;
                closestit = m;
          }
        } // end of simhit loop
        closest_simhit = closestit;
      } // end matched emtpy
      unsigned int subid = detId.subdetId();
      int layer_num = 0,ladder=0,module=0;
      if ( (subid==1)||(subid==2) ) {
        // 1 = PXB, 2 = PXF
        if ( subid ==  PixelSubdetector::PixelBarrel ) {
	  PXBDetId pxbid(detId.rawId());
	  layer_num   = pxbid.layer();
	  ladder      = pxbid.ladder();
	  module      = pxbid.module();
        } else if ( subid ==  PixelSubdetector::PixelEndcap ) {
	  PXFDetId pxfid(detId.rawId());
	  layer_num   = pxfid.disk();
        }
      int num_simhit = matched.size();
      // Fill variables pertaining to the number of RecHits in same lay,lad,mod
      int nstrk_lay = -1 , nstrk_lad= -1 , nstrk_mod = -1 ; // Stat from -1 since you should always find 1 match (itself)
      edm::Handle<SiPixelRecHitCollection> recHitColl_b;
      e.getByLabel( src_, recHitColl_b);
      for (SiPixelRecHitCollection::const_iterator iterRecHit_b = recHitColl_b->begin();iterRecHit_b!=recHitColl_b->end();iterRecHit_b++){
 	  const DetId& detId_b = iterRecHit_b->geographicalId();
	  //unsigned int subid_b = detId_b.subdetId();
	  int lay_b    = PXBDetId(detId_b.rawId()).layer(); 
	  int ladder_b = PXBDetId(detId_b.rawId()).ladder();
	  int module_b = PXBDetId(detId_b.rawId()).module();
	  if (detId.subdetId()==detId_b.subdetId()&&layer_num==lay_b){
	    nstrk_lay++;
	    if (ladder==ladder_b){
	      nstrk_lad++;
	      if (module==module_b){
		nstrk_mod++;
	      } // Same Module
	    } // Same Ladder
	  } // Same Layer
      } // End look over RecHit_b
      //
      // *****************************************************************************************
      // Fill the variables that will go into the ntuple 
      fillPRecHit(subid, iterRecHit, num_simhit, closest_simhit, geomDet,Cluster,cluster_adc,sim_stacks_struck,num_stacks_struck, nstrk_lay, nstrk_lad, nstrk_mod);
      // *****************************************************************************************
      fillEvt(e);
      pixeltree_->Fill();
      init();
      } // endif (subid==1)||(subid==2)
    } // end of rechit loop
  } // end of loop test on recHitColl size

// Now loop over SimHits and fill pixelsimhittree_
{
      Handle< std::vector<PSimHit> > pixelSimHitHandle;
      e.getByLabel("famosSimHits", "TrackerHits", pixelSimHitHandle);
      const std::vector<PSimHit>* pixelSimHits=pixelSimHitHandle.product();
            int sim_stacks_struck=0, nstrk_layer = -1, nstrk_ladder = -1, nstrk_module = -1; // Stat from -1 since you should always find 1 match (itself)
      for (std::vector<PSimHit>::const_iterator simHitIter = pixelSimHits->begin(); simHitIter<pixelSimHits->end(); simHitIter++)
      {      sim_stacks_struck =0; nstrk_layer = -1; nstrk_ladder = -1; nstrk_module = -1; // Stat from -1 since you should always find 1 match (itself)
             DetId detId(simHitIter->detUnitId());
 	     unsigned int subdetId = detId.subdetId();
             int lay    = (PXBDetId(simHitIter->detUnitId()).layer() ) ;
             int ladder = (PXBDetId(simHitIter->detUnitId()).ladder()) ;
	     int module = (PXBDetId(simHitIter->detUnitId()).module()) ;
             Handle< std::vector<PSimHit> > pixelSimHitHandle_b;
             e.getByLabel("famosSimHits", "TrackerHits", pixelSimHitHandle_b);
             const std::vector<PSimHit>* pixelSimHits_b=pixelSimHitHandle_b.product();
             for (std::vector<PSimHit>::const_iterator simHitIter_b = pixelSimHits_b->begin(); simHitIter_b<pixelSimHits_b->end(); simHitIter_b++)
             { // we need a sub loop to look at the other SimHits in the same layer...
		DetId detId_b(simHitIter_b->detUnitId());
		unsigned int subdetId_b = detId_b.subdetId();
		if (subdetId==subdetId_b){
		   int lay_b    = (PXBDetId(simHitIter_b->detUnitId()).layer() ) ;
		   int ladder_b = (PXBDetId(simHitIter_b->detUnitId()).ladder()) ;
		   int module_b = (PXBDetId(simHitIter_b->detUnitId()).module()) ;
		   if (lay==lay_b){
		      nstrk_layer++;
		      if (ladder==ladder_b){
			 nstrk_ladder++;
			 if (module==module_b){
			    nstrk_module++;
			 } // Same Module
		      } // Same Ladder
		   } // Same Layer
		} // Same DetID
             }
        // fill the ntuple
        fillSimHit(simHitIter,nstrk_layer,nstrk_ladder,nstrk_module);
        fillEvt(e);
        pixelsimhittree_->Fill();
        init();
        // end fill
      }
}
// Now loop over recotracks

  edm::Handle<View<reco::Track> >  trackCollection;
  edm::InputTag trackProducer;
  trackProducer = conf_.getParameter<edm::InputTag>("trackProducer");
  e.getByLabel(trackProducer, trackCollection);


  std::cout << " num of reco::Tracks with "
            << trackProducer.process()<<":"
            << trackProducer.label()<<":"
            << trackProducer.instance()
            << ": " << trackCollection->size() << "\n";

  int rT = 0;
  for(View<reco::Track>::size_type i=0; i<trackCollection->size(); ++i){
      ++rT;
      RefToBase<reco::Track> track(trackCollection, i);
//      std::cout << " num of hits for track " << rT << " = " << track->recHitsSize() << std::endl;
      for(trackingRecHit_iterator ih=track->recHitsBegin(); ih != track->recHitsEnd(); ++ih) {
        TrackingRecHit * hit = (*ih)->clone();
        const DetId& detId =  hit->geographicalId();
        const GeomDet* geomDet( theGeometry->idToDet(detId) );

        unsigned int subdetId = detId.subdetId();
        int layerNumber=0;
        int ringNumber = 0;
        int stereo = 0;
        std::string detname;
        if ( subdetId == StripSubdetector::TIB) {
          detname = "TIB";
          TIBDetId tibid(detId.rawId());
          layerNumber = tibid.layer();
          stereo = tibid.stereo();
        } else if ( subdetId ==  StripSubdetector::TOB ) {
          detname = "TOB";
	  TOBDetId tobid(detId.rawId());
          layerNumber = tobid.layer();
          stereo = tobid.stereo();
        } else if ( subdetId ==  StripSubdetector::TID) {
          detname = "TID";
          TIDDetId tidid(detId.rawId());
          layerNumber = tidid.wheel();
          ringNumber = tidid.ring();
          stereo = tidid.stereo();
        } else if ( subdetId ==  StripSubdetector::TEC ) {
          detname = "TEC";
          TECDetId tecid(detId.rawId());
          layerNumber = tecid.wheel();
          ringNumber = tecid.ring();
          stereo = tecid.stereo();
        } else if ( subdetId ==  PixelSubdetector::PixelBarrel ) {
          detname = "PXB";
          PXBDetId pxbid(detId.rawId());
          layerNumber = pxbid.layer();
          stereo = 1;
        } else if ( subdetId ==  PixelSubdetector::PixelEndcap ) {
          detname = "PXF";
          PXFDetId pxfid(detId.rawId());
          layerNumber = pxfid.disk();
          stereo = 1;
        }
//        std::cout << "RecHit in " << detname << " from detid " << detId.rawId()
//                  << " subdet = " << subdetId
//                  << " layer = " << layerNumber
//                  << " Stereo = " << stereo
//                  << std::endl;
        if(hit->isValid()) {
          unsigned int subid = detId.subdetId();
          if ( (subid==1)||(subid==2) ) {
            // 1 = PXB, 2 = PXF
            fillPRecHit(subid, ih, geomDet);
            fillEvt(e);
            pixeltree2_->Fill();
            init();
          }
        }
      } //end of loop on tracking rechits
  } // end of loop on recotracks

            
} // end analyze function

void ExtraHitNtuplizer::fillStubSimHit( cmsUpgrades::GlobalStub_PSimHit_Collection::const_iterator stubiter)
{
  stubSimHit_.gx = stubiter->position().x()             ;
  stubSimHit_.gy = stubiter->position().y()             ;
  stubSimHit_.gz = stubiter->position().z()  		;
  stubSimHit_.gr = stubiter->position().perp()  	;
}

void ExtraHitNtuplizer::fillStubDigiHit( cmsUpgrades::GlobalStub_PixelDigi_Collection::const_iterator stubiter)
{
  stubDigiHit_.gx = stubiter->position().x()             ;
  stubDigiHit_.gy = stubiter->position().y()             ;
  stubDigiHit_.gz = stubiter->position().z()             ;
  stubDigiHit_.gr = stubiter->position().perp()          ;
}


void ExtraHitNtuplizer::fillSimHit(std::vector<PSimHit>::const_iterator simhit_i,
                                   int nstrk_lay,int nstrk_lad,int nstrk_mod)
{          
             DetId detId(simhit_i->detUnitId());
             const GeomDet* geomDet( theGeometry->idToDet(detId) );

    float sim_x1 = (*simhit_i).entryPoint().x();
    float sim_x2 = (*simhit_i).exitPoint().x();
    simHit_.x = 0.5*(sim_x1+sim_x2);
    float sim_y1 = (*simhit_i).entryPoint().y();
    float sim_y2 = (*simhit_i).exitPoint().y();
    simHit_.y = 0.5*(sim_y1+sim_y2);
    simHit_.spreadx = (sim_x2-sim_x1);
    simHit_.spready = (sim_y2-sim_y1);
    LocalPoint lp (simHit_.x,simHit_.y);
    GlobalPoint GP = geomDet->surface().toGlobal(lp);
    simHit_.gx = GP.x();
    simHit_.gy = GP.y();
    simHit_.gz = GP.z();
    simHit_.gr = sqrt((GP.x()*GP.x())+(GP.y()*GP.y()));
    simHit_.subid  = detId.subdetId();
    simHit_.layer  = PXBDetId(simhit_i->detUnitId()).layer();
    simHit_.ladder = PXBDetId(simhit_i->detUnitId()).ladder();
    simHit_.module = PXBDetId(simhit_i->detUnitId()).module();
    simHit_.nstrk_layer  = nstrk_lay;
    simHit_.nstrk_ladder = nstrk_lad;
    simHit_.nstrk_module = nstrk_mod;
       double theta=21.0, eta=22.0, phi=23.0;
       ExtraHitNtuplizer::getthetaetaphi(double(GP.x()), double(GP.y()), double(GP.z()), theta, eta, phi);
    simHit_.theta =  float(theta);
    simHit_.eta = float(eta);
    simHit_.phi = float(phi);
//    std::cout	<<"\nExtrahit SimHit: SubId   = "		<< simHit_.subid
//             	<<"\n\tlayer, ladder , module = "		<< simHit_.layer	<<"\t, "<< simHit_.ladder	<<"\t, "<< simHit_.module
//		<<"\n\tWith N hits in the same Lay,Lad,mod "	<< simHit_.nstrk_layer	<<"\t, "<< simHit_.nstrk_ladder	<<"\t, "<< simHit_.nstrk_module
//		<<"\n\twith a Local X,Y       = "		<< simHit_.x		<<"\t, "<< simHit_.y
//		<<"\n\tLocated at X,Y,Z       = "		<< simHit_.gx		<<"\t, "<< simHit_.gy		<<"\t, "<< simHit_.gz
//		<<"\n\tLocated at R,Phi,Eta   = "		<< simHit_.gr		<<"\t, "<< simHit_.phi		<<"\t , "<< simHit_.eta
//		;
}


void ExtraHitNtuplizer::fillSRecHit(const int subid, 
                                   SiTrackerGSRecHit2DCollection::const_iterator pixeliter,
                                   const GeomDet* theGeom)
{
  LocalPoint lp = pixeliter->localPosition();
  LocalError le = pixeliter->localPositionError();

  striprecHit_.x = lp.x();
  striprecHit_.y = lp.y();
  striprecHit_.xx = le.xx();
  striprecHit_.xy = le.xy();
  striprecHit_.yy = le.yy();
  //MeasurementPoint mp = topol->measurementPosition(LocalPoint(striprecHit_.x, striprecHit_.y));
  //striprecHit_.row = mp.x();
  //striprecHit_.col = mp.y();
  GlobalPoint GP = theGeom->surface().toGlobal(pixeliter->localPosition());
  striprecHit_.gx = GP.x();
  striprecHit_.gy = GP.y();
  striprecHit_.gz = GP.z();
  striprecHit_.gr = sqrt((GP.x()*GP.x())+(GP.y()*GP.y()));
  striprecHit_.subid = subid;
}
void ExtraHitNtuplizer::fillPRecHit(const int subid, 
                                  SiPixelRecHitCollection::const_iterator pixeliter,
                                  const int num_simhit,
                                  std::vector<PSimHit>::const_iterator closest_simhit,
                                  const GeomDet* PixGeom,
                                  edmNew::DetSet<SiPixelCluster>::const_iterator cluster,
                                  float cluster_adc,
                                  const int layers_struck,
                                  const int nlayers_struck,
				  int nstrk_lay,int nstrk_lad,int nstrk_mod
                                  )
{
  LocalPoint lp = pixeliter->localPosition();
  LocalError le = pixeliter->localPositionError();

  recHit_.x = lp.x();
  recHit_.y = lp.y();
  recHit_.xx = le.xx();
  recHit_.xy = le.xy();
  recHit_.yy = le.yy();
  //MeasurementPoint mp = topol->measurementPosition(LocalPoint(recHit_.x, recHit_.y));
  //recHit_.row = mp.x();
  //recHit_.col = mp.y();
  GlobalPoint GP = PixGeom->surface().toGlobal(pixeliter->localPosition());
  recHit_.gx = GP.x();
  recHit_.gy = GP.y();
  recHit_.gz = GP.z();
  recHit_.gr = sqrt((GP.x()*GP.x())+(GP.y()*GP.y()));
  recHit_.subid = subid;
  recHit_.nsimhit = num_simhit;
  recHit_.cluster_adc = cluster_adc;
  recHit_.nstrk_layer  = nstrk_lay;
  recHit_.nstrk_ladder = nstrk_lad;
  recHit_.nstrk_module = nstrk_mod;
  double theta=21.0, eta=22.0, phi=23.0;
  ExtraHitNtuplizer::getthetaetaphi(double(GP.x()), double(GP.y()), double(GP.z()), theta, eta, phi);
  recHit_.theta = float(theta);
  recHit_.phi = float(phi);
  recHit_.eta = float(eta);
  recHit_.hnlayersstruck = nlayers_struck;
  recHit_.num_pix = cluster->size();
  recHit_.spreadx = cluster->sizeX();
  recHit_.spready = cluster->sizeY();
  // Get the physical attributes of the pixels
    float thePitchX=20.0,thePitchY=21.0,theThickness=22.0;
    const DetId& detid = pixeliter->geographicalId();
    DetId detIdObject( detid );
    const GeomDetUnit * genericDet =  theGeometry->idToDetUnit( detIdObject );
    const PixelGeomDetUnit * pixDet = dynamic_cast<const PixelGeomDetUnit*>(genericDet);
    const RectangularPixelTopology * theTopol = dynamic_cast<const RectangularPixelTopology*>( & (pixDet->specificTopology()) );
    std::pair<float,float> pitchxy = theTopol->pitch();
    thePitchX = pitchxy.first;
    thePitchY = pitchxy.second;
    theThickness = pixDet->surface().bounds().thickness();
    float xdrift,ydrift;
    getlorentz(pixDet,xdrift,ydrift);
  // Fill more things
  recHit_.layer  = PXBDetId(detid).layer() ;
  recHit_.ladder = PXBDetId(detid).ladder();
  recHit_.module = PXBDetId(detid).module();
  recHit_.pitchx = thePitchX;  
  recHit_.pitchy = thePitchY;
  recHit_.thickness = theThickness;
  if(num_simhit > 0) {
    float sim_x1 = (*closest_simhit).entryPoint().x();
    float sim_x2 = (*closest_simhit).exitPoint().x();
    recHit_.hx = 0.5*(sim_x1+sim_x2);
    float sim_y1 = (*closest_simhit).entryPoint().y();
    float sim_y2 = (*closest_simhit).exitPoint().y();
    recHit_.hy = 0.5*(sim_y1+sim_y2);
    recHit_.hspreadx = (sim_x2-sim_x1); // To be physically comparable to previous variables this needs a Lorentz Drift correction.
    recHit_.hspready = (sim_y2-sim_y1);
  // Fill some information specific to the SimHit (should be close to RecHit)
  float lpx=0.5*(sim_x1+sim_x2),lpy=0.5*(sim_y1+sim_y2),lpz=0.0;
  LocalPoint lp(lpx,lpy,lpz);
  GlobalPoint HGP = PixGeom->surface().toGlobal(lp);
  double htheta=21.0, heta=22.0, hphi=23.0;
  ExtraHitNtuplizer::getthetaetaphi(double(HGP.x()), double(HGP.y()), double(HGP.z()), htheta, heta, hphi);
  recHit_.hphi = float(hphi);
  recHit_.heta = float(heta);
  recHit_.hgx  = float(HGP.x());
  recHit_.hgy  = float(HGP.y());
  recHit_.hgz  = float(HGP.z());


    }
}
void ExtraHitNtuplizer::fillPRecHit(const int subid, 
                                  trackingRecHit_iterator ih,
                                  const GeomDet* PixGeom)
{
  TrackingRecHit * pixeliter = (*ih)->clone();
  LocalPoint lp = pixeliter->localPosition();
  LocalError le = pixeliter->localPositionError();

  recHit_.x = lp.x();
  recHit_.y = lp.y();
  recHit_.xx = le.xx();
  recHit_.xy = le.xy();
  recHit_.yy = le.yy();
  GlobalPoint GP = PixGeom->surface().toGlobal(pixeliter->localPosition());
  recHit_.gx = GP.x();
  recHit_.gy = GP.y();
  recHit_.gz = GP.z();
  recHit_.gr = sqrt((GP.x()*GP.x())+(GP.y()*GP.y()));
  recHit_.subid = subid;
}

void
ExtraHitNtuplizer::fillEvt(const edm::Event& E)
{
   evt_.run = E.id().run();
   evt_.evtnum = E.id().event();
}

void ExtraHitNtuplizer::init()
{
  evt_.init();
  recHit_.init();
  simHit_.init();
  striprecHit_.init();
}

void ExtraHitNtuplizer::evt::init()
{
  int dummy_int = 9999;
  run = dummy_int;
  evtnum = dummy_int;
}
void ExtraHitNtuplizer::getthetaetaphi(
                                  const double X, const double Y, const double Z,
                                  double& theta, double& eta, double& phi)
{
  double pi=3.1415927;
  double R=sqrt(X*X+Y*Y+Z*Z);
  if (R!=0) {
     theta=acos(Z/R);
     }
  eta = -log( tan( theta/2.0 ) );
  phi = atan2(Y,X);
  if(phi<=0) {
    phi = phi+(2.0*pi);
    }
  return;
}
void ExtraHitNtuplizer::getlorentz(const PixelGeomDetUnit* pixDet,
                               float& xdrift, float& ydrift)
// Possibly just use CondFormats/SiPixelObjects/interface/SiPixelLorentzAngle.h, look into later
{  // This function does not work yet. Basic code taken out of PixelCPEBase.cc
   // No idea when I will fix it.    It was just an FYI thing here.
   xdrift=9999;ydrift=9999;
   float theThickness = pixDet->surface().bounds().thickness();
   // from PixelCPEBase::driftDirection( GlobalVector bfield )
   GlobalVector bfield = themag->inTesla(pixDet->surface().position()) ;  // This is the input for PixelCPEBase::driftDirection
      Frame detFrame(pixDet->surface().position(), pixDet->surface().rotation());
      LocalVector Bfield = detFrame.toLocal(bfield);
/*      double langle = lorentzAngle_->getLorentzAngle(uint32_t(pixDet->geographicalId().rawId()));  
      
      float alpha2;
      alpha2 = langle*langle; // There is an bool alpha2Order (extra E*B term) option ignore for now...
        float dir_x =  ( langle * Bfield.y() + alpha2* Bfield.z()* Bfield.x() );
        float dir_y = -( langle * Bfield.x() - alpha2* Bfield.z()* Bfield.y() );
        float dir_z = -( 1 + alpha2* Bfield.z()*Bfield.z() );
        float scale = (1 + alpha2* Bfield.z()*Bfield.z() );
        LocalVector theDriftDirection = LocalVector(dir_x/scale, dir_y/scale, dir_z/scale );
        // theDriftDirection is the output of PixelCPEBase::driftDirection
   // Keep it in the same overall form as PixelCPEBase::lorentzShiftX(), hence why I am not making it simplier...
   xdrift = theDriftDirection.x()/theDriftDirection.z() * theThickness; // The max X shift in cm
        // traditionaly you have lxshift = xdrift / thePitchX / 2.;
   ydrift = theDriftDirection.y()/theDriftDirection.z() * theThickness; // The max Y shift in cm
        // traditionaly you have lyshift = ydrift / thePitchY / 2.;
  // float lxshift = xdrift / thePitchX / 2.;
  // float lyshift = ydrift / thePitchY / 2.;
*/
  return;
}

void ExtraHitNtuplizer::StubSimHit::init()
{
  float dummy_float = 9999.0;
  gx = dummy_float;
  gy = dummy_float;
  gz = dummy_float;
  gr = dummy_float;

}

void ExtraHitNtuplizer::StubDigiHit::init()
{
  float dummy_float = 9999.0;
  gx = dummy_float;
  gy = dummy_float;
  gz = dummy_float;
  gr = dummy_float;

}


void ExtraHitNtuplizer::SimHit::init()
{
  float dummy_float = 9999.0;
  x = dummy_float;
  y = dummy_float;
  gx = dummy_float;
  gy = dummy_float;
  gz = dummy_float;
  gr = dummy_float;
  subid  = 0;
  layer  = 0;
  ladder = 0;
  module = 0;
  nstrk_layer  = 0;
  nstrk_ladder = 0;
  nstrk_module = 0;
  theta = dummy_float;
  eta = dummy_float;
  phi = dummy_float;
  spreadx = dummy_float;
  spready = dummy_float;
}

void ExtraHitNtuplizer::RecHit::init()
{
  float dummy_float = 9999.0;

  x = dummy_float;
  y = dummy_float;
  xx = dummy_float;
  xy = dummy_float; 
  yy = dummy_float;
  row = dummy_float;
  col = dummy_float;
  gx = dummy_float;
  gy = dummy_float;
  gz = dummy_float;
  gr = dummy_float;
  layer = 0;
  ladder = 0;
  module = 0;
  nsimhit = 0;
  nstrk_layer  = 0;
  nstrk_ladder = 0;
  nstrk_module = 0;
  hx = dummy_float;
  hy = dummy_float;
  hgx = dummy_float;
  hgy = dummy_float;
  hgz = dummy_float;
  tx = dummy_float;
  ty = dummy_float;
  theta = dummy_float;
  phi = dummy_float;
  eta = dummy_float; 
  hphi = dummy_float;
  heta = dummy_float;
  pitchx = dummy_float;
  pitchy = dummy_float;
  num_pix = 0;
  spreadx = 0;
  spready = 0;
  hspreadx = dummy_float;
  hspready = dummy_float;
  thickness = dummy_float;
  cluster_adc = dummy_float;
  hnlayersstruck = 0;
}

//define this as a plug-in
DEFINE_ANOTHER_FWK_MODULE(ExtraHitNtuplizer);

