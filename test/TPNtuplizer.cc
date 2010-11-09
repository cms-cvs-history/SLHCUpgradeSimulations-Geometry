// File: TPNtuplizer.cc
// Description: see TPNtuplizer.h
// Authors: H. Cheung
//--------------------------------------------------------------


#include "SLHCUpgradeSimulations/Geometry/test/TPNtuplizer.h"

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// DataFormats
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/OwnVector.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"

#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h" 

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByChi2.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

// Geometry
#include "MagneticField/Engine/interface/MagneticField.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerTopology/interface/RectangularPixelTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelTopologyBuilder.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"

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

TPNtuplizer::TPNtuplizer(edm::ParameterSet const& conf) : 
  conf_(conf), 
  label_(conf.getParameter< std::vector<edm::InputTag> >("label")),
  label_tp_effic_(conf.getParameter< edm::InputTag >("label_tp_effic")),
  label_tp_fake_(conf.getParameter< edm::InputTag >("label_tp_fake")),
  UseAssociators_(conf.getParameter< bool >("UseAssociators")),
  associators_(conf.getParameter< std::vector<std::string> >("associators")),
  tfile_(0), 
  tptree_(0)
{
  tpSelector_ = TrackingParticleSelector(conf_.getParameter<double>("ptMinTP"),
                                         conf_.getParameter<double>("minRapidityTP"),
                                         conf_.getParameter<double>("maxRapidityTP"),
                                         conf_.getParameter<double>("tipTP"),
                                         conf_.getParameter<double>("lipTP"),
                                         conf_.getParameter<int>("minHitTP"),
                                         conf_.getParameter<bool>("signalOnlyTP"),
                                         conf_.getParameter<bool>("chargedOnlyTP"),
                                         conf_.getParameter<bool>("stableOnlyTP"),
                                         conf_.getParameter<std::vector<int> >("pdgIdTP"));
}

TPNtuplizer::~TPNtuplizer() { }  

void TPNtuplizer::endJob() 
{
  std::cout << " TPNtuplizer::endJob" << std::endl;
  tfile_->Write();
  tfile_->Close();
}

void TPNtuplizer::beginJob(const edm::EventSetup& es)
{
  std::cout << " TPNtuplizer::beginJob" << std::endl;
  std::string outputFile = conf_.getParameter<std::string>("OutputFile");
 
  tfile_ = new TFile ( outputFile.c_str() , "RECREATE" );
  tptree_ = new TTree("TPNtuple","Tracking Particle analyzer ntuple");
  shtree_ = new TTree("SHNtuple","Simtrack analyzer ntuple");

  int bufsize = 64000;

  //Common Branch
  tptree_->Branch("evt", &evt_, "run/I:evtnum:numtp:nseltp:nfdtp:numtk:nasstk", bufsize);
  tptree_->Branch("tpart", &tp_, 
    "tpn/I:bcross:tevt:charge:stable:status:pdgid:mathit:signal:llived:sel:gpsz:gpstat:npix:nbpix:nfpix:ndiff:nbdiff:nfdiff:pt/F:eta:tip:lip:p:e:phi:theta:rap:qual", bufsize);
  shtree_->Branch("simtrk_simhit", &simtrk_simhit_,
    "subid/I:layer:gx/F:gy:gz", bufsize);

  
  // in beginRun in MTV
  if (UseAssociators_) {
    edm::ESHandle<TrackAssociatorBase> theAssociator;
    for (unsigned int w=0;w<associators_.size();w++) {
      es.get<TrackAssociatorRecord>().get(associators_[w],theAssociator);
      associator_.push_back( theAssociator.product() );
    }
  }
}

// Functions that gets called by framework every event
void TPNtuplizer::analyze(const edm::Event& event, const edm::EventSetup& es)
{
  using namespace reco;

  // geometry setup
  edm::ESHandle<TrackerGeometry>        geometry;

  es.get<TrackerDigiGeometryRecord>().get(geometry);
  const TrackerGeometry*  theGeometry = &(*geometry);

  edm::Handle<TrackingParticleCollection>  TPCollectionHeff ;
  event.getByLabel(label_tp_effic_,TPCollectionHeff);
  const TrackingParticleCollection tPCeff = *(TPCollectionHeff.product());
  
  edm::Handle<TrackingParticleCollection>  TPCollectionHfake ;
  event.getByLabel(label_tp_fake_,TPCollectionHfake);
  const TrackingParticleCollection tPCfake = *(TPCollectionHfake.product());

  // Loop over all TP and get the hit information
  //std::cout << "number of TP in event = " << tPCeff.size() << std::endl;
  for (TrackingParticleCollection::size_type i=0; i<tPCeff.size(); i++){
    TrackingParticleRef tpr(TPCollectionHeff, i);
    TrackingParticle* tp=const_cast<TrackingParticle*>(tpr.get());
    if( (! tpSelector_(*tp))) continue;
    //int num_simhit = tp->matchedHit();
    //std::cout << "TP = " << i << " number of simhit = " << num_simhit << std::endl;
    std::vector<PSimHit> simhits = tp->trackPSimHit();
    //std::cout << "TP = " << i << " size of simhit vector = " << simhits.size() << std::endl;
    //std::cout << "TP = " << i << " end-begin vector = " << tp->pSimHit_end() - tp->pSimHit_begin()<< std::endl;
    std::string detname ;
    int layer=0;
    for (unsigned int j=0; j<simhits.size(); j++){
      DetId theDetUnitId(simhits[j].detUnitId());
      int detector = theDetUnitId.det();
      int subdetector = theDetUnitId.subdetId();
      if(detector != DetId::Tracker) continue;
      if ( subdetector ==  PixelSubdetector::PixelBarrel ) {
        detname = "PXB";
        PXBDetId pxbid(theDetUnitId.rawId());
        layer = pxbid.layer();
      } else if ( subdetector ==  PixelSubdetector::PixelEndcap ) {
        detname = "PXF";
        PXFDetId pxfid(theDetUnitId.rawId());
        layer = pxfid.disk();
      } else if ( subdetector == StripSubdetector::TIB) {
        detname = "TIB";
        TIBDetId tibid(theDetUnitId.rawId());
        layer = tibid.layer();
      } else if ( subdetector ==  StripSubdetector::TOB ) {
        detname = "TOB";
        TOBDetId tobid(theDetUnitId.rawId());
        layer = tobid.layer();
      } else if ( subdetector ==  StripSubdetector::TID) {
        detname = "TID";
        TIDDetId tidid(theDetUnitId.rawId());
        layer = tidid.wheel();
      } else if ( subdetector ==  StripSubdetector::TEC ) {
        detname = "TEC";
        TECDetId tecid(theDetUnitId.rawId());
        layer = tecid.wheel();
      }
      //std::cout << "       simhit " << j << " " << detname << " layer/disk = " << layer << std::endl;
      const GeomDet* geomDet( theGeometry->idToDet(theDetUnitId) );
      float simhitx = 0.5 * ( simhits[j].entryPoint().x() + simhits[j].exitPoint().x() );
      float simhity = 0.5 * ( simhits[j].entryPoint().y() + simhits[j].exitPoint().y() );
      LocalPoint lp(simhitx,simhity,0.);
      GlobalPoint GP = geomDet->surface().toGlobal(lp);
      float hitgx = GP.x();
      float hitgy = GP.y();
      float hitgz = GP.z();
      //std::cout << "       x " << hitgx << " y = " << hitgy << " z = " << hitgz << " r = " 
      //          << sqrt(hitgx*hitgx+hitgy*hitgy)<< std::endl;
      fill_simtrk_simhit(subdetector, layer, hitgx, hitgy, hitgz);
      shtree_->Fill();
      simtrk_simhit_.init();
    }
  }

  for (unsigned int ww=0;ww<associators_.size();ww++){
    // get some numbers for this event - very inefficient!
    int num_sel=0;
    for (TrackingParticleCollection::size_type i=0; i<tPCeff.size(); i++){
      TrackingParticleRef tpr(TPCollectionHeff, i);
      TrackingParticle* tp=const_cast<TrackingParticle*>(tpr.get());
      if( (! tpSelector_(*tp))) continue;
      ++num_sel;
    }
    edm::LogVerbatim("TPNtuplizer") << "\n# Tracking particles selected = " << num_sel << "\n";
    for (unsigned int www=0;www<label_.size();www++){
      // get track collection from event for specified collection label(s)
      edm::Handle<View<Track> >  trackCollection;
      event.getByLabel(label_[www], trackCollection);
      edm::LogVerbatim("TPNtuplizer") << "\n# of Reco tracks collection = " << trackCollection->size() << "\n";
      // do the association
      reco::RecoToSimCollection recSimColl;
      reco::SimToRecoCollection simRecColl;
      // only handle doing association in job at the mo
      if(UseAssociators_){
         recSimColl=associator_[ww]->associateRecoToSim(trackCollection, TPCollectionHfake, &event);
         simRecColl=associator_[ww]->associateSimToReco(trackCollection, TPCollectionHeff, &event);
      }
      // get number for this event and this track collection - very inefficient!
      int num_found=0;
      for (TrackingParticleCollection::size_type i=0; i<tPCeff.size(); i++){
        TrackingParticleRef tpr(TPCollectionHeff, i);
        TrackingParticle* tp=const_cast<TrackingParticle*>(tpr.get());
        if( (! tpSelector_(*tp))) continue;
        if(simRecColl.find(tpr) != simRecColl.end()) ++num_found;
      }
      int num_ass=0;
      for(View<Track>::size_type i=0; i<trackCollection->size(); ++i){
        RefToBase<Track> track(trackCollection, i);
        std::vector<std::pair<TrackingParticleRef, double> > tp;
        if(recSimColl.find(track) != recSimColl.end()){
          tp = recSimColl[track];
          if (tp.size()!=0) ++num_ass;
        }
      } 
      edm::LogVerbatim("TPNtuplizer") << "\n# Tracking particles selected and found = " << num_found << "\n";
      edm::LogVerbatim("TPNtuplizer") << "\n# Reco tracks associated = " << num_ass << "\n";

      edm::LogVerbatim("TPNtuplizer") << "\n# of TrackingParticles: " << tPCeff.size() << "\n";
      edm::LogVerbatim("TPNtuplizer") << "\n# of Reco tracks for ntuple: " << trackCollection->size() << "\n";
      int ats = 0;
      int st=0;
      for (TrackingParticleCollection::size_type i=0; i<tPCeff.size(); i++){
        TrackingParticleRef tpr(TPCollectionHeff, i);
        TrackingParticle* tp=const_cast<TrackingParticle*>(tpr.get());
        //if( (! tpSelector_(*tp))) continue;
        int selected = 0;
        if( tpSelector_(*tp)) selected = 1;
        st++;

        std::vector<std::pair<RefToBase<Track>, double> > rt;
        float quality = 0.0;
        int matched_hit = 0;
        if(simRecColl.find(tpr) != simRecColl.end()){
          rt = (std::vector<std::pair<RefToBase<Track>, double> >) simRecColl[tpr];
          if (rt.size()!=0) {
            ats++;
            edm::LogVerbatim("TPNtuplizer") << "TrackingParticle #" << i << " selected #" << st 
                                            << " with pt=" << sqrt(tp->momentum().perp2()) 
                                            << " associated with quality:" << rt.begin()->second <<"\n";
            quality = rt.begin()->second;
            matched_hit = 1;
          }
        }else{
          edm::LogVerbatim("TPNtuplizer") << "TrackingParticle #" << i << " selected #" << st
                                          << " with pt=" << sqrt(tp->momentum().perp2())
                                          << " NOT associated to any reco::Track" << "\n";
        }
        fillEvt(tPCeff.size(), num_sel, num_found, trackCollection->size(), num_ass, event);
        fillTP(st, matched_hit, quality, selected, tp);
        tptree_->Fill();
        init();
      } // end loop over tracking particles

      // next reconstructed tracks
      edm::LogVerbatim("TPNtuplizer") << "\n# of reco::Tracks = " << trackCollection->size() << "\n";
      int at=0;
      int rT=0;
      for(View<Track>::size_type i=0; i<trackCollection->size(); ++i){
        RefToBase<Track> track(trackCollection, i);
        rT++;

        std::vector<std::pair<TrackingParticleRef, double> > tp;
        if(recSimColl.find(track) != recSimColl.end()){
          tp = recSimColl[track];
          if (tp.size()!=0) {
            at++;
            edm::LogVerbatim("TPNtuplizer") << "reco::Track #" << rT << " with pt=" << track->pt() 
                                            << " associated with quality:" << tp.begin()->second <<"\n";
          }
        } else {
          edm::LogVerbatim("TPNtuplizer") << "reco::Track #" << rT << " with pt=" << track->pt()
                                          << " NOT associated to any TrackingParticle" << "\n";		  
        }
      } // end loop over reco tracks

    } // end loop on track collection label
  } // end of loop on associators

} // end analyze function

void TPNtuplizer::fillTP(const int num, const int matched_hit, const float quality, 
                         const int selected, const TrackingParticle* tp)
{
  edm::LogVerbatim("TPNtuplizer") << "Filling TP with pt= " << sqrt(tp->momentum().perp2()) << "\n";

  tp_.tpn = num;
  tp_.bcross = tp->eventId().bunchCrossing();
  tp_.tevt = tp->eventId().event();
  tp_.charge = tp->charge();
  int stable = 1;
  for( TrackingParticle::genp_iterator j = tp->genParticle_begin(); j != tp->genParticle_end(); ++j ) {
    const HepMC::GenParticle * p = j->get();
    if (p->status() != 1) {
      stable = 0; break;
    }
  }
  tp_.stable = stable;
  tp_.status = tp->status();
  tp_.pdgid = tp->pdgId();
  tp_.mathit = matched_hit;
  if (tp->eventId().bunchCrossing()== 0 && tp->eventId().event() == 0) tp_.signal = 1;
    else tp_.signal = 0;
  if(tp->longLived()) tp_.llived = 1;
    else tp_.llived = 0;
  tp_.sel = selected;
  int numgp = 0;
  for( TrackingParticle::genp_iterator j = tp->genParticle_begin(); j != tp->genParticle_end(); ++ j ) ++numgp;
  int gpstatus = -69;
  for( TrackingParticle::genp_iterator j = tp->genParticle_begin(); j != tp->genParticle_end(); ++ j ) {
    const HepMC::GenParticle * p = j->get();
    if (p->status() != 1) {
      gpstatus = p->status(); break;
    }
    gpstatus = p->status();
  }
  tp_.gpsz = numgp;
  tp_.gpstat = gpstatus;
  tp_.pt = sqrt(tp->momentum().perp2());
  tp_.eta = tp->eta();
  tp_.tip = sqrt(tp->vertex().perp2());
  tp_.lip = fabs(tp->vertex().z());
  tp_.p = tp->p();
  tp_.e = tp->energy();
  tp_.phi = tp->phi();
  tp_.theta = tp->theta();
  tp_.rap = tp->rapidity();
  tp_.qual = quality;
  // count the number of pixel hits for seeding study
  std::vector<PSimHit> simhits = tp->trackPSimHit();
  unsigned int layer=0;
  int npixhits = 0;
  int nbpixhits = 0;
  int nfpixhits = 0;
  int ndiffpixhits = 0;
  int ndiffbpixhits = 0;
  int ndifffpixhits = 0;
  // choose values to work for phase 1 also
  std::vector<bool> bpix_layerhit(4,false);
  std::vector<bool> fpix_diskhit(3,false);
  for (unsigned int j=0; j<simhits.size(); j++){
    DetId theDetUnitId(simhits[j].detUnitId());
    int detector = theDetUnitId.det();
    int subdetector = theDetUnitId.subdetId();
    if(detector != DetId::Tracker) continue;
    if ( subdetector ==  PixelSubdetector::PixelBarrel ) {
      nbpixhits++;
      PXBDetId pxbid(theDetUnitId.rawId());
      layer = pxbid.layer();
      if(layer > 0 && layer <=4) bpix_layerhit[layer-1]=true;
    } else if ( subdetector ==  PixelSubdetector::PixelEndcap ) {
      nfpixhits++;
      PXFDetId pxfid(theDetUnitId.rawId());
      layer = pxfid.disk();
      if(layer > 0 && layer <=3) fpix_diskhit[layer-1]=true;
    }
  }
  npixhits = nbpixhits + nfpixhits;
  for (unsigned int j=0; j<bpix_layerhit.size(); j++)
     if(bpix_layerhit[j]) ++ndiffbpixhits;
  for (unsigned int j=0; j<fpix_diskhit.size(); j++)
     if(fpix_diskhit[j]) ++ndifffpixhits;
  ndiffpixhits = ndiffbpixhits + ndifffpixhits;
  tp_.npix = npixhits;
  tp_.nbpix = nbpixhits;
  tp_.nfpix = nfpixhits;
  tp_.ndiff = ndiffpixhits;
  tp_.nbdiff = ndiffbpixhits;
  tp_.nfdiff = ndifffpixhits;
}

void TPNtuplizer::fillEvt(const int numtp, const int nseltp, const int nfdtp,
                          const int numtk, const int nasstk, const edm::Event& E)
{
   evt_.run = E.id().run();
   evt_.evtnum = E.id().event();
   evt_.numtp = numtp;
   evt_.nseltp = nseltp;
   evt_.nfdtp = nfdtp;
   evt_.numtk = numtk;
   evt_.nasstk = nasstk;
}

void TPNtuplizer::fill_simtrk_simhit(const int subid, const int layer, const float gx, 
              const float gy, const float gz)
{
  simtrk_simhit_.subid = subid;
  simtrk_simhit_.layer = layer;
  simtrk_simhit_.gx = gx;
  simtrk_simhit_.gy = gy;
  simtrk_simhit_.gz = gz;
}


void TPNtuplizer::init()
{
  evt_.init();
  tp_.init();
}

void TPNtuplizer::evt::init()
{
  int dummy_int = 9999;
  run = dummy_int;
  evtnum = dummy_int;
  numtp = dummy_int;
  numtk = dummy_int;
}

void TPNtuplizer::myTp::init()
{
  int dummy_int = 9999;
  float dummy_float = 9999.0;
  pt  = dummy_float;
  tpn = dummy_int;
  bcross = dummy_int;
  tevt = dummy_int;
  charge = dummy_int;
  stable = dummy_int;
  status = dummy_int;
  pdgid = dummy_int;
  mathit = dummy_int;
  signal = dummy_int;
  llived = dummy_int;
  sel = dummy_int;
  gpsz = dummy_int;
  gpstat = dummy_int;
  npix = dummy_int;
  nbpix = dummy_int;
  nfpix = dummy_int;
  ndiff = dummy_int;
  nbdiff = dummy_int;
  nfdiff = dummy_int;
  pt = dummy_float;
  eta = dummy_float;
  tip = dummy_float;
  lip = dummy_float;
  p = dummy_float;
  e = dummy_float;
  phi = dummy_float;
  theta = dummy_float;
  rap = dummy_float;
  qual = dummy_float;
}

void TPNtuplizer::mysimhit::init()
{
  int dummy_int = 9999;
  float dummy_float = 9999.0;
  subid = dummy_int;
  layer = dummy_int;
  gx = dummy_float;
  gy = dummy_float;
  gz = dummy_float;
}


//define this as a plug-in
DEFINE_ANOTHER_FWK_MODULE(TPNtuplizer);

