// -*- C++ -*-
//
// Package:    DatasetValidation/DatasetValidationTool
// Class:      DatasetValidationTool
//
/**\class DatasetValidationTool DatasetValidationTool.cc DatasetValidation/DatasetValidationTool/plugins/DatasetValidationTool.cc

 Description: Tool to validate datasets used for CMS Tracker Alignment by producing standard distributions

 Implementation:
    Creates Histograms with track variables required to validate either of the datasets- Cosmics, CDCs, Z->mumu, Upsilon->mumu,  

//
// Original Author:  Saumya Saumya
//         Created:  Tue, 06 Oct 2020 13:05:38 GMT
//
//

*/
// system include files
#include <memory>
/*
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"

#include "DataFormats/SiPixelDetId/interface/PixelBarrelName.h"
#include "DataFormats/SiPixelDetId/interface/PixelEndcapName.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"

#include "TMath.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TH2D.h"
#include "TProfile.h"
*/

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "CondFormats/DataRecord/interface/RunSummaryRcd.h"
#include "CondFormats/DataRecord/interface/SiStripCondDataRecords.h"
#include "CondFormats/RunInfo/interface/RunInfo.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/Phi.h"
#include "DataFormats/GeometryVector/interface/Theta.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Provenance/interface/RunLumiEventNumber.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackResiduals.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/TrackerRecHit2D/interface/ProjectedSiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit1D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include <TH1.h>
#include <TMath.h>
#include <TLorentzVector.h>

//
// class declaration
//

class DatasetValidationTool: public edm::one::EDAnalyzer<edm::one::WatchRuns, edm::one::WatchLuminosityBlocks, edm::one::SharedResources> 
{
   public:
      explicit DatasetValidationTool(const edm::ParameterSet&);
      ~DatasetValidationTool();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&) override; 
      virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
 
      // ----------member data ---------------------------
      edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
      edm::EDGetTokenT<reco::BeamSpot> beamspotToken_;
      edm::EDGetTokenT<reco::VertexCollection> vertexToken_;

      edm::Service<TFileService> fs;

      const bool isHit2D(const TrackingRecHit &hit);

      bool isResonance;

      int nTracks, nEvents, nTracksInEvent, nEventsInRun, nTracksInRun, nTracksInLuminosity, nEventsInLuminosity, nInvalidHits, nNotDetUnits,nZeroSubDetId;
      int nHits_2D, nHits_PIXEL, nHitsFPIX, nHits_FPIXplus, nHits_FPIXminus, nHits_TIDplus,nHits_TIDminus, nHits_TECplus, nHits_TECminus, nHits_ENDCAP, nHits_ENDCAPplus, nHits_ENDCAPminus;
      float etaMax_=3.0,M_PI_=3.14159;

      TH1D *h_magField; 

      //-----Kinematic Variables-----// 
      TH1D *h_charge;
      TH1D *h_p;
      TH1D *h_pt;
      TH1D *h_eta;
      TH1D *h_phi;
      TH1D *h_theta;
      TH1D *h_chi2;
      TH1D *h_chi2_ndf;
      TH1D *h_chi2_Prob;
      TH1D *h_d0;
      TH1D *h_dz;
      TH1D *h_dxy;
      TH1D *h_d0PV;
      TH1D *h_dzPV;
      TH1D *h_dxyPV;
      TH1D *h_d0BS;
      TH1D *h_dzBS;
      TH1D *h_dxyBS;
      
      //-----Hits-----//
      TH1I *h_nh_Valid;
      TH1I *h_nh_2D;
      TH1I *h_nh_PIXEL;
      TH1I *h_nh_BPIX;
      TH1I *h_nh_FPIX;         
      TH1I *h_nh_FPIXplus;
      TH1I *h_nh_FPIXminus;
      TH1I *h_nh_TIB;
      TH1I *h_nh_TID;
      TH1I *h_nh_TIDplus;
      TH1I *h_nh_TIDminus;
      TH1I *h_nh_TOB;
      TH1I *h_nh_TEC;
      TH1I *h_nh_TECplus;
      TH1I *h_nh_TECminus;
      TH1I *h_nh_ENDCAP;
      TH1I *h_nh_ENDCAPplus;
      TH1I *h_nh_ENDCAPminus;

      //-----Residuals-----//
      TH1D* h_Res_BPIX_xPrime; 
      TH1D* h_Res_FPIX_xPrime;
      TH1D* h_Res_TIB_xPrime;
      TH1D* h_Res_TID_xPrime;
      TH1D* h_Res_TOB_xPrime;
      TH1D* h_Res_TEC_xPrime;
      TH1D* h_Res_BPIX_yPrime;
      TH1D* h_Res_FPIX_yPrime;
      TH1D* h_Res_TIB_yPrime;
      TH1D* h_Res_TID_yPrime;
      TH1D* h_Res_TOB_yPrime;
      TH1D* h_Res_TEC_yPrime;

      //----- Resonances------//
      TH1D* h_Resonance_InvMass;
      TH1D* h_Resonance_Eta;
      TH1D* h_Resonance_Pt;
      TH1D* h_Resonance_Phi;

      //-----Multiplicity-----//
      TH1I* h_Tracks_In_Event;
      TH1I* h_Events_In_Run;
      TH1I* h_Tracks_In_Run;
      TH1I* h_Events_In_Luminosity;
      TH1I* h_Tracks_In_Luminosity;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
DatasetValidationTool::DatasetValidationTool(const edm::ParameterSet& iConfig)
 :
  tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"))),
  beamspotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("BS"))),
  vertexToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices")))
{
   //now do what ever initialization is needed
   usesResource("TFileService");
//   edm::Service<TFileService> fs;
   isResonance = iConfig.getParameter<bool>("IsResonance");

   nTracks=0; nEvents=0, nInvalidHits=0, nNotDetUnits=0, nZeroSubDetId=0;
}


DatasetValidationTool::~DatasetValidationTool()
{
}

//
// member functions
//

// ------------ method called for each event  ------------
void
DatasetValidationTool::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   using namespace edm;

   edm::Handle<reco::TrackCollection> tracksHandle_;
   iEvent.getByToken(tracksToken_,tracksHandle_);
   const reco::TrackCollection tC = *(tracksHandle_.product()); 

   // Geometry setup
   edm::ESHandle<TrackerGeometry> geometry;
   iSetup.get<TrackerDigiGeometryRecord>().get(geometry);
   const TrackerGeometry *theGeometry = &(*geometry);

   //Retrieve tracker topology from geometry
   edm::ESHandle<TrackerTopology> tTopoHandle;
   iSetup.get<TrackerTopologyRcd>().get(tTopoHandle);
   const TrackerTopology *const tTopo = tTopoHandle.product();

   //Beamspot Handle for dxy(BS)
   reco::BeamSpot beamSpot;
   edm::Handle<reco::BeamSpot> beamSpotHandle_;
   iEvent.getByToken(beamspotToken_, beamSpotHandle_);

   //Vertex Handle for dxy(PV)
   edm::Handle<reco::VertexCollection> vertexHandle_;
   iEvent.getByToken(vertexToken_,vertexHandle_);

   for(reco::TrackCollection::const_iterator track=tC.begin(); track!=tC.end(); track++)
   {
     nTracks++; nTracksInEvent++; nTracksInRun++; nTracksInLuminosity++;

     h_charge->Fill(track->charge());
     h_p->Fill(track->p());
     h_pt->Fill(track->pt());
     h_eta->Fill(track->eta());
     h_phi->Fill(track->phi());
     h_theta->Fill(track->theta());
     h_chi2->Fill(track->chi2());
     h_chi2_ndf->Fill(track->normalizedChi2());
     h_chi2_Prob->Fill(TMath::Prob(track->chi2(),track->ndof()));
     h_d0->Fill(track->d0());
     h_dz->Fill(track->dz());
     h_dxy->Fill(track->dxy());

     if (beamSpotHandle_.isValid())
     {
        beamSpot = *beamSpotHandle_;
        math::XYZPoint BSpoint(beamSpot.x0(), beamSpot.y0(), beamSpot.z0());
        double dxy = track->dxy(BSpoint);
        double dz = track->dz(BSpoint);
        h_dxyBS->Fill(dxy);
        h_d0BS->Fill(-dxy);
        h_dzBS->Fill(dz);
     }

     if(vertexHandle_.isValid())
     {
        reco::Vertex vtx;
        double min_dxy=100., dz=100.;
        for(auto vtx=vertexHandle_->cbegin();vtx!=vertexHandle_->cend();++vtx)
        {
            math::XYZPoint Vpoint(vtx->x(), vtx->y(), vtx->z());
            if(abs(min_dxy) > abs(track->dxy(Vpoint)))
            {
               min_dxy = track->dxy(Vpoint);
               dz = track->dz(Vpoint);
            }
        }
        h_dxyPV->Fill(min_dxy);
        h_d0PV->Fill(-min_dxy);
        h_dzPV->Fill(dz);
     }

     //---------------------------Resonances-----------------------------//
     if(isResonance)
     {
         double InvMass=0;
         for (reco::TrackCollection::const_iterator trk2=track+1; trk2!=tC.end(); trk2++)
         {   
             if(track->charge()*trk2->charge()>0) continue;
             TLorentzVector track1,track2,mother;
             track1.SetPtEtaPhiE(track->pt(),track->eta(),track->phi(),sqrt((track->p()*track->p())+(0.105*0.105)));
             track2.SetPtEtaPhiE(trk2->pt(),trk2->eta(),trk2->phi(),sqrt((trk2->p()*trk2->p())+(0.105*0.105)));
             mother=track1+track2; InvMass=mother.M();
             if(InvMass<0. || InvMass>5.) continue;
             h_Resonance_InvMass->Fill(InvMass);
             h_Resonance_Eta->Fill(mother.Eta());
             h_Resonance_Pt->Fill(mother.Pt());     
             h_Resonance_Phi->Fill(mother.Phi());        
         }
     }

     //-------------------- Hits and Residuals --------------------------//
     h_nh_Valid->Fill(track->numberOfValidHits());
     h_nh_BPIX->Fill(track->hitPattern().numberOfValidPixelBarrelHits());
     h_nh_FPIX->Fill(track->hitPattern().numberOfValidPixelEndcapHits());
     h_nh_TIB->Fill(track->hitPattern().numberOfValidStripTIBHits());
     h_nh_TOB->Fill(track->hitPattern().numberOfValidStripTOBHits());
     h_nh_TID->Fill(track->hitPattern().numberOfValidStripTIDHits());
     h_nh_TEC->Fill(track->hitPattern().numberOfValidStripTECHits());

     nHits_2D=0;
     nHits_PIXEL=0;
     nHits_FPIXplus=0;
     nHits_FPIXminus=0;
     nHits_TIDplus=0;
     nHits_TIDminus=0;
     nHits_TECplus=0;
     nHits_TECminus=0;
     nHits_ENDCAP=0;
     nHits_ENDCAPplus=0;
     nHits_ENDCAPminus=0;

     auto const &trajParams = track->extra()->trajParams();
     auto const &residuals = track->extra()->residuals();
     assert(trajParams.size() == track->recHitsSize());

     int h_index=0;
     for(auto iHit = track->recHitsBegin(); iHit!=track->recHitsEnd(); ++iHit,++h_index)
     {
         const DetId detId=(*iHit)->geographicalId();
         const int SubDetId = detId.subdetId();
         const GeomDet *geomDet(theGeometry->idToDet(detId));

         if (!(*iHit)->isValid()) { nInvalidHits++; continue; } 
         if (SubDetId == 0) { nZeroSubDetId++; continue;} 
         if (!(*iHit)->detUnit())  {  nNotDetUnits++; continue;  }// is it a single physical module?

         if (isHit2D(**iHit)) { ++nHits_2D; }  // 2D Hit

         float uOrientation(-999.F), vOrientation(-999.F);
        
         double resX = residuals.residualX(h_index); 
         double resY = residuals.residualY(h_index);

         // All the transformations here      
         LocalPoint lPModule(0., 0., 0.), lUDirection(1., 0., 0.), lVDirection(0., 1., 0.);
         GlobalPoint gUDirection = geomDet->surface().toGlobal(lUDirection);
         GlobalPoint gVDirection = geomDet->surface().toGlobal(lVDirection);
         GlobalPoint gPModule = geomDet->surface().toGlobal(lPModule);       

         //-------------------------------------- Filling Hits and Residuals -------------------------------------------//
         if (SubDetId == PixelSubdetector::PixelBarrel || SubDetId == StripSubdetector::TIB || SubDetId == StripSubdetector::TOB) 
         {
            uOrientation = deltaPhi(gUDirection.barePhi(), gPModule.barePhi()) >= 0. ? +1.F : -1.F;
            vOrientation = gVDirection.z() - gPModule.z() >= 0 ? +1.F : -1.F;
            switch(SubDetId)
            {
               case PixelSubdetector::PixelBarrel:
               h_Res_BPIX_xPrime->Fill(uOrientation*resX*10); //Units (mm)
               h_Res_BPIX_yPrime->Fill(vOrientation*resY*10);
               nHits_PIXEL++; 
               break;

               case StripSubdetector::TIB:
               h_Res_TIB_xPrime->Fill(uOrientation*resX*10);
               h_Res_TIB_yPrime->Fill(vOrientation*resY*10);
               break;

               case StripSubdetector::TOB:
               h_Res_TOB_xPrime->Fill(uOrientation*resX*10);
               h_Res_TOB_yPrime->Fill(vOrientation*resY*10);

               break;
            }
        }
        else if ( SubDetId == PixelSubdetector::PixelEndcap || SubDetId == StripSubdetector::TID || SubDetId == StripSubdetector::TEC)     
        {
           uOrientation = deltaPhi(gUDirection.barePhi(), gPModule.barePhi()) >= 0. ? +1.F : -1.F;
           vOrientation = gVDirection.perp() - gPModule.perp() >= 0. ? +1.F : -1.F;
           switch(SubDetId)
           {
              case PixelSubdetector::PixelEndcap:
              h_Res_FPIX_xPrime->Fill(uOrientation*resX*10);
              h_Res_FPIX_yPrime->Fill(vOrientation*resY*10);
              nHits_PIXEL++; nHits_ENDCAP++;
              if(tTopo->pxfSide(detId)==1) {  nHits_FPIXminus++; nHits_ENDCAPminus++;}
              else {  nHits_FPIXplus++; nHits_ENDCAPplus++; }
              break;

              case StripSubdetector::TID:
              h_Res_TID_xPrime->Fill(uOrientation*resX*10);
              h_Res_TID_yPrime->Fill(vOrientation*resY*10);
              nHits_ENDCAP++;
              if(tTopo->tidIsZMinusSide(detId)) {  nHits_TIDminus++; nHits_ENDCAPminus++; }
              else {  nHits_TIDplus++; nHits_ENDCAPplus++; }
              break;

              case StripSubdetector::TEC:
              h_Res_TEC_xPrime->Fill(uOrientation*resX*10);
              h_Res_TEC_yPrime->Fill(vOrientation*resY*10);
              nHits_ENDCAP++;
              if(tTopo->tecIsZMinusSide(detId)) {  nHits_TECminus++; nHits_ENDCAPminus++; }
              else {  nHits_TECplus++; nHits_ENDCAPplus++; }
              break;
           }
        } 

     }  //Hits Loop

     h_nh_2D->Fill(nHits_2D);
     h_nh_PIXEL->Fill(nHits_PIXEL);
     h_nh_FPIXplus->Fill(nHits_FPIXplus);
     h_nh_FPIXminus->Fill(nHits_FPIXminus);
     h_nh_TIDplus->Fill(nHits_TIDplus);
     h_nh_TIDminus->Fill(nHits_TIDminus);
     h_nh_TECplus->Fill(nHits_TECplus);
     h_nh_TECminus->Fill(nHits_TECminus);
     h_nh_ENDCAP->Fill(nHits_ENDCAP);
     h_nh_ENDCAPplus->Fill(nHits_ENDCAPplus);
     h_nh_ENDCAPminus->Fill(nHits_ENDCAPminus);

   } //Tracks Loop

   nEvents++;  nEventsInRun++;  nEventsInLuminosity++;
   h_Tracks_In_Event->Fill(nTracksInEvent);

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void
DatasetValidationTool::beginJob()
{
    
   /* TFileDirectory TrackDir = fs->mkdir("Track");
    TFileDirectory TrackHitsDir = fs->mkdir("TrackHits");
    TFileDirectory TrackResidualsDir = fs->mkdir("TrackResiduals");
    TFileDirectory MultiplicityDir = fs->mkdir("Multiplicity");
*/
    TH1D::SetDefaultSumw2(kTRUE);
    TH1I::SetDefaultSumw2(kTRUE);
   
    //==============Track Directory=================//
    h_charge = fs->make<TH1D>("h_charge", "Charge;Track Charge (e);Tracks (#)", 5, -2.5, 2.5);
    h_p = fs->make<TH1D>("h_p", "Momentum;Track p (GeV);Tracks (#)", 200, 0.,200.);
    h_pt = fs->make<TH1D>("h_pt", "Transverse Momentum;Track p_{T} (GeV);Tracks (#)", 100, 0., 100.);
    h_eta = fs->make<TH1D>("h_eta", "#eta Distribution;Track #eta ;Tracks (#)", 100, -etaMax_, etaMax_);
    h_phi = fs->make<TH1D>("h_phi", "#phi Distribution;Track #phi (rad);Tracks (#)", 100, -M_PI_,M_PI_);
    h_theta = fs->make<TH1D>("h_theta", "#theta Distribution ;Track #theta (rad);Tracks (#)",100,-M_PI_,M_PI_);
    h_d0 = fs->make<TH1D>("h_d0", "d_{0} ;Track d_{0} (cm);Tracks (#)", 100, -1.0,1.0);
    h_dz = fs->make<TH1D>("h_dz", "d_{z};Track d_{z} (cm);Tracks (#)", 100, -20.,20.);
    h_dxy = fs->make<TH1D>("h_dxy", "d_{xy};Track d_{xy} (cm);Tracks (#)", 100, -0.5, 0.5);
    h_chi2 = fs->make<TH1D>("h_chi2","#chi^{2} Distribution;Track #chi^{2};Tracks (#)",100,0.,100.);
    h_chi2_ndf = fs->make<TH1D>("h_chi2ndof", "#chi^{2}/NDF;Track #chi^{2}/NDF;Tracks (#)", 100, 0, 10.);
    h_chi2_Prob = fs->make<TH1D>("h_chi2Prob", "#chi^{2} Prob;Track #chi^{2} Prob;Tracks (#)", 100, 0, 10.);

    //==============Hits Directory=================//
    h_nh_Valid = fs->make<TH1I>("h_nHits", "Total Hits;Hits (#); Tracks (#)", 50, -0.5, 49.5);
    h_nh_PIXEL = fs->make<TH1I>("h_nHits_PIXEL", "Hits in PIXEL;Hits (#); Tracks (#)", 30, 0., 30.);
    h_nh_FPIX = fs->make<TH1I>("h_nHits_FPIX", "Hits in FPIX;Hits (#); Tracks (#)", 30, 0., 30.);
    h_nh_FPIXplus = fs->make<TH1I>("h_nHits_FPIXplus", "Hits in FPIX+;Hits (#); Tracks (#)", 30, 0., 30.);
    h_nh_FPIXminus = fs->make<TH1I>("h_nHits_FPIXminus", "Hits in FPIX-;Hits (#); Tracks (#)", 30, 0., 30.);
    h_nh_BPIX = fs->make<TH1I>("h_nHits_BPIX", "Hits in BPIX;Hits (#); Tracks (#)", 30, 0., 30.);
    h_nh_TIB = fs->make<TH1I>("h_nHits_TIB", "Hits in TIB;Hits (#); Tracks (#)", 30, 0., 30.);
    h_nh_TID = fs->make<TH1I>("h_nHits_TID", "Hits in TID;Hits (#); Tracks (#)", 30, 0., 30.);
    h_nh_TIDplus = fs->make<TH1I>("h_nHits_TIDplus", "Hits in TID plus;Hits [#]; Tracks [#]", 30, 0., 30.);
    h_nh_TIDminus = fs->make<TH1I>("h_nHits_TIDminus", "Hits in TID-;Hits [#]; Tracks [#]", 30, 0., 30.);
    h_nh_TOB = fs->make<TH1I>("h_nHits_TOB", "Hits in TOB;Hits (#); Tracks (#)", 30, 0., 30.);
    h_nh_TEC = fs->make<TH1I>("h_nHits_TEC", "Hits in TEC;Hits (#); Tracks (#)", 30, 0., 30.);
    h_nh_TECplus = fs->make<TH1I>("h_nHits_TECplus", "Hits in TEC+;Hits [#]; Tracks [#]", 30, 0., 30.);
    h_nh_TECminus = fs->make<TH1I>("h_nHits_TECminus", "Hits in TEC-;Hits [#]; Tracks [#]", 30, 0., 30.);
    h_nh_ENDCAP = fs->make<TH1I>("h_nHits_ENDCAP", "Hits in ENDCAP;Hits (#); Tracks (#)", 30, 0., 30.);
    h_nh_ENDCAPplus = fs->make<TH1I>("h_nHits_ENDCAPplus", "Hits ENDCAP +;Hits (#); Tracks (#)", 30, 0., 30.);
    h_nh_ENDCAPminus = fs->make<TH1I>("h_nHitsENDCAPminus", "Hits ENDCAP -;Hits (#); Tracks (#)", 30, 0., 30.);

    //==============Residuals Directory=================// 
    h_Res_BPIX_xPrime = fs->make<TH1D>("h_ResBPIXxPrime", "BPIX Track X-residuals;res_{X'} (#mm);Entries (#)", 100, -.5, .5);
    h_Res_FPIX_xPrime = fs->make<TH1D>("h_ResFPIXxPrime", "FPPIX Track X-residuals;res_{X'} (#mm);Entries (#)", 100, -.2, .2);
    h_Res_TIB_xPrime = fs->make<TH1D>("h_ResTIBxPrime", "TIB Track X-residuals;Res_{X'} (#mm);Entries (#)", 100, -2., 2.);
    h_Res_TID_xPrime = fs->make<TH1D>("h_ResTIDxPrime", "TID Track X-residuals;Res_{X'} (#mm);Entries (#)", 100, -10., 10.);
    h_Res_TOB_xPrime = fs->make<TH1D>("h_ResTOBxPrime", "TOB Track X-residuals;Res_{X'} (#mm);Entries (#)", 100, -10., 10.);
    h_Res_TEC_xPrime = fs->make<TH1D>("h_ResTECxPrime", "TEC Track X-residuals;Res_{X'} (#mm);Entries (#)", 100, -10., 10.);
    h_Res_BPIX_yPrime = fs->make<TH1D>("h_ResBPIXyPrime", "BPix Track Y-residuals;res_{Y'} (#mm);Entries (#)", 100, -2., 2.);
    h_Res_FPIX_yPrime = fs->make<TH1D>("h_ResFPIXyPrime", "FPix Track Y-residuals;res_{Y'} (#mm);Entries (#)", 100, -.5, .5);
    h_Res_TIB_yPrime = fs->make<TH1D>("h_ResTIByPrime", "TIB Track Y-residuals;Res_{Y'} (#mm);Entries (#)", 100, -40., 40.);
    h_Res_TID_yPrime = fs->make<TH1D>("h_ResTIDyPrime", "TID Track Y-residuals;Res_{Y'} (#mm);Entries (#)", 100, -40., 40.);
    h_Res_TOB_yPrime = fs->make<TH1D>("h_ResTOByPrime", "TOB Track Y-residuals;Res_{X'} (#mm);Entries (#)", 100, -40., 40.);
    h_Res_TEC_yPrime = fs->make<TH1D>("h_ResTECyPrime", "TEC Track Y-residuals;Res_{X'} (#mm);Entries (#)", 100, -40., 40.);   
}

// ------------ method called once each job just after ending the event loop  ------------
void
DatasetValidationTool::endJob()
{
   std::cout<<"Events: "<<nEvents<<std::endl;
   std::cout<<"Tracks: "<<nTracks<<std::endl;
   std::cout<<"Invalid Hits: "<<nInvalidHits<<std::endl;
   std::cout<<"SubDet0 HIts: "<<nZeroSubDetId<<std::endl;
   std::cout<<"nNotDetUnits: "<<nNotDetUnits<<std::endl;
}

// ------------ method called when starting to processes a run  -------------
void
DatasetValidationTool::beginRun(edm::Run const&, edm::EventSetup const&)
{   
   nEventsInRun=0; nTracksInRun=0;
}

// ------------ method called when ending the processing of a run  -------------
void
DatasetValidationTool::endRun(edm::Run const&, edm::EventSetup const&)
{   
   h_Events_In_Run->Fill(nEventsInRun);
   h_Tracks_In_Run->Fill(nTracksInRun);
}

// ------------ method called when starting to processes a luminosity block  ------------
void
DatasetValidationTool::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
    nEventsInLuminosity=0; nTracksInLuminosity=0;  
} 

// ------------ method called when ending the processing of a luminosity block  ------------
void
DatasetValidationTool::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
    h_Events_In_Luminosity->Fill(nEventsInLuminosity);
    h_Tracks_In_Luminosity->Fill(nTracksInLuminosity);
}

// ------------------------------ 2D Hit Function ---------------------------------------
// // WILL BREAK DOWN FOR PHASE 2 !!!!!!!!!!!!!!!!!!   
const bool DatasetValidationTool::isHit2D(const TrackingRecHit &hit)
{
     bool countStereoHitAs2D_ = true;
     // we count SiStrip stereo modules as 2D if selected via countStereoHitAs2D_
     // (since they provide theta information)
     if (!hit.isValid() || (hit.dimension() < 2 && !countStereoHitAs2D_ && !dynamic_cast<const SiStripRecHit1D *>(&hit)))
     {  return false; // real RecHit1D - but SiStripRecHit1D depends on countStereoHitAs2D_  
     }
     else 
     {
        const DetId detId(hit.geographicalId());
        if (detId.det() == DetId::Tracker)
        {
           if (detId.subdetId() == PixelSubdetector::PixelBarrel || detId.subdetId() == PixelSubdetector::PixelEndcap)
           {  return true;  // pixel is always 2D
           }
           else
           {  // should be SiStrip now
              const SiStripDetId stripId(detId);
              if(stripId.stereo())  return countStereoHitAs2D_;  // stereo modules
              else if (dynamic_cast<const SiStripRecHit1D *>(&hit) || dynamic_cast<const SiStripRecHit2D *>(&hit)) return false;  // rphi modules hit
              //the following two are not used any more since ages...
                           else if (dynamic_cast<const SiStripMatchedRecHit2D *>(&hit)) return true;  // matched is 2D       
              else if (dynamic_cast<const ProjectedSiStripRecHit2D *>(&hit))
              {
                 const ProjectedSiStripRecHit2D *pH = static_cast<const ProjectedSiStripRecHit2D *>(&hit);
                 return (countStereoHitAs2D_ && isHit2D(pH->originalHit()));  // depends on original...
              }
             else 
              { // edm::LogError("UnkownType") << "@SUB=DMRChecker::isHit2D"
                //                         << "Tracker hit not in pixel, neither SiStripRecHit[12]D nor "
                //                         << "SiStripMatchedRecHit2D nor ProjectedSiStripRecHit2D.";
                                 return false;
              }
           } 
        }
        else 
        {   // not tracker??
            // edm::LogWarning("DetectorMismatch") << "@SUB=DMRChecker::isHit2D"
            //                                  << "Hit not in tracker with 'official' dimension >=2.";
            return true;  // dimension() >= 2 so accept that...
        }
    }// never reached...
}


                
                                         
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DatasetValidationTool::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DatasetValidationTool);
