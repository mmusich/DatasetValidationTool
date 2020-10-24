// -*- C++ -*-
//
// Package:    DatasetValidation/DatasetValidationTool
// Class:      DatasetValidationTool
//
/**\class DatasetValidationTool DatasetValidationTool.cc DatasetValidation/DatasetValidationTool/plugins/DatasetValidationTool.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Saumya Saumya
//         Created:  Tue, 06 Oct 2020 13:05:38 GMT
//
//


// system include files
#include <memory>

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


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


class DatasetValidationTool : public edm::one::EDAnalyzer<edm::one::SharedResources> 
{
   public:
      explicit DatasetValidationTool(const edm::ParameterSet&);
      ~DatasetValidationTool();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<reco::TrackCollection> tracksToken_;

      edm::Service<TFileService> fs;
      int nTracks,nEvents;
      int nHits_Total,nHits_PIXEL,nHits_FPIX,nHits_FPIXplus,nHits_FPIXminus,nHits_BPIX,nHits_TIB,nHits_TID,nHits_TIDplus,nHits_TIDminus,nHits_TOB,nHits_TEC,nHits_TECplus,nHits_TECminus,nHits_ENDCAP,nHits_ENDCAPplus,nHits_ENDCAPminus;
     float etaMax_=3.0,M_PI_=3.14159;

      TH1D *hcharge;
      TH1D *hp;
      TH1D *hpt;
      TH1D *heta;
      TH1D *hphi;
      TH1D *htheta;
      TH1D *hd0;
      TH1D *hdz;
      TH1D *hdxy;
      TH1D *hchi2;
      TH1D *hchi2norm;
      TH1D *hchi2Prob;
      TH1D *hntrk;     

      TH1I *nh;
      TH1I *nh_BPIX;
      TH1I *nh_FPIX;         
      TH1I *nh_FPIXplus;
      TH1I *nh_FPIXminus;
      TH1I *nh_PIXEL;
      TH1I *nh_TIB;
      TH1I *nh_TID;
 /*     TH1I *nh_TIDplus;
      TH1I *nh_TIDminus;
   */   TH1I *nh_TOB;
      TH1I *nh_TEC;
 /*     TH1I *nh_TECplus;
   */   TH1I *nh_TECminus;
      TH1I *nh_ENDCAP;
      TH1I *nh_ENDCAPplus;
      TH1I *nh_ENDCAPminus;
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
  tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks")))
{
   //now do what ever initialization is needed
   usesResource("TFileService");

}


DatasetValidationTool::~DatasetValidationTool()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

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

   for(auto track:*tracksHandle_)
   {
     hcharge->Fill(track.charge());
     hcharge->Fill(track.charge());
     hp->Fill(track.p());
     hpt->Fill(track.pt());
     heta->Fill(track.eta());
     hphi->Fill(track.phi());
     htheta->Fill(track.theta());
     hd0->Fill(track.d0());
     hdz->Fill(track.dz());
     hdxy->Fill(track.dxy());
     hchi2->Fill(track.chi2());     
     hchi2norm->Fill(track.normalizedChi2());
     hchi2Prob->Fill(TMath::Prob(track.chi2(),track.ndof()));

     nHits_Total=0;
     nHits_PIXEL=0;
     nHits_FPIX=0;
     nHits_FPIXplus=0;
     nHits_FPIXminus=0;
     nHits_BPIX=0;
     nHits_TIB=0;
     nHits_TID=0;
//     nHits_TIDplus=0;
//     nHits_TIDminus=0;
     nHits_TOB=0;
     nHits_TEC=0;
//     nHits_TECplus=0;
//     nHits_TECminus=0;
     nHits_ENDCAP=0;
     nHits_ENDCAPplus=0;
     nHits_ENDCAPminus=0;

     int h_index=0;
     auto const &residuals = track.extra()->residuals();
     for(auto iHit = track.recHitsBegin(); iHit!=track.recHitsEnd(); ++iHit,++h_index)
     {
         const DetId detId=(*iHit)->geographicalId();
         const int SubDetId = detId.subdetId();
         if (!(*iHit)->isValid()) continue;

         nHits_Total++;
//         double resX=residuals.residualX(h_index);
//         double resY=residuals.residualY(h_index);

         //    Hit in PixelBarrel  //
         if (SubDetId == PixelSubdetector::PixelBarrel)
         {
              nHits_BPIX++; nHits_PIXEL++;
         }
         //                      Hit information in PixelEndcap                  //   
         else if (SubDetId == PixelSubdetector::PixelEndcap)
         {
              nHits_FPIX++; nHits_PIXEL++; nHits_ENDCAP++;
         }
         //                         Hit information in TIB                       //
         else if (SubDetId == SiStripDetId::TIB)
         {
              nHits_TIB++;
         }
         //                         Hit information in TID                       //
         else if (SubDetId == SiStripDetId::TID)
         {
              nHits_TID++; nHits_ENDCAP++;
         }
         //                        Hit information in TOB                       //
         else if (SubDetId == SiStripDetId::TOB)
         {
              nHits_TOB++;
         }
         //                        Hit information in TEC                       //
         else if (SubDetId == SiStripDetId::TEC)
         {
              nHits_TEC++; nHits_ENDCAP++;
         }

     }  //Hits Loop
    nh->Fill(nHits_Total);
    nh_BPIX->Fill(nHits_BPIX);
    nh_FPIX->Fill(nHits_FPIX);
    nh_FPIXplus->Fill(nHits_FPIXplus);
    nh_FPIXminus->Fill(nHits_FPIXminus);
    nh_PIXEL->Fill(nHits_PIXEL);
    nh_TIB->Fill(nHits_TIB);
    nh_TID->Fill(nHits_TID);
 //   nh_TIDplus->Fill(nHits_TIDplus);
 //   nh_TIDminus->Fill(nHits_TIDminus); 
    nh_TOB->Fill(nHits_TOB);
    nh_TEC->Fill(nHits_TEC);
 //   nh_TECplus->Fill(nHits_TECplus);
 //   nh_TECminus->Fill(nHits_TECminus);
    nh_ENDCAP->Fill(nHits_ENDCAP);
    nh_ENDCAPplus->Fill(nHits_ENDCAPplus);
    nh_ENDCAPminus->Fill(nHits_ENDCAPminus);

 
  }  //Tracks Loop

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void
DatasetValidationTool::beginJob()
{
    

    TFileDirectory EventDir = fs->mkdir("Event");
    TFileDirectory HitsDir = fs->mkdir("Hits");

    TH1D::SetDefaultSumw2(kTRUE);
   
    hcharge = EventDir.make<TH1D>("h_charge", "Charge;Charge [e];Tracks", 5, -2.5, 2.5);
    hp = EventDir.make<TH1D>("h_p", "Momentum;p [GeV];Tracks", 100, 0.,100.);
    hpt = EventDir.make<TH1D>("h_pt", "Transverse Momentum;p_{T} [GeV];Tracks", 100, 0., 100.);
    heta = EventDir.make<TH1D>("h_eta", "#eta Distribution;#eta ;Tracks", 100, -etaMax_, etaMax_);
    hphi = EventDir.make<TH1D>("h_phi", "#phi Distribution;#phi [rad];Tracks", 100, -M_PI_,M_PI_);
    htheta = EventDir.make<TH1D>("h_theta", "#theta Distribution ;#theta [rad];Tracks",100,-M_PI_,M_PI_);
    hd0 = EventDir.make<TH1D>("h_d0", "d_{0} ;Track d_{0}[cm];Tracks", 100, -1.0,1.0);
    hdz = EventDir.make<TH1D>("h_dz", "d_{z};Track d_{z} [cm];Tracks", 100, -20.,20.);
    hdxy = EventDir.make<TH1D>("h_dxy", "d_{xy};Track d_{xy} [cm];Tracks", 100, -0.5, 0.5);
    hchi2 = EventDir.make<TH1D>("h_chi2","#chi^{2} Distribution;#chi^{2};Tracks",100,0.,100.);
    hchi2norm = EventDir.make<TH1D>("h_chi2ndof", "#chi^{2}/NDF;#chi^{2}/NDF;Tracks", 100, 0, 10.);
    hchi2Prob = EventDir.make<TH1D>("h_chi2Prob", "#chi^{2} Prob;#chi^{2}/ndf;Tracks", 100, 0, 10.);

    nh = HitsDir.make<TH1I>("h_nHits", "Total Hits;Hits [#]; Tracks [#]", 50, -0.5, 49.5);
    nh_PIXEL = HitsDir.make<TH1I>("h_nHits_PIXEL", "Hits in PIXEL;Hits [#]; Tracks [#]", 30, 0., 30.);
    nh_FPIX = HitsDir.make<TH1I>("h_nHits_FPIX", "Hits in FPIX;Hits [#]; Tracks [#]", 30, 0., 30.);
    nh_FPIXplus = HitsDir.make<TH1I>("h_nHits_FPIXplus", "Hits in FPIX+;Hits [#]; Tracks [#]", 30, 0., 30.);
    nh_FPIXminus = HitsDir.make<TH1I>("h_nHits_FPIXminus", "Hits in FPIX-;Hits [#]; Tracks [#]", 30, 0., 30.);
    nh_BPIX = HitsDir.make<TH1I>("h_nHits_BPIX", "Hits in BPIX;Hits [#]; Tracks [#]", 30, 0., 30.);
    nh_TIB = HitsDir.make<TH1I>("h_nHits_TIB", "Hits in TIB;Hits [#]; Tracks [#]", 30, 0., 30.);
    nh_TID = HitsDir.make<TH1I>("h_nHits_TID", "Hits in TID;Hits [#]; Tracks [#]", 30, 0., 30.);
 //   nh_TIDplus = HitsDir.make<TH1I>("h_nHits_TIDplus", "Hits in TID plus;Hits [#]; Tracks [#]", 30, 0., 30.);
 //   nh_TIDminus = HitsDir.make<TH1I>("h_nHits_TIDminus", "Hits in TID-;Hits [#]; Tracks [#]", 30, 0., 30.);
    nh_TOB = HitsDir.make<TH1I>("h_nHits_TOB", "Hits in TOB;Hits [#]; Tracks [#]", 30, 0., 30.);
    nh_TEC = HitsDir.make<TH1I>("h_nHits_TEC", "Hits in TEC;Hits [#]; Tracks [#]", 30, 0., 30.);
//    nh_TECplus = HitsDir.make<TH1I>("h_nHits_TECplus", "Hits in TEC+;Hits [#]; Tracks [#]", 30, 0., 30.);
  //  nh_TECminus = HitsDir.make<TH1I>("h_nHits_TECminus", "Hits in TEC-;Hits [#]; Tracks [#]", 30, 0., 30.);
    nh_ENDCAP = HitsDir.make<TH1I>("h_nHits_ENDCAP", "Hits in ENDCAP;Hits [#]; Tracks [#]", 30, 0., 30.);
    nh_ENDCAPplus = HitsDir.make<TH1I>("h_nHits_ENDCAPplus", "Hits ENDCAP +;Hits [#]; Tracks [#]", 30, 0., 30.);
    nh_ENDCAPminus = HitsDir.make<TH1I>("h_nHitsENDCAPminus", "Hits ENDCAP -;Hits [#]; Tracks [#]", 30, 0., 30.);
}

// ------------ method called once each job just after ending the event loop  ------------
void
DatasetValidationTool::endJob()
{
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
