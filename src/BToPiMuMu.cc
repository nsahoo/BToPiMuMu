// -*- C++ -*-
//
// Package:    BToPiMuMu
// Class:      BToPiMuMu
// 
/**\class BToPiMuMu BToPiMuMu.cc bph-ana/BToPiMuMu/src/BToPiMuMu.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Niladribihari Sahoo <nsahoo@cern.ch>
//         Created:  Thu Oct  2 08:05:02 CEST 2014
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "MagneticField/Engine/interface/MagneticField.h"


//#include "DataFormats/Common/interface/Handle.h"
//#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"


//#include "DataFormats/TrackReco/interface/Track.h"
//#include "DataFormats/TrackReco/interface/TrackFwd.h"
//#include "DataFormats/Math/interface/Error.h"
//#include "DataFormats/Math/interface/Point3D.h"
//#include "DataFormats/VertexReco/interface/Vertex.h"
//#include "DataFormats/Math/interface/LorentzVector.h"

#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"

#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TH1.h>

using namespace std;
using namespace edm;
using namespace reco;

const int MUONMINUS_PDG_ID = 13;
const int PIONPLUS_PDG_ID = 211;
const int BPLUS_PDG_ID = 521;
const int JPSI_PDG_ID = 443;
const int PSI2S_PDG_ID = 100443;

const int ETA_PDG_ID = 221;
const int DZERO_PDG_ID = 421;
const int DSTAR2007_PDG_ID = 423;
const int DSPLUS_PDG_ID = 431;
const int DS1PLUS2460_PDG_ID = 20433;
const int SIGMACSTARPLUSPLUS_PDG_ID = 4224;
const int DELTAPLUS_PDG_ID = 2214;

const double PI = 3.141592653589793;


////// Structures //////                                                                                                                                              
struct HistArgs{
  char name[128];
  char title[128];
  int n_bins;
  double x_min;
  double x_max;
};


enum HistName{
  h_events,      
  h_mupt,
  h_mueta,
  h_mumdcabs,
  h_mumutrkr,
  h_mumutrkz,

  h_mumudca,
  h_mumuvtxcl,
  h_mumupt,
  h_mumumass,
  h_mumulxybs,

  h_mumucosalphabs,
  h_trkpt,
  h_trkdcasigbs,
  h_bvtxchisq,
  h_bvtxcl,

  h_bmass,

  kHistNameSize
};


///// Global hist args //////                                                                                                                                 

HistArgs hist_args[kHistNameSize] = {
  // name, title, n_bins, x_min, x_max                                                                                                                  
            
  {"h_events", "Processed Events", 1, 0, 1},      
  {"h_mupt", "Muon pT; [GeV]", 100, 0, 30},
  {"h_mueta", "Muon eta", 100, 0, 3},
  {"h_mumdcabs", "#mu^{-} DCA beam spot; DCA [cm]", 100, 0, 10},
  {"h_mumutrkr", "#mu^{+}#mu^{-} distance in phi-eta; [cm]", 100, 0, 50},
  {"h_mumutrkz", "#mu^{+}#mu^{-} distance in Z; [cm]", 100, 0, 100},

  {"h_mumudca", "#mu^{+}#mu^{-} DCA; [cm]", 100, 0, 20},
  {"h_mumuvtxcl", "#mu^{+}#mu^{-} vertex CL", 100, 0, 1},
  {"h_mumupt", "#mu^{+}#mu^{-} pT ; pT [GeV]", 100, 0, 50},
  {"h_mumumass", "#mu^{+}#mu^{-} invariant mass; M(#mu^{+}#mu^{-}) [GeV/c^{2}]",
   100, 2, 20},
  {"h_mumulxybs", "#mu^{+}#mu^{-} Lxy #sigma beam spot", 100, 0, 100},

  {"h_mumucosalphabs", "#mu^{+}#mu^{-} cos #alpha beam spot", 100, 0, 1},
  {"h_trkpt", "Pion track pT; pT [GeV]", 100, 0, 20},
  {"h_trkdcabssig", "Pion track DCA/#sigma beam spot; DCA/#sigma", 100, 0, 10},
  {"h_bvtxchisq", "B decay vertex chisq", 100, 0, 1000},
  {"h_bvtxcl", "B decay vertex CL", 100, 0, 1},

  {"h_bmass", "B mass; M(B) [GeV]", 100, 0, 20},

};

// Define histograms                                                                                                                                         
TH1F *histos[kHistNameSize];

//
// class declaration
//

class BToPiMuMu : public edm::EDAnalyzer {
   public:
      explicit BToPiMuMu(const edm::ParameterSet&);
      ~BToPiMuMu();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  void buildBuToPiMuMu(const edm::Event &);
  bool hasBeamSpot(const edm::Event&);
  bool hasPrimaryVertex(const edm::Event &);
  void clearVariables();



      // ----------member data ---------------------------

  // --- begin input from python file ---                                                                                                                             
  string OutputFileName_;

  // particle properties                                                                                                                                             
  ParticleMass MuonMass_;
  float MuonMassErr_;
  ParticleMass PionMass_;
  float PionMassErr_;
  double BuMass_;

  // labels                                                                                                                                                           
  edm::InputTag GenParticlesLabel_;
  edm::InputTag TriggerResultsLabel_;
  edm::InputTag BeamSpotLabel_;
  edm::InputTag VertexLabel_;
  edm::InputTag MuonLabel_;
  
  edm::InputTag TrackLabel_;
  vector<string> TriggerNames_;
  vector<string> LastFilterNames_;

  //gen particle
  bool IsMonteCarlo_;
  bool KeepGENOnly_;
  double TruthMatchMuonMaxR_;
  double TruthMatchPionMaxR_;

  //pre-selection cuts
  double MuonMinPt_;
  double MuonMaxEta_;
  double MuonMaxDcaBs_;
  double TrkMinPt_;
  double TrkMinDcaSigBs_;
  double TrkMaxR_;
  double TrkMaxZ_;
  double MuMuMaxDca_;
  double MuMuMinVtxCl_;
  double MuMuMinPt_;
  double MuMuMinInvMass_;
  double MuMuMaxInvMass_;
  double MuMuMinLxySigmaBs_;
  double MuMuMinCosAlphaBs_;

  double BMinVtxCl_;
  double BMinMass_;
  double BMaxMass_;

  // ---- end input from python file -----

  // Across the event
  map<string, string> mapTriggerToLastFilter_;
  reco::BeamSpot beamSpot_;
  edm::ESHandle<MagneticField> bFieldHandle_;
  reco::Vertex primaryVertex_;

  // ---- Root Variables ----
  TFile* fout_;
  TTree* tree_;

  unsigned int run, event, lumiblock, nprivtx;
  vector<string> *triggernames;
  vector<int> *triggerprescales;

  // dimuon
  vector<double> *mumdcabs, *mumdcabserr, *mumpx, *mumpy, *mumpz;
  vector<double> *mupdcabs, *mupdcabserr, *muppx, *muppy, *muppz;
  vector<double> *mumutrkr, *mumutrkz , *mumudca;
  vector<double> *mumuvtxcl, *mumulsbs, *mumulsbserr;
  vector<double> *mumucosalphabs, *mumucosalphabserr;
  vector<double> *mumumass, *mumumasserr;

  // soft muon variables
  vector<bool> *mumisgoodmuon, *mupisgoodmuon ;
  vector<int> *mumnpixhits, *mupnpixhits, *mumnpixlayers, *mupnpixlayers;
  vector<int> *mumntrkhits, *mupntrkhits, *mumntrklayers, *mupntrklayers;
  vector<double> *mumnormchi2, *mupnormchi2;
  vector<double> *mumdxyvtx, *mupdxyvtx, *mumdzvtx, *mupdzvtx;
  vector<string> *mumtriglastfilter, *muptriglastfilter;
  vector<double> *mumpt, *muppt, *mumeta, *mupeta;

  // Pion track
  vector<int> *trkchg; // +1 for Pi+, -1 for Pi-
  vector<double> *trkpx, *trkpy, *trkpz, *trkpt;
  vector<double> *trkdcabs, *trkdcabserr;

  // B+ and B-
  int nb;
  vector<int> *bchg; // +1 for b+, -1 for b-
  vector<double> *bpx, *bpxerr, *bpy, *bpyerr, *bpz, *bpzerr, *bmass, *bmasserr;
  vector<double> *bvtxcl, *bvtxx, *bvtxxerr, *bvtxy, *bvtxyerr, *bvtxz, *bvtxzerr;
  vector<double> *bcosalphabs, *bcosalphabserr, *blsbs, *blsbserr, *bctau, *bctauerr;
  vector<double> *bcosalphabs2d, *bcosalphabs2derr; 

  // For MC
  int genbchg; // +1 for b+, -1 for b-
  double genbpx, genbpy, genbpz;

  int gentrkchg;
  double gentrkpx, gentrkpy, gentrkpz;
  double genmumpx, genmumpy, genmumpz;
  double genmuppx, genmuppy, genmuppz;

  string decname;

  vector<bool> *istruemum, *istruemup, *istruetrk, *istruebu;

  TDatime t_begin_ , t_now_ ; 
  int n_processed_, n_selected_;

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
BToPiMuMu::BToPiMuMu(const edm::ParameterSet& iConfig):

    OutputFileName_(iConfig.getParameter<string>("OutputFileName")),

    // particle properties
    MuonMass_(iConfig.getUntrackedParameter<double>("MuonMass")),
    MuonMassErr_(iConfig.getUntrackedParameter<double>("MuonMassErr")),
    PionMass_(iConfig.getUntrackedParameter<double>("PionMass")),
    PionMassErr_(iConfig.getUntrackedParameter<double>("PionMassErr")),
    BuMass_(iConfig.getUntrackedParameter<double>("BuMass")),

    // labels
    GenParticlesLabel_(iConfig.getParameter<edm::InputTag>("GenParticlesLabel")),
    TriggerResultsLabel_(iConfig.getParameter<edm::InputTag>("TriggerResultsLabel")),
    BeamSpotLabel_(iConfig.getParameter<edm::InputTag>("BeamSpotLabel")),
    VertexLabel_(iConfig.getParameter<edm::InputTag>("VertexLabel")),
    MuonLabel_(iConfig.getParameter<edm::InputTag>("MuonLabel")),
    TrackLabel_(iConfig.getParameter<edm::InputTag>("TrackLabel")),
    TriggerNames_(iConfig.getParameter< vector<string> >("TriggerNames")),
    LastFilterNames_(iConfig.getParameter< vector<string> >("LastFilterNames")),

    // gen particle
    IsMonteCarlo_(iConfig.getUntrackedParameter<bool>("IsMonteCarlo")),
    KeepGENOnly_(iConfig.getUntrackedParameter<bool>("KeepGENOnly")),
    TruthMatchMuonMaxR_(iConfig.getUntrackedParameter<double>("TruthMatchMuonMaxR")),
    TruthMatchPionMaxR_(iConfig.getUntrackedParameter<double>("TruthMatchPionMaxR")),

    // pre-selection cuts
    MuonMinPt_(iConfig.getUntrackedParameter<double>("MuonMinPt")),
    MuonMaxEta_(iConfig.getUntrackedParameter<double>("MuonMaxEta")),
    MuonMaxDcaBs_(iConfig.getUntrackedParameter<double>("MuonMaxDcaBs")),
    TrkMinPt_(iConfig.getUntrackedParameter<double>("TrkMinPt")),
    TrkMinDcaSigBs_(iConfig.getUntrackedParameter<double>("TrkMinDcaSigBs")),
    TrkMaxR_(iConfig.getUntrackedParameter<double>("TrkMaxR")),
    TrkMaxZ_(iConfig.getUntrackedParameter<double>("TrkMaxZ")),
    MuMuMaxDca_(iConfig.getUntrackedParameter<double>("MuMuMaxDca")),
    MuMuMinVtxCl_(iConfig.getUntrackedParameter<double>("MuMuMinVtxCl")),
    MuMuMinPt_(iConfig.getUntrackedParameter<double>("MuMuMinPt")),
    MuMuMinInvMass_(iConfig.getUntrackedParameter<double>("MuMuMinInvMass")),
    MuMuMaxInvMass_(iConfig.getUntrackedParameter<double>("MuMuMaxInvMass")),
    MuMuMinLxySigmaBs_(iConfig.getUntrackedParameter<double>("MuMuMinLxySigmaBs")),
    MuMuMinCosAlphaBs_(iConfig.getUntrackedParameter<double>("MuMuMinCosAlphaBs")),
    
   
    BMinVtxCl_(iConfig.getUntrackedParameter<double>("BMinVtxCl")),
    BMinMass_(iConfig.getUntrackedParameter<double>("BMinMass")),
    BMaxMass_(iConfig.getUntrackedParameter<double>("BMaxMass")),

    tree_(0),
    triggernames(0), triggerprescales(0),
    mumdcabs(0), mumdcabserr(0), mumpx(0), mumpy(0), mumpz(0),
    mupdcabs(0), mupdcabserr(0), muppx(0), muppy(0), muppz(0),
    mumutrkr(0), mumutrkz(0), mumudca(0), mumuvtxcl(0), mumulsbs(0),
    mumulsbserr(0), mumucosalphabs(0), mumucosalphabserr(0),
    mumumass(0), mumumasserr(0),
    mumisgoodmuon(0), mupisgoodmuon(0),
    mumnpixhits(0), mupnpixhits(0), mumnpixlayers(0), mupnpixlayers(0),
    mumntrkhits(0), mupntrkhits(0), mumntrklayers(0), mupntrklayers(0),
    mumnormchi2(0), mupnormchi2(0), mumdxyvtx(0), mupdxyvtx(0),
    mumdzvtx(0), mupdzvtx(0), mumtriglastfilter(0), muptriglastfilter(0),
    mumpt(0), muppt(0), mumeta(0), mupeta(0),

    trkchg(0), trkpx(0), trkpy(0), trkpz(0), trkpt(0),
    trkdcabs(0), trkdcabserr(0),

    nb(0), bchg(0), bpx(0), bpxerr(0), bpy(0), bpyerr(0), bpz(0), bpzerr(0),
    bmass(0), bmasserr(0),
    bvtxcl(0), bvtxx(0), bvtxxerr(0), bvtxy(0), bvtxyerr(0), bvtxz(0), bvtxzerr(0),
    bcosalphabs(0), bcosalphabserr(0), blsbs(0), blsbserr(0), bctau(0), bctauerr(0),
    bcosalphabs2d(0), bcosalphabs2derr(0), 

    genbchg(0),
    genbpx(0), genbpy(0), genbpz(0),
    gentrkchg(0), gentrkpx(0), gentrkpy(0), gentrkpz(0),
    genmumpx(0), genmumpy(0), genmumpz(0),
    genmuppx(0), genmuppy(0), genmuppz(0),

    decname(""),
    istruemum(0), istruemup(0), istruetrk(0), istruebu(0)
{
   //now do what ever initialization is needed
      assert(TriggerNames_.size() == LastFilterNames_.size());
    for (size_t i = 0; i < TriggerNames_.size(); ++i)
      mapTriggerToLastFilter_[TriggerNames_[i]] = LastFilterNames_[i];
}



BToPiMuMu::~BToPiMuMu()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
BToPiMuMu::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //   using namespace edm;




}


// ------------ method called once each job just before starting event loop  ------------
void 
BToPiMuMu::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
BToPiMuMu::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
BToPiMuMu::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
BToPiMuMu::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
BToPiMuMu::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
BToPiMuMu::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
BToPiMuMu::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}





//define this as a plug-in
DEFINE_FWK_MODULE(BToPiMuMu);
