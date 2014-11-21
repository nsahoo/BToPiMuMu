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
// copyright @2014   Niladribihari Sahoo
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
//#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

//#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
//#include "DataFormats/Math/interface/Error.h"
//#include "DataFormats/Math/interface/Point3D.h"
//#include "DataFormats/VertexReco/interface/Vertex.h"
//#include "DataFormats/Math/interface/LorentzVector.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

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

  void computeLS (double, double, double, double, double, double, double,
		  double, double, double, double, double, double, double,
		  double, double, double, double, double*, double*);

  void computeCosAlpha (double, double, double, double, double,
			double,
			double, double, double, double,
			double, double, double, double,
			double, double, double, double,
			double*, double*);

  void computeCosAlpha2d(double, double, double, double, double, 
			 double, double, double, double, double,
			 double, double, double, double,
			 double, double, double, double,
			 double*, double*);

  void computeCtau(RefCountedKinematicTree, double &, double &);

  bool hasBeamSpot(const edm::Event&);

  bool hasGoodClosestApproachTracks (const reco::TransientTrack,
				     const reco::TransientTrack,
				     double&, double &, double &);

  bool hasPrimaryVertex(const edm::Event &);
  void clearVariables();
  void hltReport(const edm::Event&);
  void saveGenInfo(const edm::Event&);
  void saveTruthMatch(const edm::Event&);
  bool hasGoodPionTrack(const edm::Event&, const pat::GenericParticle, double &);
  bool matchMuonTrack (const edm::Event&, const reco::TrackRef);
  bool hasGoodMuonDcaBs (const reco::TransientTrack, double &, double &);
  bool hasGoodTrackDcaBs (const reco::TransientTrack, double &, double &);
  bool hasGoodBuVertex(const reco::TrackRef, const reco::TrackRef,
		       const reco::TrackRef,
		       double &, double &,
		       RefCountedKinematicTree &);
  bool hasGoodMuMuVertex (const reco::TransientTrack, const reco::TransientTrack,
			  reco::TransientTrack &, reco::TransientTrack &,
			  double &, double &, double &, double &, double &,
			  double &, double &, double &);
  bool hasGoodBuMass(RefCountedKinematicTree, double &);

  void saveBuToPiMuMu(RefCountedKinematicTree);
  void saveBuVertex(RefCountedKinematicTree);
  void saveBuCosAlpha(RefCountedKinematicTree);
  void saveBuCosAlpha2d(RefCountedKinematicTree); 
  void saveBuLsig(RefCountedKinematicTree);
  void saveBuCtau(RefCountedKinematicTree);


  void saveSoftMuonVariables(pat::Muon, pat::Muon, reco::TrackRef, reco::TrackRef);
  void saveDimuVariables(double, double, double, double, double, double,
			 double, double, double, double, double, double,
			 double, double);




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

  n_processed_ += 1;
  histos[h_events]->Fill(0);

  clearVariables();

  run = iEvent.id().run() ;
  event = iEvent.id().event() ;
  lumiblock = iEvent.luminosityBlock();


}


// ------------ method called once each job just before starting event loop  ------------
void 
BToPiMuMu::beginJob()
{


  t_begin_.Set();
  printf("\n ---------- Begin Job ---------- \n"); /* added */
  t_begin_.Print();


  n_processed_ = 0;
  n_selected_ = 0;


  fout_ = new TFile(OutputFileName_.c_str(), "RECREATE");
  fout_->cd();


  for(int i=0; i<kHistNameSize; i++) {
    histos[i] = new TH1F(hist_args[i].name, hist_args[i].title,
			 hist_args[i].n_bins,
			 hist_args[i].x_min, hist_args[i].x_max);
  }


  tree_ = new TTree ("tree", "BToPiMuMu");

  tree_->Branch("run", &run, "run/i");
  tree_->Branch("event", &event, "event/i");
  tree_->Branch("lumiblock", &lumiblock, "lumiblock/i");
  tree_->Branch("nprivtx", &nprivtx, "nprivtx/i");
  tree_->Branch("triggernames", &triggernames);
  tree_->Branch("triggerprescales", &triggerprescales);
  tree_->Branch("mumdcabs", &mumdcabs);
  tree_->Branch("mumdcabserr", &mumdcabserr);
  tree_->Branch("mumpx", &mumpx);
  tree_->Branch("mumpy", &mumpy);
  tree_->Branch("mumpz", &mumpz);
  tree_->Branch("mupdcabs", &mupdcabs);
  tree_->Branch("mupdcabserr", &mupdcabserr);
  tree_->Branch("muppx", &muppx);
  tree_->Branch("muppy", &muppy);
  tree_->Branch("muppz", &muppz);
  tree_->Branch("mumutrkr", &mumutrkr);
  tree_->Branch("mumutrkz", &mumutrkz);
  tree_->Branch("mumudca", &mumudca);
  tree_->Branch("mumuvtxcl", &mumuvtxcl);
  tree_->Branch("mumulsbs", &mumulsbs);
  tree_->Branch("mumulsbserr", &mumulsbserr);
  tree_->Branch("mumucosalphabs", &mumucosalphabs);
  tree_->Branch("mumucosalphabserr", &mumucosalphabserr);
  tree_->Branch("mumumass", &mumumass);
  tree_->Branch("mumumasserr", &mumumasserr);
  tree_->Branch("mumisgoodmuon", &mumisgoodmuon);
  tree_->Branch("mupisgoodmuon", &mupisgoodmuon);
  tree_->Branch("mumnpixhits", &mumnpixhits);
  tree_->Branch("mupnpixhits", &mupnpixhits);
  tree_->Branch("mumnpixlayers", &mumnpixlayers);
  tree_->Branch("mupnpixlayers", &mupnpixlayers);
  tree_->Branch("mumntrkhits", &mumntrkhits);
  tree_->Branch("mupntrkhits", &mupntrkhits);
  tree_->Branch("mumntrklayers", &mumntrklayers);
  tree_->Branch("mupntrklayers", &mupntrklayers);
  tree_->Branch("mumnormchi2", &mumnormchi2);
  tree_->Branch("mupnormchi2", &mupnormchi2);
  tree_->Branch("mumdxyvtx", &mumdxyvtx);
  tree_->Branch("mupdxyvtx", &mupdxyvtx);
  tree_->Branch("mumdzvtx", &mumdzvtx);
  tree_->Branch("mupdzvtx", &mupdzvtx);
  tree_->Branch("mumtriglastfilter", &mumtriglastfilter);
  tree_->Branch("muptriglastfilter", &muptriglastfilter);
  tree_->Branch("mumpt", &mumpt);
  tree_->Branch("muppt", &muppt);
  tree_->Branch("mumeta", &mumeta);
  tree_->Branch("mupeta", &mupeta);
  tree_->Branch("trkchg", &trkchg);
  tree_->Branch("trkpx", &trkpx);
  tree_->Branch("trkpy", &trkpy);
  tree_->Branch("trkpz", &trkpz);
  tree_->Branch("trkpt", &trkpt);
  tree_->Branch("trkdcabs", &trkdcabs);
  tree_->Branch("trkdcabserr", &trkdcabserr);

  tree_->Branch("nb", &nb, "nb/I");
  tree_->Branch("bchg", &bchg);
  tree_->Branch("bpx", &bpx);
  tree_->Branch("bpxerr", &bpxerr);
  tree_->Branch("bpy", &bpy);
  tree_->Branch("bpyerr", &bpyerr);
  tree_->Branch("bpz", &bpz);
  tree_->Branch("bpzerr", &bpzerr);
  tree_->Branch("bmass", &bmass);
  tree_->Branch("bmasserr", &bmasserr);
  tree_->Branch("bvtxcl", &bvtxcl);
  tree_->Branch("bvtxx", &bvtxx);
  tree_->Branch("bvtxxerr", &bvtxxerr);
  tree_->Branch("bvtxy", &bvtxy);
  tree_->Branch("bvtxyerr", &bvtxyerr);
  tree_->Branch("bvtxz", &bvtxz);
  tree_->Branch("bvtxzerr", &bvtxzerr);
  tree_->Branch("bcosalphabs", &bcosalphabs);
  tree_->Branch("bcosalphabserr", &bcosalphabserr);

  tree_->Branch("bcosalphabs2d", &bcosalphabs2d);
  tree_->Branch("bcosalphabs2derr", &bcosalphabs2derr);
  tree_->Branch("blsbs", &blsbs);
  tree_->Branch("blsbserr", &blsbserr);
  tree_->Branch("bctau", &bctau);
  tree_->Branch("bctauerr", &bctauerr);


}

// ------------ method called once each job just after ending the event loop  ------------
void 
BToPiMuMu::endJob() 
{

  fout_->cd();
  tree_->Write();

  for(int i = 0; i < kHistNameSize; i++) {
    histos[i]->Write();
    histos[i]->Delete();
  }

  fout_->Close();
  t_now_.Set();
  printf(" \n ---------- End Job ---------- \n" ) ; /* added */
  t_now_.Print();
  printf(" processed: %i \n selected: %i \n \
duration: %i sec \n rate: %g evts/sec\n",
	 n_processed_, n_selected_,
	 t_now_.Convert() - t_begin_.Convert(),
	 float(n_processed_)/(t_now_.Convert()-t_begin_.Convert()) );

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

void 
BToPiMuMu::clearVariables(){

  run = 0;
  event = 0;
  lumiblock = 0;
  nprivtx = 0;
  triggernames->clear();
  triggerprescales->clear();
  mumdcabs->clear(); mumdcabserr->clear(); mumpx->clear(); mumpy->clear();
  mumpz->clear();
  mupdcabs->clear(); mupdcabserr->clear(); muppx->clear(); muppy->clear();
  muppz->clear();
  mumutrkr->clear(); mumutrkz->clear();
  mumudca->clear(); mumuvtxcl->clear(); mumulsbs->clear(); mumulsbserr->clear();
  mumucosalphabs->clear(); mumucosalphabserr->clear();
  mumumass->clear(); mumumasserr->clear();
  mumisgoodmuon->clear(); mupisgoodmuon->clear();
  mumnpixhits->clear(); mupnpixhits->clear(); mumnpixlayers->clear();
  mupnpixlayers->clear();
  mumntrkhits->clear(); mupntrkhits->clear(); mumntrklayers->clear();
  mupntrklayers->clear();
  mumnormchi2->clear(); mupnormchi2->clear(); mumdxyvtx->clear(); mupdxyvtx->clear();
  mumdzvtx->clear(); mupdzvtx->clear(); mumtriglastfilter->clear();
  muptriglastfilter->clear();
  mumpt->clear(); muppt->clear();
  mumeta->clear(); mupeta->clear();
  trkchg->clear(); trkpx->clear(); trkpy->clear(); trkpz->clear();
  trkpt->clear();
  trkdcabs->clear(); trkdcabserr->clear();
  bpx->clear(); bpxerr->clear(); bpy->clear(); bpyerr->clear();
  bpz->clear(); bpzerr->clear(); nb = 0;
  bchg->clear(); bmass->clear(); bmasserr->clear();
  bvtxcl->clear(); bvtxx->clear(); bvtxxerr->clear(); bvtxy->clear();
  bvtxyerr->clear();
  bvtxz->clear(); bvtxzerr->clear(); bcosalphabs->clear(); bcosalphabserr->clear();
  blsbs->clear(); blsbserr->clear(); bctau->clear(); bctauerr->clear();
  bcosalphabs2d->clear(); bcosalphabs2derr->clear();

}

void
BToPiMuMu::hltReport(const edm::Event& iEvent)
{
  edm::Handle<edm::TriggerResults> hltTriggerResults;
  try {iEvent.getByLabel( TriggerResultsLabel_, hltTriggerResults ); }
  catch ( ... ) { edm::LogInfo("myHLT")
      << __LINE__ << " : couldn't get handle on HLT Trigger" ; }
  
  HLTConfigProvider hltConfig_;
  if (hltTriggerResults.isValid()) {
    const edm::TriggerNames& triggerNames_ = iEvent.triggerNames(*hltTriggerResults);

    for (unsigned int itrig = 0; itrig < hltTriggerResults->size(); itrig++){

      // Only consider the triggered case.
      if ((*hltTriggerResults)[itrig].accept() == 1){

        string triggername = triggerNames_.triggerName(itrig);        
        int triggerprescale = hltConfig_.prescaleValue(itrig, triggername);

        // Loop over our interested HLT trigger names to find if this event contains.
        for (unsigned int it=0; it<TriggerNames_.size(); it++){
	  if (triggername.find(TriggerNames_[it]) != string::npos) {
	    // save the no versioned case
	    triggernames->push_back(TriggerNames_[it]);
	    triggerprescales->push_back(triggerprescale);

	  }}}}}
}

bool
BToPiMuMu::hasBeamSpot(const edm::Event& iEvent)
{
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel(BeamSpotLabel_, beamSpotHandle);
  
  if ( ! beamSpotHandle.isValid() ) {
    edm::LogError("myBeam") << "No beam spot available from EventSetup" ;
    return false;
  }
  
  beamSpot_ = *beamSpotHandle;
  return true;
}

bool
BToPiMuMu::hasPrimaryVertex(const edm::Event& iEvent)
{
  edm::Handle<reco::VertexCollection> recVtxs;
  iEvent.getByLabel(VertexLabel_, recVtxs);
  nprivtx = recVtxs->size();

  for (std::vector<reco::Vertex>::const_iterator iVertex = recVtxs->begin();
       iVertex != recVtxs->end(); iVertex++) {
    primaryVertex_ = *(iVertex);
    if (primaryVertex_.isValid()) break;
  }

  if (!primaryVertex_.isValid()) return false;
 
  return true;
}

bool
BToPiMuMu::hasGoodPionTrack(const edm::Event& iEvent,
			   const pat::GenericParticle iTrack,
			   double & pion_trk_pt)
{
  reco::TrackRef theTrackRef = iTrack.track();
  if ( theTrackRef.isNull() ) return false;

  // veto muon tracks
  if ( matchMuonTrack(iEvent, theTrackRef) ) return false;
   
  // check the track kinematics
  pion_trk_pt = theTrackRef->pt();

  if ( theTrackRef->pt() < TrkMinPt_ ) return false;

  return true;
}

bool
BToPiMuMu::matchMuonTrack (const edm::Event& iEvent,
			  const reco::TrackRef theTrackRef)
{
  if ( theTrackRef.isNull() ) return false;

  edm::Handle< vector<pat::Muon> > thePATMuonHandle;
  iEvent.getByLabel(MuonLabel_, thePATMuonHandle);
  
  reco::TrackRef muTrackRef;
  for (vector<pat::Muon>::const_iterator iMuon = thePATMuonHandle->begin();
       iMuon != thePATMuonHandle->end(); iMuon++){

    muTrackRef = iMuon->innerTrack();
    if ( muTrackRef.isNull() ) continue;

    if (muTrackRef == theTrackRef) return true;
  }
  
  return false;
}

void
BToPiMuMu::buildBuToPiMuMu(const edm::Event& iEvent)
{
  // init variables
  edm::Handle< vector<pat::Muon> > patMuonHandle;
  iEvent.getByLabel(MuonLabel_, patMuonHandle);
  if( patMuonHandle->size() < 2 ) return ;

  edm::Handle< vector<pat::GenericParticle> >thePATTrackHandle;
  iEvent.getByLabel(TrackLabel_, thePATTrackHandle);

  bool passed;
  double DCAmumBS, DCAmumBSErr, DCAmupBS, DCAmupBSErr;
  double mumutrk_R, mumutrk_Z, DCAmumu;
  reco::TransientTrack refitMupTT, refitMumTT;
  double mu_mu_vtx_cl, mu_mu_pt, mu_mu_mass, mu_mu_mass_err;
  double MuMuLSBS, MuMuLSBSErr;
  double MuMuCosAlphaBS, MuMuCosAlphaBSErr;
  
  double pion_trk_pt, b_vtx_chisq, b_vtx_cl, b_mass;
  
  double DCAPionTrkBS, DCAPionTrkBSErr;                      
  RefCountedKinematicTree vertexFitTree;

  // ---------------------------------
  // loop 1: mu-
  // ---------------------------------
  for (vector<pat::Muon>::const_iterator iMuonM = patMuonHandle->begin();
       iMuonM != patMuonHandle->end(); iMuonM++){
    
    reco::TrackRef muTrackm = iMuonM->innerTrack();
    if ( muTrackm.isNull() ) continue;

    histos[h_mupt]->Fill(muTrackm->pt());
    histos[h_mueta]->Fill(muTrackm->eta());

    if ( (muTrackm->charge() != -1) ||
         (muTrackm->pt() < MuonMinPt_) ||
         (fabs(muTrackm->eta()) > MuonMaxEta_)) continue;
    
    // check mu- DCA to beam spot
    const reco::TransientTrack muTrackmTT(muTrackm, &(*bFieldHandle_));
    passed = hasGoodMuonDcaBs(muTrackmTT, DCAmumBS, DCAmumBSErr) ;
    histos[h_mumdcabs]->Fill(DCAmumBS);
    if ( ! passed ) continue;

    // ---------------------------------
    // loop 2: mu+
    // ---------------------------------
    for (vector<pat::Muon>::const_iterator iMuonP = patMuonHandle->begin();
         iMuonP != patMuonHandle->end(); iMuonP++){

      reco::TrackRef muTrackp = iMuonP->innerTrack();
      if ( muTrackp.isNull() ||
	   (muTrackp->charge() != 1) ||
	   (muTrackp->pt() < MuonMinPt_) ||
	   (fabs(muTrackp->eta()) > MuonMaxEta_)) continue;
      
      // check mu+ DCA to beam spot
      const reco::TransientTrack muTrackpTT(muTrackp, &(*bFieldHandle_));
      passed = hasGoodMuonDcaBs(muTrackpTT, DCAmupBS, DCAmupBSErr);
      if ( ! passed ) continue;
      
      if (! hasGoodClosestApproachTracks(muTrackpTT, muTrackmTT, mumutrk_R, mumutrk_Z, DCAmumu ) )  continue;  

      histos[h_mumutrkr]->Fill(mumutrk_R);
      histos[h_mumutrkz]->Fill(mumutrk_Z);
      histos[h_mumudca]->Fill(DCAmumu);


      if ( mumutrk_R > TrkMaxR_ ||
           mumutrk_Z > TrkMaxZ_ ||
           DCAmumu > MuMuMaxDca_ )   continue;     

      // check dimuon vertex
      passed = hasGoodMuMuVertex(muTrackpTT, muTrackmTT, refitMupTT, refitMumTT,
                                 mu_mu_vtx_cl, mu_mu_pt,
                                 mu_mu_mass, mu_mu_mass_err,
                                 MuMuLSBS, MuMuLSBSErr,
                                 MuMuCosAlphaBS, MuMuCosAlphaBSErr);
      
      histos[h_mumuvtxcl]->Fill(mu_mu_vtx_cl);
      histos[h_mumupt]->Fill(mu_mu_pt);
      histos[h_mumumass]->Fill(mu_mu_mass);
      histos[h_mumulxybs]->Fill(MuMuLSBS/MuMuLSBSErr);
      histos[h_mumucosalphabs]->Fill(MuMuCosAlphaBS);
      if ( !passed) continue;
 
      // ---------------------------------
      // loop 3: pion track
      // ---------------------------------
      for ( vector<pat::GenericParticle>::const_iterator iTrack
	      = thePATTrackHandle->begin();
	    iTrack != thePATTrackHandle->end(); ++iTrack ) {

	passed = hasGoodPionTrack(iEvent, *iTrack, pion_trk_pt);
	histos[h_trkpt]->Fill(pion_trk_pt);
	if (!passed) continue;
        

	// compute track DCA to beam spot
	reco::TrackRef pionTrack = iTrack->track();
	const reco::TransientTrack theTrackTT(pionTrack, &(*bFieldHandle_));

	passed = hasGoodTrackDcaBs(theTrackTT, DCAPionTrkBS, DCAPionTrkBSErr);         
	histos[h_trkdcasigbs]->Fill(DCAPionTrkBS/DCAPionTrkBSErr);
	if (!passed) continue;


	passed = hasGoodBuVertex(muTrackm, muTrackp, 
                                 pionTrack, b_vtx_chisq, b_vtx_cl,
                                 vertexFitTree);
	histos[h_bvtxchisq]->Fill(b_vtx_chisq);
	histos[h_bvtxcl]->Fill(b_vtx_cl);
	if (!passed) continue;
        
	passed = hasGoodBuMass(vertexFitTree, b_mass);
	histos[h_bmass]->Fill(b_mass);
	if (!passed) continue;

	nb++;

	// save the tree variables
	saveDimuVariables(DCAmumBS, DCAmumBSErr, DCAmupBS, DCAmupBSErr,
			  mumutrk_R, mumutrk_Z, DCAmumu, mu_mu_vtx_cl,
			  MuMuLSBS, MuMuLSBSErr,
			  MuMuCosAlphaBS, MuMuCosAlphaBSErr,
			  mu_mu_mass, mu_mu_mass_err);
        
	saveSoftMuonVariables(*iMuonM, *iMuonP, muTrackm, muTrackp);

	trkpt->push_back(pion_trk_pt);
	trkdcabs->push_back(DCAPionTrkBS);
	trkdcabserr->push_back(DCAPionTrkBSErr);

	bchg->push_back(iTrack->charge());
	bvtxcl->push_back(b_vtx_cl);
        

	saveBuToPiMuMu(vertexFitTree);
	saveBuVertex(vertexFitTree);
	saveBuCosAlpha(vertexFitTree);
	saveBuCosAlpha2d(vertexFitTree);       
	saveBuLsig(vertexFitTree);
	saveBuCtau(vertexFitTree);

      } // close kaon track loop
    } // close mu+ loop
  } // close mu- loop

  //  cout << "\n@@@ nb : " << nb << endl;  /* added cout statement */

  if ( nb > 0) 
    edm::LogInfo("myBu") << "Found " << nb << " Bu -> K+ mu mu.";    

}


void
BToPiMuMu::computeLS (double Vx, double Vy, double Vz,
		     double Wx, double Wy, double Wz,
		     double VxErr2, double VyErr2, double VzErr2,
		     double VxyCov, double VxzCov, double VyzCov,
		     double WxErr2, double WyErr2, double WzErr2,
		     double WxyCov, double WxzCov, double WyzCov,
		     double* deltaD, double* deltaDErr)
{
  *deltaD = sqrt((Vx-Wx) * (Vx-Wx) + (Vy-Wy) * (Vy-Wy) + (Vz-Wz) * (Vz-Wz));
  if (*deltaD > 0.)
    *deltaDErr = sqrt((Vx-Wx) * (Vx-Wx) * VxErr2 +
		      (Vy-Wy) * (Vy-Wy) * VyErr2 +
		      (Vz-Wz) * (Vz-Wz) * VzErr2 +
		      (Vx-Wx) * (Vy-Wy) * 2.*VxyCov +
		      (Vx-Wx) * (Vz-Wz) * 2.*VxzCov +
		      (Vy-Wy) * (Vz-Wz) * 2.*VyzCov +
		      (Vx-Wx) * (Vx-Wx) * WxErr2 +
		      (Vy-Wy) * (Vy-Wy) * WyErr2 +
		      (Vz-Wz) * (Vz-Wz) * WzErr2 +
		      (Vx-Wx) * (Vy-Wy) * 2.*WxyCov +
		      (Vx-Wx) * (Vz-Wz) * 2.*WxzCov +
		      (Vy-Wy) * (Vz-Wz) * 2.*WyzCov) / *deltaD;
  else *deltaDErr = 0.;
}


void
BToPiMuMu::computeCosAlpha (double Vx, double Vy, double Vz,
			   double Wx, double Wy, double Wz,
			   double VxErr2, double VyErr2, double VzErr2,
			   double VxyCov, double VxzCov, double VyzCov,
			   double WxErr2, double WyErr2, double WzErr2,
			   double WxyCov, double WxzCov, double WyzCov,
			   double* cosAlpha, double* cosAlphaErr)
{
  double Vnorm = sqrt(Vx*Vx + Vy*Vy + Vz*Vz);
  double Wnorm = sqrt(Wx*Wx + Wy*Wy + Wz*Wz);
  double VdotW = Vx*Wx + Vy*Wy + Vz*Wz;
  if ((Vnorm > 0.) && (Wnorm > 0.)) {
    *cosAlpha = VdotW / (Vnorm * Wnorm);
    *cosAlphaErr = sqrt( (
			  (Vx*Wnorm - VdotW*Wx) * (Vx*Wnorm - VdotW*Wx) * WxErr2 +
			  (Vy*Wnorm - VdotW*Wy) * (Vy*Wnorm - VdotW*Wy) * WyErr2 +
			  (Vz*Wnorm - VdotW*Wz) * (Vz*Wnorm - VdotW*Wz) * WzErr2 +
			  (Vx*Wnorm - VdotW*Wx) * (Vy*Wnorm - VdotW*Wy) * 2.*WxyCov +
			  (Vx*Wnorm - VdotW*Wx) * (Vz*Wnorm - VdotW*Wz) * 2.*WxzCov +
			  (Vy*Wnorm - VdotW*Wy) * (Vz*Wnorm - VdotW*Wz) * 2.*WyzCov) /
			 (Wnorm*Wnorm*Wnorm*Wnorm) +
			 ((Wx*Vnorm - VdotW*Vx) * (Wx*Vnorm - VdotW*Vx) * VxErr2 +
			  (Wy*Vnorm - VdotW*Vy) * (Wy*Vnorm - VdotW*Vy) * VyErr2 +
			  (Wz*Vnorm - VdotW*Vz) * (Wz*Vnorm - VdotW*Vz) * VzErr2 +
			  (Wx*Vnorm - VdotW*Vx) * (Wy*Vnorm - VdotW*Vy) * 2.*VxyCov +
			  (Wx*Vnorm - VdotW*Vx) * (Wz*Vnorm - VdotW*Vz) * 2.*VxzCov +
			  (Wy*Vnorm - VdotW*Vy) * (Wz*Vnorm - VdotW*Vz) * 2.*VyzCov) /
			 (Vnorm*Vnorm*Vnorm*Vnorm) ) / (Wnorm*Vnorm);
  } else {
    *cosAlpha = 0.;
    *cosAlphaErr = 0.;
  }
}


void
BToPiMuMu::computeCosAlpha2d (double Vx, double Vy, double Vz, 
			     double Wx, double Wy, double Wz,
			     double VxErr2, double VyErr2, double VzErr2,
			     double VxyCov, double VxzCov, double VyzCov,
			     double WxErr2, double WyErr2, double WzErr2,
			     double WxyCov, double WxzCov, double WyzCov,
			     double* cosAlpha2D, double* cosAlpha2DErr)
{
  double Vnorm = sqrt(Vx*Vx + Vy*Vy + Vz*Vz);
  double Wnorm = sqrt(Wx*Wx + Wy*Wy + Wz*Wz);
  double VdotW = Vx*Wx + Vy*Wy + Vz*Wz;
  if ((Vnorm > 0.) && (Wnorm > 0.)) {
    *cosAlpha2D = VdotW / (Vnorm * Wnorm);
    *cosAlpha2DErr = sqrt( (
			    (Vx*Wnorm - VdotW*Wx) * (Vx*Wnorm - VdotW*Wx) * WxErr2 +
			    (Vy*Wnorm - VdotW*Wy) * (Vy*Wnorm - VdotW*Wy) * WyErr2 +
			    (Vz*Wnorm - VdotW*Wz) * (Vz*Wnorm - VdotW*Wz) * WzErr2 +
			    (Vx*Wnorm - VdotW*Wx) * (Vy*Wnorm - VdotW*Wy) * 2.*WxyCov +
			    (Vx*Wnorm - VdotW*Wx) * (Vz*Wnorm - VdotW*Wz) * 2.*WxzCov +
			    (Vy*Wnorm - VdotW*Wy) * (Vz*Wnorm - VdotW*Wz) * 2.*WyzCov) /
			   (Wnorm*Wnorm*Wnorm*Wnorm) +
			   ((Wx*Vnorm - VdotW*Vx) * (Wx*Vnorm - VdotW*Vx) * VxErr2 +
			    (Wy*Vnorm - VdotW*Vy) * (Wy*Vnorm - VdotW*Vy) * VyErr2 +
			    (Wz*Vnorm - VdotW*Vz) * (Wz*Vnorm - VdotW*Vz) * VzErr2 +
			    (Wx*Vnorm - VdotW*Vx) * (Wy*Vnorm - VdotW*Vy) * 2.*VxyCov +
			    (Wx*Vnorm - VdotW*Vx) * (Wz*Vnorm - VdotW*Vz) * 2.*VxzCov +
			    (Wy*Vnorm - VdotW*Vy) * (Wz*Vnorm - VdotW*Vz) * 2.*VyzCov) /
			   (Vnorm*Vnorm*Vnorm*Vnorm) ) / (Wnorm*Vnorm);
  } else {
    *cosAlpha2D = 0.;
    *cosAlpha2DErr = 0.;
  }
}




bool
BToPiMuMu::hasGoodMuonDcaBs (const reco::TransientTrack muTrackTT,
			    double &muDcaBs, double &muDcaBsErr)
{
  TrajectoryStateClosestToPoint theDCAXBS =
    muTrackTT.trajectoryStateClosestToPoint(
     GlobalPoint(beamSpot_.position().x(),
		 beamSpot_.position().y(),beamSpot_.position().z()));
  
  if ( !theDCAXBS.isValid() ) return false;
  
  muDcaBs = theDCAXBS.perigeeParameters().transverseImpactParameter();
  muDcaBsErr = theDCAXBS.perigeeError().transverseImpactParameterError();
  if ( fabs(muDcaBs) > MuonMaxDcaBs_ ) return false;
  return true;
}

bool
BToPiMuMu::hasGoodTrackDcaBs (const reco::TransientTrack TrackTT,
			     double &DcaBs, double &DcaBsErr)
{
  TrajectoryStateClosestToPoint theDCAXBS =
    TrackTT.trajectoryStateClosestToPoint(
     GlobalPoint(beamSpot_.position().x(),
                 beamSpot_.position().y(),beamSpot_.position().z()));
  
  if ( !theDCAXBS.isValid() ) return false;
  
  DcaBs = theDCAXBS.perigeeParameters().transverseImpactParameter();
  DcaBsErr = theDCAXBS.perigeeError().transverseImpactParameterError();
  if ( fabs(DcaBs/DcaBsErr) < TrkMinDcaSigBs_ ) return false;    
  return true;
}

bool
BToPiMuMu::hasGoodClosestApproachTracks (const reco::TransientTrack muTrackpTT,
					const reco::TransientTrack muTrackmTT,
					double & trk_R,
					double & trk_Z,
					double & DCAmumu)
{
  ClosestApproachInRPhi ClosestApp;
  ClosestApp.calculate(muTrackpTT.initialFreeState(),
		       muTrackmTT.initialFreeState());
  if (! ClosestApp.status()) return false;
  
  GlobalPoint XingPoint = ClosestApp.crossingPoint();
  
  trk_R = sqrt(XingPoint.x()*XingPoint.x() + XingPoint.y()*XingPoint.y());
  trk_Z = fabs(XingPoint.z());

  if ((sqrt(XingPoint.x()*XingPoint.x() + XingPoint.y()*XingPoint.y()) >
       TrkMaxR_) || (fabs(XingPoint.z()) > TrkMaxZ_)) return false;

  DCAmumu = ClosestApp.distance();

  if (DCAmumu > MuMuMaxDca_) return false;
 
  return true;
}

bool
BToPiMuMu::hasGoodMuMuVertex (const reco::TransientTrack muTrackpTT,
			     const reco::TransientTrack muTrackmTT,
			     reco::TransientTrack &refitMupTT,
			     reco::TransientTrack &refitMumTT,
			     double & mu_mu_vtx_cl, double & mu_mu_pt,
			     double & mu_mu_mass, double & mu_mu_mass_err,
			     double & MuMuLSBS, double & MuMuLSBSErr,
			     double & MuMuCosAlphaBS,
			     double & MuMuCosAlphaBSErr)
{
  KinematicParticleFactoryFromTransientTrack partFactory;
  KinematicParticleVertexFitter PartVtxFitter;
  vector<RefCountedKinematicParticle> muonParticles;
  double chi = 0.;
  double ndf = 0.;
  muonParticles.push_back(partFactory.particle(muTrackmTT,
					       MuonMass_,chi,ndf,MuonMassErr_));
  muonParticles.push_back(partFactory.particle(muTrackpTT,
					       MuonMass_,chi,ndf,MuonMassErr_));
  RefCountedKinematicTree mumuVertexFitTree = PartVtxFitter.fit(muonParticles);
  if ( !mumuVertexFitTree->isValid()) return false;
  mumuVertexFitTree->movePointerToTheTop();
  RefCountedKinematicParticle mumu_KP = mumuVertexFitTree->currentParticle();
  RefCountedKinematicVertex mumu_KV = mumuVertexFitTree->currentDecayVertex();
  if ( !mumu_KV->vertexIsValid()) return false;
  mu_mu_vtx_cl = TMath::Prob((double)mumu_KV->chiSquared(),
			     int(rint(mumu_KV->degreesOfFreedom())));
  if (mu_mu_vtx_cl < MuMuMinVtxCl_) return false;
  // extract the re-fitted tracks
  mumuVertexFitTree->movePointerToTheTop();
  mumuVertexFitTree->movePointerToTheFirstChild();
  RefCountedKinematicParticle refitMum = mumuVertexFitTree->currentParticle();
  refitMumTT = refitMum->refittedTransientTrack();
  mumuVertexFitTree->movePointerToTheNextChild();
  RefCountedKinematicParticle refitMup = mumuVertexFitTree->currentParticle();
  refitMupTT = refitMup->refittedTransientTrack();
  TLorentzVector mymum, mymup, mydimu;
  mymum.SetXYZM(refitMumTT.track().momentum().x(),
		refitMumTT.track().momentum().y(),
		refitMumTT.track().momentum().z(), MuonMass_);
  mymup.SetXYZM(refitMupTT.track().momentum().x(),
		refitMupTT.track().momentum().y(),
		refitMupTT.track().momentum().z(), MuonMass_);
  mydimu = mymum + mymup;
  mu_mu_pt = mydimu.Perp();
  mu_mu_mass = mumu_KP->currentState().mass();
  mu_mu_mass_err = sqrt(mumu_KP->currentState().kinematicParametersError().
			matrix()(6,6));
  if ((mu_mu_pt < MuMuMinPt_) || (mu_mu_mass < MuMuMinInvMass_) ||
      (mu_mu_mass > MuMuMaxInvMass_)) return false;
  // compute the distance between mumu vtx and beam spot
  computeLS (mumu_KV->position().x(),mumu_KV->position().y(),0.0,
	     beamSpot_.position().x(),beamSpot_.position().y(),0.0,
	     mumu_KV->error().cxx(),mumu_KV->error().cyy(),0.0,
	     mumu_KV->error().matrix()(0,1),0.0,0.0,
	     beamSpot_.covariance()(0,0),beamSpot_.covariance()(1,1),0.0,
	     beamSpot_.covariance()(0,1),0.0,0.0,
	     &MuMuLSBS,&MuMuLSBSErr);
  if (MuMuLSBS/MuMuLSBSErr < MuMuMinLxySigmaBs_) return false;
  computeCosAlpha(mumu_KP->currentState().globalMomentum().x(),
		  mumu_KP->currentState().globalMomentum().y(),
		  0.0,
		  mumu_KV->position().x() - beamSpot_.position().x(),
		  mumu_KV->position().y() - beamSpot_.position().y(),
		  0.0,
		  mumu_KP->currentState().kinematicParametersError().matrix()(3,3),
		  mumu_KP->currentState().kinematicParametersError().matrix()(4,4),
		  0.0,
		  mumu_KP->currentState().kinematicParametersError().matrix()(3,4),
		  0.0,
		  0.0,
		  mumu_KV->error().cxx() + beamSpot_.covariance()(0,0),
		  mumu_KV->error().cyy() + beamSpot_.covariance()(1,1),
		  0.0,
		  mumu_KV->error().matrix()(0,1) + beamSpot_.covariance()(0,1),
		  0.0,
		  0.0,
		  &MuMuCosAlphaBS,&MuMuCosAlphaBSErr);

  if (MuMuCosAlphaBS < MuMuMinCosAlphaBs_) return false;
  return true;
}


bool
BToPiMuMu::hasGoodBuMass(RefCountedKinematicTree vertexFitTree,
			double & b_mass)
{
  vertexFitTree->movePointerToTheTop();
  RefCountedKinematicParticle b_KP = vertexFitTree->currentParticle();
  b_mass = b_KP->currentState().mass();
  if ( b_mass < BMinMass_ || b_mass > BMaxMass_ ) return false;
  return true;
}


bool
BToPiMuMu::hasGoodBuVertex(const reco::TrackRef mu1Track,
			  const reco::TrackRef mu2Track,
			  const reco::TrackRef pionTrack,
			  double & b_vtx_chisq, double & b_vtx_cl,
			  RefCountedKinematicTree &vertexFitTree)
{

  // do vertex fit for Bu
  KinematicParticleFactoryFromTransientTrack pFactory;
  reco::TransientTrack mu1TT(mu1Track, &(*bFieldHandle_) );
  reco::TransientTrack mu2TT(mu2Track, &(*bFieldHandle_) );
  reco::TransientTrack pionTT(pionTrack, &(*bFieldHandle_) );

  float chi = 0.;
  float ndf = 0.;
  vector<RefCountedKinematicParticle> vFitMCParticles;
  vFitMCParticles.push_back(pFactory.particle(mu1TT,MuonMass_,
					      chi,ndf,MuonMassErr_));
  vFitMCParticles.push_back(pFactory.particle(mu2TT,MuonMass_,
					      chi,ndf,MuonMassErr_));
  vFitMCParticles.push_back(pFactory.particle(pionTT, PionMass_, chi,
					      ndf, PionMassErr_));

  KinematicParticleVertexFitter fitter;
  vertexFitTree = fitter.fit(vFitMCParticles);
  if (!vertexFitTree->isValid()) return false;

  vertexFitTree->movePointerToTheTop();
  RefCountedKinematicVertex bDecayVertexMC = vertexFitTree->currentDecayVertex();
  if ( !bDecayVertexMC->vertexIsValid()) return false;
 

  b_vtx_chisq = bDecayVertexMC->chiSquared();
  if ( bDecayVertexMC->chiSquared()<0
       || bDecayVertexMC->chiSquared()>1000 ) return false;

  RefCountedKinematicVertex b_KV = vertexFitTree->currentDecayVertex();
  b_vtx_cl = ChiSquaredProbability((double)(b_KV->chiSquared()),
				   (double)(b_KV->degreesOfFreedom()));

  if ( b_vtx_cl < BMinVtxCl_ ) return false;
  
  return true;
}

void
BToPiMuMu::saveBuToPiMuMu(RefCountedKinematicTree vertexFitTree){
  vertexFitTree->movePointerToTheTop(); // B+ or B-
  RefCountedKinematicParticle b_KP = vertexFitTree->currentParticle();
  bpx->push_back(b_KP->currentState().globalMomentum().x());
  bpxerr->push_back( sqrt( b_KP->currentState().kinematicParametersError().matrix()(3,3) ) );
  bpy->push_back(b_KP->currentState().globalMomentum().y());
  bpyerr->push_back( sqrt( b_KP->currentState().kinematicParametersError().matrix()(4,4) ) );
  bpz->push_back(b_KP->currentState().globalMomentum().z());
  bpzerr->push_back( sqrt( b_KP->currentState().kinematicParametersError().matrix()(5,5) ) );
  bmass->push_back(b_KP->currentState().mass());
  bmasserr->push_back( sqrt( b_KP->currentState().kinematicParametersError().matrix()(6,6) ) );
  vertexFitTree->movePointerToTheFirstChild(); // mu1
  RefCountedKinematicParticle mu1_KP = vertexFitTree->currentParticle();
  vertexFitTree->movePointerToTheNextChild(); // mu2
  RefCountedKinematicParticle mu2_KP = vertexFitTree->currentParticle();
  RefCountedKinematicParticle mup_KP, mum_KP ;
  if ( mu1_KP->currentState().particleCharge() > 0 ) mup_KP = mu1_KP;
  if ( mu1_KP->currentState().particleCharge() < 0 ) mum_KP = mu1_KP;
  if ( mu2_KP->currentState().particleCharge() > 0 ) mup_KP = mu2_KP;
  if ( mu2_KP->currentState().particleCharge() < 0 ) mum_KP = mu2_KP;
  muppx->push_back(mup_KP->currentState().globalMomentum().x());
  muppy->push_back(mup_KP->currentState().globalMomentum().y());
  muppz->push_back(mup_KP->currentState().globalMomentum().z());
  mumpx->push_back(mum_KP->currentState().globalMomentum().x());
  mumpy->push_back(mum_KP->currentState().globalMomentum().y());
  mumpz->push_back(mum_KP->currentState().globalMomentum().z());
  vertexFitTree->movePointerToTheNextChild(); // pion track
  RefCountedKinematicParticle pion_KP = vertexFitTree->currentParticle();
  trkchg->push_back(pion_KP->currentState().particleCharge());
  trkpx->push_back(pion_KP->currentState().globalMomentum().x());
  trkpy->push_back(pion_KP->currentState().globalMomentum().y());
  trkpz->push_back(pion_KP->currentState().globalMomentum().z());
}


void
BToPiMuMu::saveBuVertex(RefCountedKinematicTree vertexFitTree){
  vertexFitTree->movePointerToTheTop();
  RefCountedKinematicVertex b_KV = vertexFitTree->currentDecayVertex();
  bvtxx->push_back((*b_KV).position().x());
  bvtxxerr->push_back(sqrt( abs(b_KV->error().cxx()) ));
  bvtxy->push_back((*b_KV).position().y());
  bvtxyerr->push_back(sqrt( abs(b_KV->error().cyy()) ));
  bvtxz->push_back((*b_KV).position().z());
  bvtxzerr->push_back(sqrt( abs(b_KV->error().czz()) ));
}


void
BToPiMuMu::saveBuCosAlpha(RefCountedKinematicTree vertexFitTree)
{
  // alpha is the angle in the transverse plane between the B0 momentum
  // and the seperation between the B0 vertex and the beamspot
  vertexFitTree->movePointerToTheTop();
  RefCountedKinematicParticle b_KP = vertexFitTree->currentParticle();
  RefCountedKinematicVertex b_KV = vertexFitTree->currentDecayVertex();
  double cosAlphaBS, cosAlphaBSErr;
  computeCosAlpha(b_KP->currentState().globalMomentum().x(),
		  b_KP->currentState().globalMomentum().y(),
		  b_KP->currentState().globalMomentum().z(),
		  b_KV->position().x() - beamSpot_.position().x(),
		  b_KV->position().y() - beamSpot_.position().y(),
		  b_KV->position().z() - beamSpot_.position().z(),
		  b_KP->currentState().kinematicParametersError().matrix()(3,3),
		  b_KP->currentState().kinematicParametersError().matrix()(4,4),
		  b_KP->currentState().kinematicParametersError().matrix()(5,5),
		  b_KP->currentState().kinematicParametersError().matrix()(3,4),
		  b_KP->currentState().kinematicParametersError().matrix()(3,5),
		  b_KP->currentState().kinematicParametersError().matrix()(4,5),
		  b_KV->error().cxx() + beamSpot_.covariance()(0,0),
		  b_KV->error().cyy() + beamSpot_.covariance()(1,1),
		  b_KV->error().czz() + beamSpot_.covariance()(2,2),
		  b_KV->error().matrix()(0,1) + beamSpot_.covariance()(0,1),
		  b_KV->error().matrix()(0,2) + beamSpot_.covariance()(0,2),
		  b_KV->error().matrix()(1,2) + beamSpot_.covariance()(1,2),
		  &cosAlphaBS,&cosAlphaBSErr);
  bcosalphabs->push_back(cosAlphaBS);
  bcosalphabserr->push_back(cosAlphaBSErr);
}

void
BToPiMuMu::saveBuCosAlpha2d(RefCountedKinematicTree vertexFitTree) 
{
  vertexFitTree->movePointerToTheTop();
  RefCountedKinematicParticle b_KP = vertexFitTree->currentParticle();
  RefCountedKinematicVertex b_KV = vertexFitTree->currentDecayVertex();
  double cosAlphaBS2D, cosAlphaBS2DErr;
  computeCosAlpha2d(b_KP->currentState().globalMomentum().x(),
		    b_KP->currentState().globalMomentum().y(),0.0,
		    b_KV->position().x() - beamSpot_.position().x(),
		    b_KV->position().y() - beamSpot_.position().y(),0.0,
		    b_KP->currentState().kinematicParametersError().matrix()(3,3),
		    b_KP->currentState().kinematicParametersError().matrix()(4,4),0.0,
		    b_KP->currentState().kinematicParametersError().matrix()(3,4),0.0,0.0,
		    b_KV->error().cxx() + beamSpot_.covariance()(0,0),
		    b_KV->error().cyy() + beamSpot_.covariance()(1,1),0.0,
		    b_KV->error().matrix()(0,1) + beamSpot_.covariance()(0,1),0.0,0.0,
		    &cosAlphaBS2D,&cosAlphaBS2DErr);
  bcosalphabs2d->push_back(cosAlphaBS2D);
  bcosalphabs2derr->push_back(cosAlphaBS2DErr);
}


void
BToPiMuMu::saveBuLsig(RefCountedKinematicTree vertexFitTree)
{
  vertexFitTree->movePointerToTheTop();
  RefCountedKinematicVertex b_KV = vertexFitTree->currentDecayVertex();
  double LSBS, LSBSErr;
  computeLS (b_KV->position().x(), b_KV->position().y(), 0.0,
	     beamSpot_.position().x(), beamSpot_.position().y(), 0.0,
	     b_KV->error().cxx(), b_KV->error().cyy(), 0.0,
	     b_KV->error().matrix()(0,1), 0.0, 0.0,
	     beamSpot_.covariance()(0,0), beamSpot_.covariance()(1,1), 0.0,
	     beamSpot_.covariance()(0,1), 0.0, 0.0,
	     &LSBS,&LSBSErr);
  blsbs->push_back(LSBS);
  blsbserr->push_back(LSBSErr);
}

void
BToPiMuMu::computeCtau(RefCountedKinematicTree vertexFitTree,
		      double &bctau, double &bctauerr)
{
  //calculate ctau = (mB*(Bvtx-Pvtx)*pB)/(|pB|**2)
  vertexFitTree->movePointerToTheTop();
  RefCountedKinematicParticle b_KP = vertexFitTree->currentParticle();
  RefCountedKinematicVertex b_KV = vertexFitTree->currentDecayVertex();
  double betagamma = (b_KP->currentState().globalMomentum().mag()/BuMass_);
  // calculate ctau error. Momentum error is negligible compared to
  // the vertex errors, so don't worry about it
  GlobalPoint BVP = GlobalPoint( b_KV->position() );
  GlobalPoint PVP = GlobalPoint( primaryVertex_.position().x(),
				 primaryVertex_.position().y(),
				 primaryVertex_.position().z() );
  GlobalVector sep3D = BVP-PVP;
  GlobalVector pBV = b_KP->currentState().globalMomentum();
  bctau = (BuMass_* (sep3D.dot(pBV)))/(pBV.dot(pBV));
  GlobalError BVE = b_KV->error();
  GlobalError PVE = GlobalError( primaryVertex_.error() );
  VertexDistance3D theVertexDistance3D;
  Measurement1D TheMeasurement = theVertexDistance3D.distance(
     VertexState(BVP, BVE), VertexState(PVP, PVE) );
  double myError = TheMeasurement.error();
  // ctau is defined by the portion of the flight distance along
  // the compoenent of the B momementum, so only consider the error
  // of that component, too, which is accomplished by scaling by
  // ((VB-VP)(dot)PB)/|VB-VP|*|PB|
  double scale = abs( (sep3D.dot(pBV))/(sep3D.mag()*pBV.mag()) );
  bctauerr = (myError*scale)/betagamma;
}


void
BToPiMuMu::saveBuCtau(RefCountedKinematicTree vertexFitTree)
{
  double bctau_temp, bctauerr_temp;
  computeCtau(vertexFitTree, bctau_temp, bctauerr_temp);
  bctau->push_back(bctau_temp);
  bctauerr->push_back(bctauerr_temp);
}

void
BToPiMuMu::saveSoftMuonVariables(pat::Muon iMuonM, pat::Muon iMuonP,
				reco::TrackRef muTrackm, reco::TrackRef muTrackp)
{
  mumisgoodmuon->push_back(muon::isGoodMuon(iMuonM, muon::TMOneStationTight));
  mupisgoodmuon->push_back(muon::isGoodMuon(iMuonP, muon::TMOneStationTight));
  mumnpixhits->push_back(muTrackm->hitPattern().numberOfValidPixelHits());
  mupnpixhits->push_back(muTrackp->hitPattern().numberOfValidPixelHits());
  mumnpixlayers->push_back(muTrackm->hitPattern().pixelLayersWithMeasurement());
  mupnpixlayers->push_back(muTrackp->hitPattern().pixelLayersWithMeasurement());
  mumntrkhits->push_back(muTrackm->hitPattern().numberOfValidTrackerHits());
  mupntrkhits->push_back(muTrackp->hitPattern().numberOfValidTrackerHits());
  mumntrklayers->push_back(muTrackm->hitPattern().trackerLayersWithMeasurement());
  mupntrklayers->push_back(muTrackp->hitPattern().trackerLayersWithMeasurement());
  mumnormchi2->push_back(muTrackm->normalizedChi2());
  mupnormchi2->push_back(muTrackp->normalizedChi2());
  mumdxyvtx->push_back(muTrackm->dxy(primaryVertex_.position()));
  mupdxyvtx->push_back(muTrackp->dxy(primaryVertex_.position()));
  mumdzvtx->push_back(muTrackm->dz(primaryVertex_.position()));
  mupdzvtx->push_back(muTrackp->dz(primaryVertex_.position()));
  mumpt->push_back(muTrackm->pt());
  muppt->push_back(muTrackp->pt());
  mumeta->push_back(muTrackm->eta());
  mupeta->push_back(muTrackp->eta());
}

void
BToPiMuMu::saveDimuVariables(double DCAmumBS, double DCAmumBSErr,
			    double DCAmupBS, double DCAmupBSErr,
			    double mumutrk_R, double mumutrk_Z,
			    double DCAmumu, double mu_mu_vtx_cl,
			    double MuMuLSBS, double MuMuLSBSErr,
			    double MuMuCosAlphaBS, double MuMuCosAlphaBSErr,
			    double mu_mu_mass, double mu_mu_mass_err)
{
  mumdcabs->push_back(DCAmumBS);
  mumdcabserr->push_back(DCAmumBSErr);
  mupdcabs->push_back(DCAmupBS);
  mupdcabserr->push_back(DCAmupBSErr);
  mumutrkr->push_back(mumutrk_R);
  mumutrkz->push_back(mumutrk_Z);
  mumudca->push_back(DCAmumu);
  mumuvtxcl->push_back(mu_mu_vtx_cl);
  mumulsbs->push_back(MuMuLSBS);
  mumulsbserr->push_back(MuMuLSBSErr);
  mumucosalphabs->push_back(MuMuCosAlphaBS);
  mumucosalphabserr->push_back(MuMuCosAlphaBSErr);
  mumumass->push_back(mu_mu_mass);
  mumumasserr->push_back(mu_mu_mass_err);
}



//define this as a plug-in
DEFINE_FWK_MODULE(BToPiMuMu);
