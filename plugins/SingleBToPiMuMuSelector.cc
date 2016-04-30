//************************************************
// SELECTOR CODE TO PRODUCE SINGLECAND NTUPLES
//  @author: N.Sahoo <nsahoo@cern.ch>    
//  @date: 2015-12-05  12:34 pm          
//************************************************
#define SingleBToPiMuMuSelector_cxx

#include <iostream>
#include <sstream>
#include <map>
#include "SingleBToPiMuMuSelector.h"
#include <TH2.h>
#include <TStyle.h>
#include <TProof.h>
#include <TLorentzVector.h>

//********************
// Global Constants
//********************
const double MUON_MASS = 0.10565837;
const double PION_MASS = 0.13957018;


//**************************
// user defined variables
//**************************
TDatime t_begin_ , t_now_ ;
int n_processed_, n_selected_; 
int n_triggers0, n_triggers1;

TTree *tree_; 

//*********************************
// Branch variables for new tree
//*********************************
int    Nb             = 0;
double Mumumass       = 0;
double Mumumasserr    = 0;
double Trkpt          = 0;
double Trkdcasigbs    = 0;

double Bmass          = 0;
double Bpt            = 0;
double Beta           = 0;
double Bphi           = 0;
int    Bchg           = 0;
double Bvtxcl         = 0;
double Blxysig        = 0;
double Bcosalphabs    = 0;
double Bcosalphabs2d  = 0;
double Bctau          = 0;

double Q2             = 0;
double DimuPt         = 0;
double DimuEta        = 0;

int    Triggers       = 0;


//*********************************
// Branches for Gen level info
//*********************************
int     trueBu       = 999; 
int     genBChg      = 999;
double  genBPt       = 0;
double  genBEta      = 0;
double  genBPhi      = 0;
double  genBVtxX     = 0;
double  genBVtxY     = 0;
double  genBVtxZ     = 0;
double  genMupPt     = 0;
double  genMupEta    = 0;
double  genMupPhi    = 0;
double  genMumPt     = 0;
double  genMumEta    = 0;
double  genMumPhi    = 0;

double genDimuPt     = 0;
double genDimuEta    = 0;
double genDimuPhi    = 0;

int     genTkChg     = 999;
double  genTkPt      = 0;
double  genTkEta     = 0;
double  genTkPhi     = 0;

double  genQ2        = 0;




void ClearEvent()
{//{{{
  Nb             = 0;
  Mumumass       = 0;
  Mumumasserr    = 0;
  Trkpt          = 0;
  Trkdcasigbs    = 0;
    
  Bmass          = 0;
  Bpt            = 0;
  Beta           = 0;
  Bphi           = 0;
  Bchg           = 0;
  Bvtxcl         = 0;
  Blxysig        = 0;
  Bcosalphabs    = 0;
  Bcosalphabs2d  = 0;
  Bctau          = 0;

  Q2             = 0;
  DimuPt         = 0;
  DimuEta        = 0;

  Triggers       = 0;
  trueBu         = 999;

  //*****
  // mc
  //*****
  genBChg = 999;
  genBPt = 0;
  genBEta= 0;
  genBPhi= 0;
  genBVtxX = 0;
  genBVtxY = 0;
  genBVtxZ = 0;
  genMupPt = 0;
  genMupEta= 0;
  genMupPhi= 0;
  genMumPt = 0;
  genMumEta= 0;
  genMumPhi= 0;

  genDimuPt = 0;
  genDimuEta = 0;
  genDimuPhi = 0;

  genTkChg = 999;
  genTkPt = 0;
  genTkEta= 0;
  genTkPhi= 0;
  genQ2 = 0;

}//}}}

void str_replace(std::string& str, const std::string& oldStr, const std::string& newStr)
{//{{{
  size_t pos = 0;
  while((pos = str.find(oldStr, pos)) != std::string::npos)
    {
      str.replace(pos, oldStr.length(), newStr);
      pos += newStr.length();
    }
}//}}}

string get_option_value(string option, string name)
{//{{{
  vector<string> args;
  istringstream f(option);
  string s;    
  while (getline(f, s, ';')) {
    args.push_back(s);
  }
  
  string value; 
  for(vector<string>::iterator it = args.begin(); it != args.end(); ++it) {
    value = *it; 
    unsigned found = value.find(name);
    if (found == 0) {
      str_replace(value, name+"=", ""); 
      break; 
    }
  }
  return value; 
}//}}}


void SingleBToPiMuMuSelector::Begin(TTree * /*tree*/)
{//{{{

  t_begin_.Set(); 
  printf("\n ---------- Begin Job ---------- \n");
  t_begin_.Print();

  n_processed_ = 0;
  n_selected_ = 0;

  n_triggers0 = 0;
  n_triggers1 = 0;

}//}}}

void SingleBToPiMuMuSelector::SlaveBegin(TTree * /*tree*/)
{//{{{

  string option = GetOption();
  tree_ = new TTree("tree", "tree"); 
  tree_->Branch("Mumumass"      , &Mumumass      , "Mumumass/D");
  tree_->Branch("Mumumasserr"   , &Mumumasserr   , "Mumumasserr/D");
  tree_->Branch("Trkpt"         , &Trkpt         , "Trkpt/D");
  tree_->Branch("Trkdcasigbs"   , &Trkdcasigbs   , "Trkdcasigbs/D");

  tree_->Branch("Bmass"         , &Bmass         , "Bmass/D");
  tree_->Branch("Bpt"           , &Bpt           , "Bpt/D");
  tree_->Branch("Beta"          , &Beta          , "Beta/D");
  tree_->Branch("Bphi"          , &Bphi          , "Bphi/D");
  tree_->Branch("Bchg"          , &Bchg          , "Bchg/I");

  tree_->Branch("Bvtxcl"        , &Bvtxcl        , "Bvtxcl/D");
  tree_->Branch("Blxysig"       , &Blxysig       , "Blxysig/D");
  tree_->Branch("Bcosalphabs"   , &Bcosalphabs   , "Bcosalphabs/D");
  tree_->Branch("Bcosalphabs2d" , &Bcosalphabs2d , "Bcosalphabs2d/D");
  tree_->Branch("Bctau"         , &Bctau         , "Bctau/D");

  tree_->Branch("Q2"            , &Q2            , "Q2/D");

  tree_->Branch("DimuPt"        , &DimuPt        , "DimuPt/D");
  tree_->Branch("DimuEta"       , &DimuEta       , "DimuEta/D");
  tree_->Branch("Triggers"      , &Triggers      , "Triggers/I");
  tree_->Branch("trueBu"        , &trueBu        , "trueBu/I");

  string datatype = get_option_value(option, "datatype");
  std::map<string,int> maptype;
  maptype.insert(std::pair<string,int>("data",1));
  maptype.insert(std::pair<string,int>("mc.nogen",2));
  maptype.insert(std::pair<string,int>("mc.lite",3));
  maptype.insert(std::pair<string,int>("mc.hlt",998));
  maptype.insert(std::pair<string,int>("mc",999));
  switch (maptype[datatype]) {
  case 1:
    break;
  case 2:
    break;
  case 3:
    tree_->Branch("genBChg"      , &genBChg      , "genBChg/I");
    tree_->Branch("genBPhi"      , &genBPhi      , "genBPhi/D");
    tree_->Branch("genMupPt"     , &genMupPt     , "genMupPt/D");
    tree_->Branch("genMupEta"    , &genMupEta    , "genMupEta/D");
    tree_->Branch("genMupPhi"    , &genMupPhi    , "genMupPhi/D");
    tree_->Branch("genMumPt"     , &genMumPt     , "genMumPt/D");
    tree_->Branch("genMumEta"    , &genMumEta    , "genMumEta/D");
    tree_->Branch("genMumPhi"    , &genMumPhi    , "genMumPhi/D");

    tree_->Branch("genDimuPt"    , &genDimuPt    , "genDimuPt/D");
    tree_->Branch("genDimuEta"   , &genDimuEta   , "genDimuEta/D");
    tree_->Branch("genDimuPhi"   , &genDimuPhi   , "genDimuPhi/D");

    tree_->Branch("genQ2"        , &genQ2        , "genQ2/D");
    break;
  case 998:
    tree_->Branch("genBChg"      , &genBChg      , "genBChg/I");
    tree_->Branch("genBPt"       , &genBPt       , "genBPt/D");
    tree_->Branch("genBEta"      , &genBEta      , "genBEta/D");
    tree_->Branch("genBPhi"      , &genBPhi      , "genBPhi/D");
    tree_->Branch("genBVtxX"     , &genBVtxX     , "genBVtxX/D");
    tree_->Branch("genBVtxY"     , &genBVtxY     , "genBVtxY/D");
    tree_->Branch("genBVtxZ"     , &genBVtxZ     , "genBVtxZ/D");
    tree_->Branch("genMupPt"     , &genMupPt     , "genMupPt/D");
    tree_->Branch("genMupEta"    , &genMupEta    , "genMupEta/D");
    tree_->Branch("genMupPhi"    , &genMupPhi    , "genMupPhi/D");
    tree_->Branch("genMumPt"     , &genMumPt     , "genMumPt/D");
    tree_->Branch("genMumEta"    , &genMumEta    , "genMumEta/D");
    tree_->Branch("genMumPhi"    , &genMumPhi    , "genMumPhi/D");

    tree_->Branch("genDimuPt"    , &genDimuPt    , "genDimuPt/D");
    tree_->Branch("genDimuEta"   , &genDimuEta   , "genDimuEta/D");
    tree_->Branch("genDimuPhi"   , &genDimuPhi   , "genDimuPhi/D");

    tree_->Branch("genQ2"        , &genQ2        , "genQ2/D");
    break;
  case 999:
    tree_->Branch("genBChg"      , &genBChg      , "genBChg/I");
    tree_->Branch("genBPt"       , &genBPt       , "genBPt/D");
    tree_->Branch("genBEta"      , &genBEta      , "genBEta/D");
    tree_->Branch("genBPhi"      , &genBPhi      , "genBPhi/D");
    tree_->Branch("genBVtxX"     , &genBVtxX     , "genBVtxX/D");
    tree_->Branch("genBVtxY"     , &genBVtxY     , "genBVtxY/D");
    tree_->Branch("genBVtxZ"     , &genBVtxZ     , "genBVtxZ/D");
    tree_->Branch("genMupPt"     , &genMupPt     , "genMupPt/D");
    tree_->Branch("genMupEta"    , &genMupEta    , "genMupEta/D");
    tree_->Branch("genMupPhi"    , &genMupPhi    , "genMupPhi/D");
    tree_->Branch("genMumPt"     , &genMumPt     , "genMumPt/D");
    tree_->Branch("genMumEta"    , &genMumEta    , "genMumEta/D");
    tree_->Branch("genMumPhi"    , &genMumPhi    , "genMumPhi/D");

    tree_->Branch("genDimuPt"    , &genDimuPt    , "genDimuPt/D");
    tree_->Branch("genDimuEta"   , &genDimuEta   , "genDimuEta/D");
    tree_->Branch("genDimuPhi"   , &genDimuPhi   , "genDimuPhi/D");

    tree_->Branch("genTkChg"     , &genTkChg     , "genTkChg/I");
    tree_->Branch("genTkPt"      , &genTkPt      , "genTkPt/D");
    tree_->Branch("genTkEta"     , &genTkEta     , "genTkEta/D");
    tree_->Branch("genTkPhi"     , &genTkPhi     , "genTkPhi/D");
    tree_->Branch("genQ2"        , &genQ2        , "genQ2/D");
    break;
  default:
    printf("No compatible datatype found. Please check use following types...\n\t\t[");
    for (std::map<string,int>::iterator iType = maptype.begin(); iType != maptype.end(); iType++){
      if (iType->second != 0) printf("%s,",iType->first.c_str());
    }
    printf("]\n");
    break;
  }

    fOutput->AddAll(gDirectory->GetList()); 

}//}}}

Bool_t SingleBToPiMuMuSelector::Process(Long64_t entry)
{//{{{

  ClearEvent();

  string option = GetOption();
  string datatype = get_option_value(option, "datatype"); 
  string cut = get_option_value(option, "cut"); 

  GetEntry(entry); 
  n_processed_ += 1; 
  Nb = nb; 

  if (triggernames->size() == 0) n_triggers0++;
  if (triggernames->size() == 1) n_triggers1++;

  //if (datatype != "data") SaveGen();
  if (datatype != "data" || datatype != "mc.nogen") SaveGen();

  int i = SelectB(cut); 
  //  if ( i != -1 && (datatype == "data" || istruebu->at(i)) ) {
  if ( i != -1 && (datatype == "data" || datatype == "mc.nogen" || datatype == "mc.lite" || datatype == "mc.hlt" || datatype == "mc") ) {
    printf("Entry#%lld, candidate#%d is selected.\n",entry,i);
    n_selected_ += 1; 
    SaveEvent(i);     
  }

  tree_->Fill();   
  return kTRUE;

}//}}}

void SingleBToPiMuMuSelector::SlaveTerminate()
{//{{{

}//}}}

void SingleBToPiMuMuSelector::Terminate()
{//{{{

  string option = GetOption();
  TString outfile = get_option_value(option, "ofile"); 
  //printf("option=%s\n",option.c_str());
  //printf("outfile=%s",outfile.Data());
    
  TFile file(outfile.Data(), "recreate"); 
  fOutput->Write();

  t_now_.Set(); 
  printf(" \n ---------- End Job ---------- \n" ) ;
  t_now_.Print();  
  printf(" processed: %i \n selected: %i \n \
            duration: %i sec \n rate: %g evts/sec\n",
	 n_processed_, n_selected_, 
	 t_now_.Convert() - t_begin_.Convert(), 
	 float(n_processed_)/(t_now_.Convert()-t_begin_.Convert()) );

  printf("TRIGGER INFO \n");
  printf("Triggers0 : %i\n", n_triggers0);
  printf("Triggers1 : %i\n", n_triggers1);

}//}}}

int SingleBToPiMuMuSelector::SelectB(string cut)
{//{{{

  int best_idx = -1; 
  double best_bvtxcl = 0.0; 

  if (cut == "cut0") {
    for (int i = 0; i < nb; i++) {

      if ( ! HasGoodDimuon(i) ) continue; 

      if (bvtxcl->at(i) > best_bvtxcl) {
	best_bvtxcl = bvtxcl->at(i); 
	best_idx = i; 
      }
    }

  }else if (cut == "nocut") {
    for (int i = 0; i < nb; i++) {
      if (bvtxcl->at(i) > best_bvtxcl) {
	best_bvtxcl = bvtxcl->at(i); 
	best_idx = i; 
      }
    }

  }else if (cut == "genonly") {
    best_idx = -1;

  }else{
    printf("WARNING: Unknown cut, apply 'genonly' by default.\n");
    best_idx = -1;
  }

  return best_idx;

}//}}}


bool SingleBToPiMuMuSelector::HasGoodDimuon(int i)
{//{{{
 
  if ( // soft muon id (cf https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#New_Version_recommended)
      mumisgoodmuon->at(i)
      && mupisgoodmuon->at(i) 
      && mumntrklayers->at(i) > 5  // 2012 Data
      && mupntrklayers->at(i) > 5  // 2012 Data 
      // && mumntrkhits->at(i) > 10 
      // && mupntrkhits->at(i) > 10 
      && mumnpixlayers->at(i) > 0  // 1,0 (old,new)
      && mupnpixlayers->at(i) > 0  // 1,0 (old,new) 
      //&& mumnormchi2->at(i) < 1.8 
      //&& mupnormchi2->at(i) < 1.8 
      && mumtrkqual->at(i) == 1
      && muptrkqual->at(i) == 1
      && fabs(mumdxyvtx->at(i)) < 0.3  // 3,0.3 (old,new)
      && fabs(mupdxyvtx->at(i)) < 0.3  // 3,0.3 (old,new)
      && fabs(mumdzvtx->at(i)) < 20   // 30,20 (old,new) 
      && fabs(mupdzvtx->at(i)) < 20   // 30,20 (old,new)
       ) return true; 
  return false; 
}//}}}


void SingleBToPiMuMuSelector::SaveEvent(int i)
{//{{{
  TLorentzVector B_4vec, Mup_4vec, Mum_4vec, Tk_4vec ;
  B_4vec.SetXYZM(bpx->at(i),bpy->at(i),bpz->at(i),bmass->at(i));
  Mup_4vec.SetXYZM(muppx->at(i),muppy->at(i),muppz->at(i),MUON_MASS);
  Mum_4vec.SetXYZM(mumpx->at(i),mumpy->at(i),mumpz->at(i),MUON_MASS);
  Tk_4vec.SetXYZM(trkpx->at(i),trkpy->at(i),trkpz->at(i),PION_MASS);

  Bmass = bmass->at(i); 
  Bchg = bchg->at(i); 
  Bvtxcl = bvtxcl->at(i); 
  Blxysig = (blsbs->at(i)/blsbserr->at(i)); 
  Bcosalphabs = bcosalphabs->at(i); 
  Bcosalphabs2d = bcosalphabs2d->at(i);
  Bctau = bctau->at(i); 
  trueBu = istruebu->at(i);

  Bpt = B_4vec.Pt(); 
  Beta = B_4vec.Eta();
  Bphi = B_4vec.Phi();
  Mumumass = mumumass->at(i); 
  Mumumasserr = mumumasserr->at(i); 
  //Trkpt = trkpt->at(i); 
  Trkpt = Tk_4vec.Pt();
  Trkdcasigbs = fabs( trkdcabs->at(i)/trkdcabserr->at(i) ); 
  Q2 = pow(mumumass->at(i),2);

  DimuPt  = (Mup_4vec+Mum_4vec).Pt();
  DimuEta = (Mup_4vec+Mum_4vec).Eta();

  Triggers = triggernames->size();

}//}}}

void SingleBToPiMuMuSelector::SaveGen()
{//{{{
  TLorentzVector genB_4vec, genMup_4vec, genMum_4vec, genTk_4vec ;
  genB_4vec.SetXYZM(genbpx,genbpy,genbpz,5.279);
  genMup_4vec.SetXYZM(genmuppx,genmuppy,genmuppz,MUON_MASS);
  genMum_4vec.SetXYZM(genmumpx,genmumpy,genmumpz,MUON_MASS);
  genTk_4vec.SetXYZM(gentrkpx,gentrkpy,gentrkpz,PION_MASS);


  genBChg      = genbchg;
  genBPt       = genB_4vec.Pt();
  genBEta      = genB_4vec.Eta();
  genBPhi      = genB_4vec.Phi();
  genBVtxX     = 0;//Should be at PV?
  genBVtxY     = 0;
  genBVtxZ     = 0;
  genMupPt     = genMup_4vec.Pt();
  genMupEta    = genMup_4vec.Eta();
  genMupPhi    = genMup_4vec.Phi();
  genMumPt     = genMum_4vec.Pt();
  genMumEta    = genMum_4vec.Eta();
  genMumPhi    = genMum_4vec.Phi();

  genDimuPt    = (genMup_4vec+genMum_4vec).Pt();
  genDimuEta   = (genMup_4vec+genMum_4vec).Eta();
  genDimuPhi   = (genMup_4vec+genMum_4vec).Phi();

  genTkChg     = gentrkchg;
  genTkPt      = genTk_4vec.Pt();
  genTkEta     = genTk_4vec.Eta();
  genTkPhi     = genTk_4vec.Phi();
  genQ2        = (genMup_4vec+genMum_4vec).Mag2();
    
 
}//}}}


#ifndef __CINT__ 
#include <algorithm>

char* get_option(char ** begin, char ** end, const std::string & option)
{//{{{
  char ** itr = std::find(begin, end, option);
  if (itr != end && ++itr != end)  return *itr;
  return 0;
}//}}}

bool option_exists(char** begin, char** end, const std::string& option)
{//{{{
  return std::find(begin, end, option) != end;
}//}}}

void print_usage()
{//{{{
  cerr << "Usage: SingleBuToKstarMuMuSelector datatype cut infile outfile [-n] [-s] [-j] [-h]\n"
       << "  datatype: data, mc, mc.nogen, mc.lite, mc.hlt\n"
       << "  cut     : cut0, nocut, genonly.\n"
       << "Options: \n" 
       << "  -h \t\tPrint this info.\n"
       << "  -n \t\tNumber of entries.\n" 
       << "  -s \t\tStarting run number.\n"
       << "  -j \t\tNumber of workers.\n" 
       << endl; 
}//}}}

int main(int argc, char** argv) {
  if ( (argc < 3) or option_exists(argv, argv+argc, "-h") ){
    print_usage() ;  
    return -1; 
  }

  TString datatype = argv[1]; 
  TString cut = argv[2]; 
  TString infile = argv[3]; 
  TString outfile = argv[4]; 

  Printf("datatype: '%s'", datatype.Data());
  Printf("cut: '%s'", cut.Data());
  Printf("input file: '%s'", infile.Data());
  Printf("output file: '%s'", outfile.Data());

  TChain *ch = new TChain("tree"); 
  ch->Add(infile.Data()); 

  char *j = get_option(argv, argv+argc, "-j");
  if (j) {
    TProof::Open(Form("workers=%s", j));
    ch->SetProof(); 
  }

  Long64_t nentries = 1000000000; 
  char * n = get_option(argv, argv+argc, "-n");  
  if (n){
    nentries = atoi(n);
  }
    
  int     iStart = 0;
  char *s = get_option(argv, argv+argc, "-s");
  if (s) {
    iStart = atoi(s);
    if (iStart > ch->GetEntries()){
      printf("ERROR: Number of entries is %lld.\n",ch->GetEntries());
      return -1;
    }
  }

  TString option; 
  option.Form("datatype=%s;cut=%s;ofile=sel_%s_%s_%s_s%d.root", datatype.Data(), cut.Data(), outfile.Data(), datatype.Data(), cut.Data(), iStart);

    
  // It's not allowed to run with fat trees!
  if (datatype.Data() == "mc" && (!(s) || !(n))){
    printf("WARNING: You must specify #entries(-n) and start run(-s) for datatype '%s'.\n",datatype.Data());
    return -1;
  }
    
  ch->Process("SingleBToPiMuMuSelector.cc+", option, nentries, iStart); 

  gSystem->Exit(0);

  return 0 ;
}

#endif
