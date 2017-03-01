//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jul  2 12:00:13 2015 by ROOT version 6.02/05
// from TChain tupel/MuonTree/
//////////////////////////////////////////////////////////
//#include "Tupel/Tupel/interface/correction.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TChain.h>
#include <iomanip>
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TH3.h"
#include "TRandom.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "math.h"
#include <fstream>
#include <string>
#include <iostream>
#include <TStyle.h>
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF2.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TPostScript.h"
#include "TLine.h"
#include "TEllipse.h"
#include "TMath.h"
#include "TLatex.h"
#include <vector>
#include <TSystem.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>
#include <algorithm>
#include "vector"
using namespace std;
//TTBarSolver ttsolver;

   UInt_t          Run=0;
   ULong64_t       LumiSection=0;
   ULong64_t       Event=0;
   Bool_t          RecoSuccess=0;
   Float_t         MassTT=0;
   TBranch        *b_Run;   //!
   TBranch        *b_LumiSection;   //!
   TBranch        *b_Event;   //!
   TBranch        *b_RecoSuccess;   //!
   TBranch        *b_MassTT;   //!

  void branchAdd(TTree *tree){
   tree->SetBranchAddress("Run", &Run, &b_Run);
   tree->SetBranchAddress("LumiSection", &LumiSection, &b_LumiSection);
   tree->SetBranchAddress("Event", &Event, &b_Event);
   tree->SetBranchAddress("RecoSuccess", &RecoSuccess, &b_RecoSuccess);
   tree->SetBranchAddress("MassTT", &MassTT, &b_MassTT);
  }

