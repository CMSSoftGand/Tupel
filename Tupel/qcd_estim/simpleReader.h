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
#include "Tupel/Tupel/plugins/NeutrinoSolver.cc"
#include "Tupel/Tupel/interface/NeutrinoSolver.h"
#include "Tupel/Tupel/interface/BtagUncertaintyComputer.h"
#include "Tupel/Tupel/plugins/BtagUncertaintyComputer.cc"
#include "Tupel/Tupel/plugins/TTBarSolver.cc"
#include "Tupel/Tupel/interface/BTagCalibration.h"
#include "Tupel/Tupel/interface/TTBarSolver.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
//#include "../interface/standalone_LumiReWeighting.h"
#include "Tupel/Tupel/interface/Permutation.h"
#include "Tupel/Tupel/plugins/Permutation.cc"
#include "Tupel/Tupel/interface/Hypotheses.h"
#include "Tupel/Tupel/plugins/Hypotheses.cc"
#include "vector"
using namespace std;
//TTBarSolver ttsolver;

   // 
   // Declaration of leaf types
   vector<double>  *Uncorec_METPt;
   vector<double>  *Uncorec_METPhi;
   vector<double>  *METRawPt=0;
   vector<double>  *METPt=0;
   vector<double>  *METPx=0;
   vector<double>  *METPy=0;
   vector<double>  *METPz=0;
   vector<double>  *METE=0;
   vector<double>  *METsigx2=0;
   vector<double>  *METsigxy=0;
   vector<double>  *METsigy2=0;
   vector<double>  *METsig=0;
   UInt_t          event=0;
   Int_t           realdata=0;
   UInt_t          run=0;
   UInt_t          lumi_=0;
   Int_t           bxnumber=0;
   Double_t        EvtInfo_NumVtx=0;
   Double_t        first_PV=0;
   Double_t        PU_npT=0;
   Double_t        PU_npIT=0;
   vector<double>  *St03Pt=0;
   vector<double>  *St03Eta=0;
   vector<double>  *St03Phi=0;
   vector<double>  *St03E=0;
   vector<double>  *St03M=0;
   vector<double>  *St03Id=0;
   vector<double>  *St03Status=0;
   vector<double>  *St03MotherId=0;
   vector<double>  *St03NumberMom=0;
   vector<bool>    *matchGjet=0;
   vector<double>  *GjDoughterPt_=0;
   vector<double>  *GjDoughterEta_=0;
   vector<double>  *GjDoughterPhi_=0;
   vector<double>  *GjDoughterE_=0;
   vector<double>  *MGjPt=0;
   vector<double>  *MGjeta=0;
   vector<double>  *MGjphi=0;
   vector<double>  *MGjE=0;
   Double_t        HLT_IsoMu24=0;
   Double_t        HLT_IsoMu22=0;
   Double_t        HLT_IsoTkMu24=0;
   Double_t        HLT_IsoTkMu22=0;
   Double_t        HLT_Mu17_Mu8=0;
   Double_t        HLT_Mu17_TkMu8=0;
   Double_t        HLT_Elec17_Elec8=0;
   Double_t        HLT_IsoMu24_eta2p1=0;
   Double_t        HLT_Ele32_eta2p1_WPTight_Gsf=0;
   Double_t        Flag_HBHENoiseFilter=0;
   Double_t        Flag_HBHENoiseIsoFilter=0;
   Double_t        Flag_globalTightHalo2016Filter=0;
   Double_t        Flag_EcalDeadCellTriggerPrimitiveFilter=0;
   Double_t        Flag_goodVertices=0;
   Double_t        Flag_eeBadScFilter=0;
   Double_t        Flag_muonBadTrackFilter=0;
   Double_t        Flag_chargedHadronTrackResolutionFilter=0;

   vector<double>  *patMuonPt_=0;
   vector<double>  *patMuonEta_=0;
   vector<double>  *patMuonPhi_=0;
   vector<double>  *patMuonVtxZ_=0;
   vector<double>  *patMuonEn_=0;
   vector<double>  *patMuonCharge_=0;
   vector<double>  *patMuonDxy_=0;
   vector<double>  *patMuonCombId_=0;
   vector<double>  *patMuonLooseId_=0;
   vector<double>  *patMuonMediumId_=0;
   vector<double>  *patMuonTightId_=0;
   vector<double>  *patMuonDetIsoRho_=0;
   vector<double>  *patMuonPfIsoDbeta_=0;
   vector<double>  *patMuonM_=0;
   vector<double>  *patMuonPx_=0;
   vector<double>  *patMuonPy_=0;
   vector<double>  *patMuonPz_=0;
   vector<double>  *patMuonGlobalType_=0;
   vector<double>  *patMuonTrackerType_=0;
   vector<double>  *patMuonPFType_=0;
   vector<double>  *patMuonIsoSumPt_=0;
   vector<double>  *patMuonIsoRelative_=0;
   vector<double>  *patMuonIsoCalComb_=0;
   vector<double>  *patMuonIsoDY_=0;
   vector<double>  *patMuonChi2Ndoff_=0;
   vector<double>  *patMuonNhits_=0;
   vector<double>  *patMuonNMatches_=0;
   vector<double>  *patMuonDz_=0;
   vector<double>  *patMuonPhits_=0;
   vector<double>  *patMuonTkLayers_=0;
   vector<double>  *patMuon_PF_IsoSumChargedHadronPt_=0;
   vector<double>  *patMuon_PF_IsoSumNeutralHadronEt_=0;
   vector<double>  *patMuon_PF_IsoDY_=0;
   vector<unsigned>  *patElecIdveto_=0;
   vector<unsigned>  *patElecIdloose_=0;
   vector<unsigned>  *patElecIdmedium_=0;
   vector<unsigned>  *patElecIdtight_=0;
   vector<unsigned>  *patElecIdnonTrig80_=0;
   vector<unsigned>  *patElecIdnonTrig90_=0;
   vector<unsigned>  *patElecIdTrig80_=0;
   vector<unsigned>  *patElecIdTrig90_=0;
   vector<double>  *patElecdEtaInSeed_=0;
   vector<double>  *patElecdPhiIn_=0;
   vector<double>  *patElechOverE_=0;
   vector<double>  *patElecsigmaIetaIeta_=0;
   vector<double>  *patElecfull5x5_sigmaIetaIeta_=0;
   vector<double>  *patElecooEmooP_=0;
   vector<double>  *patElecd0_=0;
   vector<double>  *patElecdz_=0;
   vector<int>     *patElecexpectedMissingInnerHits_=0;
   vector<int>     *patElecpassConversionVeto_=0;
   vector<double>  *patElecDz_=0;
   vector<double>  *patElecPt_=0;
   vector<double>  *patElecEta_=0;
   vector<double>  *patElecScEta_=0;
   vector<double>  *patElecPhi_=0;
   vector<double>  *patElecEnergy_=0;
   vector<double>  *patElecCharge_=0;
   vector<double>  *patElecchIso03_=0;
   vector<double>  *patElecnhIso03_=0;
   vector<double>  *patElecphIso03_=0;
   vector<double>  *patElecpuChIso03_=0;
   vector<double>  *patElecPfIso_=0;
   vector<double>  *patElecPfIsodb_=0;
   vector<double>  *patElecPfIsoRho_=0;
   Double_t        rhoPrime=0;
   vector<double>  *charged_=0;
   Double_t        AEff=0;
   vector<double>  *patElec_mva_presel_=0;
   vector<double>  *patJetPfAk04En_=0;
   vector<double>  *patJetPfAk04Pt_=0;
   vector<double>  *patJetPfAk04PtJERSmear=0;
   vector<double>  *patJetPfAk04PtJERSmearUp=0;
   vector<double>  *patJetPfAk04PtJERSmearDn=0;
   vector<double>  *patJetPfAk04EnJERSmear=0;
   vector<double>  *patJetPfAk04EnJERSmearUp=0;
   vector<double>  *patJetPfAk04EnJERSmearDn=0;
   vector<double>  *patJetPfAk04Eta_=0;
   vector<double>  *patJetPfAk04Phi_=0;
   vector<double>  *patJetPfAk04LooseId_=0;
   vector<double>  *patJetPfAk04Et_=0;
   vector<double>  *patJetPfAk04PartonFlavour_=0;
   vector<double>  *patJetPfAk04RawPt_=0;
   vector<double>  *patJetPfAk04RawEn_=0;
   vector<double>  *patJetPfAk04HadEHF_=0;
   vector<double>  *patJetPfAk04EmEHF_=0;
   vector<double>  *patJetPfAk04chf_=0;
   vector<double>  *patJetPfAk04nhf_=0;
   vector<double>  *patJetPfAk04cemf_=0;
   vector<double>  *patJetPfAk04nemf_=0;
   vector<double>  *patJetPfAk04cmult_=0;
   vector<double>  *patJetPfAk04nconst_=0;
   vector<double>  *patJetPfAk04BDiscdeepFudsg_;
   vector<double>  *patJetPfAk04BDiscdeepFb_;
   vector<double>  *patJetPfAk04BDiscdeepFbb_;
   vector<double>  *patJetPfAk04BDiscdeepFc_;
   vector<double>  *patJetPfAk04BDiscdeepFcc_;
   vector<double>  *patJetPfAk04BDiscCSVv2_=0;
   vector<double>  *patJetPfAk04BDiscpfCMVA_=0;
   vector<double>  *patJetPfAk04BDiscCSVV1_=0;
   vector<double>  *patJetPfAk04BDiscCSVSLV1_=0;
   vector<double>  *unc_=0;
   vector<int>     *patJetPfAk04DoughterId_=0;
   vector<double>  *patJetPfAk04DoughterPt_=0;
   vector<double>  *patJetPfAk04DoughterEta_=0;
   vector<double>  *patJetPfAk04DoughterPhi_=0;
   vector<double>  *patJetPfAk04DoughterE_=0;
   vector<double>  *patJetPfAk04PtUp_=0;
   vector<double>  *patJetPfAk04PtDn_=0;
   vector<double>  *id1_pdfInfo_=0;
   vector<double>  *id2_pdfInfo_=0;
   vector<double>  *x1_pdfInfo_=0;
   vector<double>  *x2_pdfInfo_=0;
   vector<double>  *scalePDF_pdfInfo_=0;
   Double_t        ptHat_=0;
   Double_t        mcWeight_=0;
   vector<double>  *mcWeights_=0;
   Double_t        nup=0;

   // List of branches
   TBranch        *b_Uncorec_METPt;   //!
   TBranch        *b_Uncorec_METPhi;   //!
   TBranch        *b_METRawPt;   //!
   TBranch        *b_METPt;   //!
   TBranch        *b_METPx;   //!
   TBranch        *b_METPy;   //!
   TBranch        *b_METPz;   //!
   TBranch        *b_METE;   //!
   TBranch        *b_METsigx2;   //!
   TBranch        *b_METsigxy;   //!
   TBranch        *b_METsigy2;   //!
   TBranch        *b_METsig;   //!
   TBranch        *b_event;   //!
   TBranch        *b_realdata;   //!
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_bxnumber;   //!
   TBranch        *b_EvtInfo_NumVtx;   //!
   TBranch        *b_first_PV;   //!
   TBranch        *b_PU_npT;   //!
   TBranch        *b_PU_npIT;   //!
   TBranch        *b_Dr01LepPt;   //!
   TBranch        *b_Dr01LepEta;   //!
   TBranch        *b_Dr01LepPhi;   //!
   TBranch        *b_Dr01LepE;   //!
   TBranch        *b_Dr01LepM;   //!
   TBranch        *b_Dr01LepId;   //!
   TBranch        *b_Dr01LepStatus;   //!
   TBranch        *b_Dr01LepMomId;   //!
   TBranch        *b_Bare01LepPt;   //!
   TBranch        *b_Bare01LepEta;   //!
   TBranch        *b_Bare01LepPhi;   //!
   TBranch        *b_Bare01LepE;   //!
   TBranch        *b_Bare01LepM;   //!
   TBranch        *b_Bare01LepId;   //!
   TBranch        *b_Bare01LepStatus;   //!
   TBranch        *b_Bare01LepMomId;   //!
   TBranch        *b_St03Pt;   //!
   TBranch        *b_St03Eta;   //!
   TBranch        *b_St03Phi;   //!
   TBranch        *b_St03E;   //!
   TBranch        *b_St03M;   //!
   TBranch        *b_St03Id;   //!
   TBranch        *b_St03Status;   //!
   TBranch        *b_St03MotherId;   //!
   TBranch        *b_St03NumberMom;   //!
   TBranch        *b_St01PhotonPt;   //!
   TBranch        *b_St01PhotonEta;   //!
   TBranch        *b_St01PhotonPhi;   //!
   TBranch        *b_St01PhotonE;   //!
   TBranch        *b_St01PhotonM;   //!
   TBranch        *b_St01PhotonId;   //!
   TBranch        *b_St01PhotonMomId;   //!
   TBranch        *b_St01PhotonNumberMom;   //!
   TBranch        *b_St01PhotonStatus;   //!
   TBranch        *b_matchGjet;   //!
   TBranch        *b_MGjPt;   //!
   TBranch        *b_MGjeta;   //!
   TBranch        *b_MGjphi;   //!
   TBranch        *b_MGjE;   //!
   TBranch        *b_HLT_IsoMu24;   //!
   TBranch        *b_HLT_IsoMu22;   //!
   TBranch        *b_HLT_IsoTkMu24;   //!
   TBranch        *b_HLT_IsoTkMu22;   //!
   TBranch        *b_HLT_Mu17_Mu8;   //!
   TBranch        *b_HLT_Mu17_TkMu8;   //!
   TBranch        *b_HLT_Elec17_Elec8;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1;
   TBranch        *b_HLT_Ele32_eta2p1_WPTight_Gsf;
   TBranch        *b_Flag_HBHENoiseFilter;
   TBranch        *b_Flag_HBHENoiseIsoFilter;
   TBranch        *b_Flag_globalTightHalo2016Filter;
   TBranch        *b_Flag_EcalDeadCellTriggerPrimitiveFilter;
   TBranch        *b_Flag_goodVertices;
   TBranch        *b_Flag_eeBadScFilter;
   TBranch        *b_Flag_chargedHadronTrackResolutionFilter;
   TBranch        *b_Flag_muonBadTrackFilter;
   TBranch        *b_patMuonPt_;   //!
   TBranch        *b_patMuonEta_;   //!
   TBranch        *b_patMuonPhi_;   //!
   TBranch        *b_patMuonVtxZ_;   //!
   TBranch        *b_patMuonEn_;   //!
   TBranch        *b_patMuonCharge_;   //!
   TBranch        *b_patMuonDxy_;   //!
   TBranch        *b_patMuonCombId_;   //!
   TBranch        *b_patMuonLooseId_;   //!
   TBranch        *b_patMuonMediumId_;   //!
   TBranch        *b_patMuonTightId_;   //!
   TBranch        *b_patMuonDetIsoRho_;   //!
   TBranch        *b_patMuonPfIsoDbeta_;   //!
   TBranch        *b_patMuonM_;   //!
   TBranch        *b_patMuonPx_;   //!
   TBranch        *b_patMuonPy_;   //!
   TBranch        *b_patMuonPz_;   //!
   TBranch        *b_patMuonGlobalType_;   //!
   TBranch        *b_patMuonTrackerType_;   //!
   TBranch        *b_patMuonPFType_;   //!
   TBranch        *b_patMuonIsoSumPt_;   //!
   TBranch        *b_patMuonIsoRelative_;   //!
   TBranch        *b_patMuonIsoCalComb_;   //!
   TBranch        *b_patMuonIsoDY_;   //!
   TBranch        *b_patMuonChi2Ndoff_;   //!
   TBranch        *b_patMuonNhits_;   //!
   TBranch        *b_patMuonNMatches_;   //!
   TBranch        *b_patMuonDz_;   //!
   TBranch        *b_patMuonPhits_;   //!
   TBranch        *b_patMuonTkLayers_;   //!
   TBranch        *b_patMuon_PF_IsoSumChargedHadronPt_;   //!
   TBranch        *b_patMuon_PF_IsoSumNeutralHadronEt_;   //!
   TBranch        *b_patMuon_PF_IsoDY_;   //!
   TBranch        *b_patElecdEtaInSeed_;   //!
   TBranch        *b_patElecIdveto_;   //!
   TBranch        *b_patElecIdloose_;   //!
   TBranch        *b_patElecIdmedium_;   //!
   TBranch        *b_patElecIdtight_;   //!
   TBranch        *b_patElecIdnonTrig80_;   //!
   TBranch        *b_patElecIdnonTrig90_;   //!
   TBranch        *b_patElecIdTrig80_;   //!
   TBranch        *b_patElecIdTrig90_;   //!
   TBranch        *b_patElecdPhiIn_;   //!
   TBranch        *b_patElechOverE_;   //!
   TBranch        *b_patElecsigmaIetaIeta_;   //!
   TBranch        *b_patElecfull5x5_sigmaIetaIeta_;   //!
   TBranch        *b_patElecooEmooP_;   //!
   TBranch        *b_patElecd0_;   //!
   TBranch        *b_patElecdz_;   //!
   TBranch        *b_patElecexpectedMissingInnerHits_;   //!
   TBranch        *b_patElecpassConversionVeto_;   //!
   TBranch        *b_patElecDz_;   //!
   TBranch        *b_patElecPt_;   //!
   TBranch        *b_patElecEta_;   //!
   TBranch        *b_patElecScEta_;   //!
   TBranch        *b_patElecPhi_;   //!
   TBranch        *b_patElecEnergy_;   //!
   TBranch        *b_patElecCharge_;   //!
   TBranch        *b_patElecchIso03_;   //!
   TBranch        *b_patElecnhIso03_;   //!
   TBranch        *b_patElecphIso03_;   //!
   TBranch        *b_patElecpuChIso03_;   //!
   TBranch        *b_patElecPfIso_;   //!
   TBranch        *b_patElecPfIsodb_;   //!
   TBranch        *b_patElecPfIsoRho_;   //!
   TBranch        *b_rhoPrime;   //!
   TBranch        *b_charged_;   //!
   TBranch        *b_AEff;   //!
   TBranch        *b_patElec_mva_presel_;   //!
   TBranch        *b_patJetPfAk04En_;   //!
   TBranch        *b_patJetPfAk04Pt_;   //!
   TBranch        *b_patJetPfAk04PtJERSmear;   //!
   TBranch        *b_patJetPfAk04PtJERSmearUp;   //!
   TBranch        *b_patJetPfAk04PtJERSmearDn;   //!
   TBranch        *b_patJetPfAk04EnJERSmear;   //!
   TBranch        *b_patJetPfAk04EnJERSmearUp;   //!
   TBranch        *b_patJetPfAk04EnJERSmearDn;   //!
   TBranch        *b_patJetPfAk04Eta_;   //!
   TBranch        *b_patJetPfAk04Phi_;   //!
   TBranch        *b_patJetPfAk04LooseId_;   //!
   TBranch        *b_patJetPfAk04Et_;   //!
   TBranch        *b_patJetPfAk04PartonFlavour_;   //!
   TBranch        *b_patJetPfAk04RawPt_;   //!
   TBranch        *b_patJetPfAk04RawEn_;   //!
   TBranch        *b_patJetPfAk04HadEHF_;   //!
   TBranch        *b_patJetPfAk04EmEHF_;   //!
   TBranch        *b_patJetPfAk04chf_;   //!
   TBranch        *b_patJetPfAk04nhf_;   //!
   TBranch        *b_patJetPfAk04cemf_;   //!
   TBranch        *b_patJetPfAk04nemf_;   //!
   TBranch        *b_patJetPfAk04cmult_;   //!
   TBranch        *b_patJetPfAk04nconst_;   //!
   TBranch        *b_patJetPfAk04BDiscdeepFudsg_;   //!
   TBranch        *b_patJetPfAk04BDiscdeepFb_;   //!
   TBranch        *b_patJetPfAk04BDiscdeepFbb_;   //!
   TBranch        *b_patJetPfAk04BDiscdeepFc_;   //!
   TBranch        *b_patJetPfAk04BDiscdeepFcc_;   //!
   TBranch        *b_patJetPfAk04BDiscCSVv2_;   //!
   TBranch        *b_patJetPfAk04BDiscpfCMVA_;   //!
   TBranch        *b_patJetPfAk04BDiscCSVV1_;   //!
   TBranch        *b_patJetPfAk04BDiscCSVSLV1_;   //!
   TBranch        *b_unc_;   //!
   TBranch        *b_patJetPfAk04DoughterId_;   //!
   TBranch        *b_patJetPfAk04DoughterPt_;   //!
   TBranch        *b_patJetPfAk04DoughterEta_;   //!
   TBranch        *b_patJetPfAk04DoughterPhi_;   //!
   TBranch        *b_patJetPfAk04DoughterE_;   //!
   TBranch        *b_patJetPfAk04PtUp_;   //!
   TBranch        *b_patJetPfAk04PtDn_;   //!
   TBranch        *b_caloJetPt_;   //!
   TBranch        *b_caloJetRawPt_;   //!
   TBranch        *b_caloJetEn_;   //!
   TBranch        *b_caloJetEta_;   //!
   TBranch        *b_caloJetPhi_;   //!
   TBranch        *b_caloJetHadEHF_;   //!
   TBranch        *b_caloJetEmEHF_;   //!
   TBranch        *b_caloJetEmFrac_;   //!
   TBranch        *b_caloJetn90_;   //!
   TBranch        *b_PhotonPt;   //!
   TBranch        *b_PhotonEta;   //!
   TBranch        *b_PhotonPhi;   //!
   TBranch        *b_PhotonIsoEcal;   //!
   TBranch        *b_PhotonIsoHcal;   //!
   TBranch        *b_PhotonPfIsoChargdH;   //!
   TBranch        *b_PhotonPfIsoNeutralH;   //!
   TBranch        *b_PhotonPfIsoPhoton;   //!
   TBranch        *b_PhotonPfIsoPuChargedH;   //!
   TBranch        *b_PhotonPfIsoEcalCluster;   //!
   TBranch        *b_PhotonPfIsoHcalCluster;   //!
   TBranch        *b_PhotonE3x3;   //!
   TBranch        *b_PhotonSigmaIetaIeta;   //!
   TBranch        *b_PhotonId;   //!
   TBranch        *b_id1_pdfInfo_;   //!
   TBranch        *b_id2_pdfInfo_;   //!
   TBranch        *b_x1_pdfInfo_;   //!
   TBranch        *b_x2_pdfInfo_;   //!
   TBranch        *b_scalePDF_pdfInfo_;   //!
   TBranch        *b_ptHat_;   //!
   TBranch        *b_mcWeight_;   //!
   TBranch        *b_mcWeights_;   //!
   TBranch        *b_nup;   //!

  void branchAdd(TTree *tree){
   tree->SetBranchAddress("Uncorec_METPt", &Uncorec_METPt, &b_Uncorec_METPt);
   tree->SetBranchAddress("Uncorec_METPhi", &Uncorec_METPhi, &b_Uncorec_METPhi);
   tree->SetBranchAddress("METPt", &METPt, &b_METPt);
   tree->SetBranchAddress("METPx", &METPx, &b_METPx);
   tree->SetBranchAddress("METPy", &METPy, &b_METPy);
   tree->SetBranchAddress("METPz", &METPz, &b_METPz);
   tree->SetBranchAddress("METE", &METE, &b_METE);
   tree->SetBranchAddress("METsigx2", &METsigx2, &b_METsigx2);
   tree->SetBranchAddress("METsigxy", &METsigxy, &b_METsigxy);
   tree->SetBranchAddress("METsigy2", &METsigy2, &b_METsigy2);
   tree->SetBranchAddress("METsig", &METsig, &b_METsig);
   tree->SetBranchAddress("event", &event, &b_event);
   tree->SetBranchAddress("realdata", &realdata, &b_realdata);
   tree->SetBranchAddress("run", &run, &b_run);
   tree->SetBranchAddress("lumi", &lumi_, &b_lumi);
   tree->SetBranchAddress("bxnumber", &bxnumber, &b_bxnumber);
   tree->SetBranchAddress("EvtInfo_NumVtx", &EvtInfo_NumVtx, &b_EvtInfo_NumVtx);
   tree->SetBranchAddress("first_PV", &first_PV, &b_first_PV);
   tree->SetBranchAddress("PU_npT", &PU_npT, &b_PU_npT);
   tree->SetBranchAddress("PU_npIT", &PU_npIT, &b_PU_npIT);
   tree->SetBranchAddress("St03Pt", &St03Pt, &b_St03Pt);
   tree->SetBranchAddress("St03Eta", &St03Eta, &b_St03Eta);
   tree->SetBranchAddress("St03Phi", &St03Phi, &b_St03Phi);
   tree->SetBranchAddress("St03E", &St03E, &b_St03E);
   tree->SetBranchAddress("St03M", &St03M, &b_St03M);
   tree->SetBranchAddress("St03Id", &St03Id, &b_St03Id);
   tree->SetBranchAddress("St03Status", &St03Status, &b_St03Status);
   tree->SetBranchAddress("St03MotherId", &St03MotherId, &b_St03MotherId);
   tree->SetBranchAddress("St03NumberMom", &St03NumberMom, &b_St03NumberMom);
   tree->SetBranchAddress("matchGjet", &matchGjet, &b_matchGjet);
   tree->SetBranchAddress("MGjPt", &MGjPt, &b_MGjPt);
   tree->SetBranchAddress("MGjeta", &MGjeta, &b_MGjeta);
   tree->SetBranchAddress("MGjphi", &MGjphi, &b_MGjphi);
   tree->SetBranchAddress("MGjE", &MGjE, &b_MGjE);
   tree->SetBranchAddress("HLT_IsoMu24", &HLT_IsoMu24, &b_HLT_IsoMu24);
   tree->SetBranchAddress("HLT_IsoMu22_v", &HLT_IsoMu22, &b_HLT_IsoMu22);
   tree->SetBranchAddress("HLT_IsoTkMu24", &HLT_IsoTkMu24, &b_HLT_IsoTkMu24);
   tree->SetBranchAddress("HLT_IsoTkMu22_v", &HLT_IsoTkMu22, &b_HLT_IsoTkMu22);
   tree->SetBranchAddress("HLT_Ele32_eta2p1_WPTight_Gsf", &HLT_Ele32_eta2p1_WPTight_Gsf, &b_HLT_Ele32_eta2p1_WPTight_Gsf);

   tree->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, &b_Flag_HBHENoiseFilter);
   tree->SetBranchAddress("Flag_HBHENoiseIsoFilter", &Flag_HBHENoiseIsoFilter, &b_Flag_HBHENoiseIsoFilter);
   tree->SetBranchAddress("Flag_globalTightHalo2016Filter", &Flag_globalTightHalo2016Filter, &b_Flag_globalTightHalo2016Filter);
   tree->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, &b_Flag_EcalDeadCellTriggerPrimitiveFilter);
   tree->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices, &b_Flag_goodVertices);
   tree->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter, &b_Flag_eeBadScFilter);
   tree->SetBranchAddress("Flag_chargedHadronTrackResolutionFilter", &Flag_chargedHadronTrackResolutionFilter, &b_Flag_chargedHadronTrackResolutionFilter);
   tree->SetBranchAddress("Flag_muonBadTrackFilter", &Flag_muonBadTrackFilter, &b_Flag_muonBadTrackFilter);
   tree->SetBranchAddress("patMuonPt_", &patMuonPt_, &b_patMuonPt_);
   tree->SetBranchAddress("patMuonEta_", &patMuonEta_, &b_patMuonEta_);
   tree->SetBranchAddress("patMuonPhi_", &patMuonPhi_, &b_patMuonPhi_);
   tree->SetBranchAddress("patMuonVtxZ_", &patMuonVtxZ_, &b_patMuonVtxZ_);
   tree->SetBranchAddress("patMuonEn_", &patMuonEn_, &b_patMuonEn_);
   tree->SetBranchAddress("patMuonCharge_", &patMuonCharge_, &b_patMuonCharge_);
   tree->SetBranchAddress("patMuonDxy_", &patMuonDxy_, &b_patMuonDxy_);
   tree->SetBranchAddress("patMuonCombId_", &patMuonCombId_, &b_patMuonCombId_);
   tree->SetBranchAddress("patMuonLooseId_", &patMuonLooseId_, &b_patMuonLooseId_);
   tree->SetBranchAddress("patMuonMediumId_", &patMuonMediumId_, &b_patMuonMediumId_);
   tree->SetBranchAddress("patMuonTightId_", &patMuonTightId_, &b_patMuonTightId_);
   tree->SetBranchAddress("patMuonDetIsoRho_", &patMuonDetIsoRho_, &b_patMuonDetIsoRho_);
   tree->SetBranchAddress("patMuonPfIsoDbeta_", &patMuonPfIsoDbeta_, &b_patMuonPfIsoDbeta_);
   tree->SetBranchAddress("patMuonM_", &patMuonM_, &b_patMuonM_);
   tree->SetBranchAddress("patMuonPx_", &patMuonPx_, &b_patMuonPx_);
   tree->SetBranchAddress("patMuonPy_", &patMuonPy_, &b_patMuonPy_);
   tree->SetBranchAddress("patMuonPz_", &patMuonPz_, &b_patMuonPz_);
   tree->SetBranchAddress("patMuonGlobalType_", &patMuonGlobalType_, &b_patMuonGlobalType_);
   tree->SetBranchAddress("patMuonTrackerType_", &patMuonTrackerType_, &b_patMuonTrackerType_);
   tree->SetBranchAddress("patMuonPFType_", &patMuonPFType_, &b_patMuonPFType_);
   tree->SetBranchAddress("patMuonIsoSumPt_", &patMuonIsoSumPt_, &b_patMuonIsoSumPt_);
   tree->SetBranchAddress("patMuonIsoRelative_", &patMuonIsoRelative_, &b_patMuonIsoRelative_);
   tree->SetBranchAddress("patMuonIsoCalComb_", &patMuonIsoCalComb_, &b_patMuonIsoCalComb_);
   tree->SetBranchAddress("patMuonIsoDY_", &patMuonIsoDY_, &b_patMuonIsoDY_);
   tree->SetBranchAddress("patMuonChi2Ndoff_", &patMuonChi2Ndoff_, &b_patMuonChi2Ndoff_);
   tree->SetBranchAddress("patMuonNhits_", &patMuonNhits_, &b_patMuonNhits_);
   tree->SetBranchAddress("patMuonNMatches_", &patMuonNMatches_, &b_patMuonNMatches_);
   tree->SetBranchAddress("patMuonDz_", &patMuonDz_, &b_patMuonDz_);
   tree->SetBranchAddress("patMuonPhits_", &patMuonPhits_, &b_patMuonPhits_);
   tree->SetBranchAddress("patMuonTkLayers_", &patMuonTkLayers_, &b_patMuonTkLayers_);
   tree->SetBranchAddress("patMuon_PF_IsoSumChargedHadronPt_", &patMuon_PF_IsoSumChargedHadronPt_, &b_patMuon_PF_IsoSumChargedHadronPt_);
   tree->SetBranchAddress("patMuon_PF_IsoSumNeutralHadronEt_", &patMuon_PF_IsoSumNeutralHadronEt_, &b_patMuon_PF_IsoSumNeutralHadronEt_);
   tree->SetBranchAddress("patMuon_PF_IsoDY_", &patMuon_PF_IsoDY_, &b_patMuon_PF_IsoDY_);
   tree->SetBranchAddress("patElecdEtaInSeed_", &patElecdEtaInSeed_, &b_patElecdEtaInSeed_);
   tree->SetBranchAddress("patElecIdveto_", &patElecIdveto_, &b_patElecIdveto_);
   tree->SetBranchAddress("patElecIdloose_", &patElecIdloose_, &b_patElecIdloose_);
   tree->SetBranchAddress("patElecIdmedium_", &patElecIdmedium_, &b_patElecIdmedium_);
   tree->SetBranchAddress("patElecIdtight_", &patElecIdtight_, &b_patElecIdtight_);
   tree->SetBranchAddress("patElecdPhiIn_", &patElecdPhiIn_, &b_patElecdPhiIn_);
   tree->SetBranchAddress("patElechOverE_", &patElechOverE_, &b_patElechOverE_);
   tree->SetBranchAddress("patElecsigmaIetaIeta_", &patElecsigmaIetaIeta_, &b_patElecsigmaIetaIeta_);
   tree->SetBranchAddress("patElecfull5x5_sigmaIetaIeta_", &patElecfull5x5_sigmaIetaIeta_, &b_patElecfull5x5_sigmaIetaIeta_);
   tree->SetBranchAddress("patElecooEmooP_", &patElecooEmooP_, &b_patElecooEmooP_);
   tree->SetBranchAddress("patElecd0_", &patElecd0_, &b_patElecd0_);
   tree->SetBranchAddress("patElecdz_", &patElecdz_, &b_patElecdz_);
   tree->SetBranchAddress("patElecexpectedMissingInnerHits_", &patElecexpectedMissingInnerHits_, &b_patElecexpectedMissingInnerHits_);
   tree->SetBranchAddress("patElecpassConversionVeto_", &patElecpassConversionVeto_, &b_patElecpassConversionVeto_);
   tree->SetBranchAddress("patElecDz_", &patElecDz_, &b_patElecDz_);
   tree->SetBranchAddress("patElecPt_", &patElecPt_, &b_patElecPt_);
   tree->SetBranchAddress("patElecEta_", &patElecEta_, &b_patElecEta_);
   tree->SetBranchAddress("patElecScEta_", &patElecScEta_, &b_patElecScEta_);
   tree->SetBranchAddress("patElecPhi_", &patElecPhi_, &b_patElecPhi_);
   tree->SetBranchAddress("patElecEnergy_", &patElecEnergy_, &b_patElecEnergy_);
   tree->SetBranchAddress("patElecCharge_", &patElecCharge_, &b_patElecCharge_);
   tree->SetBranchAddress("patElecchIso03_", &patElecchIso03_, &b_patElecchIso03_);
   tree->SetBranchAddress("patElecnhIso03_", &patElecnhIso03_, &b_patElecnhIso03_);
   tree->SetBranchAddress("patElecphIso03_", &patElecphIso03_, &b_patElecphIso03_);
   tree->SetBranchAddress("patElecpuChIso03_", &patElecpuChIso03_, &b_patElecpuChIso03_);
   tree->SetBranchAddress("patElecPfIso_", &patElecPfIso_, &b_patElecPfIso_);
   tree->SetBranchAddress("patElecPfIsodb_", &patElecPfIsodb_, &b_patElecPfIsodb_);
   tree->SetBranchAddress("patElecPfIsoRho_", &patElecPfIsoRho_, &b_patElecPfIsoRho_);
   tree->SetBranchAddress("rhoPrime", &rhoPrime, &b_rhoPrime);
   tree->SetBranchAddress("charged_", &charged_, &b_charged_);
   tree->SetBranchAddress("AEff", &AEff, &b_AEff);
   tree->SetBranchAddress("patElec_mva_presel_", &patElec_mva_presel_, &b_patElec_mva_presel_);
   tree->SetBranchAddress("patJetPfAk04En_", &patJetPfAk04En_, &b_patJetPfAk04En_);
   tree->SetBranchAddress("patJetPfAk04Pt_", &patJetPfAk04Pt_, &b_patJetPfAk04Pt_);
   tree->SetBranchAddress("patJetPfAk04PtJERSmear", &patJetPfAk04PtJERSmear, &b_patJetPfAk04PtJERSmear);
   tree->SetBranchAddress("patJetPfAk04PtJERSmearUp", &patJetPfAk04PtJERSmearUp, &b_patJetPfAk04PtJERSmearUp);
   tree->SetBranchAddress("patJetPfAk04PtJERSmearDn", &patJetPfAk04PtJERSmearDn, &b_patJetPfAk04PtJERSmearDn);
   tree->SetBranchAddress("patJetPfAk04EnJERSmear", &patJetPfAk04EnJERSmear, &b_patJetPfAk04EnJERSmear);
   tree->SetBranchAddress("patJetPfAk04EnJERSmearUp", &patJetPfAk04EnJERSmearUp, &b_patJetPfAk04EnJERSmearUp);
   tree->SetBranchAddress("patJetPfAk04EnJERSmearDn", &patJetPfAk04EnJERSmearDn, &b_patJetPfAk04EnJERSmearDn);
   tree->SetBranchAddress("patJetPfAk04Eta_", &patJetPfAk04Eta_, &b_patJetPfAk04Eta_);
   tree->SetBranchAddress("patJetPfAk04Phi_", &patJetPfAk04Phi_, &b_patJetPfAk04Phi_);
   tree->SetBranchAddress("patJetPfAk04LooseId_", &patJetPfAk04LooseId_, &b_patJetPfAk04LooseId_);
   tree->SetBranchAddress("patJetPfAk04Et_", &patJetPfAk04Et_, &b_patJetPfAk04Et_);
   tree->SetBranchAddress("patJetPfAk04PartonFlavour_", &patJetPfAk04PartonFlavour_, &b_patJetPfAk04PartonFlavour_);
   tree->SetBranchAddress("patJetPfAk04RawPt_", &patJetPfAk04RawPt_, &b_patJetPfAk04RawPt_);
   tree->SetBranchAddress("patJetPfAk04RawEn_", &patJetPfAk04RawEn_, &b_patJetPfAk04RawEn_);
   tree->SetBranchAddress("patJetPfAk04HadEHF_", &patJetPfAk04HadEHF_, &b_patJetPfAk04HadEHF_);
   tree->SetBranchAddress("patJetPfAk04EmEHF_", &patJetPfAk04EmEHF_, &b_patJetPfAk04EmEHF_);
   tree->SetBranchAddress("patJetPfAk04chf_", &patJetPfAk04chf_, &b_patJetPfAk04chf_);
   tree->SetBranchAddress("patJetPfAk04nhf_", &patJetPfAk04nhf_, &b_patJetPfAk04nhf_);
   tree->SetBranchAddress("patJetPfAk04cemf_", &patJetPfAk04cemf_, &b_patJetPfAk04cemf_);
   tree->SetBranchAddress("patJetPfAk04nemf_", &patJetPfAk04nemf_, &b_patJetPfAk04nemf_);
   tree->SetBranchAddress("patJetPfAk04cmult_", &patJetPfAk04cmult_, &b_patJetPfAk04cmult_);
   tree->SetBranchAddress("patJetPfAk04nconst_", &patJetPfAk04nconst_, &b_patJetPfAk04nconst_);
   tree->SetBranchAddress("patJetPfAk04BDiscdeepFudsg_", &patJetPfAk04BDiscdeepFudsg_, &b_patJetPfAk04BDiscdeepFudsg_);
   tree->SetBranchAddress("patJetPfAk04BDiscdeepFb_", &patJetPfAk04BDiscdeepFb_, &b_patJetPfAk04BDiscdeepFb_);
   tree->SetBranchAddress("patJetPfAk04BDiscdeepFbb_", &patJetPfAk04BDiscdeepFbb_, &b_patJetPfAk04BDiscdeepFbb_);
   tree->SetBranchAddress("patJetPfAk04BDiscdeepFc_", &patJetPfAk04BDiscdeepFc_, &b_patJetPfAk04BDiscdeepFc_);
   tree->SetBranchAddress("patJetPfAk04BDiscdeepFcc_", &patJetPfAk04BDiscdeepFcc_, &b_patJetPfAk04BDiscdeepFcc_);
   tree->SetBranchAddress("patJetPfAk04BDiscCSVv2_", &patJetPfAk04BDiscCSVv2_, &b_patJetPfAk04BDiscCSVv2_);
   tree->SetBranchAddress("patJetPfAk04BDiscpfCMVA_", &patJetPfAk04BDiscpfCMVA_, &b_patJetPfAk04BDiscpfCMVA_);
   tree->SetBranchAddress("unc_", &unc_, &b_unc_);
   tree->SetBranchAddress("patJetPfAk04PtUp_", &patJetPfAk04PtUp_, &b_patJetPfAk04PtUp_);
   tree->SetBranchAddress("patJetPfAk04PtDn_", &patJetPfAk04PtDn_, &b_patJetPfAk04PtDn_);
   tree->SetBranchAddress("id1_pdfInfo_", &id1_pdfInfo_, &b_id1_pdfInfo_);
   tree->SetBranchAddress("id2_pdfInfo_", &id2_pdfInfo_, &b_id2_pdfInfo_);
   tree->SetBranchAddress("x1_pdfInfo_", &x1_pdfInfo_, &b_x1_pdfInfo_);
   tree->SetBranchAddress("x2_pdfInfo_", &x2_pdfInfo_, &b_x2_pdfInfo_);
   tree->SetBranchAddress("scalePDF_pdfInfo_", &scalePDF_pdfInfo_, &b_scalePDF_pdfInfo_);
   tree->SetBranchAddress("ptHat_", &ptHat_, &b_ptHat_);
   tree->SetBranchAddress("mcWeight_", &mcWeight_, &b_mcWeight_);
   tree->SetBranchAddress("mcWeights_", &mcWeights_, &b_mcWeights_);
   tree->SetBranchAddress("nup", &nup, &b_nup);
  }

    double pi = 3.1415926535897932384626433832795028841971693;
    double DeltaR(double eta1, double eta2, double phi1, double phi2)
        {
        double deta = eta2 - eta1;
        double dphi = phi2 - phi1;
        if (fabs(dphi) > pi) dphi = 6.28 - fabs(dphi);
        double DELTAR = sqrt(pow(dphi,2)+pow(deta,2))*1.0;
        return DELTAR;
        }
    double DeltaPhi(double phi1, double phi2)
        {
        double dphi = phi2 - phi1;
        if (fabs(dphi) > pi) dphi = 6.28 - fabs(dphi);
        return dphi;
        }

