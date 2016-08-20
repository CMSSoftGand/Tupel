#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
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
#include "TChain.h"
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
<<<<<<< HEAD
#include "data/era2016/NeutrinoSolver.cc"
#include "data/era2016/NeutrinoSolver.h"
//#include "/user/mgul/Higgs_tottbar/anlyzer765/CMSSW_7_6_5/src/backgourd/ttbar/BTagCalibrationStandalone.h"
//#include "/user/mgul/Higgs_tottbar/anlyzer765/CMSSW_7_6_5/src/backgourd/ttbar/BTagCalibrationStandalone.cc"
#include "/user/mgul/Higgs_tottbar/anlyzer765/CMSSW_7_6_5/src/backgourd/ttbar/BtagUncertaintyComputer.h"
#include "/user/mgul/Higgs_tottbar/anlyzer765/CMSSW_7_6_5/src/backgourd/ttbar/BtagUncertaintyComputer.cc"
#include "BTagCalibration.h"
=======
#include "plugins/NeutrinoSolver.cc"
#include "interface/NeutrinoSolver.h"
//#include "/user/mgul/Higgs_tottbar/anlyzer765/CMSSW_7_6_5/src/backgourd/ttbar/BTagCalibrationStandalone.h"
//#include "/user/mgul/Higgs_tottbar/anlyzer765/CMSSW_7_6_5/src/backgourd/ttbar/BTagCalibrationStandalone.cc"
#include "interface/BtagUncertaintyComputer.h"
#include "plugins/BtagUncertaintyComputer.cc"
#include "interface/BTagCalibration.h"
>>>>>>> 22e0d1f5b734148c99f2f8f885198e7750e4829b
//#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"

#include "simpleReader.h"
<<<<<<< HEAD
#include "/user/mgul/Higgs_tottbar/anlyzer765/CMSSW_7_6_5/src/backgourd/ttbar/ttbar0/standalone_LumiReWeighting.h"
=======
#include "interface/standalone_LumiReWeighting.h"
>>>>>>> 22e0d1f5b734148c99f2f8f885198e7750e4829b
//#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
//#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
//#include "JetMETCorrections/Modules/interface/JetResolution.h"
using namespace std;

/*std::vector<float> getJetResolutionScales(float pt, float eta, float genjpt)
{
    std::vector<float> res(3,1.0);
      float ptSF(1.0), ptSF_err(0.0);
        if(TMath::Abs(eta)<0.5)       { ptSF=1.095; ptSF_err = 0.018; }
          else if(TMath::Abs(eta)<0.8)  { ptSF=1.120; ptSF_err = 0.028; }
            else if(TMath::Abs(eta)<1.1)  { ptSF=1.097; ptSF_err = 0.017; }
              else if(TMath::Abs(eta)<1.3)  { ptSF=1.103; ptSF_err = 0.033; }
                else if(TMath::Abs(eta)<1.7)  { ptSF=1.118; ptSF_err = 0.014; }
                  else if(TMath::Abs(eta)<1.9)  { ptSF=1.100; ptSF_err = 0.033; }
                    else if(TMath::Abs(eta)<2.1)  { ptSF=1.162; ptSF_err = 0.044; }
                      else if(TMath::Abs(eta)<2.3)  { ptSF=1.160; ptSF_err = 0.048; }
                        else if(TMath::Abs(eta)<2.5)  { ptSF=1.161; ptSF_err = 0.060; }
                          else if(TMath::Abs(eta)<2.8)  { ptSF=1.209; ptSF_err = 0.059; }
                            else if(TMath::Abs(eta)<3.0)  { ptSF=1.564; ptSF_err = 0.321; }
                              else if(TMath::Abs(eta)<3.2)  { ptSF=1.384; ptSF_err = 0.033; }
                                else if(TMath::Abs(eta)<5.0)  { ptSF=1.216; ptSF_err = 0.050; }
                                  res[0] = TMath::Max((Float_t)0.,(Float_t)(genjpt+(ptSF)*(pt-genjpt)))/pt;
                                    res[1] = TMath::Max((Float_t)0.,(Float_t)(genjpt+(ptSF-ptSF_err)*(pt-genjpt)))/pt;
                                      res[2] = TMath::Max((Float_t)0.,(Float_t)(genjpt+(ptSF+ptSF_err)*(pt-genjpt)))/pt;
                                        return res;
}
*/
  void simpleReader(TString fin,TString fout){
  TFile *f = TFile::Open(fin);
  if (f->IsZombie()) {
     printf("Input root files doesn't open, please have a look:\n");
     return;
     }
  cout<<"this is ff:  "<<fin<<"   ;  "<<fout<<endl;
  TTree *t = (TTree*)f->Get("tupel/MuonTree");
  TFile *theFile = new TFile (fout+".root","RECREATE");
  theFile->cd();  
  //  bool is_mu (true);
  bool is_mu (false);
  bool run_sys (false);
<<<<<<< HEAD
  TString era("/user/mgul/Higgs_tottbar/anlyzer808/CMSSW_8_0_11/src/Tupel/Tupel/data/era2016/");
=======
  TString era("/user/mgul/Higgs_tottbar/Anz_8011/CMSSW_8_0_11/src/Tupel/Tupel/data/era2016/");
>>>>>>> 22e0d1f5b734148c99f2f8f885198e7750e4829b
  ////  standalone_LumiReWeighting puWeight(201525,0), puUp(201525,1), puDown(201525,-1);
          TString btagUncUrl(era+"btagSFactors.csv");
          gSystem->ExpandPathName(btagUncUrl);
          std::vector<BTagCalibrationReader *> sfbReaders, sflReaders;
          TString btagEffExpUrl(era+"expTageff.root");
          gSystem->ExpandPathName(btagEffExpUrl);
          std::map<TString, TGraphAsymmErrors *> expBtagEff, expBtagEffPy8;
          BTagSFUtil myBTagSFUtil;
          if (!realdata){
            BTagCalibration btvcalib("csvv2", btagUncUrl.Data());
            sfbReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "mujets", "central") );
            sfbReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "mujets", "down") );
            sfbReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "mujets", "up") );
            sflReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "incl", "central") );
            sflReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "incl", "down") );
            sflReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "incl", "up") );
            TFile *beffIn=TFile::Open(btagEffExpUrl);
            expBtagEff["b"]=(TGraphAsymmErrors *)beffIn->Get("b");
            expBtagEff["c"]=(TGraphAsymmErrors *)beffIn->Get("c");
            expBtagEff["udsg"]=(TGraphAsymmErrors *)beffIn->Get("udsg");
            beffIn->Close();
            }
//     experimental syst Uncer
     std::vector<TString> expSysts;
     expSysts.push_back("UncJER");
     //expSysts.push_back("UncJEC");
     expSysts.push_back("Pileup");
     expSysts.push_back("Trigger");
     if (is_mu)expSysts.push_back("MuEfficiency");
     if (!is_mu)expSysts.push_back("EleEfficiency");
     expSysts.push_back("BtagEff");
     expSysts.push_back("CtagEff");
     expSysts.push_back("LtagEff");
     Int_t nExpSysts=expSysts.size();

     if (is_mu){TDirectory *mujets_2btag = theFile->mkdir("mujets_2btag");mujets_2btag->cd();}
     if (!is_mu){TDirectory *eljets_2btag = theFile->mkdir("eljets_2btag");eljets_2btag->cd();}
     TH1::SetDefaultSumw2();
     TH2::SetDefaultSumw2();
     branchAdd(t);

//  TH1D* h_jets_pt = new TH1D("MZ","M(Z)#rightarrow #mu#mu",40, 71,111.);
     TH1D* h_events_eachCut = new TH1D("events_eachCut","events_eachCut",10, 0.0,10.0);
     TH1D* h_PtJets_afterCut = new TH1D("jets_Pt_afterCut","jets_Pt_afterCut",10, 0.0,1200.);
     TH1D* h_massJets_afterCut = new TH1D("jets_mass_afterCut","jets_mass_afterCut",10, 0.0,200.);
     TH1D* h_Ejets_afterCut = new TH1D("jets_energy_afterCut","jets_energy_afterCut",10, 0.0,2500.);
     TH1D* h_sum_jets_lepPt_afterCut = new TH1D("jets_lep_Pt_afterCut","jets_lep_Pt_afterCut",10, 0.0,1600.);	
//------------------------------MET--------------------------------------------------//
     TH1D* h_met_Px = new TH1D("MET_Px","MET_Px",50, -200.0,200.);
     TH1D* h_met_Py = new TH1D("MET_Py","MET_Py",50, -200.0,200.);
     TH1D* h_met_Pz = new TH1D("MET_Pz","MET_Pz",50, -200.0,200.);
     TH1D* h_met_E = new TH1D("MET_E","MET_E",100, 0.0,200.);
//------------------------------Muon-------------------------------------------------//
     TH1D* h_patMuonPfIsoDbeta_after = new TH1D("patMuonPfIsoDbeta_after","patMuonPfIsoDbeta_after",50, 0.0,0.17);
     TH1D* h_patMuonPfIsoDbeta = new TH1D("patMuonPfIsoDbeta","patMuonPfIsoDbeta",50, 0.0,2.);
     TH1D* h_Pt_lep = new TH1D("Pt_elec","Pt_elec",50, 20.0,200.);
     TH1D* h_Eta_lep = new TH1D("Eta_elec","Eta_elec",50, -2.6,2.6);
     TH1D* h_Phi_lep = new TH1D("Phi_elec","Phi_elec",50,-3.4,3.4);
     TH1D* h_E_lep = new TH1D("E_elec","E_elec",50, 0.0,400.);
//-----------------------------Lep W  ------------------------------------------------//
     TH1D* h_Pt_muon_met = new TH1D("Pt_LepW","Pt_LepW",50, 0.0,400.);
     TH1D* h_Eta_muon_met = new TH1D("Eta_LepW","Eta_LepW",50, -2.6,2.6);
     TH1D* h_Phi_muon_met = new TH1D("Phi_LepW","Phi_LepW",50, -3.4,3.4);
     TH1D* h_M_muon_met = new TH1D("M_LepW","M_LepW",50, 0.0,200.0);
//------------------------------ Jets ------------------------------------------------//
  TH1D* h_no_Jets = new TH1D("No_ofJets","No_ofJets",7, 0.5,7.5);
    TH1D* h_PtJet1 = new TH1D("PtJet1","PtJet1",50, 0.0,400.);
      TH1D* h_EtaJet1 = new TH1D("EtaJet1","Eta_Jet1",50, -2.5,2.5);
        TH1D* h_PhiJet1 = new TH1D("PhiJet1","PhiJet1",50, -3.4,3.4);
          TH1D* h_EJet1 = new TH1D("EJet1","EJet1",50, 0.0,700.);
            TH1D* h_MJet1 = new TH1D("MJet1","MJet1",50, 0.0,100.);

  TH1D* h_PtJet2 = new TH1D("PtJet2","PtJet2",50, 0.0,300.);
      TH1D* h_EtaJet2 = new TH1D("EtaJet2","Eta_Jet2",50, -2.5,2.5);
        TH1D* h_PhiJet2 = new TH1D("PhiJet2","PhiJet2",50, -3.4,3.4);
          TH1D* h_EJet2 = new TH1D("EJet2","EJet2",50, 0.0,500.);
            TH1D* h_MJet2 = new TH1D("MJet2","MJet2",50, 0.0,80.);

  TH1D* h_PtJet3 = new TH1D("PtJet3","PtJet3",50, 0.0,200.);
    TH1D* h_EtaJet3 = new TH1D("EtaJet3","Eta_Jet3",50, -2.5,2.5);
      TH1D* h_PhiJet3 = new TH1D("PhiJet3","PhiJet3",50, -3.4,3.4);
        TH1D* h_EJet3 = new TH1D("EJet3","EJet3",50, 0.0,400.);
          TH1D* h_MJet3 = new TH1D("MJet3","MJet3",50, 0.0,60.);

    TH1D* h_PtJet4 = new TH1D("PtJet4","PtJet4",50, 0.0,150.);
      TH1D* h_EtaJet4 = new TH1D("EtaJet4","Eta_Jet4",50, -2.5,2.5);
        TH1D* h_PhiJet4 = new TH1D("PhiJet4","PhiJet4",50, -3.4,3.4);
          TH1D* h_EJet4 = new TH1D("EJet4","EJet4",50, 0.0,300.);
            TH1D* h_MJet4 = new TH1D("MJet4","MJet4",50, 0.0,50.);
//------------------------------bJets-------------------------------------------------//
  TH1D* h_n_bjets = new TH1D("n_bjets","n_bjets",5, 0.5,5.5);
  TH1D* h_n_bjets_aft = new TH1D("n_bjets_aft","n_bjets_aft",5, 1,6);
  TH1D* h_n_cjets_aft = new TH1D("n_cjets_aft","n_cjets_aft",5, 1,6);
    TH1D* h_Pt_bJet1 = new TH1D("Pt_bJet1","Pt_bJet1",50, 0.0,300.);
      TH1D* h_Eta_bJet1 = new TH1D("Eta_bJet1","Eta_bJet1",50, -2.5,2.5);
        TH1D* h_Phi_bJet1 = new TH1D("Phi_bJet1","Phi_bJet1",50, -3.4,3.4);
          TH1D* h_E_bJet1 = new TH1D("E_bJet1","E_bJet1",50, 0.0,600.0);
            TH1D* h_M_bJet1 = new TH1D("M_bJet1","M_bJet1",50, 0.0,80.0);

  TH1D* h_Pt_bJet2 = new TH1D("Pt_bJet2","Pt_bJet2",50, 0.0,200.);
    TH1D* h_Eta_bJet2 = new TH1D("Eta_bJet2","Eta_bJet2",50, -2.5,2.5);
      TH1D* h_Phi_bJet2 = new TH1D("Phi_bJet2","Phi_bJet2",50, -3.4,3.4);
        TH1D* h_E_bJet2 = new TH1D("E_bJet2","E_bJet2",50, 0.0,400.0);
          TH1D* h_M_bJet2 = new TH1D("M_bJet2","M_bJet2",50, 0.0,60.0);

  TH1D* h_DPhi_bj12 = new TH1D("DPhi_bj12","DPhi_bj12",50, 3.14,3.14);
      TH1D* h_DR_bj12 = new TH1D("DR_bj12","DR_bj12",50, 0.0,10);
//------------------------------Light jets-------------------------------------------//
  TH1D* h_n_ljets = new TH1D("n_ljets","n_ljets",7, 0.5,7.5);
  TH1D* h_n_ljets_aft = new TH1D("n_ljets_aft","n_ljets_aft",5, 1,6);
    TH1D* h_Pt_lJet1 = new TH1D("Pt_lJet1","Pt_lJet1",50, 0.0,300.);
      TH1D* h_Eta_lJet1 = new TH1D("Eta_lJet1","Eta_lJet1",50, -2.5,2.5);
        TH1D* h_Phi_lJet1 = new TH1D("Phi_lJet1","Phi_lJet1",50, -3.4,3.4);
          TH1D* h_E_lJet1 = new TH1D("E_lJet1","E_lJet1",50, 0.0,400.0);
            TH1D* h_M_lJet1 = new TH1D("M_lJet1","M_lJet1",50, 0.0,80.0);

  TH1D* h_Pt_lJet2 = new TH1D("Pt_lJet2","Pt_lJet2",50, 0.0,200.);
    TH1D* h_Eta_lJet2 = new TH1D("Eta_lJet2","Eta_lJet2",50, -2.5,2.5);
      TH1D* h_Phi_lJet2 = new TH1D("Phi_lJet2","Phi_lJet2",50, -3.4,3.4);
        TH1D* h_E_lJet2 = new TH1D("E_lJet2","E_lJet2",50, 0.0,300.0);
          TH1D* h_M_lJet2 = new TH1D("M_lJet2","M_lJet2",50, 0.0,60.0);
//------------------------------Had W -----------------------------------------------//
  TH1D* h_Pt_hadW = new TH1D("Pt_hadW","Pt_hadW",50, 0.0,400.);
    TH1D* h_Eta_hadW = new TH1D("Eta_hadW","Eta_hadW",50, -3.5,3.5);
      TH1D* h_Phi_hadW = new TH1D("Phi_hadW","Phi_hadW",50, -3.4,3.4);
        TH1D* h_E_hadW = new TH1D("E_hadW","E_hadW",50, 0.0,800.0);
          TH1D* h_M_hadW = new TH1D("M_hadW","M_hadW",50, 30.0,200.0);
//------------------------------Leptonic Top ----------------------------------------//
  TH1D* h_M_lepTop  = new TH1D("M_lepTop","M_lepTop",50, 0.0,400.0);
    TH1D* h_Pt_lepTop = new TH1D("Pt_lepTop","Pt_lepTop",50, 0.0,500.0);
//------------------------------Hadronic Top----------------------------------------//
  TH1D* h_M_hadTop  = new TH1D("M_hadTop","M_hadTop",50, 0.0,400.0);
    TH1D* h_Pt_hadTop  = new TH1D("Pt_hadTop","Pt_hadTop",50, 0.0,500.0);
//------------------------------Heavy Higgs----------------------------------------//
  TH1D* h_Heavy_H   = new TH1D("Mass_H","Mass_H",100, 200.0,1200.0);
    TH1D* hPt_Heavy_H   = new TH1D("Pt_H","Pt_H",50, 0.0,500.0);
//      TH1D* hEta_Heavy_H   = new TH1D("Pt_H","Pt_H",50, 0.0,500.0);
//------------------------------------------------------------------------------------//
  TH1D* hDPhi_j1_l = new TH1D("DPhi_j1_l","DPhi j1 l",15, -3.5, 3.5);
  TH1D* hDPhi_j2_l = new TH1D("DPhi_j2_l","DPhi j2 l",15, -3.5, 3.5);
  TH1D* hdR_jet1_lep = new TH1D("dR_j1_l","dR j1 l",15, 0., 5);
  TH1D* hdR_jet2_lep = new TH1D("dR_j2_l","dR j2 l",15, 0., 5);
  TH1D* hDPhi_j1_l_aftCut = new TH1D("DPhi_j1_l_aftCut","DPhi j1 l aftCut",15, -3.5, 3.5);
  TH1D* hDPhi_j2_l_aftCut = new TH1D("DPhi_j2_l_aftCut","DPhi j2 l aftCut",15, -3.5, 3.5);
  TH1D* hDPhi_nu_j = new TH1D("DPhi_nu_j","DPhi nu j",15, -3.5, 3.5);
  TH1D* hDPhi_nu_bj = new TH1D("DPhi_nu_bj","DPhi nu bj",15, -3.5, 3.5);
  TH1D* hDPhi_j_bj = new TH1D("DPhi_j_bj","DPhi j bj",15, -3.5, 3.5);

    TH1D* hDPhiM600_j_l = new TH1D("DPhiM600_j_l","DPhiM600 j l",15, -3.5, 3.5);
    TH1D* hDPhiM600_nu_bj = new TH1D("DPhiM600_nu_bj","DPhiM600 nu bj",15, -3.5, 3.5);
    TH1D* hDPhiM600_j_bj = new TH1D("DPhiM600_j_bj","DPhiM600 j bj",15, -3.5, 3.5);

      TH1D* hDPhiM800_j_l = new TH1D("DPhiM800_j_l","DPhiM800 j l",15, -3.5, 3.5);
      TH1D* hDPhiM800_nu_bj = new TH1D("DPhiM800_nu_bj","DPhiM800 nu bj",15, -3.5, 3.5);
      TH1D* hDPhiM800_j_bj = new TH1D("DPhiM800_j_bj","DPhiM800 j bj",15, -3.5, 3.5);

TH2D* hdphi_mass = new TH2D("dphi_mass","dphi_mass",15,-3.15,3.15,15,0,1000);
  TH1D* Phi_jet = new TH1D("Phi_jet","Phi jet",50, -3.5, 3.5);
    TH1D* hDPhi_bjets_lep = new TH1D("DPhi_bjets_lep","DPhi bjets lep",50, -3.5, 3.5);
      TH1D* hDPhi_Lepbjet_Hadbjet = new TH1D("DPhi_Lepbjet_Hadbjet","DPhi Lepbjet Hadbjet",50, -3.5, 3.5);
        TH1D* h_EvtInfo_NumVtx  = new TH1D("EvtInfo_NumVtx","EvtInfo_NumVtx",40, 0.0,40.0);
        TH1D* h_EvtInfo_NumVtx_w  = new TH1D("EvtInfo_NumVtx_w","EvtInfo_NumVtx_w",40, 0.0,40.0);
          TH1D* h_PU_npT  = new TH1D("PU_npT","PU_npT",50, 0.0,40.0);
                    TH1D* h_bDisct_CSVv2 = new TH1D("bDiscCSVv2","bDiscCSVv2",50, 0.0,1.0);
                      TH1D* h_bDisct_CSVv2aft = new TH1D("CSV_v2","CSV_v2",50, 0.0,1.0);
                        TH1D* h_pfJet_cmult = new TH1D("Jet_cMult","Jet_cMult",10, 0.0,40.0);
    std::map<TString, TH2 *> plots2d;
    plots2d["metptshapes_exp"]= new TH2F("metptshapes_exp", ";Missing transverse energy [GeV];Events" , 100,0.,200., 2*nExpSysts,0,2*nExpSysts);
    plots2d["nbjetshapes_exp"]= new TH2F("nbjetshapes_exp", "; bJets Multiplicity ;Events" , 5,1,6, 2*nExpSysts,0,2*nExpSysts);
    plots2d["bjetptshapes_exp"]= new TH2F("bjetptshapes_exp", "; bJets pt ;Events" , 50,0.,300., 2*nExpSysts,0,2*nExpSysts);
    plots2d["nljetshapes_exp"]= new TH2F("nljetshapes_exp", "; light Jets Multiplicity ;Events" , 5,1,6, 2*nExpSysts,0,2*nExpSysts);
    plots2d["ljetptshapes_exp"]= new TH2F("ljetptshapes_exp", "; light Jets pt ;Events" , 50,0,300, 2*nExpSysts,0,2*nExpSysts);
    plots2d["ncjetshapes_exp"]= new TH2F("ncjetshapes_exp", "; c Jets Multiplicity ;Events" , 5,1,6, 2*nExpSysts,0,2*nExpSysts);
    for(Int_t isyst=0; isyst<nExpSysts; isyst++){
      for(Int_t ivar=0; ivar<2; ivar++){
        TString label(expSysts[isyst] + (ivar==0 ? "Down" : "Up"));
        plots2d["metptshapes_exp"]->GetYaxis()->SetBinLabel(2*isyst+ivar+1, label);
        plots2d["nbjetshapes_exp"]->GetYaxis()->SetBinLabel(2*isyst+ivar+1, label);
        plots2d["bjetptshapes_exp"]->GetYaxis()->SetBinLabel(2*isyst+ivar+1, label);
        plots2d["nljetshapes_exp"]->GetYaxis()->SetBinLabel(2*isyst+ivar+1, label);
        plots2d["ljetptshapes_exp"]->GetYaxis()->SetBinLabel(2*isyst+ivar+1, label);
        plots2d["ncjetshapes_exp"]->GetYaxis()->SetBinLabel(2*isyst+ivar+1, label);
      }
    }

        double real_count=0;
        double complex_count=0;
        double solved_complex_count=0;

//        vector<float> weight_;
        vector<float> e_met;
        vector<float> px_met;
        vector<float> py_met;
        vector<float> dphi_nu_j;
        vector<float> dphi_nu_bj;
        vector<float> dphi_j_bj;
        vector<float> muon_pt;
        vector<float> muon_eta;
        vector<float> dphi_j_l;
        vector<float> n_ofjets;
        vector<float> n_ofbjets;
        vector<float> n_ofljets;
        vector<float> jet1_pt;
        vector<float> bjet1_pt;
        vector<float> bjet2_eta;
        vector<float> ht;
        vector<float> jets_e;
        vector<float> jets_m;
        vector<float> jets_l_pt;
        vector<float> hadW_pt;
        vector<float> hadW_rap;
        vector<float> lepW_pt;
        vector<float> lepW_rap;
        vector<float> top_m;
        vector<float> top_pt;
        vector<float> top_rap;
        vector<float> atop_m;
        vector<float> atop_pt;
        vector<float> atop_rap;
        vector<float> higgs_pt;
        vector<float> higgs_m;
        vector<float> higgs_rap;

        TLorentzVector v_elec;
        TLorentzVector v_muon;
        TLorentzVector v_muon9;
        TLorentzVector v_met;
        TLorentzVector v_lepW;
        TLorentzVector v_jets;
        TLorentzVector v_jetAll[1000];
        TLorentzVector v_jet[100];
        TLorentzVector v_bjets;
        TLorentzVector v_bjets9;
        TLorentzVector v_Hadbjet;
        TLorentzVector v_Lepbjet;


        TLorentzVector v_ljet_temp;
        TLorentzVector v_ljets;
        TLorentzVector hadWT_;

        TLorentzVector v_Heavy_H;
        TTree *weight_tree;
        float met_E=0.;

        Int_t nentries(t->GetEntriesFast());
          nentries=5000;
        for (int jentry=0; jentry < nentries; jentry++)
        {
        if(jentry%1000==0)cout<<" << "<<jentry<<"/"<<nentries<<endl;
        double w=1;
<<<<<<< HEAD
        std::vector<TGraph *>puWgtGr;
        std::vector<float> puWeight(3,1.0);
        TString puWgtUrl("/user/mgul/Higgs_tottbar/anlyzer765/CMSSW_7_6_5/src/backgourd/ttbar/pileupWgts.root");
=======
        float iSecret;
        srand (time(NULL));
        iSecret = rand() % 100 + 1;
        bool run_273158_274093 (iSecret <= 5.);
        bool run_274094_276097 (5. < iSecret && iSecret < 63.);
        bool elrun_below_273726 (iSecret <= 4.);
        bool elrun_above_273726 (4. < iSecret && iSecret < 42.);
        std::vector<TGraph *>puWgtGr;
        std::vector<float> puWeight(3,1.0);
        TString puWgtUrl(era+"pileupWgts.root");
>>>>>>> 22e0d1f5b734148c99f2f8f885198e7750e4829b
        if (!realdata){
          TFile *fIn=TFile::Open(puWgtUrl);
          if(fIn){
          puWgtGr.push_back( (TGraph *)fIn->Get("puwgts_nom") );
          puWgtGr.push_back( (TGraph *)fIn->Get("puwgts_down") );
          puWgtGr.push_back( (TGraph *)fIn->Get("puwgts_up") );
          fIn->Close();
            }
          }
          if(!realdata){
          if(puWgtGr.size())
          {
          puWeight[0]=puWgtGr[0]->Eval(int(PU_npT));
          puWeight[1]=puWgtGr[1]->Eval(int(PU_npT));
          puWeight[2]=puWgtGr[2]->Eval(int(PU_npT));
          }
          }
          w=puWeight[0];
        t->GetEntry(jentry);
        h_events_eachCut->Fill(0);
        if (is_mu) if (HLT_IsoMu20!=1 && HLT_IsoTkMu20!=1)continue;
        if (!is_mu)if (HLT_Ele23_WPLoose_Gsf_v!=1)continue;
        h_events_eachCut->Fill(1);
        //noise
          if(Flag_HBHENoiseFilter!=1 || Flag_HBHENoiseIsoFilter!=1 || Flag_CSCTightHalo2015Filter!=1)continue;
          if (Flag_EcalDeadCellTriggerPrimitiveFilter!=1 || Flag_goodVertices!=1 || Flag_eeBadScFilter!=1)continue;
      int n_pat_elec = 0,n_elec15=0;
      std::vector<TLorentzVector>lep_vector;
      TLorentzVector v_lep;
      vector<unsigned int> n_lep_v;
      vector<int> lep_charge;
      vector<unsigned int> n_lep15_v;

      std::vector<TLorentzVector>el_vector;
      vector<int>el_charge;
      vector<unsigned int> n_elec_v;
      vector<unsigned int> n_elec15_v;

      for (unsigned int elec =0; elec < patElecPt_->size(); ++elec){
        if(patElecPt_->at(elec)>20 && fabs(patElecScEta_->at(elec))<2.5 && patElecIdveto_->at(elec) !=0){
            n_elec15_v.push_back(n_elec15);
            n_elec15++;
            }
      if (patElecPt_->at(elec) > 25. && fabs(patElecScEta_->at(elec)) <2.5 && patElecIdtight_->at(elec)>0 ){


     if (1.4442 < fabs( patElecScEta_->at(elec)) && fabs(patElecScEta_->at(elec))< 1.5660)continue;
      n_pat_elec++;
      n_elec_v.push_back(n_pat_elec);
      v_elec.SetPtEtaPhiE(patElecPt_->at(elec),patElecEta_->at(elec),patElecPhi_->at(elec),patElecEnergy_->at(elec));
      el_vector.push_back(v_elec);
      el_charge.push_back(patElecCharge_->at(elec));
      }
      }
      //////-----------------------------------muon------------------------------------------------------------//
      std::vector<TLorentzVector>mu_vector;
      vector<int>mu_charge;
      vector<unsigned int> n_muon_v;
      vector<unsigned int> n_muon15_v;
      vector<double> mu_eta_v;
      float Pt_muon=0,Eta_muon=0,Phi_muon=0,E_muon=0,muon_TIso=0;
      int n_muon=0,n_muon15=0;
      for(unsigned int mu=0; mu<patMuonPt_->size();mu++){


          muon_TIso=patMuonPfIsoDbeta_->at(mu);
          h_patMuonPfIsoDbeta->Fill(muon_TIso,w);
          if(patMuonPt_->at(mu) > 22. && fabs(patMuonEta_->at(mu))<2.4 &&
             patMuonTightId_->at(mu)>0 && patMuonPfIsoDbeta_->at(mu) < 0.15)
            {

            double mu_eta = patMuonEta_->at(mu);      
            mu_eta_v.push_back(mu_eta);
            n_muon++;
            n_muon_v.push_back(n_muon);
            v_muon.SetPtEtaPhiE(patMuonPt_->at(mu),patMuonEta_->at(mu),
            patMuonPhi_->at(mu),patMuonEn_->at(mu));
            mu_vector.push_back(v_muon);
            mu_charge.push_back(patMuonCharge_->at(mu));
            muon_TIso=patMuonPfIsoDbeta_->at(mu);
            }
            if(patMuonPt_->at(mu)>10 && fabs(patMuonEta_->at(mu))<2.4 && patMuonLooseId_->at(mu)>0 && fabs(patMuonPfIsoDbeta_->at(mu))<0.25){
            n_muon15_v.push_back(n_muon15);
            n_muon15++;
            }
            }
        if (!is_mu){n_lep_v = n_elec_v; n_lep15_v = n_elec15_v; lep_vector = el_vector; lep_charge = el_charge; v_lep=v_elec;}
        if (is_mu){n_lep_v = n_muon_v; n_lep15_v = n_muon15_v; lep_vector = mu_vector; lep_charge = mu_charge; v_lep=v_muon;}
            //-------------MET--------------------------------------------------//
            for (unsigned int nu =0; nu < METPt->size(); nu++)
                {
                v_met.SetPxPyPzE(METPx->at(0),METPy->at(0),METPz->at(0),METE->at(0));
                }
            //----------------------------------PF Jets------------------------------------------------//

            int n_pat_bjets=0,ncjets=0;
            int n_ljets=0;
            float jet_cMult=0;
            int n_pat_jets=0;
            double DR_mu_j=9999;
            vector<unsigned int> no_jets;
            vector<unsigned int> n_bjets_v;
            vector<TLorentzVector> bjets_v;
            vector<TLorentzVector> ljets_v;
            vector<TLorentzVector> vec_bjets;
            vector<unsigned int> n_ljets_v;
            vector<float> csv_v;
            double csv=-99, jet_pt=0.,nom=1., denom=1.,ind_jet_const=0.;
            vector <TLorentzVector> temp1;
            vector <TLorentzVector> v_ljetsAll;
            TLorentzVector lorenz_v; 
            vector<TLorentzVector>jet_vector;
            double JER_Uncer = 0.;
            TLorentzVector metp4; 
            for (unsigned int pf=0;pf < patJetPfAk04PtJERSmear->size();++pf){
              double Up = patJetPfAk04PtJERSmearUp->at(pf);//Fix Me I will be after or before the cut
              double Dn = patJetPfAk04PtJERSmearDn->at(pf);
              double Cn = patJetPfAk04PtJERSmear->at(pf);
              double sig_Up=abs(Up-Cn)/Cn;
              double sig_Dn=abs(Dn-Cn)/Cn;
              JER_Uncer = max(sig_Up,sig_Dn);
              
              if(patJetPfAk04PtJERSmear->at(pf) > 30.&& fabs(patJetPfAk04Eta_->at(pf)) < 2.4 && patJetPfAk04LooseId_->at(pf)>0){
              n_pat_jets++;
              TLorentzVector v_jetsTemp;
              v_jetsTemp.SetPtEtaPhiE(patJetPfAk04PtJERSmear->at(pf),patJetPfAk04Eta_->at(pf),patJetPfAk04Phi_->at(pf),patJetPfAk04En_->at(pf));

              if (n_lep_v.size()>0)DR_mu_j= DeltaR(lep_vector[0].Eta(), v_jetsTemp.Eta(), lep_vector[0].Phi(), v_jetsTemp.Phi());
              if(DR_mu_j<0.4)continue;
              v_jets.SetPtEtaPhiE(patJetPfAk04PtJERSmear->at(pf),patJetPfAk04Eta_->at(pf),patJetPfAk04Phi_->at(pf),patJetPfAk04En_->at(pf));
              jet_vector.push_back(v_jets);
              h_bDisct_CSVv2->Fill(patJetPfAk04BDiscCSVv2_->at(pf),w);
              if(patJetPfAk04BDiscCSVv2_->at(pf) > csv) csv = patJetPfAk04BDiscCSVv2_->at(pf);
              csv_v.push_back(patJetPfAk04BDiscCSVv2_->at(pf));
              jet_cMult=patJetPfAk04cmult_->at(pf);
              h_no_Jets->Fill(n_pat_jets,w);
              no_jets.push_back(n_pat_jets);
              jet_pt=patJetPfAk04PtJERSmear->at(pf); 
              int jflav( abs(patJetPfAk04PartonFlavour_->at(pf)) );
              bool isBTagged(csv>0.800);

              if(!realdata){
                float jptForBtag(v_jets.Pt()>1000. ? 999. : v_jets.Pt()), jetaForBtag(fabs(v_jets.Eta()));
                float expEff(1.0), jetBtagSF(1.0);
                if(abs(jflav)==4){
                  ncjets++;
                  expEff    = expBtagEff["c"]->Eval(jptForBtag);
                  jetBtagSF = sfbReaders[0]->eval( BTagEntry::FLAV_C, jetaForBtag, jptForBtag);
                  jetBtagSF *= expEff>0 ? 1 : 0. ;
                  }
                else if(abs(jflav)==5){
                  n_pat_bjets++;
                  expEff    = expBtagEff["b"]->Eval(jptForBtag);
                  jetBtagSF = sfbReaders[0]->eval( BTagEntry::FLAV_B, jetaForBtag, jptForBtag);
                  jetBtagSF *= expEff>0 ? 1 : 0. ;
                  }
                else{
                  n_ljets++;
                  expEff    = expBtagEff["udsg"]->Eval(jptForBtag);
                  jetBtagSF = sflReaders[0]->eval( BTagEntry::FLAV_UDSG, jetaForBtag, jptForBtag);
                  jetBtagSF *= expEff>0 ? 1 : 0. ;
                  }
                myBTagSFUtil.modifyBTagsWithSF(isBTagged,    jetBtagSF,     expEff);
              }
              if(isBTagged){
                n_bjets_v.push_back(n_pat_jets); 
                bjets_v.push_back(v_jets);
                v_bjets=v_jets;
                TLorentzVector *lepT, *bjetsT;    
                lepT= &v_lep;
                bjetsT= &v_bjets;
                double test;
                NeutrinoSolver NS(lepT, bjetsT);
            metp4 = TLorentzVector(NS.GetBest(v_met.X(), v_met.Y(), 1., 1., 0.,test));
            if (test == -1)continue;
                }
              if(!isBTagged){//Fix me, I have some some events with isBTagged==0 but csv>0.800
                temp1.push_back(v_jets);
                n_ljets_v.push_back(n_pat_jets);
                ljets_v.push_back(v_jets);
                v_ljets=v_jets;
                }
              }
            }
              float PtJets=0,massJets=0,EJets=0;
              for (unsigned int i=0;i<jet_vector.size();i++)
              {
              if (jet_vector.size()>0){
              v_jetAll[i]=jet_vector[i];
              PtJets+=v_jetAll[i].Pt();
              massJets+=v_jetAll[i].M();
              EJets+=v_jetAll[i].E();
                }
              }
              if (is_mu)if (n_elec_v.size() !=0 || n_elec15_v.size() !=0 )continue;
              if (!is_mu)if (n_muon_v.size() !=0 || n_muon15_v.size() !=0)continue;
              h_events_eachCut->Fill(2);
              if(n_lep_v.size() != 1 || n_lep15_v.size() !=1)continue;

              h_events_eachCut->Fill(3);
double trans_m_w=sqrt(pow(lep_vector[0].Pt() + v_met.Pt(), 2) - pow(lep_vector[0].Px() + v_met.Px(), 2) - pow(lep_vector[0].Py() + v_met.Py(), 2));
              if (trans_m_w < 50.)continue;
              if(no_jets.size() < 4)continue;
              h_events_eachCut->Fill(4);
              h_events_eachCut->Fill(5);
              h_bDisct_CSVv2aft->Fill(csv);
              if (n_bjets_v.size() < 2)continue;
              h_events_eachCut->Fill(6);
              if(n_ljets_v.size() < 2)continue;
              h_events_eachCut->Fill(7);
              h_EvtInfo_NumVtx->Fill(EvtInfo_NumVtx);
              h_EvtInfo_NumVtx_w->Fill(EvtInfo_NumVtx,w);
              h_PU_npT->Fill(PU_npT,w);
              h_n_ljets_aft->Fill(ljets_v.size(),w);
              h_n_bjets_aft->Fill(n_pat_bjets,w);
              h_n_cjets_aft->Fill(ncjets,w);

      // Trig SF
      int mu_id=-999;
      int e_id=-999;
      for (unsigned int i=0; i<St03Id->size(); ++i){ 
        if (fabs(St03Id->at(i)) == 13) mu_id=St03Id->at(i);
        if (fabs(St03Id->at(i)) == 11) e_id=St03Id->at(i);
         }
<<<<<<< HEAD
      TString trigSFurl("/user/mgul/Higgs_tottbar/anlyzer765/CMSSW_7_6_5/src/backgourd/ttbar/TriggerSF_v3.root");
      TString mutrigSFurl("/user/mgul/Higgs_tottbar/anlyzer765/CMSSW_7_6_5/src/backgourd/ttbar/SingleMuonTrigger_Z_RunCD_Reco76X_Feb15.root");
      TString eselurl("/user/mgul/Higgs_tottbar/anlyzer765/CMSSW_7_6_5/src/backgourd/ttbar/CutBasedID_TightWP_76X_18Feb.txt_SF2D.root");
//      TString trigSFurl("/user/mgul/Higgs_tottbar/anlyzer765/CMSSW_7_6_5/src/data/singleMuon/singleMuonC/leptonEfficiencies.root");
=======
      TString eltrigSFurl(era+"TriggerSF_v3.root");
      TString mutrigSFurl(era+"SingleMuonTrigger_Z_RunBCD_prompt80X_7p65.root");
      TString eselurl;
      if (elrun_below_273726)eselurl=era+"egammaEffiID_SF2D_below_273726.root";
      if (elrun_above_273726)eselurl=era+"egammaEffiID_SF2D_above_273726.root";
>>>>>>> 22e0d1f5b734148c99f2f8f885198e7750e4829b
      std::map<TString,TH2 *> lepEffH;
      std::vector<float> lepTriggerSF(3,1.0),lepSelSF(3,1.0);
      if(!realdata)
        {
        float trigSF(1.0), trigSFUnc(0.03);
        TFile *eselfIn=TFile::Open(eselurl);
        lepEffH["e_sel"]=(TH2 *)eselfIn->Get("EGamma_SF2D");
        for(auto& it : lepEffH) it.second->SetDirectory(0);
        eselfIn->Close();
        TString prefix("m");
        if (!is_mu)prefix="e";
        if(lepEffH.find(prefix+"_sel")!=lepEffH.end())
        {
          for(UInt_t il=0; il<TMath::Min((UInt_t)1,(UInt_t)n_lep_v.size()); il++)
          {
            float minEtaForEff( lepEffH[prefix+"_sel"]->GetXaxis()->GetXmin() ), maxEtaForEff( lepEffH[prefix+"_sel"]->GetXaxis()->GetXmax()-0.01 );
            float etaForEff=TMath::Max(TMath::Min(float(fabs(v_lep.Eta())),maxEtaForEff),minEtaForEff);
            Int_t etaBinForEff=lepEffH[prefix+"_sel"]->GetXaxis()->FindBin(etaForEff);
            
            float minPtForEff( lepEffH[prefix+"_sel"]->GetYaxis()->GetXmin() ), maxPtForEff( lepEffH[prefix+"_sel"]->GetYaxis()->GetXmax()-0.01 );
            float ptForEff=TMath::Max(TMath::Min(float(v_lep.Pt()),maxPtForEff),minPtForEff);
            Int_t ptBinForEff=lepEffH[prefix+"_sel"]->GetYaxis()->FindBin(ptForEff);

            float selSF(lepEffH[prefix+"_sel"]->GetBinContent(etaBinForEff,ptBinForEff));
            float selSFUnc(lepEffH[prefix+"_sel"]->GetBinError(etaBinForEff,ptBinForEff));
            lepSelSF[0]*=selSF;      lepSelSF[1]*=(selSF-selSFUnc);       lepSelSF[2]*=(selSF+selSFUnc);
<<<<<<< HEAD
//cout<<"minEtaForEff,  "<<minEtaForEff<<",  etaBinForEff, "<<etaBinForEff<<",   minPtForEff,   "<<minPtForEff<<",  selSF,  "<<selSF<<endl;
          }
        }
    if (!is_mu){
        TFile *fIn=TFile::Open(trigSFurl);
=======
//        cout<<"lepSelSF[0]:  "<<lepSelSF[0]<<",  
          }
        }
    if (!is_mu){
        TFile *fIn=TFile::Open(eltrigSFurl);
>>>>>>> 22e0d1f5b734148c99f2f8f885198e7750e4829b
        lepEffH["Ele23_WPLoose_Gsf"]=(TH2 *)fIn->Get("Ele23_WPLoose_Gsf");
        for(auto& it : lepEffH) it.second->SetDirectory(0);
        fIn->Close();
        if(lepEffH.find("Ele23_WPLoose_Gsf")!=lepEffH.end()){
          for(UInt_t il=0; il<TMath::Min((UInt_t)1,(UInt_t)n_lep_v.size()); il++)
            {
//            Int_t ilIdx=n_lep_v[il];
            float minEtaForEff(lepEffH["Ele23_WPLoose_Gsf"]->GetYaxis()->GetXmin()), maxEtaForEff( lepEffH["Ele23_WPLoose_Gsf"]->GetYaxis()->GetXmax()-0.01 );
            float etaForEff=TMath::Max(TMath::Min(float(fabs(v_lep.Eta())),maxEtaForEff),minEtaForEff);
            Int_t etaBinForEff=lepEffH["Ele23_WPLoose_Gsf"]->GetYaxis()->FindBin(etaForEff);

            float minPtForEff( lepEffH["Ele23_WPLoose_Gsf"]->GetXaxis()->GetXmin() ), maxPtForEff( lepEffH["Ele23_WPLoose_Gsf"]->GetXaxis()->GetXmax()-0.01 );
            float ptForEff=TMath::Max(TMath::Min(float(v_lep.Pt()),maxPtForEff),minPtForEff);
            Int_t ptBinForEff=lepEffH["Ele23_WPLoose_Gsf"]->GetXaxis()->FindBin(ptForEff);

            trigSF=(lepEffH["Ele23_WPLoose_Gsf"]->GetBinContent(ptBinForEff, etaBinForEff));
            trigSFUnc=(lepEffH["Ele23_WPLoose_Gsf"]->GetBinError(etaBinForEff,ptBinForEff));
//            float trigSF(1.0), trigSFUnc(0.03);
            lepTriggerSF[0]*=trigSF; lepTriggerSF[1]*=(trigSF-trigSFUnc); lepTriggerSF[2]*=(trigSF+trigSFUnc);
              }
            }
          }
      if (is_mu){
         TFile *mufIn=TFile::Open(mutrigSFurl);
<<<<<<< HEAD
         mufIn->cd("runD_IsoMu20_OR_IsoTkMu20_HLTv4p3_PtEtaBins");
         gDirectory->Get("runD_IsoMu20_OR_IsoTkMu20_HLTv4p3_PtEtaBins");
=======
         if (run_273158_274093)mufIn->cd("IsoMu20_OR_IsoTkMu20_PtEtaBins_Run273158_to_274093");
         if (run_274094_276097)mufIn->cd("IsoMu20_OR_IsoTkMu20_PtEtaBins_Run274094_to_276097");
         gDirectory->Get("efficienciesDATA");
>>>>>>> 22e0d1f5b734148c99f2f8f885198e7750e4829b
         lepEffH["pt_abseta_ratio"]=(TH2 *)gDirectory->Get("pt_abseta_ratio");
         for(auto& it : lepEffH) it.second->SetDirectory(0);
         mufIn->Close();
          if(lepEffH.find("pt_abseta_ratio")!=lepEffH.end()){
            for(UInt_t il=0; il<TMath::Min((UInt_t)1,(UInt_t)n_lep_v.size()); il++)
            {
              float minEtaForEff(lepEffH["pt_abseta_ratio"]->GetYaxis()->GetXmin()), maxEtaForEff( lepEffH["pt_abseta_ratio"]->GetYaxis()->GetXmax()-0.01 );
              float etaForEff=TMath::Max(TMath::Min(float(fabs(v_lep.Eta())),maxEtaForEff),minEtaForEff);
              Int_t etaBinForEff=lepEffH["pt_abseta_ratio"]->GetYaxis()->FindBin(etaForEff);

              float minPtForEff( lepEffH["pt_abseta_ratio"]->GetXaxis()->GetXmin() ), maxPtForEff( lepEffH["pt_abseta_ratio"]->GetXaxis()->GetXmax()-0.01 );
              float ptForEff=TMath::Max(TMath::Min(float(v_lep.Pt()),maxPtForEff),minPtForEff);
              Int_t ptBinForEff=lepEffH["pt_abseta_ratio"]->GetXaxis()->FindBin(ptForEff);

              trigSF=(lepEffH["pt_abseta_ratio"]->GetBinContent(ptBinForEff, etaBinForEff));
              trigSFUnc=(lepEffH["pt_abseta_ratio"]->GetBinError(etaBinForEff,ptBinForEff));
              lepTriggerSF[0]*=trigSF; lepTriggerSF[1]*=(trigSF-trigSFUnc); lepTriggerSF[2]*=(trigSF+trigSFUnc);

            }
            }
          }
         }
    if(!realdata )w*=lepTriggerSF[0]*lepSelSF[0];

//  HadW reconstruction
    TLorentzVector hadWT_; 
    double mmin=9999;
    unsigned int ljet_ind1=99, ljet_ind2=99;
    if (csv < 0.800){
    for (unsigned int i=0; i<n_ljets_v.size(); i++){
      for (unsigned int j=i+1; j<n_ljets_v.size(); j++){
      if (i==j )continue;
        double mmin_temp = (ljets_v[i] + ljets_v[j]).M(); 
        if (fabs(mmin_temp-80.4) < fabs (mmin-80.4)){
          mmin=mmin_temp; ljet_ind1=i;ljet_ind2=j;
        hadWT_=ljets_v[i] + ljets_v[j];
//        cout<<"this is csv:  "<<csv<<endl;
      cout<<"csv: "<<csv<<", i: "<<i<<", J:  "<<j<<", n_ljets_v.size:  "<<n_ljets_v.size()<<",  ljets_v[i].M: "<<ljets_v[i].M()<<",  ljets_v[j]).M(), "<<ljets_v[j].M()<<endl;
        }
      }
    }
  }
    cout<<"ljet_ind1: "<<ljet_ind1<<",  ljet_ind2:  "<<ljet_ind2<<",  hadWT_.M:  "<<hadWT_.M()<<endl;



//----------------------------Muon parameters--------------------------------------------//
 h_met_Pz->Fill(metp4.Pz(),w);
float        Pt_lep = v_lep.Pt();
        h_Pt_lep->Fill(Pt_lep,w);
//        Eta_lep = v_lep.Eta();
        h_Eta_lep->Fill(v_lep.Eta(),w);
//        Phi_lep = v_lep.Phi();
        h_Phi_lep->Fill(v_lep.Phi(),w);
//        E_lep = v_lep.E();
        h_E_lep->Fill(v_lep.E(),w);
        h_patMuonPfIsoDbeta_after->Fill(muon_TIso,w);
//-----------------------------Jets param------------------------------------------------//
//--------------------All jets-------------------------
        h_pfJet_cmult->Fill(jet_cMult,w);
        h_PtJets_afterCut->Fill(PtJets,w);
        h_massJets_afterCut->Fill(massJets,w);
        h_Ejets_afterCut->Fill(EJets,w);
//-------------------------------
        float PtJet1 = v_jetAll[0].Pt();
        h_PtJet1->Fill(PtJet1,w);
        float EtaJet1 = v_jetAll[0].Eta();
        h_EtaJet1->Fill(EtaJet1,w);
        float PhiJet1 = v_jetAll[0].Phi();
        h_PhiJet1->Fill(PhiJet1,w);
        float EJet1 = v_jetAll[0].E();
        h_EJet1->Fill(EJet1,w);
        float MJet1 = v_jetAll[0].M();
        h_MJet1->Fill(MJet1,w);

//------------------------------------------
        float PtJet2 = v_jetAll[1].Pt();
        h_PtJet2->Fill(PtJet2,w);
        float EtaJet2 = v_jetAll[1].Eta();
        h_EtaJet2->Fill(EtaJet2,w);
        float PhiJet2 = v_jetAll[1].Phi();
        h_PhiJet2->Fill(PhiJet2,w);
        float EJet2 = v_jetAll[1].E();
        h_EJet2->Fill(EJet2,w);
        float MJet2 = v_jetAll[1].M();
        h_MJet2->Fill(MJet2,w);
//------------------------------------------
        float PtJet3 = v_jetAll[2].Pt();
        h_PtJet3->Fill(PtJet3,w);
        float EtaJet3 = v_jetAll[2].Eta();
        h_EtaJet3->Fill(EtaJet3,w);
        float PhiJet3 = v_jetAll[2].Phi();
        h_PhiJet3->Fill(PhiJet3,w);
        float EJet3 = v_jetAll[2].E();
        h_EJet3->Fill(EJet3,w);
        float MJet3 = v_jetAll[2].M();
        h_MJet3->Fill(MJet3,w);
//------------------------------------------
        float PtJet4 = v_jetAll[3].Pt();
        h_PtJet4->Fill(PtJet4,w);
        float EtaJet4 = v_jetAll[3].Eta();
        h_EtaJet4->Fill(EtaJet4,w);
        float PhiJet4 = v_jetAll[3].Phi();
        h_PhiJet4->Fill(PhiJet4,w);
        float EJet4 = v_jetAll[3].E();
        h_EJet4->Fill(EJet4,w);
        float MJet4 = v_jetAll[3].M();
        h_MJet4->Fill(MJet4,w);
//-------------------------bjets------------------------------------------//
//        h_n_bjets->Fill(n_bjets_v.size(),w);
        float Pt_bJet1 = bjets_v[0].Pt();
        h_Pt_bJet1->Fill(Pt_bJet1,w);
        float Eta_bJet1 = bjets_v[0].Eta();
        h_Eta_bJet1->Fill(Eta_bJet1,w);
        float Phi_bJet1 = bjets_v[0].Phi();
        h_Phi_bJet1->Fill(Phi_bJet1,w);
        float E_bJet1 = bjets_v[0].E();
        h_E_bJet1->Fill(E_bJet1,w);
        float M_bJet1 = bjets_v[0].M();
        h_M_bJet1->Fill(M_bJet1,w);
//------------------------------------
//---------------------------Light jets----------------------------------//
//        h_n_ljets->Fill(n_ljets_v.size(),w);
        float Pt_lJet1 = temp1[0].Pt();
        h_Pt_lJet1->Fill(Pt_lJet1,w);
        float Eta_lJet1 =  temp1[0].Eta();
        h_Eta_lJet1->Fill(Eta_lJet1,w);
        float Phi_lJet1 =  temp1[0].Phi();
        h_Phi_lJet1->Fill(Phi_lJet1,w);
        float E_lJet1 =  temp1[0].E();
        h_E_lJet1->Fill(E_lJet1,w);
        float M_lJet1 =  temp1[0].M();
        h_M_lJet1->Fill(M_lJet1,w);

if (metp4.Pt()>0)h_met_E->Fill(metp4.Pt(),w);
//  exp systematics

    if (run_sys)
    {
    for(size_t ivar=0; ivar<expSysts.size(); ivar++)
    {
     TString varName=expSysts[ivar]; 
     bool updateBtag(varName.Contains("tagEff"));
     bool updateJER(varName.Contains("UncJER"));
     for(int isign=0; isign<2; isign++)
      {
      float newWgt(w);
      if(varName=="Pileup" && puWeight[0]!=0)
        newWgt*=(isign==0 ? puWeight[1]/puWeight[0] : puWeight[2]/puWeight[0]);
      if(varName=="Trigger")
        newWgt *= (isign==0 ? lepTriggerSF[1]/lepTriggerSF[0] : lepTriggerSF[2]/lepTriggerSF[0]);
      if(varName=="MuEfficiency" || "EleEfficiency")
        newWgt *= (isign==0 ? lepSelSF[1]/lepSelSF[0] : lepSelSF[2]/lepSelSF[0]);
      //jets

      std::vector<TLorentzVector> varBJets,varLightJets;
      TLorentzVector jetDiff(0,0,0,0),jetSum(0,0,0,0);
      if( updateBtag  ){   //Fix Me
        varBJets=bjets_v;
        varLightJets=ljets_v;
      }
      else
      {
      for (unsigned int ij=0; ij<patJetPfAk04PtJERSmear->size();ij++){
      int k =no_jets[ij];
      int jflav( abs(patJetPfAk04PartonFlavour_->at(ij)) );
      TLorentzVector jp4;
      jp4.SetPtEtaPhiE(patJetPfAk04PtJERSmear->at(ij),patJetPfAk04Eta_->at(ij),patJetPfAk04Phi_->at(ij),patJetPfAk04En_->at(ij));
      jetDiff -= jp4;

        jetDiff -= jp4;
        jetSum += jp4;
        if(jp4.Pt()<30) continue;
        if(fabs(jp4.Eta()) > 2.4) continue;//Eta value Fix Me 

              if (is_mu)if (n_elec_v.size() !=0 || n_elec15_v.size() !=0 )continue;
              if (!is_mu)if (n_muon_v.size() !=0 || n_muon15_v.size() !=0)continue;
              if(n_lep_v.size() != 1 || n_lep15_v.size() !=1)continue;

double trans_m_w=sqrt(pow(lep_vector[0].Pt() + v_met.Pt(), 2) - pow(lep_vector[0].Px() + v_met.Px(), 2) - pow(lep_vector[0].Py() + v_met.Py(), 2));
              if (trans_m_w < 50.)continue;
              if(no_jets.size() < 4)continue;
              if (n_bjets_v.size() < 1)continue;
              if(n_ljets_v.size() < 2)continue;

        float csv_new = csv_v[0];  // Fix me it is 0 or k
        bool isBTagged(csv_new>0.800);
        if(!realdata){
          float jptForBtag(bjets_v[0].Pt()>1000. ? 999. : bjets_v[0].Pt()), jetaForBtag(fabs(bjets_v[0].Eta()));
          float expEff(1.0), jetBtagSF(1.0);
          if(jflav==4){
            expEff        = expBtagEff["c"]->Eval(jptForBtag);
            int idx(0);
            if(varName=="CtagEff")   idx=(isign==0 ? 1 : 2);
            jetBtagSF  = sfbReaders[idx]->eval( BTagEntry::FLAV_C, jetaForBtag, jptForBtag);
            jetBtagSF *= expEff>0 ? 1 : 0. ;
            }
          else if(jflav==5){
            expEff=expBtagEff["b"]->Eval(jptForBtag);
            int idx(0);
            if(varName=="BtagEff")   idx=(isign==0 ? 1 : 2);
            jetBtagSF=sfbReaders[idx]->eval( BTagEntry::FLAV_B, jetaForBtag, jptForBtag);
            jetBtagSF *= expEff>0 ? 1 : 0. ;
            }
          else{
            expEff=expBtagEff["udsg"]->Eval(jptForBtag);
            int idx(0);
            if(varName=="LtagEff")   idx=(isign==0 ? 1 : 2);
            jetBtagSF=sflReaders[idx]->eval( BTagEntry::FLAV_UDSG, jetaForBtag, jptForBtag);
            jetBtagSF *= expEff>0 ? 1 : 0. ;
          }
          myBTagSFUtil.modifyBTagsWithSF(isBTagged,    jetBtagSF,     expEff);
          }
          if(isBTagged) varBJets.push_back(bjets_v[0]);
          else          varLightJets.push_back(bjets_v[0]);
          }
        }
        if (metp4.Pt()>0)plots2d["metptshapes_exp"]->Fill(metp4.Pt(),2*ivar+isign,newWgt);
        if (varBJets.size() > 0.)plots2d["nbjetshapes_exp"]->Fill(varBJets.size(),2*ivar+isign,newWgt);
        if (varBJets.size() > 0.)plots2d["bjetptshapes_exp"]->Fill(varBJets[0].Pt(),2*ivar+isign,newWgt);
        if (varLightJets.size() > 0.)plots2d["nljetshapes_exp"]->Fill(varLightJets.size(),2*ivar+isign,newWgt);
        if (varLightJets.size() > 0.)plots2d["ljetptshapes_exp"]->Fill(varLightJets[0].Pt(),2*ivar+isign,newWgt);
        if (ncjets > 0.)plots2d["ncjetshapes_exp"]->Fill(ncjets,2*ivar+isign,newWgt);
      }
    }
  }


 // TTree weight_tree;
  weight_tree = new TTree("tree","tree");
  weight_tree->Fill();
        }//entries loop
        theFile->Write();
        theFile->Close();
        }//function loop
