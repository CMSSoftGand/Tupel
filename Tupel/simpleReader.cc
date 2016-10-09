#include "simpleReader.h"
#include "atest_class.h"
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
//    bool is_mu (true);
  bool is_mu (false);
  bool run_sys (false);
  TString era("/user/mgul/Higgs_tottbar/Anz_8011/CMSSW_8_0_11/src/Tupel/Tupel/data/era2016/");
  ////  standalone_LumiReWeighting puWeight(201525,0), puUp(201525,1), puDown(201525,-1);
//https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80X
          TString btagUncUrl(era+"cMVAv2_ichep.csv");
          cout<<"string is : "<<btagUncUrl<<endl;
          gSystem->ExpandPathName(btagUncUrl);
          std::vector<BTagCalibrationReader *> sfbReaders, sflReaders;
          TString btagEffExpUrl(era+"expTageff.root");
          gSystem->ExpandPathName(btagEffExpUrl);
          std::map<TString, TGraphAsymmErrors *> expBtagEff, expBtagEffPy8;
          BTagSFUtil myBTagSFUtil;
          if (!realdata){
            BTagCalibration btvcalib("cMVAv2", btagUncUrl.Data());
            sfbReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "ttbar", "central") );
            sfbReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "ttbar", "down") );
            sfbReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "ttbar", "up") );
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
     expSysts.push_back("UncJES");
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
     TH1D* h_Pt_lep = new TH1D("Pt_lep","Pt_lep",50, 20.0,200.);
     TH1D* h_Eta_lep = new TH1D("Eta_lep","Eta_lep",50, -2.6,2.6);
     TH1D* h_Phi_lep = new TH1D("Phi_lep","Phi_lep",50,-3.4,3.4);
     TH1D* h_E_lep = new TH1D("E_lep","E_lep",50, 0.0,400.);
//-----------------------------Lep W  ------------------------------------------------//
     TH1D* h_Pt_lepW = new TH1D("Pt_LepW","Pt_LepW",50, 0.0,400.);
     TH1D* h_M_lepW = new TH1D("M_LepW","M_LepW",50, 0.0,200.0);
//------------------------------ Jets ------------------------------------------------//
     TH1D* h_no_Jets = new TH1D("No_ofJets","No_ofJets",7, 0.5,7.5);
     TH1F** h_PtJet=new TH1F*[5];TH1F** h_EtaJet=new TH1F*[5];TH1F** h_PhiJet=new TH1F*[5];TH1F** h_EJet=new TH1F*[5];TH1F** h_MJet=new TH1F*[5];
     char nam_pt[100];char nam_eta[100];char nam_phi[100];char nam_energy[100];char nam_mass[100];
     for (int i=1; i<5; i++){
       sprintf(nam_pt,"PtJet%i",i);
       sprintf(nam_eta,"EtaJet%i",i);
       sprintf(nam_phi,"PhiJet%i",i);
       sprintf(nam_energy,"EJet%i",i);
       sprintf(nam_mass,"MJet%i",i);
     h_PtJet[i] = new TH1F(nam_pt,nam_pt,50, 0.0,300.);
     h_EtaJet[i] = new TH1F(nam_eta,nam_eta,50, -2.5,2.5);
     h_PhiJet[i] = new TH1F(nam_phi,nam_phi,50, -3.4,3.4);
     h_EJet[i] = new TH1F(nam_energy,nam_energy,50, 0.0,700.);
     h_MJet[i] = new TH1F(nam_mass,nam_mass,50, 0.0,100.);
      }
//------------------------------bJets-------------------------------------------------//
  TH1D* h_n_bjets = new TH1D("n_bjets","n_bjets",5, 0.5,5.5);
  TH1D* h_n_bjets_aft = new TH1D("n_bjets_aft","n_bjets_aft",5, 1,6);
  TH1D* h_n_cjets_aft = new TH1D("n_cjets_aft","n_cjets_aft",5, 1,6);
  TH1F** h_Pt_bJet=new TH1F*[3];TH1F** h_Eta_bJet=new TH1F*[3];TH1F** h_Phi_bJet=new TH1F*[3];TH1F** h_E_bJet=new TH1F*[3];TH1F** h_M_bJet=new TH1F*[3];
  for (int i=1; i<3; i++){
    sprintf(nam_pt,"Pt_bJet%i",i);
    sprintf(nam_eta,"Eta_bJet%i",i);
    sprintf(nam_phi,"Phi_bJet%i",i);
    sprintf(nam_energy,"E_bJet%i",i);
    sprintf(nam_mass,"M_bJet%i",i);
    h_Pt_bJet[i] = new TH1F(nam_pt,nam_pt,50, 0.0,250.);
    h_Eta_bJet[i] = new TH1F(nam_eta,nam_eta,50, -2.5,2.5);
    h_Phi_bJet[i] = new TH1F(nam_phi,nam_phi,50, -3.4,3.4);
    h_E_bJet[i] = new TH1F(nam_energy,nam_energy,50, 0.0,500.);
    h_M_bJet[i] = new TH1F(nam_mass,nam_mass,50, 0.0,70.);
    }
    TH1D* h_DPhi_bj12 = new TH1D("DPhi_bj12","DPhi_bj12",50, 3.14,3.14);
    TH1D* h_DR_bj12 = new TH1D("DR_bj12","DR_bj12",50, 0.0,10);
//------------------------------Light jets-------------------------------------------//
  TH1D* h_n_ljets = new TH1D("n_ljets","n_ljets",7, 0.5,7.5);
  TH1D* h_n_ljets_aft = new TH1D("n_ljets_aft","n_ljets_aft",5, 1,6);
  TH1F** h_Pt_lJet=new TH1F*[3];TH1F** h_Eta_lJet=new TH1F*[3];TH1F** h_Phi_lJet=new TH1F*[3];TH1F** h_E_lJet=new TH1F*[3];TH1F** h_M_lJet=new TH1F*[3];
  for (int i=1; i<3; i++){
    sprintf(nam_pt,"Pt_lJet%i",i);
    sprintf(nam_eta,"Eta_lJet%i",i);
    sprintf(nam_phi,"Phi_lJet%i",i);
    sprintf(nam_energy,"E_lJet%i",i);
    sprintf(nam_mass,"M_lJet%i",i);
    h_Pt_lJet[i] = new TH1F(nam_pt,nam_pt,50, 0.0,250.);
    h_Eta_lJet[i] = new TH1F(nam_eta,nam_eta,50, -2.5,2.5);
    h_Phi_lJet[i] = new TH1F(nam_phi,nam_phi,50, -3.4,3.4);
    h_E_lJet[i] = new TH1F(nam_energy,nam_energy,50, 0.0,350.);
    h_M_lJet[i] = new TH1F(nam_mass,nam_mass,50, 0.0,70.);
    }
//------------------------------Had W -----------------------------------------------//
  TH1D* h_Pt_hadW = new TH1D("Pt_hadW","Pt_hadW",50, 0.0,400.);
          TH1D* h_M_hadW = new TH1D("M_hadW","M_hadW",50, 0.0,400.0);
//------------------------------Leptonic Top ----------------------------------------//
  TH1D* h_M_lepTop  = new TH1D("M_lepTop","M_lepTop",50, 0.0,500.0);
    TH1D* h_Pt_lepTop = new TH1D("Pt_lepTop","Pt_lepTop",50, 0.0,600.0);
//------------------------------Hadronic Top----------------------------------------//
  TH1D* h_M_hadTop  = new TH1D("M_hadTop","M_hadTop",50, 0.0,500.0);
    TH1D* h_Pt_hadTop  = new TH1D("Pt_hadTop","Pt_hadTop",50, 0.0,600.0);
//------------------------------Heavy Higgs----------------------------------------//
  TH1D* h_M_tt   = new TH1D("Mass_H","Mass_H",100, 200.0,1000.0);
    TH1D* h_Pt_tt   = new TH1D("Pt_H","Pt_H",50, 0.0,500.0);
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
                    TH1D* h_bDisct_cMVAv2 = new TH1D("bDiscMVAv2","bDiscMVAv2",50, 0.0,1.0);
                      TH1D* h_bDisct_cMVAv2aft = new TH1D("cMVA_v2","cMVA_v2",50, 0.0,1.0);
                        TH1D* h_pfJet_cmult = new TH1D("Jet_cMult","Jet_cMult",10, 0.0,40.0);
                        TH1D* h_NuDiscr = new TH1D("NuDiscr","NuDiscr()",25, 0.0,10.0);
                        TH1D* h_MassDiscr = new TH1D("MassDiscr","MassDiscr()",25, 0.0,20.0);

    std::map<TString, TH2 *> plots2d;
    plots2d["metptshapes_exp"]= new TH2F("metptshapes_exp", ";Missing transverse energy [GeV];Events" , 100,0.,200., 2*nExpSysts,0,2*nExpSysts);
    plots2d["nbjetshapes_exp"]= new TH2F("nbjetshapes_exp", "; bJets Multiplicity ;Events" , 5,1,6, 2*nExpSysts,0,2*nExpSysts);
    plots2d["bjetptshapes_exp"]= new TH2F("bjetptshapes_exp", "; bJets pt ;Events" , 50,0.,250., 2*nExpSysts,0,2*nExpSysts);
    plots2d["nljetshapes_exp"]= new TH2F("nljetshapes_exp", "; light Jets Multiplicity ;Events" , 5,1,6, 2*nExpSysts,0,2*nExpSysts);
    plots2d["ljetptshapes_exp"]= new TH2F("ljetptshapes_exp", "; light Jets pt ;Events" , 50,0,300, 2*nExpSysts,0,2*nExpSysts);
    plots2d["ncjetshapes_exp"]= new TH2F("ncjetshapes_exp", "; c Jets Multiplicity ;Events" , 5,1,6, 2*nExpSysts,0,2*nExpSysts);
    plots2d["jetptshapes_exp"]= new TH2F("jetptshapes_exp", "; Jets pt ;Events" , 50,0.,300., 2*nExpSysts,0,2*nExpSysts);
    for(Int_t isyst=0; isyst<nExpSysts; isyst++){
      for(Int_t ivar=0; ivar<2; ivar++){
        TString label(expSysts[isyst] + (ivar==0 ? "Down" : "Up"));
        plots2d["metptshapes_exp"]->GetYaxis()->SetBinLabel(2*isyst+ivar+1, label);
        plots2d["nbjetshapes_exp"]->GetYaxis()->SetBinLabel(2*isyst+ivar+1, label);
        plots2d["bjetptshapes_exp"]->GetYaxis()->SetBinLabel(2*isyst+ivar+1, label);
        plots2d["nljetshapes_exp"]->GetYaxis()->SetBinLabel(2*isyst+ivar+1, label);
        plots2d["ljetptshapes_exp"]->GetYaxis()->SetBinLabel(2*isyst+ivar+1, label);
        plots2d["ncjetshapes_exp"]->GetYaxis()->SetBinLabel(2*isyst+ivar+1, label);
        plots2d["jetptshapes_exp"]->GetYaxis()->SetBinLabel(2*isyst+ivar+1, label);
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
	Int_t nentries(t->GetEntriesFast());
//        nentries=100000;
        for (int jentry=0; jentry < nentries; jentry++)
        {
        t->GetEntry(jentry);
        if(jentry%10000==0)cout<<" << "<<jentry<<"/"<<nentries<<endl;
        double w=1;
/*        float iSecret;
        srand (time(NULL));
        iSecret = rand() % 100 + 1;
        bool run_273158_274093 (iSecret <= 5.);
        bool run_274094_276097 (5. < iSecret);// && iSecret < 63.);//Fix me when available the rest
        bool elrun_below_273726 (iSecret <= 4.);
        bool elrun_above_273726 (4. < iSecret );//&& iSecret < 42.);//Fix me when available the rest
*/
        std::vector<TGraph *>puWgtGr;
        std::vector<float> puWeight(3,1.0);
        TString puWgtUrl(era+"pileupWgts.root");
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
        h_events_eachCut->Fill(0);
	if (first_PV!=1)continue;
        if (is_mu) if (HLT_IsoMu22!=1 && HLT_IsoTkMu22!=1)continue;
        if (!is_mu)if (HLT_Ele32_eta2p1_WPTight_Gsf!=1)continue;
        h_events_eachCut->Fill(1);
        //noise
//cout<<"these are filters: "<<Flag_HBHENoiseFilter<<",  "<<Flag_HBHENoiseIsoFilter<<",  "<<Flag_CSCTightHalo2015Filter<<endl;
//          if(Flag_HBHENoiseFilter!=1 || Flag_HBHENoiseIsoFilter!=1 || Flag_CSCTightHalo2015Filter!=1)continue;
//          if(Flag_HBHENoiseFilter!=1)continue;

          if( Flag_HBHENoiseIsoFilter!=1)continue;
          if( Flag_globalTightHalo2016Filter!=1)continue;
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
      if (patElecPt_->at(elec) > 35. && fabs(patElecScEta_->at(elec)) <2.1 && patElecIdtight_->at(elec)>0 ){


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
          if(patMuonPt_->at(mu) > 26. && fabs(patMuonEta_->at(mu))<2.4 &&
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
            vector<float> cMVA_v;
            double cMVA=-99,cMVAv2=-999, jet_pt=0.,nom=1., denom=1.,ind_jet_const=0.;
            vector <TLorentzVector> temp1;
            vector <TLorentzVector> v_ljetsAll;
            TLorentzVector lorenz_v; 
	    TLorentzVector lorenz_v_raw;
            vector<TLorentzVector>jet_vector;
            vector<double> JER_Uncer;
            TLorentzVector metp4; 
	    int looseId=0;
            for (unsigned int pf=0;pf < patJetPfAk04PtJERSmear->size();++pf){
lorenz_v_raw.SetPtEtaPhiE(patJetPfAk04RawPt_->at(pf),patJetPfAk04Eta_->at(pf),patJetPfAk04Phi_->at(pf),patJetPfAk04RawEn_->at(pf));
              double Up = patJetPfAk04PtJERSmearUp->at(pf);//Fix Me I will be after or before the cut
              double Dn = patJetPfAk04PtJERSmearDn->at(pf);
              double Cn = patJetPfAk04PtJERSmear->at(pf);
              double sig_Up=abs(Up-Cn)/Cn;
              double sig_Dn=abs(Dn-Cn)/Cn;
              JER_Uncer.push_back(max(sig_Up,sig_Dn));
//cout<<"tis is jet loose id:  "<<patJetPfAk04LooseId_->at(pf)<<endl;
//else looseId=0;
              if(patJetPfAk04PtJERSmear->at(pf) > 30.&& fabs(patJetPfAk04Eta_->at(pf)) < 2.4 && patJetPfAk04LooseId_->at(pf)>0){
              n_pat_jets++;
              TLorentzVector v_jetsTemp;
              v_jetsTemp.SetPtEtaPhiE(patJetPfAk04PtJERSmear->at(pf),patJetPfAk04Eta_->at(pf),patJetPfAk04Phi_->at(pf),patJetPfAk04En_->at(pf));

              if (n_lep_v.size()>0)DR_mu_j= DeltaR(lep_vector[0].Eta(), v_jetsTemp.Eta(), lep_vector[0].Phi(), v_jetsTemp.Phi());
              if(DR_mu_j<0.4)continue;
              v_jets.SetPtEtaPhiE(patJetPfAk04PtJERSmear->at(pf),patJetPfAk04Eta_->at(pf),patJetPfAk04Phi_->at(pf),patJetPfAk04En_->at(pf));
              jet_vector.push_back(v_jets);
              h_bDisct_cMVAv2->Fill(patJetPfAk04BDiscpfCMVA_->at(pf),w);
              if(patJetPfAk04BDiscpfCMVA_->at(pf) > cMVA) cMVA = patJetPfAk04BDiscpfCMVA_->at(pf);
	            cMVAv2=patJetPfAk04BDiscpfCMVA_->at(pf);
              cMVA_v.push_back(patJetPfAk04BDiscpfCMVA_->at(pf));
              jet_cMult=patJetPfAk04cmult_->at(pf);
              h_no_Jets->Fill(n_pat_jets,w);
              no_jets.push_back(n_pat_jets);
              jet_pt=patJetPfAk04PtJERSmear->at(pf); 
              int jflav( abs(patJetPfAk04PartonFlavour_->at(pf)) );
              bool isBTagged(cMVAv2>0.185);
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
                myBTagSFUtil.modifyBTagsWithSF(isBTagged,    jetBtagSF,     expEff);//btag SF works, certified
              }
              if(isBTagged){
                n_bjets_v.push_back(n_pat_jets); 
                bjets_v.push_back(v_jets);
                v_bjets=v_jets;
                TLorentzVector *lepT, *bjetsT;    
                lepT= &v_lep;
                bjetsT= &v_bjets;
                double test;
                bool info(true);
                NeutrinoSolver NS(lepT, bjetsT);
                metp4 = TLorentzVector(NS.GetBest(v_met.X(), v_met.Y(), 1., 1., 0.,test));
                //if (test == -1 && (metp4+lep_vector[0]).M()==0.1057)continue;
                if (test == -1 )continue;
              TLorentzVector *tt1;
              tt1=&v_jets;
//              TTBarSolver ttsolver;
//              ttsolver.Solve(tt1, tt1, tt1, tt1, tt1, &metp4);
//              cout<<"this is tt1 Pt:  "<<tt1->Pt()<<",  "<<metp4.Pt()<<endl;
                }

              if(!isBTagged){//Fix me, I have some some events with isBTagged==0 but cMVA>0.185
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
              h_events_eachCut->Fill(4);
              if(no_jets.size() < 4)continue;
              h_events_eachCut->Fill(5);
              h_bDisct_cMVAv2aft->Fill(cMVA);
              if (n_bjets_v.size() < 2)continue;
              h_events_eachCut->Fill(6);
              if(n_ljets_v.size() < 2)continue;
              h_events_eachCut->Fill(7);
              h_PU_npT->Fill(PU_npT,w);
              h_n_ljets_aft->Fill(ljets_v.size(),w);
              h_n_bjets_aft->Fill(n_pat_bjets,w);
              h_n_cjets_aft->Fill(ncjets,w);
              if (n_bjets_v.size() < 2)continue;
              h_EvtInfo_NumVtx->Fill(EvtInfo_NumVtx);
              h_EvtInfo_NumVtx_w->Fill(EvtInfo_NumVtx,w);
//----------------Rochester algorithm----------------
    int lepChar=lep_charge[0];
    double nuDisc, massDiscr, nuchi2;
    TTBarSolver solver;
    TLorentzVector *wja;
    TLorentzVector wja_new;
    TLorentzVector *wjb;
    TLorentzVector wjb_new;
    TLorentzVector *Blep_temp;
    TLorentzVector *Bhad_temp;
    TLorentzVector Bhad;
    TLorentzVector Blep;
    TLorentzVector Whad;
    TLorentzVector Thad;
    TLorentzVector Wlep;
    TLorentzVector Tlep;
    Permutation test;
    Permutation best_permutation;

    for (unsigned int i=0; i<ljets_v.size(); i++){
    for (unsigned int j=i+1;j<ljets_v.size(); j++){
      if (i==j)continue;
    for (unsigned int k=0; k<bjets_v.size(); k++){
    for (unsigned int l=k+1; l<bjets_v.size(); l++){
      if (k==l)continue;
    Permutation permutation(&ljets_v[i], &ljets_v[j], &bjets_v[l], &bjets_v[k], &v_lep, &v_met,lepChar);
    solver.Solve(permutation);
    if(permutation.Prob()  < best_permutation.Prob())
    {
      best_permutation = permutation;
      Bhad=best_permutation.bhadt();
      Blep=best_permutation.blept();
      Whad=best_permutation.WHad();
      Thad=best_permutation.THad();
      Wlep=best_permutation.WLep();
      nuDisc=best_permutation.NuDiscr();
      nuchi2=best_permutation.NuChisq();
      massDiscr=best_permutation.MassDiscr();
    }
    }}}}
    if (Thad.M()!=0){
    Tlep=Blep+metp4+lep_vector[0];
    h_M_hadW->Fill(Whad.M(),w);
    h_Pt_hadW->Fill(Whad.Pt(),w);
    h_M_hadTop->Fill(Thad.M(),w);
    h_Pt_hadTop->Fill(Thad.Pt(),w);
    h_M_lepTop->Fill((Blep+metp4+lep_vector[0]).M(),w);
    h_Pt_lepTop->Fill(Tlep.Pt(),w);
    h_M_lepW->Fill(Wlep.M(),w);
    h_Pt_lepW->Fill(Wlep.Pt(),w);
    h_M_tt->Fill((Tlep+Thad).M(),w);
    h_Pt_tt->Fill((Tlep+Thad).Pt(),w);
    }

//-----------------------------------------------
      // Trig SF
bool is_sf(true);
if (is_sf){
      int mu_id=-999;
      int e_id=-999;
      for (unsigned int i=0; i<St03Id->size(); ++i){ 
        if (fabs(St03Id->at(i)) == 13) mu_id=St03Id->at(i);
        if (fabs(St03Id->at(i)) == 11) e_id=St03Id->at(i);
         }
      TString eltrigSFurl(era+"SF_HLT_Ele32_eta2p1_WPTight_Gsf.root");
      TString mutrigSFurl(era+"SingleMuonTrigger_Z_RunBCD_prompt80X_7p65.root");
      TString eIDurl;
// for 12.9/pb https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaIDRecipesRun2#Cut_based_electron_identificatio
//for 7.6/fb https://twiki.cern.ch/twiki/bin/view/CMS/MuonWorkInProgressAndPagResults#Results_on_7_6_fb
        eIDurl=era+"egammaEff_SF2D.root";
      TString muIDurl(era+"MuonID_Z_RunBCD_prompt80X_7p65.root");
      TString muIsourl(era+"MuonIso_Z_RunBCD_prompt80X_7p65.root");
      std::map<TString,TH2 *> lepEffH;
      std::vector<float> lepTriggerSF(3,1.0);
      std::vector<float> lepidSF(3,1.0),lepisoSF(3,1.0);
      float trigSF(1.0), trigSFUnc(0.03);
      float idSF(1.0), idSFUnc(0.03);
      float isoSF(1.0), isoSFUnc(0.03);
      TString Trig_hist("muon_trigg");
      if (!is_mu)Trig_hist="Ele32_eta2p1_WPTight_Gsf_EffData";
      TString id_hist("muon_id");
      if(!is_mu)id_hist="elec_id";

   if (!realdata){
   if (is_mu ){
      TFile *mufIn=TFile::Open(mutrigSFurl);
      mufIn->cd("IsoMu22_OR_IsoTkMu22_PtEtaBins_Run274094_to_276097");
      gDirectory->Get("efficienciesDATA");
         lepEffH[Trig_hist]=(TH2 *)gDirectory->Get("efficienciesDATA/pt_abseta_DATA");
         for(auto& it : lepEffH) it.second->SetDirectory(0);
         mufIn->Close();

      TFile *muIDfIn=TFile::Open(muIDurl);
      muIDfIn->cd("MC_NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1");
      gDirectory->Get("efficienciesDATA");
         lepEffH[id_hist]=(TH2 *)gDirectory->Get("efficienciesDATA/abseta_pt_DATA");
         for(auto& it : lepEffH) it.second->SetDirectory(0);
         muIDfIn->Close();

      TFile *muIsofIn=TFile::Open(muIsourl);
      muIsofIn->cd("MC_NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1");
      gDirectory->Get("efficienciesDATA");
         lepEffH["muon_iso"]=(TH2 *)gDirectory->Get("efficienciesDATA/pt_abseta_DATA");
         for(auto& it : lepEffH) it.second->SetDirectory(0);
         muIsofIn->Close();
      }
   if (!is_mu){
      TFile *fIn=TFile::Open(eltrigSFurl);
      lepEffH[Trig_hist]=(TH2 *)fIn->Get("Ele32_eta2p1_WPTight_Gsf__EffData");
      for(auto& it : lepEffH) it.second->SetDirectory(0);
              fIn->Close();

      TFile *eidfIn=TFile::Open(eIDurl);
      lepEffH[id_hist]=(TH2 *)eidfIn->Get("EGamma_SF2D");
      for(auto& it : lepEffH) it.second->SetDirectory(0);
              eidfIn->Close();
      }
          if(lepEffH.find(Trig_hist)!=lepEffH.end()){
            for(UInt_t il=0; il<TMath::Min((UInt_t)1,(UInt_t)n_lep_v.size()); il++)
            {
              float minEtaForEff(lepEffH[Trig_hist]->GetYaxis()->GetXmin()), maxEtaForEff( lepEffH[Trig_hist]->GetYaxis()->GetXmax()-0.01 );
              float etaForEff=TMath::Max(TMath::Min(float(fabs(v_lep.Eta())),maxEtaForEff),minEtaForEff);
              Int_t etaBinForEff=lepEffH[Trig_hist]->GetYaxis()->FindBin(etaForEff);

              float minPtForEff( lepEffH[Trig_hist]->GetXaxis()->GetXmin() ), maxPtForEff( lepEffH[Trig_hist]->GetXaxis()->GetXmax()-0.01 );
              float ptForEff=TMath::Max(TMath::Min(float(v_lep.Pt()),maxPtForEff),minPtForEff);
              Int_t ptBinForEff=lepEffH[Trig_hist]->GetXaxis()->FindBin(ptForEff);
              trigSF=(lepEffH[Trig_hist]->GetBinContent(ptBinForEff, etaBinForEff));
              trigSFUnc=(lepEffH[Trig_hist]->GetBinError(ptBinForEff,etaBinForEff));
              lepTriggerSF[0]*=trigSF; lepTriggerSF[1]*=(trigSF-trigSFUnc); lepTriggerSF[2]*=(trigSF+trigSFUnc);
            }
      }
          if(lepEffH.find(id_hist)!=lepEffH.end()){
            for(UInt_t il=0; il<TMath::Min((UInt_t)1,(UInt_t)n_lep_v.size()); il++)
            {
              float minPtForEff(lepEffH[id_hist]->GetYaxis()->GetXmin()), maxPtForEff( lepEffH[id_hist]->GetYaxis()->GetXmax()-0.01 );
              float ptForEff=TMath::Max(TMath::Min(float(fabs(v_lep.Pt())),maxPtForEff),minPtForEff);
              Int_t ptBinForEff=lepEffH[id_hist]->GetYaxis()->FindBin(ptForEff);

              float minEtaForEff( lepEffH[id_hist]->GetXaxis()->GetXmin() ), maxEtaForEff( lepEffH[id_hist]->GetXaxis()->GetXmax()-0.01 );
              float etaForEff=TMath::Max(TMath::Min(float(v_lep.Eta()),maxEtaForEff),minEtaForEff);
              Int_t etaBinForEff=lepEffH[id_hist]->GetXaxis()->FindBin(etaForEff);
              idSF=(lepEffH[id_hist]->GetBinContent(etaBinForEff, ptBinForEff));
              idSFUnc=(lepEffH[id_hist]->GetBinError(etaBinForEff,ptBinForEff));
              lepidSF[0]*=idSF; lepidSF[1]*=(idSF-idSFUnc); lepidSF[2]*=(idSF+idSFUnc);
            }
      }
          if(lepEffH.find("muon_iso")!=lepEffH.end()){
            for(UInt_t il=0; il<TMath::Min((UInt_t)1,(UInt_t)n_lep_v.size()); il++)
            {
              float minEtaForEff(lepEffH["muon_iso"]->GetYaxis()->GetXmin()), maxEtaForEff( lepEffH["muon_iso"]->GetYaxis()->GetXmax()-0.01 );
              float etaForEff=TMath::Max(TMath::Min(float(fabs(v_lep.Eta())),maxEtaForEff),minEtaForEff);
              Int_t etaBinForEff=lepEffH["muon_iso"]->GetYaxis()->FindBin(etaForEff);

              float minPtForEff( lepEffH["muon_iso"]->GetXaxis()->GetXmin() ), maxPtForEff( lepEffH["muon_iso"]->GetXaxis()->GetXmax()-0.01 );
              float ptForEff=TMath::Max(TMath::Min(float(v_lep.Pt()),maxPtForEff),minPtForEff);
              Int_t ptBinForEff=lepEffH["muon_iso"]->GetXaxis()->FindBin(ptForEff);
              isoSF=(lepEffH["muon_iso"]->GetBinContent(ptBinForEff, etaBinForEff));
              isoSFUnc=(lepEffH["muon_iso"]->GetBinError(ptBinForEff,etaBinForEff));
              lepisoSF[0]*=isoSF; lepisoSF[1]*=(isoSF-isoSFUnc); lepisoSF[2]*=(isoSF+isoSFUnc);
            }
      }
  }
     w*=lepTriggerSF[0]*lepidSF[0]*lepisoSF[0];
}
/*     cout<<"this is only w:  "<<w<<endl;
     w*=lepTriggerSF[0];
     cout<<"this is w1:  "<<w<<endl;
     w*=lepTriggerSF[0]*lepidSF[0];
     cout<<"this is  w2:  "<<w<<endl;
     w*=lepTriggerSF[0]*lepidSF[0]*lepisoSF[0];
     cout<<"this is  w3:  "<<w<<endl;
*/



//h_M_hadW->Fill(best_permutation.WHad().M(),w);
// h_M_hadTop->Fill(had_tMass.at(0),w);
/*h_Pt_hadW->Fill(best_permutation.THad().M(),w);
h_M_lepTop->Fill(best_permutation.TLep().M(),w);
h_Pt_lepTop->Fill(best_permutation.TLep().Pt(),w);*/
/*
//  cout<<"out side,  hadWT_;  "<<hadWT_.M()<<endl;
//  lepW reconstruction
 TLorentzVector lepWT_ = mu_vector[0]+metp4;
//     cout<<"this is lep W; mu_vector.size:  "<<mu_vector.size()<<", metp4 , "<<metp4.Pt()<<",  lepWT_.M   "<<lepWT_.M()<<",  "<<lepWT_.Pt()<<endl;

*/
//----------------------------Muon parameters--------------------------------------------//
 h_met_Pz->Fill(metp4.Pz(),w);
float        Pt_lep = v_lep.Pt();
        h_Pt_lep->Fill(Pt_lep,w);
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
        for (int i=1;i<5;i++){
        h_PtJet[i]->Fill(v_jetAll[i-1].Pt(),w);
        h_EtaJet[i]->Fill(v_jetAll[i-1].Eta(),w);
        h_PhiJet[i]->Fill(v_jetAll[i-1].Phi(),w);
        h_EJet[i]->Fill(v_jetAll[i-1].E(),w);
        h_MJet[i]->Fill(v_jetAll[i-1].M(),w);
        }

//-------------------------bjets------------------------------------------//
        for (int i=1;i<3;i++){
        h_Pt_bJet[i]->Fill(bjets_v[i-1].Pt(),w);
        h_Eta_bJet[i]->Fill(bjets_v[i-1].Eta(),w);
        h_Phi_bJet[i]->Fill(bjets_v[i-1].Phi(),w);
        h_E_bJet[i]->Fill(bjets_v[i-1].E(),w);
        h_M_bJet[i]->Fill(bjets_v[i-1].M(),w);
        }
//------------------------------------
//---------------------------Light jets----------------------------------//
        for (int i=1;i<3;i++){
        h_Pt_lJet[i]->Fill(temp1[i-1].Pt(),w);
        h_Eta_lJet[i]->Fill(temp1[i-1].Eta(),w);
        h_Phi_lJet[i]->Fill(temp1[i-1].Phi(),w);
        h_E_lJet[i]->Fill(temp1[i-1].E(),w);
        h_M_lJet[i]->Fill(temp1[i-1].M(),w);
        }

if (metp4.Pt()>0)h_met_E->Fill(metp4.Pt(),w);

//  exp systematics

    if (!run_sys)continue;
   /* 
    for(size_t ivar=0; ivar<expSysts.size(); ivar++)
    {
     TString varName=expSysts[ivar]; 
     bool updateBtag(varName.Contains("tagEff"));
     bool updateJER(varName.Contains("UncJER"));
     for(int isign=0; isign<2; isign++)
      {
      float newWgt(w);
      if(varName=="Pileup" && puWeight[0]!=0){
        newWgt*=(isign==0 ? puWeight[1]/puWeight[0] : puWeight[2]/puWeight[0]);
      }
      if(varName=="Trigger")
        newWgt *= (isign==0 ? lepTriggerSF[1]/lepTriggerSF[0] : lepTriggerSF[2]/lepTriggerSF[0]);
      if(varName=="MuEfficiency" || "EleEfficiency")
        newWgt *= (isign==0 ? lepidSF[1]/lepidSF[0] : lepidSF[2]/lepidSF[0]);
      //jets
      std::vector<TLorentzVector> varBJets,varLightJets;
      TLorentzVector jetDiff(0,0,0,0),jetSum(0,0,0,0);
      vector<TLorentzVector> jp44;
      if( updateBtag  ){   //Fix Me
        varBJets=bjets_v;
        varLightJets=ljets_v;
      }
      else
      {

      TLorentzVector jp4;
      TLorentzVector jp4Up;
      TLorentzVector jp4Dn;
      int jflav (0);
      for (unsigned int ij=0; ij<no_jets.size();ij++){
      jflav= abs(patJetPfAk04PartonFlavour_->at(ij)) ;
      jp4.SetPtEtaPhiE(patJetPfAk04PtJERSmear->at(ij),patJetPfAk04Eta_->at(ij),patJetPfAk04Phi_->at(ij),patJetPfAk04En_->at(ij));
      jp4Up.SetPtEtaPhiE(patJetPfAk04PtJERSmearUp->at(ij),patJetPfAk04Eta_->at(ij),patJetPfAk04Phi_->at(ij),patJetPfAk04En_->at(ij));
      jp4Dn.SetPtEtaPhiE(patJetPfAk04PtJERSmear->at(ij),patJetPfAk04Eta_->at(ij),patJetPfAk04Phi_->at(ij),patJetPfAk04En_->at(ij));
      jp44.push_back(jp4);
      }
        jetDiff -= jp4;
        jetSum += jp4;
      if(varName!="UncJES"){
        jp4 *=(1.0+(isign==0?-1.:1.)*unc_->at(ivar));
       }
      if(varName=="UncJER"){
        jp4 = (isign==0? jp4Dn:jp4Up); 
        }
        float cMVA_new = cMVA_v[0];  // Fix me it is 0 or k
        bool isBTagged(cMVA_new>0.800);
        if(!realdata){
          float jptForBtag(jp4.Pt()>1000. ? 999. : jp4.Pt()), jetaForBtag(fabs(jp4.Eta()));
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
        if (metp4.Pt()>0)plots2d["metptshapes_exp"]->Fill(metp4.Pt(),2*ivar+isign,newWgt);
        if (varBJets.size() > 0.)plots2d["nbjetshapes_exp"]->Fill(varBJets.size(),2*ivar+isign,newWgt);
        if (varBJets.size() > 0.)plots2d["bjetptshapes_exp"]->Fill(varBJets[0].Pt(),2*ivar+isign,newWgt);
        if (varLightJets.size() > 0.)plots2d["nljetshapes_exp"]->Fill(varLightJets.size(),2*ivar+isign,newWgt);
        if (varLightJets.size() > 0.)plots2d["ljetptshapes_exp"]->Fill(varLightJets[0].Pt(),2*ivar+isign,newWgt);
        if (ncjets > 0.)plots2d["ncjetshapes_exp"]->Fill(ncjets,2*ivar+isign,newWgt);
        if (jp44.size()>0)plots2d["jetptshapes_exp"]->Fill(jp44[0].Pt(),2*ivar+isign,newWgt);
      }
    }*/

 // TTree weight_tree;
  weight_tree = new TTree("tree","tree");
  weight_tree->Fill();
        }//entries loop
        theFile->Write();
        theFile->Close();
        }//function loop
