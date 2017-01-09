#include "simpleReader.h"
using namespace std;

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
    cout<<"Running Muon Channel: = "<<is_mu<<endl;
// xsec order:    {ttbar, st_tch, st_tw, wjets,    DY,    WW,    WZ,    ZZ,   ,ttW,    ttZ,    qcd pt 15-20,       qcd pt 20-30,      qcd pt 30-50,       qcd pt 50-80,        qcd pt 80-120, qcd pt 120-170, qcd pt 170-300, qcd pt 300-470, qcd pt 470-600,     qcd pt 600-800, qcd pt 800-1000, qcd pt 1000-Inf
    float xsec[] ={831.76, 70.70, 19.55, 61526.7, 5765.4, 118.7, 47.13 ,16.523,0.2043, 0.2529, 1273190000.*0.003, 558528000.*0.0053, 139803000.*0.01182, 19222500.*0.02276, 2758420.*0.03844, 469797.*0.05362, 117989.*0.07335,7820.25*0.10196, 645.528* 0.12242, 187.109*0.13412, 32.3486*0.14552, 10.4305*0.15544}; 
//                     qcd pt20-30,        qcd pt 30-50,     , qcd pt 50-80,       qcd pt 80-120,  qcd pt 120-170,   qcd pt 170-300, qcd pt 300-Inf
    float emqcd_xsec[]={558528000.*0.0053, 139803000.*0.01182, 19222500.*0.02276, 2758420.*0.03844, 469797.*0.05362, 117989.*0.07335, 7820.25*0.10196};
  TString era("/afs/cern.ch/work/m/mgul/public/Hto_ttbar/deepCSV8024/CMSSW_8_0_24/src/Tupel/Tupel/data/moriond17/");
  ////  standalone_LumiReWeighting puWeight(201525,0), puUp(201525,1), puDown(201525,-1);
//https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80X
          TString btagUncUrl(era+"cMVAv2_ichep.csv");
          cout<<"string is : "<<btagUncUrl<<endl;
          gSystem->ExpandPathName(btagUncUrl);
          std::vector<BTagCalibrationReader *> sfbReaders,sfcReaders, sflReaders;
          TString btagEffExpUrl(era+"expTageff.root");
          gSystem->ExpandPathName(btagEffExpUrl);
          std::map<TString, TGraphAsymmErrors *> expBtagEff, expBtagEffPy8;
          BTagSFUtil myBTagSFUtil;
          if (!realdata){
vector<string> bsystypes;
bsystypes.push_back("central");bsystypes.push_back("down");bsystypes.push_back("up");
            BTagCalibration btvcalib("cMVAv2", btagUncUrl.Data());
for (unsigned int i=0; i<bsystypes.size();i++){
            sfbReaders.push_back( new BTagCalibrationReader(BTagEntry::OP_MEDIUM, bsystypes.at(i).c_str()) );
	    sfcReaders.push_back( new BTagCalibrationReader(BTagEntry::OP_MEDIUM, bsystypes.at(i).c_str()) );
            sflReaders.push_back( new BTagCalibrationReader(BTagEntry::OP_MEDIUM, bsystypes.at(i).c_str()) );
sfbReaders[i]->load(btvcalib,BTagEntry::FLAV_B,"ttbar");
sfcReaders[i]->load(btvcalib,BTagEntry::FLAV_C,"ttbar");
sflReaders[i]->load(btvcalib,BTagEntry::FLAV_UDSG,"incl");
}
            TFile *beffIn=TFile::Open(btagEffExpUrl);
            expBtagEff["b"]=(TGraphAsymmErrors *)beffIn->Get("b");
            expBtagEff["c"]=(TGraphAsymmErrors *)beffIn->Get("c");
            expBtagEff["udsg"]=(TGraphAsymmErrors *)beffIn->Get("udsg");
            beffIn->Close();
            }

      // Trig SF
      int mu_id=-999;
      int e_id=-999;
      TString eltrigSFurl(era+"SF_HLT_Ele32_eta2p1_WPTight_Gsf.root");
      TString mutrigSFurl(era+"TriggerSF_v1_IsoMu24.root");
      TString eIDurl;
// for 12.9/pb https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaIDRecipesRun2#Cut_based_electron_identificatio
//for 7.6/fb https://twiki.cern.ch/twiki/bin/view/CMS/MuonWorkInProgressAndPagResults#Results_on_7_6_fb
        eIDurl=era+"egammaEff_SF2D.root";
      TString muIDurl(era+"MuonID.root");
      TString muIsourl(era+"MuonIso.root");
      std::map<TString,TH2 *> lepEffH;
      std::vector<float> lepTriggerSF(3,1.0);
      std::vector<float> lepidSF(3,1.0),lepisoSF(3,1.0);
      float trigSF(1.0), trigSFUnc(0.03);
      float idSF(1.0), idSFUnc(0.03);
      float isoSF(1.0), isoSFUnc(0.03);
      TString Trig_hist("muon_trigg");
      if (!is_mu)Trig_hist="Ele32_eta2p1_WPTight_Gsf";
      TString id_hist("muon_id");
      if(!is_mu)id_hist="elec_id";
   if (!realdata){
   if (is_mu ){
      TFile *mufIn=TFile::Open(mutrigSFurl);
//      mufIn->cd("IsoMu22_OR_IsoTkMu22_PtEtaBins_Run274094_to_276097");
//      gDirectory->Get("efficienciesDATA");
//         lepEffH[Trig_hist]=(TH2 *)gDirectory->Get("efficienciesDATA/pt_abseta_DATA");
         lepEffH[Trig_hist]=(TH2 *)mufIn->Get("IsoMu24_OR_IsoTkMu24");
         for(auto& it : lepEffH) it.second->SetDirectory(0);
         mufIn->Close();

      TFile *muIDfIn=TFile::Open(muIDurl);
         lepEffH[id_hist]=(TH2 *)muIDfIn->Get("histo2D_MC");
         for(auto& it : lepEffH) it.second->SetDirectory(0);
         muIDfIn->Close();

      TFile *muIsofIn=TFile::Open(muIsourl);
         lepEffH["muon_iso"]=(TH2 *)muIsofIn->Get("histo2D_MC");
         for(auto& it : lepEffH) it.second->SetDirectory(0);
         muIsofIn->Close();
      }
   if (!is_mu){
      TFile *fIn=TFile::Open(eltrigSFurl);
      lepEffH[Trig_hist]=(TH2 *)fIn->Get("Ele32_eta2p1_WPTight_Gsf");
      for(auto& it : lepEffH) it.second->SetDirectory(0);
              fIn->Close();

      TFile *eidfIn=TFile::Open(eIDurl);
      lepEffH[id_hist]=(TH2 *)eidfIn->Get("EGamma_SF2D");
      for(auto& it : lepEffH) it.second->SetDirectory(0);
              eidfIn->Close();
      }
     }
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
     if (is_mu){TDirectory *mujets_2btag = theFile->mkdir("mujets_2btag");mujets_2btag->cd();}
     if (!is_mu){TDirectory *eljets_2btag = theFile->mkdir("eljets_2btag");eljets_2btag->cd();}
     TH1::SetDefaultSumw2();
     TH2::SetDefaultSumw2();
     branchAdd(t);
	   Int_t nentries(t->GetEntriesFast());
     float scale=-999999.0;
     float lumi = 12877.401701;
     if (!is_mu)lumi = 12883.8601;
        if (fin.Contains( "TTJets"))scale=(xsec[0]*lumi)/nentries;
        if (fin.Contains( "SingleT_t.root"))scale=(xsec[1]*lumi)/nentries;
        if (fin.Contains( "_tW_nohad"))scale=(xsec[2]*lumi)/nentries;
        if (fin.Contains( "WJets"))scale=(xsec[3]*lumi)/nentries;
        if (fin.Contains( "DYJets"))scale=(xsec[4]*lumi)/nentries;
        if (fin.Contains( "_WWToLNuQQ"))scale=(xsec[5]*lumi)/nentries;
        if (fin.Contains( "_WZ"))scale=(xsec[6]*lumi)/nentries;
        if (fin.Contains( "_ZZ"))scale=(xsec[7]*lumi)/nentries;
        if (fin.Contains( "TTWJetsToLNu"))scale=(xsec[8]*lumi)/nentries;
        if (fin.Contains( "TTZToLLNuNu_M-10"))scale=(xsec[9]*lumi)/nentries;
        if (fin.Contains( "_QCD_Pt_15to20"))scale=(xsec[10]*lumi)/nentries;
        if (fin.Contains( "_QCD_Pt_20to30"))scale=(xsec[11]*lumi)/nentries;
        if (fin.Contains( "_QCD_Pt_30to50"))scale=(xsec[12]*lumi)/nentries;
        if (fin.Contains( "_QCD_Pt_50to80"))scale=(xsec[13]*lumi)/nentries;
        if (fin.Contains( "_QCD_Pt_80to120"))scale=(xsec[14]*lumi)/nentries;
        if (fin.Contains( "_QCD_Pt_120to170"))scale=(xsec[15]*lumi)/nentries;
        if (fin.Contains( "_QCD_Pt_170to300"))scale=(xsec[16]*lumi)/nentries;
        if (fin.Contains( "_QCD_Pt_300to470"))scale=(xsec[17]*lumi)/nentries;
        if (fin.Contains( "_QCD_Pt_470to600"))scale=(xsec[18]*lumi)/nentries;
        if (fin.Contains( "_QCD_Pt_600to800"))scale=(xsec[19]*lumi)/nentries;
        if (fin.Contains( "_QCD_Pt_800to1000"))scale=(xsec[20]*lumi)/nentries;
        if (fin.Contains( "_QCD_Pt_1000toInf"))scale=(xsec[21]*lumi)/nentries;
        
        if (fin.Contains( "_QCDEM_Pt_20to30"))scale=(emqcd_xsec[0]*lumi)/nentries;
        if (fin.Contains( "_QCDEM_Pt_30to50"))scale=(emqcd_xsec[1]*lumi)/nentries;
        if (fin.Contains( "_QCDEM_Pt_50to80"))scale=(emqcd_xsec[2]*lumi)/nentries;
        if (fin.Contains( "_QCDEM_Pt_80to120"))scale=(emqcd_xsec[3]*lumi)/nentries;
        if (fin.Contains( "_QCDEM_Pt_120to170"))scale=(emqcd_xsec[4]*lumi)/nentries;
        if (fin.Contains( "_QCDEM_Pt_170to300"))scale=(emqcd_xsec[5]*lumi)/nentries;
        if (fin.Contains( "_QCDEM_Pt_300toInf"))scale=(emqcd_xsec[6]*lumi)/nentries;

        cout<<"this is scale:   "<<scale<<",  nentries:  "<<nentries<<endl;
//  TH1D* h_jets_pt = new TH1D("MZ","M(Z)#rightarrow #mu#mu",40, 71,111.);
     TH1D* h_events_counter = new TH1D("events_counter","events_counter",8, 0.0,8.0);
     TH1D* h_events_eachCut = new TH1D("events_eachCut","events_eachCut",6, 0.0,6.0);
     TH1D* h_events_eachCutp15_40 = new TH1D("events_eachCutp15_40","events_eachCutp15_40",6, 0.0,6.0);
     TH1D* h_events_eachCutp40_70 = new TH1D("events_eachCutp40_70","events_eachCutp40_70",6, 0.0,6.0);
     TH1D* h_events_eachCutp70_10 = new TH1D("events_eachCutp70_10","events_eachCutp70_10",6, 0.0,6.0);
     TH1D* h_events_eachCutp10_50 = new TH1D("events_eachCutp10_50","events_eachCutp10_50",6, 0.0,6.0);
     TH1D* h_events_eachCutp50 = new TH1D("events_eachCutp50","events_eachCutp50",6, 0.0,6.0);
     TH1D* h_events_eachCut_el = new TH1D("events_eachCut_el","events_eachCut_el",3, 0.0,3.0);
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

     TH1D* h_M_t = new TH1D("M_T","M_T",50, 0.0,500.0);
TH1D** h_M_tQCD=new TH1D*[6]; char nam_mw[100];
for (int i=1; i<6; i++){
  sprintf(nam_mw,"nam_mW%i",i);
     h_M_tQCD[i] = new TH1D(nam_mw,nam_mw,50, 0.0,500.0);
}
     TH1D* h_M_t_aft = new TH1D("M_T_aft","M_T_aft",50, 0.0,500.0);
//------------------------------ Jets ------------------------------------------------//
     TH1D* h_no_Jets = new TH1D("n_Jets","n_Jets",11, -0.5,10.5);
     TH1D* h_no_Jets20 = new TH1D("n_Jets20","n_Jets20",11, -0.5,10.5);
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
  TH1D* h_n_bjets = new TH1D("n_bjets","n_bjets",6, -0.5,5.5);
  TH1D* h_n_bjets_corr = new TH1D("n_bjets_corr","n_bjets",6, -0.5,5.5);
  TH1D* h_n_bjets_aft = new TH1D("n_bjets_aft","n_bjets",6, -0.5,5.5);
  TH1D* h_n_cjets = new TH1D("n_cjets","n cjets",6, -0.5,5.5);
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
  Float_t bins[] = { 250.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 610.0, 640.0, 680.0, 730.0, 800.0, 920.0, 1200.0};
  Int_t  binnum = sizeof(bins)/sizeof(Float_t) - 1;
  TH1D* h_M_tt_binned   = new TH1D("Mass_H_binned","t#bar{t} mass [GeV]",binnum,bins);
  TH1D* h_M_tt   = new TH1D("Mass_H","t#bar{t} mass [GeV]",100, 250.0,1200.0);
    TH1D* h_Pt_tt   = new TH1D("Pt_H","t#bar{t} [GeV]",50, 0.0,500.0);
    TH1D* leptop_cos_ttrest   = new TH1D("cos_theta","cos theta in t#bar{t} RestF ",25,-1., 1.0);

char nam_massb[7]; char nam_massq[7]; char theta[7];
TH1D** h_M_tt_binnedQCD=new TH1D*[5];TH1D** h_M_ttQCD=new TH1D*[5];TH1D** leptop_cos_ttrestQCD=new TH1D*[5];
for (int i=1; i<4; i++){
  sprintf(nam_massb,"Mass_H_binned%i",i);sprintf(nam_massq,"Mass_H%i",i);sprintf(theta,"cos_theta%i",i);
  h_M_tt_binnedQCD[i]   = new TH1D(nam_massb,nam_massb,binnum,bins);
  h_M_ttQCD[i]   = new TH1D(nam_massq,nam_massq,100, 250.0,1200.0);
  leptop_cos_ttrestQCD[i]   = new TH1D(theta,theta,15,-1., 1.0);
}
//------------------------------------------------------------------------------------//

        TH1D* h_EvtInfo_NumVtx  = new TH1D("EvtInfo_NumVtx","EvtInfo_NumVtx",40, 0.0,40.0);
        TH1D* h_EvtInfo_NumVtx_w  = new TH1D("EvtInfo_NumVtx_w","EvtInfo_NumVtx_w",40, 0.0,40.0);
        TH1D* h_PU_npT  = new TH1D("PU_npT","PU_npT",50, 0.0,40.0);
        TH1D* h_bDisct_cMVAv2 = new TH1D("bDiscMVAv2","bDiscMVAv2",50, 0.0,1.0);
        TH1D* h_bDisct_cMVAv2aft = new TH1D("cMVA_v2","cMVA_v2",50, 0.0,1.0);
        TH1D* h_pfJet_cmult = new TH1D("Jet_cMult","Jet_cMult",10, 0.0,40.0);
        TH1D* h_NuDiscr = new TH1D("NuDiscr","NuDiscr()",25, 0.0,10.0);
        TH1D* h_MassDiscr = new TH1D("MassDiscr","MassDiscr()",25, 0.0,20.0);
        TH1D* h_lheSigEvn = new TH1D("lheSigEvn","lheSigEvn",10, 0.0,10.0);

        TLorentzVector v_elec;
        TLorentzVector v_muon;
        TLorentzVector v_muon9;
        TLorentzVector v_met;
        TLorentzVector v_jets;
        TLorentzVector v_jetAll[1000];
        TLorentzVector v_jet[100];
        TLorentzVector v_bjets;
        TLorentzVector v_bjets9;
        TLorentzVector v_ljets;
//        TTree *weight_tree;
        nentries=5000;
        for (int jentry=0; jentry < nentries; jentry++)
        {
        t->GetEntry(jentry);
        if(jentry%1==0)cout<<" << "<<jentry<<"/"<<nentries<<endl;
/*        float iSecret;
        srand (time(NULL));
        iSecret = rand() % 100 + 1;
        bool run_273158_274093 (iSecret <= 5.);
        bool run_274094_276097 (5. < iSecret);// && iSecret < 63.);//Fix me when available the rest
        bool elrun_below_273726 (iSecret <= 4.);
        bool elrun_above_273726 (4. < iSecret );//&& iSecret < 42.);//Fix me when available the rest
*/
          if(!realdata){
          if(puWgtGr.size())
          {
          puWeight[0]=puWgtGr[0]->Eval(int(PU_npT));
          puWeight[1]=puWgtGr[1]->Eval(int(PU_npT));
          puWeight[2]=puWgtGr[2]->Eval(int(PU_npT));
          }
          }
        double w=1;
        w=puWeight[0];
        if (!realdata)w*=mcWeight_;
        h_events_counter->Fill(0.,w);
        bool no_trig_samp(false);
        if (fin.Contains( "qcd_EMEnriched") || fin.Contains( "WWToLNuQQ" ) || fin.Contains( "tW_nohad" ) || fin.Contains( "TeV_WZ" )|| fin.Contains( "TeV_ZZ" ) || fin.Contains( "V_W1Jets") || fin.Contains( "V_W2Jets") || fin.Contains( "V_W3Jets") || fin.Contains( "V_W4Jets") )no_trig_samp=true;
//        if (lheSigEvn!=1 )continue;
        if (is_mu) if (HLT_IsoMu24!=1 && HLT_IsoTkMu24!=1)continue;
        if (!is_mu)if (HLT_Ele32_eta2p1_WPTight_Gsf!=1)continue;
        h_events_eachCut->Fill(0.,w);
        h_events_eachCutp15_40->Fill(0.,w);
        h_events_eachCutp40_70->Fill(0.,w);
        h_events_eachCutp70_10->Fill(0.,w);
        h_events_eachCutp10_50->Fill(0.,w);
        h_events_eachCutp50->Fill(0.,w);
//        if (fin.Contains( "HToTT-semilep"))if (lheSigEvn!=3)continue;//Resonance: 5%=1, 10%=3, 25%=5, 50%=7;;; Interference: 5%=2, 10%=4, 25%=6, 50%=8;
        //noise
	      if (first_PV!=1)continue;
          if (Flag_HBHENoiseFilter!=1 || Flag_HBHENoiseIsoFilter!=1 || Flag_globalTightHalo2016Filter!=1)continue;
          if (Flag_EcalDeadCellTriggerPrimitiveFilter!=1 || Flag_goodVertices!=1 || Flag_eeBadScFilter!=1)continue;

      int n_pat_elec = 0,n_elec20=0;
      std::vector<TLorentzVector>lep_vector, loose_lep_vector;
      TLorentzVector v_lep;
      vector<unsigned int> n_lep_v;
      vector<int> lep_charge;
      vector<unsigned int> n_lep10_v;
      std::vector<TLorentzVector>el_vector, loose_el_vector;
      vector<int>el_charge;
      vector<unsigned int> n_elec_v;
      vector<unsigned int> n_elec20_v;
      for (unsigned int elec =0; elec < patElecPt_->size(); ++elec){
        if(patElecPt_->at(elec)>20 && fabs(patElecScEta_->at(elec))<2.5 && patElecIdveto_->at(elec) !=0){
            TLorentzVector loose_el_vector_;
            loose_el_vector_.SetPtEtaPhiE(patElecPt_->at(elec),patElecEta_->at(elec),patElecPhi_->at(elec),patElecEnergy_->at(elec));
            loose_el_vector.push_back(loose_el_vector_);
            n_elec20_v.push_back(n_elec20);
            n_elec20++;
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
      bool loose_mu_cut(false);
      std::vector<TLorentzVector>mu_vector,loose_mu_vector;
      vector<int>mu_charge;
      vector<unsigned int> n_muon_v;
      vector<unsigned int> n_muon10_v;
      int n_muon=0,n_muon10=0;
      float muon_TIso=0;

      float mu_pt=0,mu_eta=0.,mu_Tid=0., mu_iso=0., mu_phi=0.,mu_E=0., mu_Lid=0.;
      for(unsigned int mu=0; mu<patMuonPt_->size();mu++){
           mu_pt=patMuonPt_->at(mu); mu_eta=patMuonEta_->at(mu); mu_Tid=patMuonTightId_->at(mu); mu_iso=patMuonPfIsoDbeta_->at(mu);
           mu_phi=patMuonPhi_->at(mu); mu_E=patMuonEn_->at(mu); mu_Lid=patMuonLooseId_->at(mu);
          h_patMuonPfIsoDbeta->Fill(patMuonPfIsoDbeta_->at(mu),w);
          if (mu_pt>10 && fabs(mu_eta)<2.4 && mu_Lid>0 && fabs(mu_iso) <0.25)loose_mu_cut=true;
          if(mu_pt > 26. && fabs(mu_eta)<2.4 && mu_Tid>0 && fabs(mu_iso) < 0.15)
            {n_muon++;
            n_muon_v.push_back(n_muon);
            v_muon.SetPtEtaPhiE(mu_pt,mu_eta,mu_phi,mu_E);
            mu_vector.push_back(v_muon);
            mu_charge.push_back(patMuonCharge_->at(mu));
            muon_TIso=patMuonPfIsoDbeta_->at(mu);
            }
        if(loose_mu_cut){
          n_muon10_v.push_back(n_muon10);
          TLorentzVector v_muon_;
          v_muon_.SetPtEtaPhiE(mu_pt,mu_eta,mu_phi,mu_E);
          loose_mu_vector.push_back(v_muon_);
          n_muon10++;
          }
            }

        if (!is_mu){n_lep_v = n_elec_v; n_lep10_v = n_elec20_v; lep_vector = el_vector;loose_lep_vector = loose_el_vector; lep_charge = el_charge; v_lep=v_elec;}
        if (is_mu){n_lep_v = n_muon_v; n_lep10_v = n_muon10_v; lep_vector = mu_vector;loose_lep_vector = loose_mu_vector; lep_charge = mu_charge; v_lep=v_muon;}
            //-------------MET--------------------------------------------------//
            for (unsigned int nu =0; nu < METPt->size(); nu++)
                {
                v_met.SetPxPyPzE(METPx->at(0),METPy->at(0),METPz->at(0),METE->at(0));
                }
            //----------------------------------PF Jets------------------------------------------------//
int test_njets20=0;
            int n_pat_bjets=0,ncjets=0;
            int n_ljets=0;
            float jet_cMult=0;
            int n_pat_jets=0;
            double DR_mu_j=9999,Up=0., Cn=0.,Dn=0.,sig_Up=0.,sig_Dn=0.,jesUnc=0.;
            vector<unsigned int> no_jets;
            vector<unsigned int> n_bjets_v;
            vector<TLorentzVector> bjets_v;
            vector<TLorentzVector> ljets_v;
            vector<TLorentzVector> vec_bjets;
            vector<unsigned int> n_ljets_v;
            vector<unsigned int> n_bjets_vector;
            vector<float> cMVA_v;
            double cMVA=-99,cMVAv2=-999, nom=1., denom=1.,ind_jet_const=0.;
            vector <TLorentzVector> temp1;
            vector <TLorentzVector> v_ljetsAll;
            TLorentzVector jp4_Up;
            TLorentzVector jp4_Dn; 
            vector<TLorentzVector>jet_vector;
            vector<double> JER_Uncer;
            TLorentzVector metp4;
            int n_bjets_t=0; 
            for (unsigned int pf=0;pf < patJetPfAk04PtJERSmear->size();++pf){
              if(patJetPfAk04PtJERSmear->at(pf) > 20.&& fabs(patJetPfAk04Eta_->at(pf)) < 2.4 && patJetPfAk04LooseId_->at(pf)>0){
                TLorentzVector v_jetsTemp20;
       v_jetsTemp20.SetPtEtaPhiE(patJetPfAk04PtJERSmear->at(pf),patJetPfAk04Eta_->at(pf),patJetPfAk04Phi_->at(pf),patJetPfAk04En_->at(pf));
       if (n_lep_v.size()>0)DR_mu_j= DeltaR(loose_lep_vector[0].Eta(), v_jetsTemp20.Eta(), loose_lep_vector[0].Phi(), v_jetsTemp20.Phi());
                if(DR_mu_j>=0.4){
                  test_njets20++;
              }}

              if(patJetPfAk04PtJERSmear->at(pf) > 30.&& fabs(patJetPfAk04Eta_->at(pf)) < 2.4 && patJetPfAk04LooseId_->at(pf)>0){
                TLorentzVector v_jetsTemp;
                v_jetsTemp.SetPtEtaPhiE(patJetPfAk04PtJERSmear->at(pf),patJetPfAk04Eta_->at(pf),patJetPfAk04Phi_->at(pf),patJetPfAk04En_->at(pf));
                if (lep_vector.size()>0)DR_mu_j= DeltaR(lep_vector[0].Eta(), v_jetsTemp.Eta(), lep_vector[0].Phi(), v_jetsTemp.Phi());
                if(DR_mu_j<0.4)continue;
              n_pat_jets++;
              Up = patJetPfAk04PtJERSmearUp->at(pf);//Fix Me I will be after or before the cut
              Dn = patJetPfAk04PtJERSmearDn->at(pf);
              Cn = patJetPfAk04PtJERSmear->at(pf);
              sig_Up=abs(Up-Cn)/Cn;
              sig_Dn=abs(Dn-Cn)/Cn;
              JER_Uncer.push_back(max(sig_Up,sig_Dn));
              jesUnc= unc_->at(pf);
              v_jets.SetPtEtaPhiE(patJetPfAk04PtJERSmear->at(pf),patJetPfAk04Eta_->at(pf),patJetPfAk04Phi_->at(pf),patJetPfAk04En_->at(pf));
              jp4_Up.SetPtEtaPhiE(patJetPfAk04PtJERSmearUp->at(pf),patJetPfAk04Eta_->at(pf),patJetPfAk04Phi_->at(pf),patJetPfAk04En_->at(pf));
              jp4_Dn.SetPtEtaPhiE(patJetPfAk04PtJERSmearDn->at(pf),patJetPfAk04Eta_->at(pf),patJetPfAk04Phi_->at(pf),patJetPfAk04En_->at(pf));
              jet_vector.push_back(v_jets);
              h_bDisct_cMVAv2->Fill(patJetPfAk04BDiscpfCMVA_->at(pf),w);
              if(patJetPfAk04BDiscpfCMVA_->at(pf) > cMVA) cMVA = patJetPfAk04BDiscpfCMVA_->at(pf);
	            cMVAv2=patJetPfAk04BDiscpfCMVA_->at(pf);
              cMVA_v.push_back(patJetPfAk04BDiscpfCMVA_->at(pf));
              jet_cMult=patJetPfAk04cmult_->at(pf);
              no_jets.push_back(n_pat_jets);
              int jflav( abs(patJetPfAk04PartonFlavour_->at(pf)) );
              if (cMVAv2>0.185){
                n_bjets_t++;
              n_bjets_vector.push_back(n_bjets_t);
              }
              bool isBTagged(cMVAv2>0.185);
              if(!realdata){
                float jptForBtag(v_jets.Pt()>1000. ? 999. : v_jets.Pt()), jetaForBtag(fabs(v_jets.Eta()));
                float expEff(1.0), jetBtagSF(1.0);
                if(abs(jflav)==4){
                  ncjets++;
                  expEff    = expBtagEff["c"]->Eval(jptForBtag);
                  jetBtagSF = sfcReaders[0]->eval( BTagEntry::FLAV_C, jetaForBtag, jptForBtag);
                //  cout<<"jetBtagSF for c:  "<<jetBtagSF<<endl;
                  jetBtagSF *= expEff>0 ? expBtagEff["c"]->Eval(jptForBtag)/expBtagEff["c"]->Eval(jptForBtag) : 0.;
                  }
                else if(abs(jflav)==5){
                  n_pat_bjets++;
                  expEff    = expBtagEff["b"]->Eval(jptForBtag);
                  jetBtagSF = sfbReaders[0]->eval( BTagEntry::FLAV_B, jetaForBtag, jptForBtag);
                  jetBtagSF *= expEff>0 ? expBtagEff["b"]->Eval(jptForBtag)/expBtagEff["b"]->Eval(jptForBtag) : 0.;
                  }
                else{
                  n_ljets++;
                  expEff    = expBtagEff["udsg"]->Eval(jptForBtag);
                  jetBtagSF = sflReaders[0]->eval( BTagEntry::FLAV_UDSG, jetaForBtag, jptForBtag);
                  jetBtagSF *= expEff> 0 ? expBtagEff["udsg"]->Eval(jptForBtag)/expBtagEff["udsg"]->Eval(jptForBtag) : 0.;
                  }
                myBTagSFUtil.modifyBTagsWithSF(isBTagged,    jetBtagSF,     expEff);//btag SF works, certified
              }

//              cout<<"this is isBTagged:  "<<isBTagged<<endl;
              if(isBTagged){
                n_bjets_v.push_back(n_pat_bjets); 
                bjets_v.push_back(v_jets);
                v_bjets=v_jets;
                TLorentzVector *lepT, *bjetsT;    
                lepT= &v_lep;
                bjetsT= &v_bjets;
                double test;
                NeutrinoSolver NS(lepT, bjetsT);
                metp4 = TLorentzVector(NS.GetBest(v_met.X(), v_met.Y(), 1., 1., 0.,test));
                if (test == -1 )continue;
                }

              if(!isBTagged){//Fix me, I have some some events with isBTagged==0 but cMVA>0.185
                temp1.push_back(v_jets);
                n_ljets_v.push_back(n_pat_jets);
                ljets_v.push_back(v_jets);
                v_ljets=v_jets;
                }
              }
            }

              h_n_bjets->Fill(n_bjets_t,w);
              h_n_bjets_corr->Fill(n_pat_bjets,w);

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


//-----------------------------------------------
      // Trig SF
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
     if (!realdata)w*=lepTriggerSF[0]*lepidSF[0]*lepisoSF[0];


     cout<<"this is only w:  "<<w<<endl;
     w*=lepTriggerSF[0];
     cout<<"this is w1:  "<<w<<endl;
     w*=lepTriggerSF[0]*lepidSF[0];
     cout<<"this is  w2:  "<<w<<endl;
     w*=lepTriggerSF[0]*lepidSF[0]*lepisoSF[0];
     cout<<"this is  w3:  "<<w<<endl;
     w*=lepTriggerSF[0]*lepidSF[0]*lepisoSF[0];
cout<<"this is final w: :::::::::"<<w<<endl;


   bool lep0_cut(false), lep1_cut(false), trans_mW_cut(false),jets_cut(false), bjets_cut(false), ljets_cut(false);
              if (is_mu)if (n_elec_v.size() ==0 && n_elec20_v.size() ==0 )lep0_cut=true;
              if (!is_mu)if (n_muon_v.size() ==0 && n_muon10_v.size() ==0)lep0_cut=true;
              if(n_lep_v.size() == 1 && n_lep10_v.size() ==1)lep1_cut=true;
double trans_m_w=0.;
if (lep_vector.size()>0)trans_m_w=sqrt(pow(lep_vector[0].Pt() + v_met.Pt(), 2) - pow(lep_vector[0].Px() + v_met.Px(), 2) - pow(lep_vector[0].Py() + v_met.Py(), 2));
              if (trans_m_w > 50.)trans_mW_cut=true;
              if(n_pat_jets > 3)jets_cut=true;
              if (bjets_v.size() > 1)bjets_cut=true;
              if(n_ljets_v.size() > 1)ljets_cut=true;
//----------------Rochester algorithm----------------
if (lep0_cut ){h_events_eachCut->Fill(1.,w);
    if (lep1_cut){h_events_eachCut->Fill(2.,w);
    h_M_t->Fill(trans_m_w,w);
       if (trans_mW_cut){h_events_eachCut->Fill(3,w);h_events_eachCut_el->Fill(0.,w);
          h_no_Jets->Fill(n_pat_jets,w);
          h_no_Jets20->Fill(test_njets20,w);
          if (jets_cut){ h_events_eachCut->Fill(4.,w);h_events_eachCut_el->Fill(1.,w);
              if (bjets_cut){h_events_eachCut->Fill(5.,w);h_events_eachCut_el->Fill(2.,w);
                  if (ljets_cut){h_PU_npT->Fill(PU_npT,w);
                    h_n_ljets_aft->Fill(ljets_v.size(),w);
                    h_n_bjets_aft->Fill(n_bjets_v.size(),w);
                    h_n_cjets->Fill(ncjets,w);
                    h_EvtInfo_NumVtx->Fill(EvtInfo_NumVtx);
                    h_EvtInfo_NumVtx_w->Fill(EvtInfo_NumVtx,w);
                    h_M_t_aft->Fill(trans_m_w,w);


    int lepChar=lep_charge[0];
    double nuDisc, massDiscr, nuchi2;
    TTBarSolver solver;
    TLorentzVector Bhad;
    TLorentzVector Blep;
    TLorentzVector Whad;
    TLorentzVector Thad;
    TLorentzVector Wlep;
    TLorentzVector Tlep;
    double costheta=999.;
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
      Tlep=best_permutation.TLep();
      Wlep=best_permutation.WLep();
      nuDisc=best_permutation.NuDiscr();
      nuchi2=best_permutation.NuChisq();
      massDiscr=best_permutation.MassDiscr();

      hyp::TTbar ttang(best_permutation);
      auto ttcm = ttang.to_CM();
      auto tlepcm = ttang.tlep().to_CM();
      auto thadcm = ttang.thad().to_CM();
      costheta = ttang.unit3D().Dot(ttcm.tlep().unit3D());
    }
    }}}}

//----------------------------Muon parameters--------------------------------------------//
        h_lheSigEvn->Fill(lheSigEvn,w);
        h_met_Pz->Fill(metp4.Pz(),w);
        float Pt_lep = v_lep.Pt();
        h_Pt_lep->Fill(Pt_lep,w);
        h_Eta_lep->Fill(v_lep.Eta(),w);
//        Phi_lep = v_lep.Phi();
        h_Phi_lep->Fill(v_lep.Phi(),w);
//        E_lep = v_lep.E();
        h_E_lep->Fill(v_lep.E(),w);
        h_patMuonPfIsoDbeta_after->Fill(muon_TIso,w);

//-----------------------------Jets param------------------------------------------------//
//--------------------All jets-------------------------
        h_bDisct_cMVAv2aft->Fill(cMVA,w);
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
    if (Thad.M()!=0){
    h_M_hadW->Fill(Whad.M(),w);
    h_Pt_hadW->Fill(Whad.Pt(),w);
    h_M_hadTop->Fill(Thad.M(),w);
    h_Pt_hadTop->Fill(Thad.Pt(),w);
    h_M_lepTop->Fill(Tlep.M(),w);
    h_Pt_lepTop->Fill(Tlep.Pt(),w);
    h_M_lepW->Fill(Wlep.M(),w);
    h_Pt_lepW->Fill(Wlep.Pt(),w);
    h_M_tt->Fill((Tlep+Thad).M(),w);
    h_M_tt_binned->Fill((Tlep+Thad).M(),w);
    h_Pt_tt->Fill((Tlep+Thad).Pt(),w);
    if (costheta!=999.)leptop_cos_ttrest->Fill(costheta,w);
    }
}}}}}}
 // TTree weight_tree;
//  weight_tree = new TTree("tree","tree");
//  weight_tree->Fill();
        }//entries loop
        theFile->Write();
        theFile->Close();
        }//function loop
