    #include "simpleReader.h"
    void simpleReader(TString fin,TString fout){
      TFile *f = TFile::Open(fin);
      if (f->IsZombie()) {
      printf("Input root files doesn't open, please have a look:\n");
      return;}
    cout<<"this is fileIn:  "<<fin<<"   ;  "<<fout<<endl;
    TTree *t = (TTree*)f->Get("tupel/MuonTree");
    TFile *theFile = new TFile (fout+".root","RECREATE");
    theFile->cd();  
    bool is_mu (true);
//  bool is_mu (false);
    cout<<"Running Muon Channel: = "<<is_mu<<endl;
//-----------------------------------------------------------------------//


//-----------------------------------------------------------------------//
// xsec order:    ttbar, st_tch, st_tw, wjets,    DY,    WW,    WZ,    ZZ,   ,ttW,    ttZ,    qcd pt 15-20,       qcd pt 20-30,      qcd pt 30-50,       qcd pt 50-80,        qcd pt 80-120, qcd pt 120-170, qcd pt 170-300, qcd pt 300-470, qcd pt 470-600,     qcd pt 600-800, qcd pt 800-1000, qcd pt 1000-Inf
    float xsec[] ={831.76, 70.70, 19.55, 61526.7, 5765.4, 118.7, 47.13 ,16.523,0.2043, 0.2529, 1273190000.*0.003, 558528000.*0.0053, 139803000.*0.01182, 19222500.*0.02276, 2758420.*0.03844, 469797.*0.05362, 117989.*0.07335,7820.25*0.10196, 645.528* 0.12242, 187.109*0.13412, 32.3486*0.14552, 10.4305*0.15544}; 
//                     qcd pt20-30,        qcd pt 30-50,     , qcd pt 50-80,       qcd pt 80-120,  qcd pt 120-170,   qcd pt 170-300, qcd pt 300-Inf
    float emqcd_xsec[]={558528000.*0.0053, 139803000.*0.01182, 19222500.*0.02276, 2758420.*0.03844, 469797.*0.05362, 117989.*0.07335, 7820.25*0.10196};
    TString era("data/moriond17/");
//https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco
          TString btagUncUrl(era+"cMVAv2_Moriond17_B_H.csv");
          cout<<"string is : "<<btagUncUrl<<endl;
          gSystem->ExpandPathName(btagUncUrl);
          std::vector<BTagCalibrationReader *> sfbReaders,sfcReaders, sflReaders;
          TString btagEffExpUrl(era+"expTageff.root");
          gSystem->ExpandPathName(btagEffExpUrl);
          std::map<TString, TGraphAsymmErrors *> expBtagEff, expBtagEffPy8;
          BTagSFUtil myBTagSFUtil;
cout<<"this is realdata: "<<realdata<<endl;
BTagCalibration btvcalib("cMVAv2", btagUncUrl.Data());
std::map<BTagEntry::JetFlavor,BTagCalibrationReader *> btvCalibReaders;
BTagCalibrationReader reader(BTagEntry::OP_MEDIUM,"central", {"up", "down"});
reader.load(btvcalib,BTagEntry::FLAV_B,"ttbar"); 
reader.load(btvcalib,BTagEntry::FLAV_C,"ttbar"); 
reader.load(btvcalib,BTagEntry::FLAV_UDSG,"incl"); 


          if (!realdata){
            TFile *beffIn=TFile::Open(btagEffExpUrl);
            expBtagEff["b"]=(TGraphAsymmErrors *)beffIn->Get("b");
            expBtagEff["c"]=(TGraphAsymmErrors *)beffIn->Get("c");
            expBtagEff["udsg"]=(TGraphAsymmErrors *)beffIn->Get("udsg");
            beffIn->Close();
            }

      // Trig SF
      int mu_id=-999;
      int e_id=-999;
      TString eltrigSFurl(era+"eleTriggerSFGsf25_27_Run2016All_v1.root");
// new files for egamma id and efficiencies
//https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2
      TString elTightIDSFurl(era+"egammaEffi_TightCutID_txt_EGM2D.root");
      TString elEffSFurl(era+"egammaEffi_Recotxt_EGM2D.root");
      TString mutrigSFurl_BCDEF(era+"MuTrigg_EfficienciesAndSF_RunBtoF.root");
      TString mutrigSFurl_GH(era+"MuTrigg_EfficienciesAndSF_GH.root");
// for 12.9/pb https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaIDRecipesRun2#Cut_based_electron_identificatio
//for 7.6/fb https://twiki.cern.ch/twiki/bin/view/CMS/MuonWorkInProgressAndPagResults#Results_on_7_6_fb
      TString lepIDurl_BCDEF(era+"MuonID_SF_BCDEF.root");
      TString lepIDurl_GH(era+"MuonID_SF_GH.root");
      
      TString muIsourl_BCDEF(era+"MuonIso_SF_BCDEF.root");
      TString muIsourl_GH(era+"MuonIso_SF_GH.root");
      
      std::map<TString,TH2 *> lepEffH;
      TString Trig_hist("muon_trigg");
      TString Trig_histB_F("muon_trigg_bcdef");
      TString Trig_histG_H("muon_trigg_gh");

      if (!is_mu)Trig_hist="Ele32_eta2p1_WPTight_Gsf";
      TString eid_hist("elec_id"),esel_hist("elec_sel");
      TString lepidBCDEF_hist("muon_idbcdef"), lepisoBCDEF_hist("muon_isobcdef");
      TString lepidGH_hist("muon_idgh"), lepisoGH_hist("muon_isogh");
      TString lepid_hist("muon_id"), lepiso_hist("muon_iso");
      TString id_hist("muon_id");
      if(!is_mu)id_hist="elec_id";
      if (!realdata){
      if (is_mu ){
      TFile *mufInB_F=TFile::Open(mutrigSFurl_BCDEF);
         lepEffH[Trig_histB_F]=(TH2 *)mufInB_F->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/pt_abseta_ratio")->Clone();
         lepEffH[Trig_histB_F]->SetDirectory(0);
         mufInB_F->Close();
      TFile *mufInG_H=TFile::Open(mutrigSFurl_GH);
         lepEffH[Trig_histG_H]=(TH2 *)mufInG_H->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/pt_abseta_ratio")->Clone();
         lepEffH[Trig_histG_H]->SetDirectory(0);
         mufInG_H->Close();

      TFile *muIDfInBDCEF=TFile::Open(lepIDurl_BCDEF);
        lepEffH[lepidBCDEF_hist]=(TH2 *)muIDfInBDCEF->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio")->Clone();
        lepEffH[lepidBCDEF_hist]->SetDirectory(0);
        muIDfInBDCEF->Close();
      TFile *muIDfInGH=TFile::Open(lepIDurl_GH);
        lepEffH[lepidGH_hist]=(TH2 *)muIDfInGH->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio")->Clone();
        lepEffH[lepidGH_hist]->SetDirectory(0);
        muIDfInGH->Close();

      TFile *muIsofInBCDEF=TFile::Open(muIsourl_BCDEF);
         lepEffH[lepisoBCDEF_hist]=(TH2 *)muIsofInBCDEF->Get("TightISO_TightID_pt_eta/pt_abseta_ratio")->Clone();
         lepEffH[lepisoBCDEF_hist]->SetDirectory(0);
         muIsofInBCDEF->Close();
      TFile *muIsofInGH=TFile::Open(muIsourl_GH);
         lepEffH[lepisoGH_hist]=(TH2 *)muIsofInGH->Get("TightISO_TightID_pt_eta/pt_abseta_ratio")->Clone();
         lepEffH[lepisoGH_hist]->SetDirectory(0);
         muIsofInGH->Close();
        }
   if (!is_mu){
      TFile *eIDfIn=TFile::Open(elTightIDSFurl);
      lepEffH[eid_hist]=(TH2 *)eIDfIn->Get("EGamma_SF2D");
      lepEffH[eid_hist]->SetDirectory(0);
        eIDfIn->Close();
      TFile *eSelfIn=TFile::Open(elEffSFurl);
      lepEffH[esel_hist]=(TH2 *)eSelfIn->Get("EGamma_SF2D");
      lepEffH[esel_hist]->SetDirectory(0);
        eSelfIn->Close();
      TFile *fIn=TFile::Open(eltrigSFurl);
      lepEffH[Trig_hist]=(TH2 *)fIn->Get("Ele25_eta2p1_WPTight_Gsf");
      for(auto& it : lepEffH) it.second->SetDirectory(0);
              fIn->Close();

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
     float luminosity = 12877.401701;
     if (!is_mu)luminosity = 12883.8601;
        if (fin.Contains( "TTJets"))scale=(xsec[0]*luminosity)/nentries;
        if (fin.Contains( "SingleT_t.root"))scale=(xsec[1]*luminosity)/nentries;
        if (fin.Contains( "_tW_nohad"))scale=(xsec[2]*luminosity)/nentries;
        if (fin.Contains( "WJets"))scale=(xsec[3]*luminosity)/nentries;
        if (fin.Contains( "DYJets"))scale=(xsec[4]*luminosity)/nentries;
        if (fin.Contains( "_WWToLNuQQ"))scale=(xsec[5]*luminosity)/nentries;
        if (fin.Contains( "_WZ"))scale=(xsec[6]*luminosity)/nentries;
        if (fin.Contains( "_ZZ"))scale=(xsec[7]*luminosity)/nentries;
        if (fin.Contains( "TTWJetsToLNu"))scale=(xsec[8]*luminosity)/nentries;
        if (fin.Contains( "TTZToLLNuNu_M-10"))scale=(xsec[9]*luminosity)/nentries;
        if (fin.Contains( "_QCD_Pt_15to20"))scale=(xsec[10]*luminosity)/nentries;
        if (fin.Contains( "_QCD_Pt_20to30"))scale=(xsec[11]*luminosity)/nentries;
        if (fin.Contains( "_QCD_Pt_30to50"))scale=(xsec[12]*luminosity)/nentries;
        if (fin.Contains( "_QCD_Pt_50to80"))scale=(xsec[13]*luminosity)/nentries;
        if (fin.Contains( "_QCD_Pt_80to120"))scale=(xsec[14]*luminosity)/nentries;
        if (fin.Contains( "_QCD_Pt_120to170"))scale=(xsec[15]*luminosity)/nentries;
        if (fin.Contains( "_QCD_Pt_170to300"))scale=(xsec[16]*luminosity)/nentries;
        if (fin.Contains( "_QCD_Pt_300to470"))scale=(xsec[17]*luminosity)/nentries;
        if (fin.Contains( "_QCD_Pt_470to600"))scale=(xsec[18]*luminosity)/nentries;
        if (fin.Contains( "_QCD_Pt_600to800"))scale=(xsec[19]*luminosity)/nentries;
        if (fin.Contains( "_QCD_Pt_800to1000"))scale=(xsec[20]*luminosity)/nentries;
        if (fin.Contains( "_QCD_Pt_1000toInf"))scale=(xsec[21]*luminosity)/nentries;
        
        if (fin.Contains( "_QCDEM_Pt_20to30"))scale=(emqcd_xsec[0]*luminosity)/nentries;
        if (fin.Contains( "_QCDEM_Pt_30to50"))scale=(emqcd_xsec[1]*luminosity)/nentries;
        if (fin.Contains( "_QCDEM_Pt_50to80"))scale=(emqcd_xsec[2]*luminosity)/nentries;
        if (fin.Contains( "_QCDEM_Pt_80to120"))scale=(emqcd_xsec[3]*luminosity)/nentries;
        if (fin.Contains( "_QCDEM_Pt_120to170"))scale=(emqcd_xsec[4]*luminosity)/nentries;
        if (fin.Contains( "_QCDEM_Pt_170to300"))scale=(emqcd_xsec[5]*luminosity)/nentries;
        if (fin.Contains( "_QCDEM_Pt_300toInf"))scale=(emqcd_xsec[6]*luminosity)/nentries;
        
//        nloSF=((double)(SumofPosWeights+SumofNegWeights)/(SumofPosWeights-SumofNegWeights))*mcWeight_;
        double nloSF=1.0;
        if (fin.Contains( "ST_sch"))nloSF=(4390570.0 + 1019900.0)/(4390570.0 - 1019900.0);
        if (fin.Contains( "TTWToLNu"))nloSF=(1076350.0 + 344472)/(1076350.0 - 344472);     
        if (fin.Contains( "TTWToQQ"))nloSF=(836059.0 + 266635.0)/(836059.0 - 266635.0);   
        if (fin.Contains( "TTZToLLNuNu"))nloSF=(790230.0 + 288031.0)/(790230.0 - 288031.0);
        if (fin.Contains( "TTZToQQ"))nloSF=(621075.0 + 224734.0)/(621075.0 - 224734.0);
cout<<"this is nloSF:  "<<nloSF<<endl;
//  TH1D* h_jets_pt = new TH1D("MZ","M(Z)#rightarrow #mu#mu",40, 71,111.);
     TH1D* h_mcWeight = new TH1D("mcWeight","mcWeight",2, 0.0,2.0);
     TH1D* h_events_counter = new TH1D("events_counter","events_counter",8, 0.0,8.0);
     TH1D* h_events_eachCut = new TH1D("events_eachCut","events_eachCut",6, 0.0,6.0);
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
//------------------------------Lepton-------------------------------------------------//
       std::map<TString,TH1 *> histos;
       for(int ireg=0; ireg<2; ireg++)
             {
             TString pf(ireg==0 ? "_ee" : "_eb");
             histos["sigmaietaieta"+pf]=new TH1D("sigmaietaieta"+pf,";#sigma(i#eta,i#eta);Electrons",25,0,0.04);
             histos["detain"+pf]=new TH1D("detain"+pf,";#Delta#eta(in);Electrons",25,-0.015,0.15);
             histos["dphiin"+pf]=new TH1D("dephiin"+pf,";#Delta#phi(in);Electrons",25,-0.05,0.05);
             histos["hovere"+pf]=new TH1D("hovere"+pf,";H/E;Electrons",25,0,0.05);
             histos["eoverpinv"+pf]=new TH1D("eoverpinv"+pf,";1/E-1/p [GeV^{-1}];Electrons",25,0,0.1);
             histos["d0"+pf]=new TH1D("d0"+pf,";d_{0} [cm];Electrons",25,-0.05,0.05);
             histos["dz"+pf]=new TH1D("dz"+pf,";d_{z} [cm];Electrons",25,-0.5,0.5);
             histos["misshits"+pf]=new TH1D("misshits"+pf,";Missing hits;Electrons",4,0,4);
             histos["convveto"+pf]=new TH1D("convveto"+pf,";Conversion veto flag;Electrons",4,0,4);
             histos["reliso"+pf]=new TH1D("reliso"+pf,";Relative isolation;Electrons",25,0,0.1);
             }
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
  TH1D* h_n_bjets_aft = new TH1D("n_bjets_aft","n_bjets",6, -0.5,5.5);
  TH1D* h_n_cjets = new TH1D("n_cjets","n cjets",6, -0.5,5.5);
  histos["deepcsv_btag"]=new TH1D("deepcsv_btag","DeepCSV btag",25,0.2,1.2);
  histos["deepcsv_ctag"]=new TH1D("deepcsv_ctag","DeepCSV ctag",25,0.2,1.2);
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
  Float_t bins[] = { 
          250.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 
          560.0, 580.0, 610.0, 640.0, 680.0, 730.0, 800.0, 920.0, 1200.0};
  Float_t cos_bins[]={0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
  Int_t cos_binnum=sizeof(cos_bins)/sizeof(Float_t) - 1;
  Int_t  binnum = sizeof(bins)/sizeof(Float_t) - 1;
  TH1D* h_M_tt_binned   = new TH1D("Mass_H_binned","t#bar{t} mass [GeV]",binnum,bins);
  TH1D* h_M_tt   = new TH1D("Mass_H","t#bar{t} mass [GeV]",100, 250.0,1200.0);
    TH1D* h_Pt_tt   = new TH1D("Pt_H","t#bar{t} [GeV]",50, 0.0,500.0);
    TH1D* leptop_cos_ttrest   = new TH1D("cos_theta","cos#theta in t#bar{t} RestF ",25,-1., 1.0);
    TH2D* h_M_tt_cos   = new TH2D("Mass_H_costheta1","t#bar{t} mass vs cosine [GeV]",binnum,bins,cos_binnum, cos_bins);
//------------------------------------------------------------------------------------//

        TH1D* h_EvtInfo_NumVtx  = new TH1D("EvtInfo_NumVtx","EvtInfo_NumVtx",40, 0.0,40.0);
        TH1D* h_EvtInfo_NumVtx1  = new TH1D("EvtInfo_NumVtx1","EvtInfo_NumVtx1",40, 0.0,40.0);
        TH1D* h_EvtInfo_NumVtx2  = new TH1D("EvtInfo_NumVtx2","EvtInfo_NumVtx2",40, 0.0,40.0);
        TH1D* h_EvtInfo_NumVtx3  = new TH1D("EvtInfo_NumVtx3","EvtInfo_NumVtx3",40, 0.0,40.0);
        TH1D* h_EvtInfo_NumVtx4  = new TH1D("EvtInfo_NumVtx4","EvtInfo_NumVtx4",40, 0.0,40.0);
        TH1D* h_EvtInfo_NumVtx_w  = new TH1D("EvtInfo_NumVtx_w","EvtInfo_NumVtx_w",40, 0.0,40.0);
        TH1D* h_PU_npT  = new TH1D("PU_npT","PU_npT",50, 0.0,40.0);
        TH1D* h_bDisct_cMVAv2 = new TH1D("bDiscMVAv2","bDiscMVAv2",50, 0.0,1.0);
        TH1D* h_bDisct_cMVAv2aft = new TH1D("cMVA_v2","cMVA_v2",50, 0.0,1.0);
        TH1D* h_pfJet_cmult = new TH1D("Jet_cMult","Jet_cMult",10, 0.0,40.0);
        TH1D* h_NuDiscr = new TH1D("NuDiscr","NuDiscr()",25, 0.0,10.0);
        TH1D* h_MassDiscr = new TH1D("MassDiscr","MassDiscr()",25, 0.0,20.0);

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
    TFile *synchroFile = new TFile (fout+"synchro.root","RECREATE");
    synchroFile->cd();  
  TTree *tree = new TTree("sync","sync");
 //       nentries=4200;

  string oFile = "_mu";
  if (!is_mu)oFile = "_el";
    ofstream outputFile_Objects, outputFile_HLT,outputFile_0lep,outputFile_1lep,
    outputFile_W_T, outputFile_4jets, outputFile_2bjets,outputFile_ljets, outputFile_RecoSuccess,outputFile_RecoFail,outputFile_RecoSuccessAll, outputFile_TTObjects;
    outputFile_Objects.open("testing/Object.txt",std::ios::trunc);
    outputFile_HLT.open("testing/events_HLT.txt",std::ios::trunc);
    outputFile_0lep.open("testing/events_0lep.txt",std::ios::trunc);
    outputFile_1lep.open("testing/events_1lep.txt",std::ios::trunc);
    outputFile_W_T.open("testing/events_W_T.txt",std::ios::trunc);
    outputFile_4jets.open("testing/events_4jets.txt",std::ios::trunc);
    outputFile_2bjets.open("testing/events_2bjets.txt",std::ios::trunc);
    outputFile_ljets.open("testing/events_ljets.txt",std::ios::trunc);
    outputFile_RecoSuccess.open("testing/events_RecoSuccess.txt",std::ios::trunc);
    outputFile_RecoFail.open("testing/events_RecoFail.txt",std::ios::trunc);
    outputFile_RecoSuccessAll.open("testing/events_RecoSuccessAll.txt",std::ios::trunc);
    outputFile_TTObjects.open("testing/events_TTObjects.txt",std::ios::trunc);

    outputFile_Objects<< "*************   Event   "<<run<<":"<< lumi_ << ":"<< event <<"   **************\n"<<endl; 
    outputFile_HLT<< "run:lumi:event"<<endl;
    outputFile_0lep<< "run:lumi:event"<<endl;
    outputFile_1lep<< "run:lumi:event"<<endl;
    outputFile_W_T<< "run:lumi:event"<<endl;
    outputFile_4jets<< "run:lumi:event"<<endl;
    outputFile_2bjets<< "run:lumi:event"<<endl;
    outputFile_ljets<< "run:lumi:event"<<endl;
    outputFile_RecoSuccess<< "run:lumi:event"<<endl;
    outputFile_RecoFail<< "run:lumi:event"<<endl;
    outputFile_RecoSuccessAll<< "run:lumi:event"<<endl;

        for (int jentry=0; jentry < nentries; jentry++)
        {
        t->GetEntry(jentry);
        if(jentry%1000==0)cout<<" << "<<jentry<<"/"<<nentries<<endl;

unsigned int Run=run;
ULong64_t LumiSection = lumi_;
ULong64_t Event = event;
Bool_t         RecoSuccess=0;
float        MassTT;
          if(!realdata){
          if(puWgtGr.size())
          {
          puWeight[0]=puWgtGr[0]->Eval(int(PU_npT));
          puWeight[1]=puWgtGr[1]->Eval(int(PU_npT));
          puWeight[2]=puWgtGr[2]->Eval(int(PU_npT));
          }
          }
        if (is_mu){
        float iSecret;
        srand (time(NULL));
        iSecret = rand() % 100 + 1;
        if (iSecret <= 55.){
          lepid_hist=lepidBCDEF_hist;
          lepiso_hist=lepisoBCDEF_hist;
          Trig_hist = Trig_histB_F;
          }
        else {
          lepid_hist=lepidGH_hist;
          lepiso_hist=lepisoGH_hist;
          Trig_hist = Trig_histG_H;
          }
        }
        double w=1;
        if (!realdata)w*=puWeight[0]*nloSF*mcWeight_;
        if (mcWeight_>0)h_mcWeight->Fill(0.,mcWeight_);
        else h_mcWeight->Fill(1.0,mcWeight_);
        h_events_counter->Fill(0.,w);
        bool triggerName=false;
        triggerName = is_mu ? HLT_IsoMu24==1 || HLT_IsoTkMu24==1 : HLT_Ele32_eta2p1_WPTight_Gsf==1; 
        if (!triggerName)continue;
    outputFile_HLT<< run<<":"<< lumi_ << ":"<< event <<","<<endl;
        h_events_eachCut->Fill(0.,w);
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
          bool passTightId = false;
                              if( (fabs(patElecScEta_->at(elec))<=1.479
                                && fabs(patElecfull5x5_sigmaIetaIeta_->at(elec)) < 0.00998
                                && fabs(patElecdEtaInSeed_->at(elec)) < 0.00308
                                && fabs(patElecdPhiIn_->at(elec)) < 0.0816
                                && fabs(patElechOverE_->at(elec)) < 0.0414
                                && fabs(patElecooEmooP_->at(elec))< 0.0129   
                                && fabs(patElecd0_->at(elec))     < 0.05   
                                && fabs(patElecdz_->at(elec))     < 0.10
                                && fabs(patElecexpectedMissingInnerHits_->at(elec)) <=1
                                && fabs(patElecpassConversionVeto_->at(elec)) 
                                )
                              ||
                               (fabs(patElecScEta_->at(elec))>1.479
                                && fabs(patElecfull5x5_sigmaIetaIeta_->at(elec)) < 0.0292
                                && fabs(patElecdEtaInSeed_->at(elec))            < 0.00605
                                && fabs(patElecdPhiIn_->at(elec))                < 0.0394
                                && fabs(patElechOverE_->at(elec))                < 0.0641 
                                && fabs(patElecooEmooP_->at(elec))               < 0.0129   
                                && fabs(patElecd0_->at(elec))                    < 0.10 
                                && fabs(patElecdz_->at(elec))                    < 0.20
                                && fabs(patElecexpectedMissingInnerHits_->at(elec))<=1
                                && fabs(patElecpassConversionVeto_->at(elec))
                                ))passTightId=true;
            bool passVetoId = false; 
                              if((fabs(patElecScEta_->at(elec))<=1.479
                                && fabs(patElecfull5x5_sigmaIetaIeta_->at(elec))< 0.0115 
                                && fabs(patElecdEtaInSeed_->at(elec)) < 0.00749 
                                && fabs(patElecdPhiIn_->at(elec)) < 0.228
                                && fabs(patElechOverE_->at(elec)) < 0.356 
                                && fabs(patElecooEmooP_->at(elec)) < 0.299
                               // && fabs(patElecd0_->at(elec)) < 0.05
                               // && fabs(patElecdz_->at(elec)) < 0.10 
                                && fabs(patElecexpectedMissingInnerHits_->at(elec)) <= 2 
                                && fabs(patElecpassConversionVeto_->at(elec))
                                )
                              ||
                               (fabs(patElecScEta_->at(elec))>1.479
                                && fabs(patElecfull5x5_sigmaIetaIeta_->at(elec)) < 0.037
                                && fabs(patElecdEtaInSeed_->at(elec))            < 0.00895
                                && fabs(patElecdPhiIn_->at(elec))                < 0.213 
                                && fabs(patElechOverE_->at(elec))                < 0.211
                                && fabs(patElecooEmooP_->at(elec))               < 0.15
                               // && fabs(patElecd0_->at(elec))                    < 0.10 
                               // && fabs(patElecdz_->at(elec))                    < 0.20 
                                && fabs(patElecexpectedMissingInnerHits_->at(elec)) <= 3 
                                && fabs(patElecpassConversionVeto_->at(elec))
                                ))passVetoId = true;
                                

      bool ip_cut=false;
	if (fabs(patElecScEta_->at(elec))<=1.479){if (fabs(patElecd0_->at(elec)) < 0.05 && fabs(patElecdz_->at(elec)) < 0.10)ip_cut=true;}
	if (fabs(patElecScEta_->at(elec))>1.479 ){if (fabs(patElecd0_->at(elec)) < 0.10 && fabs(patElecdz_->at(elec)) < 0.20)ip_cut=true;}
        if(patElecPt_->at(elec)>20 && fabs(patElecScEta_->at(elec))<2.5 && patElecIdveto_->at(elec) !=0){
            TLorentzVector loose_el_vector_;
            loose_el_vector_.SetPtEtaPhiE(patElecPt_->at(elec),patElecEta_->at(elec),patElecPhi_->at(elec),patElecEnergy_->at(elec));
            loose_el_vector.push_back(loose_el_vector_);
            n_elec20_v.push_back(n_elec20);
            n_elec20++;
            }
      if (patElecPt_->at(elec) > 35. && fabs(patElecScEta_->at(elec)) <2.1 && patElecIdtight_->at(elec)>0 && ip_cut ){
     if (1.4442 < fabs( patElecScEta_->at(elec)) && fabs(patElecScEta_->at(elec))< 1.5660)continue;
      n_pat_elec++;
      n_elec_v.push_back(n_pat_elec);
      v_elec.SetPtEtaPhiE(patElecPt_->at(elec),patElecEta_->at(elec),patElecPhi_->at(elec),patElecEnergy_->at(elec));
      el_vector.push_back(v_elec);
      el_charge.push_back(patElecCharge_->at(elec));
  outputFile_Objects<<"Electrons: "<<endl;
  outputFile_Objects<<"#"<<elec+1<<endl;
  outputFile_Objects<<"Momentum (pt, eta, phi): "<< patElecPt_->at(elec)<<" : "<<patElecEta_->at(elec)<< " : "<<patElecPhi_->at(elec)<<endl;
  outputFile_Objects<<"Pass trigger preselection:   " <<patElec_mva_presel_->at(elec)<<endl;
  outputFile_Objects<<"Pass cut-based ID (veto to tight): "<< patElecIdveto_->at(elec)<<" , "<<patElecIdloose_->at(elec)
                    <<" , "<<patElecIdmedium_->at(elec)<<" , "<<patElecIdtight_->at(elec)<<endl;

      TString pf(fabs(patElecScEta_->at(elec)) > 1.479 ? "_ee" : "_eb");
      histos["sigmaietaieta"+pf]->Fill(patElecfull5x5_sigmaIetaIeta_->at(elec),w);
      histos["detain"+pf]->Fill(patElecdEtaInSeed_->at(elec),w);
      histos["dphiin"+pf]->Fill(patElecdPhiIn_->at(elec),w);
      histos["hovere"+pf]->Fill(patElechOverE_->at(elec),w);
      histos["eoverpinv"+pf]->Fill(patElecooEmooP_->at(elec),w);
      histos["d0"+pf]->Fill(patElecd0_->at(elec),w);
      histos["dz"+pf]->Fill(patElecdz_->at(elec),w);
      histos["misshits"+pf]->Fill(patElecexpectedMissingInnerHits_->at(elec),w);
      histos["convveto"+pf]->Fill(patElecpassConversionVeto_->at(elec),w);
      histos["reliso"+pf]->Fill(patElecPfIsoRho_->at(elec),w);
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

      for(unsigned int mu=0; mu<patMuonPt_->size();mu++){
          h_patMuonPfIsoDbeta->Fill(patMuonPfIsoDbeta_->at(mu),w);
          if(patMuonPt_->at(mu) > 26. && fabs(patMuonEta_->at(mu))<2.4 && patMuonTightId_->at(mu)>0 && fabs(patMuonPfIsoDbeta_->at(mu)) < 0.15)
            {n_muon++;
            n_muon_v.push_back(n_muon);
            v_muon.SetPtEtaPhiE(patMuonPt_->at(mu),patMuonEta_->at(mu),patMuonPhi_->at(mu),patMuonEn_->at(mu));
            mu_vector.push_back(v_muon);
            mu_charge.push_back(patMuonCharge_->at(mu));
            muon_TIso=patMuonPfIsoDbeta_->at(mu);

    outputFile_Objects<<"\n"<<"Muons:"<<endl;
    outputFile_Objects<<"#"<<mu+1<<endl;
    outputFile_Objects<<"Momentum (pt, eta, phi): "<<patMuonPt_->at(mu)<< " : "<<patMuonEta_->at(mu)<<" : "<<patMuonPhi_->at(mu)<<endl;
    outputFile_Objects<<"Isolation: "<<patMuonPfIsoDbeta_->at(mu)<<endl;
    outputFile_Objects<<"Pass ID (loose to tight): "<<patMuonLooseId_->at(mu)<<"   "<<patMuonMediumId_->at(mu)<<"   "<<patMuonTightId_->at(mu)<<endl;

            }
        if(patMuonPt_->at(mu)>10 && fabs(patMuonEta_->at(mu))<2.4 && patMuonLooseId_->at(mu)>0 && fabs(patMuonPfIsoDbeta_->at(mu)) <0.25){
          n_muon10++;
          n_muon10_v.push_back(n_muon10);
          TLorentzVector v_muon_;
          v_muon_.SetPtEtaPhiE(patMuonPt_->at(mu),patMuonEta_->at(mu),patMuonPhi_->at(mu),patMuonEn_->at(mu));
          loose_mu_vector.push_back(v_muon_);
          }
            }

        if (!is_mu){n_lep_v = n_elec_v; n_lep10_v = n_elec20_v; lep_vector = el_vector;loose_lep_vector = loose_el_vector; lep_charge = el_charge; v_lep=v_elec;}
        if (is_mu){n_lep_v = n_muon_v; n_lep10_v = n_muon10_v; lep_vector = mu_vector;loose_lep_vector = loose_mu_vector; lep_charge = mu_charge; v_lep=v_muon;}
            //-------------MET--------------------------------------------------//
        float uncorec_metpt=0.,uncorec_metphi=0.;
            for (unsigned int nu =0; nu < METPt->size(); nu++)
                {
                uncorec_metpt=Uncorec_METPt->at(nu);
                uncorec_metphi=Uncorec_METPhi->at(nu);
                v_met.SetPxPyPzE(METPx->at(0),METPy->at(0),METPz->at(0),METE->at(0));
                }
            //----------------------------------PF Jets------------------------------------------------//
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
            TLorentzVector metp4,lorenz_v_raw;
                TLorentzVector *lepT, *bjetsT;    
                TLorentzVector *lepT_, *bjetsT_;   
                double test;
                double nstest_   = numeric_limits<double>::max();
                double nsbest_   = numeric_limits<double>::max();
 
            int n_bjets_t=0; 
            for (unsigned int pf=0;pf < patJetPfAk04PtJERSmear->size();++pf){
        //    for (unsigned int pf=0;pf < patJetPfAk04Pt_->size();++pf)

              /*v_jets.SetPtEtaPhiE(patJetPfAk04Pt_->at(pf),patJetPfAk04Eta_->at(pf),patJetPfAk04Phi_->at(pf),patJetPfAk04En_->at(pf));
              float genJet_pt(0);
              if (MGjPt->size()>0)genJet_pt = MGjPt->at(0);
              if (!realdata && v_jets.Pt()>0 && genJet_pt>0)
              
                float jerSmear=getJetResolutionScales(v_jets.Pt(),v_jets.Eta(),genJet_pt)[0];
                v_jets *= jerSmear;
              */

              Up = patJetPfAk04PtJERSmearUp->at(pf);//Fix Me I will be after or before the cut
              Dn = patJetPfAk04PtJERSmearDn->at(pf);
              Cn = patJetPfAk04PtJERSmear->at(pf);
              sig_Up=abs(Up-Cn)/Cn;
              sig_Dn=abs(Dn-Cn)/Cn;
              JER_Uncer.push_back(max(sig_Up,sig_Dn));
              jesUnc= unc_->at(pf);
    lorenz_v_raw.SetPtEtaPhiE(patJetPfAk04RawPt_->at(pf),patJetPfAk04Eta_->at(pf),patJetPfAk04Phi_->at(pf),patJetPfAk04RawEn_->at(pf));
     
              if(patJetPfAk04PtJERSmear->at(pf) > 20.&& fabs(patJetPfAk04Eta_->at(pf)) < 2.4 && patJetPfAk04LooseId_->at(pf)>0){
    
    outputFile_Objects<<"\n"<<"Jets: "<<endl;
    outputFile_Objects<<"#"<<pf+1<<endl;
    outputFile_Objects<<"Raw Momentum (pt, eta, phi, m): "<<patJetPfAk04RawPt_->at(pf)<<" : "<<patJetPfAk04Eta_->at(pf)
                      <<" : "<<patJetPfAk04Phi_->at(pf)<<" : "<<lorenz_v_raw.M()<<endl;
    outputFile_Objects<<"Fully Corrected pt:  "<<patJetPfAk04PtJERSmear->at(pf)<<endl;
    outputFile_Objects<<"CSV: "<<patJetPfAk04BDiscCSVv2_->at(pf)<<" ,  cMVA:  "<<patJetPfAk04BDiscpfCMVA_->at(pf)<< endl;
    outputFile_Objects<<"deepCSVbtag:  "<<patJetPfAk04BDiscdeepFb_->at(pf)+patJetPfAk04BDiscdeepFbb_->at(pf)<<", deepCSVctag: "
    <<patJetPfAk04BDiscdeepFc_->at(pf)+patJetPfAk04BDiscdeepFcc_->at(pf)/(1-(patJetPfAk04BDiscdeepFb_->at(pf)+patJetPfAk04BDiscdeepFbb_->at(pf)))<<endl;

          //    if(v_jets.Pt() > 20.&& fabs(v_jets.Eta()) < 2.4 && patJetPfAk04LooseId_->at(pf)>0){}
                TLorentzVector v_jetsTemp;
                v_jetsTemp.SetPtEtaPhiE(patJetPfAk04PtJERSmear->at(pf),patJetPfAk04Eta_->at(pf),patJetPfAk04Phi_->at(pf),patJetPfAk04EnJERSmear->at(pf));
          //      v_jetsTemp=v_jets;
                if (loose_lep_vector.size()>0)DR_mu_j= DeltaR(loose_lep_vector[0].Eta(), v_jetsTemp.Eta(), loose_lep_vector[0].Phi(), v_jetsTemp.Phi());
                if(DR_mu_j<0.4)continue;
    outputFile_Objects<<" DR(mu, j):  "<<DR_mu_j<<endl;
              n_pat_jets++;
              v_jets.SetPtEtaPhiE(patJetPfAk04PtJERSmear->at(pf),patJetPfAk04Eta_->at(pf),patJetPfAk04Phi_->at(pf),patJetPfAk04EnJERSmear->at(pf));
/*              v_jets.SetPtEtaPhiE(patJetPfAk04Pt_->at(pf),patJetPfAk04Eta_->at(pf),patJetPfAk04Phi_->at(pf),patJetPfAk04En_->at(pf));
          float genJet_pt(0);
          if (MGjPt->size()>0)genJet_pt = MGjPt->at(0);
              if (!realdata && v_jets.Pt()>0 && genJet_pt>0)
              {
                float jerSmear=getJetResolutionScales(v_jets.Pt(),v_jets.Eta(),genJet_pt)[0];
                v_jets *= jerSmear;
              }
 */             jp4_Up.SetPtEtaPhiE(patJetPfAk04PtJERSmearUp->at(pf),patJetPfAk04Eta_->at(pf),patJetPfAk04Phi_->at(pf),patJetPfAk04EnJERSmearUp->at(pf));
              jp4_Dn.SetPtEtaPhiE(patJetPfAk04PtJERSmearDn->at(pf),patJetPfAk04Eta_->at(pf),patJetPfAk04Phi_->at(pf),patJetPfAk04EnJERSmearDn->at(pf));
              jet_vector.push_back(v_jets);
              h_bDisct_cMVAv2->Fill(patJetPfAk04BDiscpfCMVA_->at(pf),w);
              if(patJetPfAk04BDiscpfCMVA_->at(pf) > cMVA) cMVA = patJetPfAk04BDiscpfCMVA_->at(pf);
	            cMVAv2=patJetPfAk04BDiscpfCMVA_->at(pf);
              cMVA_v.push_back(patJetPfAk04BDiscpfCMVA_->at(pf));
              jet_cMult=patJetPfAk04cmult_->at(pf);
              no_jets.push_back(n_pat_jets);
              int jflav( abs(patJetPfAk04PartonFlavour_->at(pf)) );
              histos["deepcsv_btag"]->Fill(patJetPfAk04BDiscdeepFb_->at(pf)+patJetPfAk04BDiscdeepFbb_->at(pf),w);
              histos["deepcsv_ctag"]->Fill(patJetPfAk04BDiscdeepFc_->at(pf)+patJetPfAk04BDiscdeepFcc_->at(pf)
              /(1-(patJetPfAk04BDiscdeepFb_->at(pf)+patJetPfAk04BDiscdeepFbb_->at(pf))),w);
              if (cMVAv2>0.4432){
                n_bjets_t++;
              n_bjets_vector.push_back(n_bjets_t);
              }
              bool isBTagged(cMVAv2>0.4432);
              if(!realdata){
                float jptForBtag(v_jets.Pt()>1000. ? 999. : v_jets.Pt()), jetaForBtag(fabs(v_jets.Eta()));
                float expEff(1.0), jetBtagSF(1.0);
                if(abs(jflav)==5){
                expEff    = expBtagEff["b"]->Eval(jptForBtag);
                jetBtagSF    = reader.eval_auto_bounds("central", BTagEntry::FLAV_B,jetaForBtag,jptForBtag);
                jetBtagSF *= expEff>0 ? 1 : 0.;
                }
                else if(abs(jflav)==4){
                expEff    = expBtagEff["c"]->Eval(jptForBtag);
                jetBtagSF    = reader.eval_auto_bounds("central", BTagEntry::FLAV_C,jetaForBtag,jptForBtag);
                jetBtagSF *= expEff>0 ? 1 : 0.;
                }
                else{
                expEff    = expBtagEff["udsg"]->Eval(jptForBtag);
                jetBtagSF    = reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG,jetaForBtag,jptForBtag);
                myBTagSFUtil.modifyBTagsWithSF(isBTagged, jetBtagSF, expEff);
                jetBtagSF *= expEff> 0 ? 1 : 0.;
                }
                myBTagSFUtil.modifyBTagsWithSF(isBTagged,    jetBtagSF,     expEff);
                }

              if(isBTagged){
                n_bjets_v.push_back(n_pat_jets); 
                bjets_v.push_back(v_jets);
                v_bjets=v_jets;
          //      TLorentzVector *lepT, *bjetsT;    
       /*         lepT= &v_lep;
                bjetsT= &v_bjets;
       //       cout<<"n_pat_jets:  "<<n_pat_jets<<", cMVAv2:  "<<cMVAv2<<endl;
                NeutrinoSolver NS(lepT, bjetsT,80.,173.);
                metp4 = TLorentzVector(NS.GetBest(v_met.X(), v_met.Y(), 1., 1., 0.,test));
                if(test > 0. && test < 10000. && N_right_->GetXaxis()->GetXmin()>=0. && N_right_->GetXaxis()->GetXmax()<=150.)
                {
                  //nstest_ = -1.*Log(N_right_->Interpolate(Sqrt(test)));
                  nstest_ = -1. *Log(N_right_->Interpolate(Sqrt(test)));
                 }
                 if (nstest_ < nstestbest_){nstestbest_=nstest_; bjetsT_=&v_bjets;cout<<"I am here:  "<<n_pat_jets<<endl;}

        cout<<",  n_pat_jets:  "<<n_pat_jets<<",   test:  "<<nstestbest_<<",  nstest_:  "<<nstest_<<", bjetsT:  "<<bjetsT_->Pt() <<",  "<<bjetsT->Pt() <<endl;
                */}

              if(!isBTagged){//Fix me, I have some some events with isBTagged==0 but cMVA>0.185
                temp1.push_back(v_jets);
                n_ljets_v.push_back(n_pat_jets);
                ljets_v.push_back(v_jets);
                v_ljets=v_jets;
                }
              }
            }


    outputFile_Objects<<"\n"<< "MET:  "<<endl;
    outputFile_Objects<<" Raw MET (pt, phi):  "<<uncorec_metpt<<",   "<<uncorec_metphi<<endl;
    outputFile_Objects<<" Corrected MET (pt, phi):  "<<v_met.Pt()<<", "<<v_met.Phi()<<endl;

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
          std::vector<float> lepTriggerSF(3,1.0);
          std::vector<float> lepidSF(3,1.0),lepisoSF(3,1.0);
          float trigSF(1.0), trigSFUnc(0.03);
          float idSF(1.0), idSFUnc(0.03);
          float isoSF(1.0), isoSFUnc(0.03);
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
          if(lepEffH.find(lepid_hist)!=lepEffH.end()){
            for(UInt_t il=0; il<TMath::Min((UInt_t)1,(UInt_t)n_lep_v.size()); il++)
            {
              float minPtForEff(lepEffH[lepid_hist]->GetXaxis()->GetXmin()), maxPtForEff( lepEffH[lepid_hist]->GetXaxis()->GetXmax()-0.01 );
              float ptForEff=TMath::Max(TMath::Min(float(fabs(v_lep.Pt())),maxPtForEff),minPtForEff);
              Int_t ptBinForEff=lepEffH[lepid_hist]->GetXaxis()->FindBin(ptForEff);
              float minEtaForEff( lepEffH[lepid_hist]->GetYaxis()->GetXmin() ), maxEtaForEff( lepEffH[lepid_hist]->GetYaxis()->GetXmax()-0.01 );
              float etaForEff=TMath::Max(TMath::Min(float(v_lep.Eta()),maxEtaForEff),minEtaForEff);
              Int_t etaBinForEff=lepEffH[lepid_hist]->GetYaxis()->FindBin(etaForEff);
              idSF=(lepEffH[lepid_hist]->GetBinContent(ptBinForEff,etaBinForEff));
              idSFUnc=(lepEffH[lepid_hist]->GetBinError(ptBinForEff,etaBinForEff));
              lepidSF[0]*=idSF; lepidSF[1]*=(idSF-idSFUnc); lepidSF[2]*=(idSF+idSFUnc);

//              cout<<"this is muID: ptForEff: "<<ptForEff<<",  etaForEff:  "<<etaForEff<<",v_lep.Eta()  "<<v_lep.Eta()<<",  SF: "<<idSF<<endl;
             }
      }
       if(lepEffH.find(eid_hist)!=lepEffH.end()){//pt on Y-axis and Eta on X, 7 changes needed to reverse the order, special look at GetBinContent(ptBinForEff,etaBinForEff)) line.
            for(UInt_t il=0; il<TMath::Min((UInt_t)1,(UInt_t)n_lep_v.size()); il++)
            {
              float minPtForEff(lepEffH[eid_hist]->GetYaxis()->GetXmin()), maxPtForEff( lepEffH[eid_hist]->GetYaxis()->GetXmax()-0.01 );
              float ptForEff=TMath::Max(TMath::Min(float(fabs(v_lep.Pt())),maxPtForEff),minPtForEff);
              Int_t ptBinForEff=lepEffH[eid_hist]->GetYaxis()->FindBin(ptForEff);
              float minEtaForEff( lepEffH[eid_hist]->GetXaxis()->GetXmin() ), maxEtaForEff( lepEffH[eid_hist]->GetXaxis()->GetXmax()-0.01 );
              float etaForEff=TMath::Max(TMath::Min(float(v_lep.Eta()),maxEtaForEff),minEtaForEff);
              Int_t etaBinForEff=lepEffH[eid_hist]->GetXaxis()->FindBin(etaForEff);
              idSF=(lepEffH[eid_hist]->GetBinContent(etaBinForEff, ptBinForEff));
              idSFUnc=(lepEffH[eid_hist]->GetBinError(etaBinForEff, ptBinForEff));
              lepidSF[0]*=idSF; lepidSF[1]*=(idSF-idSFUnc); lepidSF[2]*=(idSF+idSFUnc);
              }
       }
       if(lepEffH.find(esel_hist)!=lepEffH.end()){//Pt on Y-axis and eta on X-axis
             for(UInt_t il=0; il<TMath::Min((UInt_t)1,(UInt_t)n_lep_v.size()); il++)
             {
               float minPtForEff(lepEffH[esel_hist]->GetYaxis()->GetXmin()), maxPtForEff( lepEffH[esel_hist]->GetYaxis()->GetXmax()-0.01 );
               float ptForEff=TMath::Max(TMath::Min(float(fabs(v_lep.Pt())),maxPtForEff),minPtForEff);
               Int_t ptBinForEff=lepEffH[esel_hist]->GetYaxis()->FindBin(ptForEff);
               float minEtaForEff( lepEffH[esel_hist]->GetXaxis()->GetXmin() ), maxEtaForEff( lepEffH[esel_hist]->GetXaxis()->GetXmax()-0.01 );
               float etaForEff=TMath::Max(TMath::Min(float(v_lep.Eta()),maxEtaForEff),minEtaForEff);
               Int_t etaBinForEff=lepEffH[esel_hist]->GetXaxis()->FindBin(etaForEff);
               isoSF=(lepEffH[esel_hist]->GetBinContent(etaBinForEff, ptBinForEff));
               isoSFUnc=(lepEffH[esel_hist]->GetBinError(etaBinForEff, ptBinForEff));
               lepisoSF[0]*=isoSF; lepisoSF[1]*=(isoSF-isoSFUnc); lepisoSF[2]*=(isoSF+isoSFUnc);
              }
      }
       if(lepEffH.find(lepiso_hist)!=lepEffH.end()){
            for(UInt_t il=0; il<TMath::Min((UInt_t)1,(UInt_t)n_lep_v.size()); il++)
            {
              float minEtaForEff(lepEffH[lepiso_hist]->GetYaxis()->GetXmin()), maxEtaForEff( lepEffH[lepiso_hist]->GetYaxis()->GetXmax()-0.01 );
              float etaForEff=TMath::Max(TMath::Min(float(fabs(v_lep.Eta())),maxEtaForEff),minEtaForEff);
              Int_t etaBinForEff=lepEffH[lepiso_hist]->GetYaxis()->FindBin(etaForEff);

              float minPtForEff( lepEffH[lepiso_hist]->GetXaxis()->GetXmin() ), maxPtForEff( lepEffH[lepiso_hist]->GetXaxis()->GetXmax()-0.01 );
              float ptForEff=TMath::Max(TMath::Min(float(v_lep.Pt()),maxPtForEff),minPtForEff);
              Int_t ptBinForEff=lepEffH[lepiso_hist]->GetXaxis()->FindBin(ptForEff);
              isoSF=(lepEffH[lepiso_hist]->GetBinContent(ptBinForEff, etaBinForEff));
              isoSFUnc=(lepEffH[lepiso_hist]->GetBinError(ptBinForEff,etaBinForEff));
              lepisoSF[0]*=isoSF; lepisoSF[1]*=(isoSF-isoSFUnc); lepisoSF[2]*=(isoSF+isoSFUnc);
//              cout<<"this is mu Iso: ptForEff: "<<ptForEff<<",  etaForEff:  "<<etaForEff<<",v_lep.Eta()  "<<v_lep.Eta()<<",  SF: "<<isoSF<<endl;
            }
      }
     if (!realdata)w*=lepTriggerSF[0]*lepidSF[0]*lepisoSF[0];
//     cout<<"this is only lepidSF[0]:  "<<lepidSF[0]<<endl;
/*     w*=lepTriggerSF[0];
     cout<<"this is w1:  "<<w<<endl;
     w*=lepTriggerSF[0]*lepidSF[0];
     cout<<"this is  w2:  "<<w<<endl;
     w*=lepTriggerSF[0]*lepidSF[0]*lepisoSF[0];
     cout<<"this is  w3:  "<<w<<endl;
     w*=lepTriggerSF[0]*lepidSF[0]*lepisoSF[0];
cout<<"this is final w: :::::::::"<<lepisoSF[0]<<endl;
*/
   bool jets2_cut(false);
   if (n_pat_jets > 1)jets2_cut=true;
   bool lep0_cut(false), lep1_cut(false), trans_mW_cut(false),jets_cut(false), bjets_cut(false), ljets_cut(false);
              if (is_mu)if (n_elec_v.size() ==0 && n_elec20_v.size() ==0 )lep0_cut=true;
              if (!is_mu)if (n_muon_v.size() ==0 && n_muon10_v.size() ==0)lep0_cut=true;
              if(n_lep_v.size() == 1 && n_lep10_v.size() ==1)lep1_cut=true;
   double trans_m_w=0.;
   if (lep_vector.size()>0)trans_m_w=sqrt(pow(lep_vector[0].Pt() + v_met.Pt(), 2) - pow(lep_vector[0].Px() + v_met.Px(), 2) - pow(lep_vector[0].Py() + v_met.Py(), 2));
              if (trans_m_w > 50.)trans_mW_cut=true;
              if(n_pat_jets > 3)jets_cut=true;
              if (n_bjets_v.size() > 1)bjets_cut=true;
              if(n_ljets_v.size() > 1)ljets_cut=true;
if (run==283307 && lumi_==213 && event==392668382)
   {
        cout<< run<<":"<<lumi_<<":"<<event<<"met: "<<v_met.Pt()<<", lep Pt(): "<<lep_vector[0].Pt()<<",  mass of W_T: "<<trans_m_w<<endl;
         }


//cout<<"#jets:  "<<n_pat_jets<<", #bjets:  "<<bjets_v.size()<<", #ljets:  "<<n_ljets_v.size()<<endl;
//----------------Rochester algorithm----------------

//cout<<"this is puWeight:  "<<puWeight[0]<<",  lepidSF[0]: "<<lepidSF[0]<<",  lepisoSF[0]: "<<lepisoSF[0]<<", nloSF*mcWeight_:  "<<w<<endl;
    TLorentzVector nuT;
    TLorentzVector Bhad;
    TLorentzVector Blep;
    TLorentzVector Whad;
    TLorentzVector Thad;
    TLorentzVector Wlep;
    TLorentzVector Tlep;
    double costheta=999.;
    if (lep0_cut ){h_events_eachCut->Fill(1.,w);
outputFile_0lep<< run<<":"<< lumi_ << ":"<< event <<","<<endl;
    if (lep1_cut){h_events_eachCut->Fill(2.,w);
outputFile_1lep<< run<<":"<< lumi_ << ":"<< event <<","<<endl;
    h_M_t->Fill(trans_m_w,w);
       if (trans_mW_cut){h_events_eachCut->Fill(3,w);h_events_eachCut_el->Fill(0.,w);
outputFile_W_T<< run<<":"<< lumi_ << ":"<< event <<","<<endl;
          h_no_Jets->Fill(n_pat_jets,w);
          if (jets_cut){ h_events_eachCut->Fill(4.,w);h_events_eachCut_el->Fill(1.,w);
outputFile_4jets<< run<<":"<< lumi_ << ":"<< event <<","<<endl;
              if (bjets_cut){h_events_eachCut->Fill(5.,w);h_events_eachCut_el->Fill(2.,w);
outputFile_2bjets<< run<<":"<< lumi_ << ":"<< event <<","<<endl;
                    h_PU_npT->Fill(PU_npT,w);
                    h_n_ljets_aft->Fill(ljets_v.size(),w);
                    h_n_bjets_aft->Fill(n_bjets_v.size(),w);
                    h_n_cjets->Fill(ncjets,w);
                    h_EvtInfo_NumVtx->Fill(EvtInfo_NumVtx);
                    h_EvtInfo_NumVtx1->Fill(EvtInfo_NumVtx,puWeight[0]);
                    h_EvtInfo_NumVtx2->Fill(EvtInfo_NumVtx,w);
                    h_EvtInfo_NumVtx3->Fill(EvtInfo_NumVtx,lepidSF[0]);
                    h_EvtInfo_NumVtx4->Fill(EvtInfo_NumVtx,lepisoSF[0]);
                    h_EvtInfo_NumVtx_w->Fill(EvtInfo_NumVtx,w);
                    h_M_t_aft->Fill(trans_m_w,w);
    
    int lepChar=lep_charge[0];
    double nuDisc, massDiscr, nuchi2;
    TTBarSolver solver;
    Permutation best_permutation;
outputFile_ljets<< run<<":"<< lumi_ << ":"<< event <<","<<endl;
outputFile_ljets<< "Number of non-bjets:  "<<n_ljets_v.size()<<","<<endl;
    outputFile_TTObjects<< "*************   Event   "<<run<<":"<< lumi_ << ":"<< event <<"   **************\n"<<endl; 
//-----------------------------------------------------------------------//
  double nschi    = numeric_limits<double>::max();
  double res      = numeric_limits<double>::max();
  double nstest   = numeric_limits<double>::max();
  double masstest = numeric_limits<double>::max();
  double btagtest = numeric_limits<double>::max();
  double ptratios = numeric_limits<double>::max();
  double qgtest   = numeric_limits<double>::max();
  bool IsComp=0;
  bool IsWHadComp=0;
  bool IsTHadComp=0;
  bool IsTLepComp=0;
//cout<<"this is metp4.Px():  "<<v_met.Px()<<",,   py:  "<<v_met.Py()<<", pz:  "<<v_met.Pz()<<endl;                          
//cout<<"this is metp4.Px():  "<<metp4.Px()<<",,   py:  "<<metp4.Py()<<", pz:  "<<metp4.Pz()<<endl;                          
//cout<<"this is lep from analyzer: "<<v_bjets.Pt()<<endl;

//-----------------------------------------------------------------------//
 vector<TLorentzVector>AllTypesJets;
     for (unsigned int i=0; i<bjets_v.size(); i++){
          AllTypesJets.push_back(bjets_v[i]);
          }
     for (unsigned int i=0; i<ljets_v.size(); i++){
          AllTypesJets.push_back(ljets_v[i]);
          }

    for (unsigned int k=0; k<bjets_v.size(); k++){
    for (unsigned int l=0; l<bjets_v.size(); l++){
            if (l==k) continue; 

     for (unsigned int i=0; i<AllTypesJets.size(); i++){
            if ( i==k || i==l)continue;
      for (unsigned int j=i+1;j<AllTypesJets.size(); j++){
          if (j==i || j==l || j==k)continue;

//------------------------------

outputFile_TTObjects<<",  bjetsize: "<<bjets_v.size()<<" all jets size: "<<AllTypesJets.size()<<",  k, l, i, j:  "<<k<<", "<<l<<", "<<i<<", "<<j<<" , lep bjets Pt:  "<<bjets_v[k].Pt()<<", had bjets Pt: "<<bjets_v[l].Pt()<<", jet1 Pt: "<<AllTypesJets[i].Pt()<<", jet2: "<<AllTypesJets[j].Pt()<<endl;
outputFile_TTObjects<<" lep Pt:  "<<v_lep.Pt()<<", MET:  "<<v_met.Pt()<<endl;
//------------------------------
//if ((run==283885 && lumi_==643 && event==1178801269)||(run==283885 && lumi_==1036 && event==1876957143)||
//(run==283885 && lumi_==580 && event==1062700251)||(run==283885 && lumi_==1035 && event==1874322805)||
//if(run==283270&& lumi_==199 && event==245933524)
//{
//cout<< run<<":"<<lumi_<<":"<<event<<endl;
//cout<<"  bjetsize: "<<bjets_v.size()<<" all jets size: "<<AllTypesJets.size()<<",  k, l, i, j:  "<<k<<", "<<l<<", "<<i<<", "<<j<<" , lep bjets Pt:  "<<bjets_v[k].Pt()<<", had bjets Pt: "<<bjets_v[l].Pt()<<", jet1 Pt: "<<AllTypesJets[i].Pt()<<", jet2: "<<AllTypesJets[j].Pt()<<",  lep pt: "<<v_lep.Pt()<<", met:  "<<v_met.Pt()<<endl;

        Permutation permutation(&AllTypesJets[i], &AllTypesJets[j], &bjets_v[k], &bjets_v[l], &v_lep, &v_met,lepChar);
    solver.Solve(permutation);
    if(permutation.Prob()  < best_permutation.Prob())
    {
      best_permutation = permutation;
      nuT=best_permutation.Nu();
      IsWHadComp=best_permutation.IsWHadComplete();
      IsTHadComp=best_permutation.IsTHadComplete();
      IsTLepComp=best_permutation.IsTLepComplete();
      IsComp=best_permutation.IsComplete();
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


outputFile_TTObjects<<"IsWHadComp:  "<<IsWHadComp<<", IsTHadComp:  "<<IsTHadComp<<",  IsTLepComp:  "<<IsTLepComp<<"IsComp:  "<<IsComp<<"mass of TTbar:  "<< (Thad+Tlep).M()<< ", nuT: "<<nuT.Pt()<<", Bhad, Blep:  "<<Bhad.Pt()<<", "<<Blep.Pt()<<", Wlep, Whad:  "<<Wlep.Pt()<<", "<<Whad.P() <<endl;
outputFile_TTObjects<<"nuDisc:  "<< nuDisc <<", nuchi2:  "<<nuchi2<<", massDiscr: "<<massDiscr<< endl;
    }
//cout<<" nuchi2:  "<<nuchi2<<",  nuDisc:  "<< nuDisc <<", massDiscr: "<<massDiscr<< endl;
//    cout<<"mass of ttbar:  "<<(Thad+Tlep).M()<<endl;
//}

    }}}}
 
//---------------------------Light jets----------------------------------//
        for (int i=1;i<3;i++){
        h_Pt_lJet[i]->Fill(temp1[i-1].Pt(),w);
        h_Eta_lJet[i]->Fill(temp1[i-1].Eta(),w);
        h_Phi_lJet[i]->Fill(temp1[i-1].Phi(),w);
        h_E_lJet[i]->Fill(temp1[i-1].E(),w);
        h_M_lJet[i]->Fill(temp1[i-1].M(),w);
        }
//cout<<"this is per: nu pz:  "<<nuT.Pz()<<endl;
// Add a new section for testing categories
//----------------------------Muon parameters--------------------------------------------//
        h_met_Pz->Fill(metp4.Pz(),w);
        float Pt_lep = v_lep.Pt();
        h_Pt_lep->Fill(Pt_lep,w);
        h_Eta_lep->Fill(v_lep.Eta(),w);
        h_Phi_lep->Fill(v_lep.Phi(),w);
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
        if (metp4.Pt()>0)
         { h_met_E->Fill(metp4.Pt(),w);
           h_met_Px->Fill(metp4.Px(),w);
           h_met_Py->Fill(metp4.Py(),w);
           h_met_Pz->Fill(metp4.Pz(),w);
         }
//cout<<"Whad.M() after full selec:  "<<Whad.M()<<",  Thad.M():  "<<Thad.M()<<endl;
    if (Thad.M()!=0){
outputFile_RecoSuccessAll<< run<<":"<< lumi_ << ":"<< event <<","<<endl;
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
  RecoSuccess = 1;
    if (costheta!=999.){
      leptop_cos_ttrest->Fill(costheta,w);
      h_M_tt_cos->Fill((Tlep+Thad).M(),fabs(costheta),w);
      }
    }
    if (Thad.M()==0){
outputFile_RecoFail<< run<<":"<< lumi_ << ":"<< event <<","<<endl;
    }
outputFile_RecoSuccess<< run<<":"<< lumi_ << ":"<< event <<","<<endl;
MassTT=(Thad+Tlep).M();
      tree->Branch("Run",&Run,"Run/i");
      tree->Branch("LumiSection",&LumiSection,"LumiSection/l");
      tree->Branch("Event",&Event,"Event/l");
      tree->Branch("RecoSuccess",&RecoSuccess,"RecoSuccess/O");
      tree->Branch("MassTT",&MassTT,"MassTT/F");

tree->Fill();//-- works

}}}}}

//tree = new TTree("tree","tree");
//    tree = new TTree("tree","synchro");
//  tree->Write();
//  tree->Fill();
        }//entries loop
  //tree->Write();
        theFile->Write();
        theFile->Close();
        synchroFile->Write();
        synchroFile->Close();
        }//function loop
