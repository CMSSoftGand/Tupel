#include "simpleReader.h"
using namespace std;

  void template_27_10_16(TString fin,TString fout){
  TFile *f = TFile::Open(fin);
  if (f->IsZombie()) {
     printf("Input root files doesn't open, please have a look:\n");
     return;
     }
  cout<<"this is ff:  "<<fin<<"   ;  "<<fout<<endl;
  TTree *t = (TTree*)f->Get("tupel/MuonTree");
  TFile *theFile = new TFile (fout+".root","RECREATE");
  theFile->cd();  
    bool is_mu (true);
//  bool is_mu (false);
    cout<<"Running Muon Channel: = "<<is_mu<<endl;
  bool run_sys (true);
  TString era("/user/mgul/Higgs_tottbar/Anz_8011/CMSSW_8_0_11/src/Tupel/Tupel/data/era2016/");
  bool no_trig_samp(false);
  if (fin.Contains( "qcd_EMEnriched") || fin.Contains( "WWToLNuQQ" ) || fin.Contains( "tW_nohad" ) || fin.Contains( "TeV_WZ" )|| fin.Contains( "TeV_ZZ" ))no_trig_samp=true;
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
/*
      int mu_id=-999;
      int e_id=-999;
      for (unsigned int i=0; i<St03Id->size(); ++i){ 
        if (fabs(St03Id->at(i)) == 13) mu_id=St03Id->at(i);
        if (fabs(St03Id->at(i)) == 11) e_id=St03Id->at(i);
         }
         */
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
   if (is_mu ){
      TFile *mufIn=TFile::Open(mutrigSFurl);
//      mufIn->cd("IsoMu22_OR_IsoTkMu22_PtEtaBins_Run274094_to_276097");
//      gDirectory->Get("efficienciesDATA");
//         lepEffH[Trig_hist]=(TH2 *)gDirectory->Get("efficienciesDATA/pt_abseta_DATA");
         if (!no_trig_samp)lepEffH[Trig_hist]=(TH2 *)mufIn->Get("IsoMu24_OR_IsoTkMu24");
         else lepEffH[Trig_hist]=(TH2 *)mufIn->Get("IsoMu24_OR_IsoTkMu24__EffData");
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
      if (!no_trig_samp)lepEffH[Trig_hist]=(TH2 *)fIn->Get("Ele32_eta2p1_WPTight_Gsf");
      else lepEffH[Trig_hist]=(TH2 *)fIn->Get("Ele32_eta2p1_WPTight_Gsf__EffData");
      for(auto& it : lepEffH) it.second->SetDirectory(0);
              fIn->Close();

      TFile *eidfIn=TFile::Open(eIDurl);
      lepEffH[id_hist]=(TH2 *)eidfIn->Get("EGamma_SF2D");
      for(auto& it : lepEffH) it.second->SetDirectory(0);
              eidfIn->Close();
      }

//     experimental syst Uncer
     std::vector<TString> expSysts;
     expSysts.push_back("Pileup");
     expSysts.push_back("Trigger");
     if (is_mu)expSysts.push_back("MuEfficiency");
     if (!is_mu)expSysts.push_back("EleEfficiency");
     expSysts.push_back("UncJES");
     expSysts.push_back("UncJER");
     expSysts.push_back("BtagEff");
     expSysts.push_back("CtagEff");
     expSysts.push_back("LtagEff");
     Int_t nExpSysts=expSysts.size();

     if (is_mu){TDirectory *mujets_2btag = theFile->mkdir("mujets_2btag");mujets_2btag->cd();}
     if (!is_mu){TDirectory *eljets_2btag = theFile->mkdir("eljets_2btag");eljets_2btag->cd();}
     TH1::SetDefaultSumw2();
     TH2::SetDefaultSumw2();
     branchAdd(t);

//------------------------------Heavy Higgs----------------------------------------//
  TString sample_name;
  if (fin.Contains( "TTJets_powheg"))sample_name="TT";
  if (fin.Contains( "WJets"))sample_name="WJets";
  if (fin.Contains( "DYJetsToLL"))sample_name="ZJets";
  if (fin.Contains( "V_WZ") || fin.Contains( "V_ZZ"))sample_name="VV";
  if (fin.Contains( "SingleT_t_"))sample_name="tChannel";
  if (fin.Contains( "SingleT_tW"))sample_name="tWChannel";
  if (fin.Contains( "TTZToLLNuNu") || fin.Contains( "TTWJets") )sample_name="TTV";
  if (fin.Contains( "QCD_Pt_") && !fin.Contains( "qcd_EMEnriched"))sample_name="QCDmujets";
  if (fin.Contains( "qcd_EMEnriched"))sample_name="QDCejets";
  std::vector<TString>samples;
  samples.push_back(sample_name);
  std::vector<TString> Up_Dn;
  Up_Dn.push_back("Down");Up_Dn.push_back("Up");
  TString chan, mc_stat;
  chan = is_mu ? "m" : "e";
  mc_stat = is_mu ? "mujets" : "ejets";
  std::vector<TString> sysNames; 
  sysNames.push_back("pileup");
  sysNames.push_back("CMS_scale_j_13TeV");       sysNames.push_back("CMS_res_j_13TeV");
  sysNames.push_back("CMS_METunclustered_13TeV");sysNames.push_back("CMS_eff_b_13TeV");
  sysNames.push_back("CMS_fake_b_13TeV");        sysNames.push_back("CMS_eff_"+chan);
  sysNames.push_back("CMS_eff_trigger_"+chan);   sysNames.push_back("QCDscale_ttbar");
  if (fin.Contains( "TTJets_powheg")){
  char itest[20];
  for (int i =1; i <20; i++)
    {
    sprintf(itest,"%i",i);sysNames.push_back("TT_CMS_httbar_"+mc_stat+"_MCstatBin"+itest);
    }}
  cout<<"this is sample name:  "<<samples[0]<<endl;
  Float_t bins[] = { 250.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 610.0, 640.0, 680.0, 730.0, 800.0, 920.0, 1200.0};
  Int_t  binnum = sizeof(bins)/sizeof(Float_t) - 1;
  std::map<TString, TH1 *> plots1d;
  TString sysName;
  for(size_t n=0; n<sysNames.size(); n++){
    sysName=sysNames[n];
    TString sname  =samples[0];
    for(Int_t ivar=0; ivar<2; ivar++){
    TString nlabel(sname+"_"+sysName + (ivar==0 ? "Up" : "Down"));
    plots1d[nlabel]   = new TH1D(nlabel,"t#bar{t} mass [GeV]",binnum,bins);
    }}
    TH1D* H_mass   = new TH1D("mass_ttbar","t#bar{t} mass [GeV]",binnum,bins);
    TH1D* h_ttbar_avg   = new TH1D("mass_ttbar_avg","t#bar{t} mass [GeV]",binnum,bins);
    TH1D* h_ttbar_pdf   = new TH1D("h_ttbar_pdf","t#bar{t} mass [GeV]",binnum,bins);
    TH1D* leptop_cos_ttrest   = new TH1D("cos_theta","cos theta in t#bar{t} RestF ",25,-1., 1.0);

    std::map<int, TH1 *> h_ttMpdfs;
    char nam_ttbar[110];
    for (int i=0; i<=110; i++){
      sprintf(nam_ttbar,"ttbar_mass%i",i);
      h_ttMpdfs[i] = new TH1F(nam_ttbar,nam_ttbar,binnum,bins);
        }

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
	Int_t nentries(t->GetEntriesFast());
        nentries=10000;
        for (int jentry=0; jentry < nentries; jentry++)
        {
        t->GetEntry(jentry);
        if(jentry%10==0)cout<<" << "<<jentry<<"/"<<nentries<<endl;
        double w=1;
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
          w=puWeight[0];
        
        if (fin.Contains("HToTT-semilep"))if (lheSigEvn!=1 )continue;
        if(!no_trig_samp)if (is_mu) if (HLT_IsoMu24!=1 && HLT_IsoTkMu24!=1)continue;
        if(!no_trig_samp)if (!is_mu)if (HLT_Ele32_eta2p1_WPTight_Gsf!=1)continue;
//        if (fin.Contains( "HToTT-semilep"))if (lheSigEvn!=3)continue;//Resonance: 5%=1, 10%=3, 25%=5, 50%=7;;; Interference: 5%=2, 10%=4, 25%=6, 50%=8;
        if (!realdata)w*=mcWeight_;
        //noise
	      if (first_PV!=1)continue;
          if(Flag_HBHENoiseFilter!=1 || Flag_HBHENoiseIsoFilter!=1)continue;
          if (!no_trig_samp)if (Flag_globalTightHalo2016Filter!=1)continue;
          if (Flag_EcalDeadCellTriggerPrimitiveFilter!=1 || Flag_goodVertices!=1 || Flag_eeBadScFilter!=1)continue;
      int n_pat_elec = 0,n_elec15=0;
      std::vector<TLorentzVector>lep_vector, loose_lep_vector;
      TLorentzVector v_lep;
      vector<unsigned int> n_lep_v;
      vector<int> lep_charge;
      vector<unsigned int> n_lep20_v;
      std::vector<TLorentzVector>el_vector,loose_el_vector;
      vector<int>el_charge;
      vector<unsigned int> n_elec_v;
      vector<unsigned int> n_elec20_v;
      for (unsigned int elec =0; elec < patElecPt_->size(); ++elec){
        if(patElecPt_->at(elec)>20 && fabs(patElecScEta_->at(elec))<2.5 && patElecIdveto_->at(elec) !=0){
            TLorentzVector loose_el_vector_;
            loose_el_vector_.SetPtEtaPhiE(patElecPt_->at(elec),patElecEta_->at(elec),patElecPhi_->at(elec),patElecEnergy_->at(elec));
            loose_el_vector.push_back(loose_el_vector_);
            n_elec20_v.push_back(n_elec15);
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
      std::vector<TLorentzVector>mu_vector, loose_mu_vector;
      vector<int>mu_charge;
      vector<unsigned int> n_muon_v;
      vector<unsigned int> n_muon10_v;
      vector<double> mu_eta_v;
      float Pt_muon=0,Eta_muon=0,Phi_muon=0,E_muon=0,muon_TIso=0;
      int n_muon=0,n_muon15=0;
      for(unsigned int mu=0; mu<patMuonPt_->size();mu++){
          muon_TIso=patMuonPfIsoDbeta_->at(mu);
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
            TLorentzVector v_muon_;
            v_muon_.SetPtEtaPhiE(patMuonPt_->at(mu),patMuonEta_->at(mu), patMuonPhi_->at(mu),patMuonEn_->at(mu));
            loose_mu_vector.push_back(v_muon_);
            n_muon10_v.push_back(n_muon15);
            n_muon15++;
            }
            }
        if (!is_mu){n_lep_v = n_elec_v; n_lep20_v = n_elec20_v; lep_vector = el_vector;loose_lep_vector = loose_el_vector; lep_charge = el_charge; v_lep=v_elec;}
        if (is_mu){n_lep_v = n_muon_v; n_lep20_v = n_muon10_v; lep_vector = mu_vector;loose_lep_vector = loose_mu_vector; lep_charge = mu_charge; v_lep=v_muon;}
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
            vector<TLorentzVector> bjets_v, bjets_vUp,bjets_vDn;
            vector<TLorentzVector> ljets_v, ljets_vUp,ljets_vDn;
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
                if (n_lep_v.size()>0)DR_mu_j= DeltaR(loose_lep_vector[0].Eta(), v_jetsTemp.Eta(), loose_lep_vector[0].Phi(), v_jetsTemp.Phi());
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
              TString UpDn;
              jet_vector.push_back(v_jets);
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
                  jetBtagSF = sfbReaders[0]->eval( BTagEntry::FLAV_C, jetaForBtag, jptForBtag);
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
              if(isBTagged){
                n_bjets_v.push_back(n_pat_bjets); 
                bjets_v.push_back(v_jets);
                bjets_vUp.push_back(jp4_Up);
                bjets_vDn.push_back(jp4_Dn);
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
                ljets_vUp.push_back(jp4_Up);
                ljets_vDn.push_back(jp4_Dn);
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
              if (is_mu)if (n_elec_v.size() !=0 || n_elec20_v.size() !=0 )continue;
              if (!is_mu)if (n_muon_v.size() !=0 || n_muon10_v.size() !=0)continue;
              if(n_lep_v.size() != 1 || n_lep20_v.size() !=1)continue;
double trans_m_w=sqrt(pow(lep_vector[0].Pt() + v_met.Pt(), 2) - pow(lep_vector[0].Px() + v_met.Px(), 2) - pow(lep_vector[0].Py() + v_met.Py(), 2));
              if(no_jets.size() < 4)continue;
              if (bjets_v.size() < 2)continue;
              if(n_ljets_v.size() < 1)continue;

//cout<<"test all the cuts  "<<endl;
//----------------Rochester algorithm----------------
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
//    TLorentzVector ttcm;
    //auto ttcm=0., tlepcm=0., thadcm=0.;
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
   //   auto tlepcm = ttang.tlep().to_CM();
   //   auto thadcm = ttang.thad().to_CM();
     costheta = ttang.unit3D().Dot(ttcm.tlep().unit3D());      
    }
    }}}}
    TLorentzVector ThadDn;
    TLorentzVector TlepDn;
    Permutation best_permutationDn;
    for (unsigned int i=0; i<ljets_vDn.size(); i++){
    for (unsigned int j=i+1;j<ljets_vDn.size(); j++){
      if (i==j)continue;
    for (unsigned int k=0; k<bjets_vDn.size(); k++){
    for (unsigned int l=k+1; l<bjets_vDn.size(); l++){
      if (k==l)continue;
    Permutation permutationDn(&ljets_vDn[i], &ljets_vDn[j], &bjets_vDn[l], &bjets_vDn[k], &v_lep, &v_met,lepChar);
    solver.Solve(permutationDn);
    if(permutationDn.Prob()  < best_permutationDn.Prob())
    {
      best_permutationDn = permutationDn;
      ThadDn=best_permutationDn.THad();
      TlepDn=best_permutationDn.TLep();
    }
    }}}}

    TLorentzVector ThadUp;
    TLorentzVector TlepUp;
    Permutation best_permutationUp;
    for (unsigned int i=0; i<ljets_vUp.size(); i++){
    for (unsigned int j=i+1;j<ljets_vUp.size(); j++){
      if (i==j)continue;
    for (unsigned int k=0; k<bjets_vUp.size(); k++){
    for (unsigned int l=k+1; l<bjets_vUp.size(); l++){
      if (k==l)continue;
    Permutation permutationUp(&ljets_vUp[i], &ljets_vUp[j], &bjets_vUp[l], &bjets_vUp[k], &v_lep, &v_met,lepChar);
    solver.Solve(permutationUp);
    if(permutationUp.Prob()  < best_permutationUp.Prob())
    {
      best_permutationUp = permutationUp;
      ThadUp=best_permutationUp.THad();
      TlepUp=best_permutationUp.TLep();
    }
    }}}}

//-----------------------------------------------
      // Trig SF

/*      int mu_id=-999;
      int e_id=-999;
      for (unsigned int i=0; i<St03Id->size(); ++i){ 
        if (fabs(St03Id->at(i)) == 13) mu_id=St03Id->at(i);
        if (fabs(St03Id->at(i)) == 11) e_id=St03Id->at(i);
         }
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
*/
   if (!realdata){
/*   if (is_mu ){
      TFile *mufIn=TFile::Open(mutrigSFurl);
//      mufIn->cd("IsoMu22_OR_IsoTkMu22_PtEtaBins_Run274094_to_276097");
//      gDirectory->Get("efficienciesDATA");
//         lepEffH[Trig_hist]=(TH2 *)gDirectory->Get("efficienciesDATA/pt_abseta_DATA");
         if (!no_trig_samp)lepEffH[Trig_hist]=(TH2 *)mufIn->Get("IsoMu24_OR_IsoTkMu24");
         else lepEffH[Trig_hist]=(TH2 *)mufIn->Get("IsoMu24_OR_IsoTkMu24__EffData");
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
      if (!no_trig_samp)lepEffH[Trig_hist]=(TH2 *)fIn->Get("Ele32_eta2p1_WPTight_Gsf");
      else lepEffH[Trig_hist]=(TH2 *)fIn->Get("Ele32_eta2p1_WPTight_Gsf__EffData");
      for(auto& it : lepEffH) it.second->SetDirectory(0);
              fIn->Close();

      TFile *eidfIn=TFile::Open(eIDurl);
      lepEffH[id_hist]=(TH2 *)eidfIn->Get("EGamma_SF2D");
      for(auto& it : lepEffH) it.second->SetDirectory(0);
              eidfIn->Close();
      }
  */    
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
     if (!realdata)w*=lepTriggerSF[0]*lepidSF[0]*lepisoSF[0];
/*
     cout<<"this is only w:  "<<w<<endl;
     w*=lepTriggerSF[0];
     cout<<"this is w1:  "<<w<<endl;
     w*=lepTriggerSF[0]*lepidSF[0];
     cout<<"this is  w2:  "<<w<<endl;
     w*=lepTriggerSF[0]*lepidSF[0]*lepisoSF[0];
     cout<<"this is  w3:  "<<w<<endl;
     w*=lepTriggerSF[0]*lepidSF[0]*lepisoSF[0];
cout<<"this is final w: :::::::::"<<w<<endl;
*/  
    
        int pdfsets = 110;
        double newWgtpdf(w);
        double totalWgt=0., alphas_109=0., alphas_110=0.;
        for (int i=0; i<=pdfsets; i++){
          totalWgt+=mcWeights_->at(9+i);
          alphas_109=mcWeights_->at(109);alphas_110=mcWeights_->at(110);
        }
        double wgtAvg=totalWgt/100.;
        float m_tt=(Tlep+Thad).M();
        for (int i=0; i<=pdfsets; i++){
          if(Thad.M()!=0)h_ttMpdfs[i]  ->Fill(m_tt,w*mcWeights_->at(i));
        }

    if (Thad.M()!=0){
      float m_ttbar=(Tlep+Thad).M();
      for(size_t n=0; n<sysNames.size(); n++){
         TString sysName=sysNames[n];
         TString sname  =samples[0];
         for(Int_t isign=0; isign<2; isign++){
           float newWgt(w);
           if (sysName=="pileup"){
           newWgt*=(isign==0 ? puWeight[1]/puWeight[0] : puWeight[2]/puWeight[0]);
           }
           if (sysName=="CMS_scale_j_13TeV"){
           newWgt*=(1.0+(isign==0?-1.:1.)*jesUnc);
           }
         TLorentzVector jp4, jetDiff(0,0,0,0),jetSum(0,0,0,0), varlp4(v_lep);
        if (sysName=="CMS_eff_b_13TeV" ||sysName=="CMS_fake_b_13TeV"){
//        varlp4 = (1.0+(isign==0?-1.:1.)*0.02 ) *lp4;

        float cMVA_new ;  
        std::vector<TLorentzVector> varBJets,varLightJets;
        for (unsigned int ij=0; ij<patJetPfAk04PtJERSmear->size();ij++)
          {
          if (patJetPfAk04PtJERSmear->at(ij)<30. || fabs(patJetPfAk04Eta_->at(ij))>2.4 || patJetPfAk04LooseId_->at(ij)==0)continue;
          TLorentzVector v_jetsTemp;
          v_jetsTemp.SetPtEtaPhiE(patJetPfAk04PtJERSmear->at(ij),patJetPfAk04Eta_->at(ij),patJetPfAk04Phi_->at(ij),patJetPfAk04En_->at(ij));
          if (n_lep_v.size()>0)DR_mu_j= DeltaR(lep_vector[0].Eta(), v_jetsTemp.Eta(), lep_vector[0].Phi(), v_jetsTemp.Phi());
          if(DR_mu_j<0.4)continue;
          int jflav( abs(patJetPfAk04PartonFlavour_->at(ij))) ;
          bool isBTagged(patJetPfAk04BDiscpfCMVA_->at(ij)>0.185);
          jp4.SetPtEtaPhiE(patJetPfAk04PtJERSmear->at(ij),patJetPfAk04Eta_->at(ij),patJetPfAk04Phi_->at(ij),patJetPfAk04En_->at(ij));
          jetDiff -= jp4;
          jetSum += jp4;

          if(!realdata)
          {
          float jptForBtag(jp4.Pt()>1000. ? 999. : jp4.Pt()), jetaForBtag(fabs(jp4.Eta()));
          float expEff(1.0), jetBtagSF(1.0);
          if(jflav==4){
            expEff        = expBtagEff["c"]->Eval(jptForBtag);
            int idx(0);
            if(sysName=="CMS_eff_b_13TeV"){ 
            idx=(isign==0 ? 1 : 2);
            jetBtagSF  = sfbReaders[idx]->eval( BTagEntry::FLAV_C, jetaForBtag, jptForBtag);
            jetBtagSF *= expEff>0 ? expBtagEff["c"]->Eval(jptForBtag)/expBtagEff["c"]->Eval(jptForBtag) : 0.;
            newWgt*=jetBtagSF;         
            }}
          else if(jflav==5){
            expEff=expBtagEff["b"]->Eval(jptForBtag);
            int idx(0);
            if(sysName=="CMS_eff_b_13TeV") 
            {idx=(isign==0 ? 1 : 2);
            jetBtagSF=sfbReaders[idx]->eval( BTagEntry::FLAV_B, jetaForBtag, jptForBtag);
            jetBtagSF *= expEff>0 ? expBtagEff["b"]->Eval(jptForBtag)/expBtagEff["b"]->Eval(jptForBtag) : 0.;
            newWgt*=jetBtagSF;
//            cout<<"ij, idx:  "<<ij<<", "<<idx<<",  expEff:  "<<expEff<<", ->Eval(jptForBtag)  "<<(expBtagEff["b"]->Eval(jptForBtag))/(expBtagEff["b"]->Eval(jptForBtag))<<",   jetBtagSF:  "<<jetBtagSF<<",  newWgt:  "<<newWgt<<endl;
            }
            }
          else{
            expEff=expBtagEff["udsg"]->Eval(jptForBtag);
            int idx(0);
            if(sysName=="CMS_fake_b_13TeV"){ 
            idx=(isign==0 ? 1 : 2);
            jetBtagSF=sflReaders[idx]->eval( BTagEntry::FLAV_UDSG, jetaForBtag, jptForBtag);
            jetBtagSF *= expEff>0 ? expBtagEff["udsg"]->Eval(jptForBtag)/expBtagEff["udsg"]->Eval(jptForBtag) : 0.;
            newWgt*=jetBtagSF;
          }}
          myBTagSFUtil.modifyBTagsWithSF(isBTagged,    jetBtagSF,     expEff);
          }
          if(isBTagged){ varBJets.push_back(jp4);}
          else          varLightJets.push_back(jp4);
          }
         }

          TLorentzVector varMet(metp4);
          varMet += jetDiff;
          varMet += (varlp4-v_lep);
         
          if (sysName=="CMS_eff_m" || sysName=="CMS_eff_e"){ 
          newWgt *= (isign==0 ? lepidSF[1]/lepidSF[0] : lepidSF[2]/lepidSF[0]);
          newWgt *= (isign==0 ? lepisoSF[1]/lepisoSF[0]:lepisoSF[2]/lepisoSF[0]);
          }
          if (sysName=="CMS_eff_trigger_m" || sysName=="CMS_eff_trigger_e"){
            newWgt *= (isign==0 ? lepTriggerSF[1]/lepTriggerSF[0] : lepTriggerSF[2]/lepTriggerSF[0]);
          }
          if (sysName.Contains("MCstatBin"))newWgt=w;
          if (sysName=="QCDscale_ttbar"){
            double newWgt = w;
            for (int b=0; b <= binnum; b++){
            vector<double> pdf_array;
            double binCont = 0.0,aph_s1=0., aph_s2=0.,pdfUnc =0.0,pdfUnc_aphaDn=0.,pdfUnc_aphaUp=0.;
//https://twiki.cern.ch/twiki/bin/viewauth/CMS/LHEReaderCMSSW?rev=7
//https://arxiv.org/pdf/1510.03865v2.pdf ,See equation 20... and use quantile 
            for (int i=9; i<109; i++){
            binCont = h_ttMpdfs[i]->GetBinContent(b);
            pdf_array.push_back(binCont);
            }
            sort(begin(pdf_array), end(pdf_array));
            pdfUnc = (pdf_array.at(83)-pdf_array.at(15))/2;
            pdfUnc_aphaDn=TMath::Sqrt(pdfUnc*pdfUnc + (h_ttMpdfs[109]->GetBinContent(b))*(h_ttMpdfs[109]->GetBinContent(b)));
            pdfUnc_aphaUp=TMath::Sqrt(pdfUnc*pdfUnc + (h_ttMpdfs[110]->GetBinContent(b))*(h_ttMpdfs[110]->GetBinContent(b)));
            newWgt = (isign==0 ? pdfUnc_aphaDn : pdfUnc_aphaUp);
            TString nlabel(sname+"_"+sysName + (isign==0 ? "Down" : "Up"));
            plots1d[nlabel]   ->Fill(m_ttbar,newWgt);
            }
          }
          if (sysName=="CMS_res_j_13TeV"){
             newWgt = w;
             m_ttbar = isign==0 ? (TlepDn+ThadDn).M(): (TlepUp+ThadUp).M();
             TString nlabel(sname+"_"+sysName + (isign==0 ? "Down" : "Up"));
             plots1d[nlabel]->Fill(m_ttbar,newWgt);
            }
         TString nlabel(sname+"_"+sysName + (isign==0 ? "Down" : "Up"));
         if (sysName!="QCDscale_ttbar" && sysName!="CMS_res_j_13TeV")plots1d[nlabel]->Fill(m_ttbar,newWgt);
           
                                  }}

      H_mass ->Fill((Tlep+Thad).M(),w);
      if (costheta!=999.)leptop_cos_ttrest->Fill(costheta,w);
//    h_M_tt_binned->Fill((Tlep+Thad).M(),w);
    }

    

/**///////////////////////////////////////////////////////////////////////////////////
//  exp systematics
/*    if (run_sys){cout<<"sys tests  "<<endl;
    for(int isign=0; isign<2; isign++)
    {float newWgt(w);
    if(puWeight[0]!=0){
      newWgt*=(isign==0 ? puWeight[1]/puWeight[0] : puWeight[2]/puWeight[0]);
      if (Thad.M()>0 &&isign==1)h_PuUp->Fill((Tlep+Thad).M(),newWgt);
      if (Thad.M()>0 &&isign==0)h_PuDn->Fill((Tlep+Thad).M(),newWgt);
      if (Thad.M()>0 &&isign==0)h_PuNn->Fill((Tlep+Thad).M(),    w );
    }}
    for(int isign=0; isign<2; isign++)
    {float newWgt(w);
    newWgt *= (isign==0 ? lepTriggerSF[1]/lepTriggerSF[0] : lepTriggerSF[2]/lepTriggerSF[0]);
    if (Thad.M()>0 &&isign==1)h_TrigUp->Fill((Tlep+Thad).M(),newWgt);
    if (Thad.M()>0 &&isign==0)h_TrigDn->Fill((Tlep+Thad).M(),newWgt);
    if (Thad.M()>0 &&isign==0)h_TrigNn->Fill((Tlep+Thad).M(),     w);
    }
    for(int isign=0; isign<2; isign++)
    {if (Thad.M()>0){
      float newWgt(w);
    newWgt *= (isign==0 ? lepidSF[1]/lepidSF[0] : lepidSF[2]/lepidSF[0]);
    if (isign==1)h_LepEffUp->Fill((Tlep+Thad).M(),newWgt);
    if (isign==0)h_LepEffDn->Fill((Tlep+Thad).M(),newWgt);
    if (isign==0)h_LepEffNn->Fill((Tlep+Thad).M(),     w);
    }}
    for(int isign=0; isign<2; isign++)
    {if (metp4.Pt()>0 && sig_Dn>0){
      float newWgt(w);
    newWgt *= (isign==0 ? Dn/Cn : Up/Cn);
    if (isign==1)h_JERUp->Fill(metp4.Pt(),newWgt);
    if (isign==0)h_JERDn->Fill(metp4.Pt(),newWgt);
    if (isign==0)h_JERNn->Fill(metp4.Pt(),     w);
    }}
    
    for(int isign=0; isign<3; isign++)
    {if (v_jets.Pt()>0){
      float newWgt(w);
    newWgt*=(1.0+(isign==0?-1.:1.)*jesUnc);
    if (isign==1)h_JESUp->Fill(v_jets.Pt(),newWgt);
    if (isign==0)h_JESDn->Fill(v_jets.Pt(),newWgt);
    if (isign==2)h_JESNn->Fill(v_jets.Pt(), w);
    }}
  
        float cMVA_new = cMVA_v[0];  
        bool isBTagged(cMVA_new>0.185);
        TLorentzVector jp4;
        int n_bj=0,n_cj=0,n_lj=0;
        std::vector<TLorentzVector> varBJets,varLightJets;
        for(int isign=0; isign<3; isign++)
        {float newWgt(w);
        if(!realdata)
        {
        for (unsigned int ij=0; ij<patJetPfAk04PtJERSmear->size();ij++)
           {
          if (patJetPfAk04PtJERSmear->at(ij)<30. || fabs(patJetPfAk04Eta_->at(ij))>2.4)continue;
          int jflav( abs(patJetPfAk04PartonFlavour_->at(ij))) ;
          jp4.SetPtEtaPhiE(patJetPfAk04PtJERSmear->at(ij),patJetPfAk04Eta_->at(ij),patJetPfAk04Phi_->at(ij),patJetPfAk04En_->at(ij));
          float jptForBtag(jp4.Pt()>1000. ? 999. : jp4.Pt()), jetaForBtag(fabs(jp4.Eta()));
          float expEff(1.0), jetBtagSF(1.0);
          if(jflav==4){
            expEff        = expBtagEff["c"]->Eval(jptForBtag);
            int idx(0);
            if (isign==0)idx=0;if (isign==1)idx=1;if (isign==2)idx=2;
            jetBtagSF  = sfbReaders[idx]->eval( BTagEntry::FLAV_C, jetaForBtag, jptForBtag);
            jetBtagSF *= expEff>0 ? expBtagEff["c"]->Eval(jptForBtag)/expBtagEff["c"]->Eval(jptForBtag) : 0.;
         
            }
          else if(jflav==5){
            n_bj++;
            expEff=expBtagEff["b"]->Eval(jptForBtag);
            int idx(0);
            if (isign==0)idx=0;if (isign==1)idx=1;if (isign==2)idx=2;
            jetBtagSF=sfbReaders[idx]->eval( BTagEntry::FLAV_B, jetaForBtag, jptForBtag);
            jetBtagSF *= expEff>0 ? expBtagEff["b"]->Eval(jptForBtag)/expBtagEff["b"]->Eval(jptForBtag) : 0.;
            newWgt*=jetBtagSF;
            }
          else{
            expEff=expBtagEff["udsg"]->Eval(jptForBtag);
            int idx(0);
            if (isign==0)idx=0;if (isign==1)idx=1;if (isign==2)idx=2;
            jetBtagSF=sflReaders[idx]->eval( BTagEntry::FLAV_UDSG, jetaForBtag, jptForBtag);
            jetBtagSF *= expEff>0 ? expBtagEff["udsg"]->Eval(jptForBtag)/expBtagEff["udsg"]->Eval(jptForBtag) : 0.;
          }
          myBTagSFUtil.modifyBTagsWithSF(isBTagged,    jetBtagSF,     expEff);
          }
          if(isBTagged) varBJets.push_back(jp4);
          else          varLightJets.push_back(jp4);
          }
          if (Thad.M()>0){
          if (isign==2){h_btagUp->Fill((Tlep+Thad).M(),newWgt);}
          if (isign==0){h_btagNn->Fill((Tlep+Thad).M(),newWgt);}
          if (isign==1){h_btagDn->Fill((Tlep+Thad).M(),newWgt);}
          }
      }
    }
*/


/*    if (!run_sys)continue;
    for(size_t ivar=0; ivar<expSysts.size(); ivar++)
    { 
     TString varName=expSysts[ivar]; 
     bool updateBtag(varName.Contains("tagEff"));
     bool updateJES(varName.Contains("UncJES"));
     bool updateJER(varName.Contains("UncJER"));
     for(int isign=0; isign<2; isign++)
      {
      float newWgt(w);
      if(varName=="Pileup" && puWeight[0]!=0){
        newWgt*=(isign==0 ? puWeight[1]/puWeight[0] : puWeight[2]/puWeight[0]);
      }
      if(varName=="Trigger"){
        newWgt *= (isign==0 ? lepTriggerSF[1]/lepTriggerSF[0] : lepTriggerSF[2]/lepTriggerSF[0]);
      }
      if(varName=="MuEfficiency" || "EleEfficiency"){
        newWgt *= (isign==0 ? lepidSF[1]/lepidSF[0] : lepidSF[2]/lepidSF[0]);
      }
      //jets
      std::vector<TLorentzVector> varBJets,varLightJets;
      TLorentzVector jetDiff(0,0,0,0),jetSum(0,0,0,0);
      vector<TLorentzVector> jp44;
      vector<TLorentzVector> jp41;
      if(! updateBtag){   //Fix Me
        varBJets=bjets_v;
        varLightJets=ljets_v;
//        cout<<"no tag is here: "<<varName<<endl;
      }

      TLorentzVector jp4;
      TLorentzVector jp4Up;
      TLorentzVector jp4Dn;
      int jflav (0),n_bj(0);
      float jecUnc = 0.;
      for (unsigned int ij=0; ij<patJetPfAk04PtJERSmear->size();ij++)
      {
      if (patJetPfAk04PtJERSmear->at(ij)<30. || fabs(patJetPfAk04Eta_->at(ij))>2.4)continue;
      jecUnc = unc_->at(ij);
      jflav= abs(patJetPfAk04PartonFlavour_->at(ij)) ;
      jp4.SetPtEtaPhiE(patJetPfAk04PtJERSmear->at(ij),patJetPfAk04Eta_->at(ij),patJetPfAk04Phi_->at(ij),patJetPfAk04En_->at(ij));
      jp4Up.SetPtEtaPhiE(patJetPfAk04PtJERSmearUp->at(ij),patJetPfAk04Eta_->at(ij),patJetPfAk04Phi_->at(ij),patJetPfAk04En_->at(ij));
      jp4Dn.SetPtEtaPhiE(patJetPfAk04PtJERSmearDn->at(ij),patJetPfAk04Eta_->at(ij),patJetPfAk04Phi_->at(ij),patJetPfAk04En_->at(ij));
      jp41.push_back(jp4);
      jp44.push_back(jp4);
        jetDiff -= jp4;
        jetSum += jp4;
      if(varName == "UncJES"){
        jp44[0] *=(1.0+(isign==0?-1.:1.)*jecUnc);
       }
      if(varName == "UncJER"){
        jp44[0] = (isign==0? jp4Dn:jp4Up); 
        }
      if (updateBtag)
        {
        float cMVA_new = cMVA_v[0];  // Fix me it is 0 or k
        bool isBTagged(cMVA_new>0.185);
        if(!realdata){
          float jptForBtag(jp4.Pt()>1000. ? 999. : jp4.Pt()), jetaForBtag(fabs(jp4.Eta()));
          float expEff(1.0), jetBtagSF(1.0);
          if(jflav==4){
            expEff        = expBtagEff["c"]->Eval(jptForBtag);
            int idx(0);
            if(varName=="CtagEff")   idx=(isign==0 ? 1 : 2);
            jetBtagSF  = sfbReaders[idx]->eval( BTagEntry::FLAV_C, jetaForBtag, jptForBtag);
            jetBtagSF *= expEff>0 ? expBtagEff["c"]->Eval(jptForBtag)/expBtagEff["c"]->Eval(jptForBtag) : 0.;
            }
          else if(jflav==5){
            n_bj++;
            expEff=expBtagEff["b"]->Eval(jptForBtag);
            int idx(0);
            if(varName=="BtagEff")   idx=(isign==0 ? 1 : 2);
            jetBtagSF=sfbReaders[idx]->eval( BTagEntry::FLAV_B, jetaForBtag, jptForBtag);
            jetBtagSF *= expEff>0 ? expBtagEff["b"]->Eval(jptForBtag)/expBtagEff["b"]->Eval(jptForBtag) : 0.;
            }
          else{
            expEff=expBtagEff["udsg"]->Eval(jptForBtag);
            int idx(0);
            if(varName=="LtagEff")   idx=(isign==0 ? 1 : 2);
            jetBtagSF=sflReaders[idx]->eval( BTagEntry::FLAV_UDSG, jetaForBtag, jptForBtag);
            jetBtagSF *= expEff>0 ? expBtagEff["udsg"]->Eval(jptForBtag)/expBtagEff["udsg"]->Eval(jptForBtag) : 0.;
          }
          myBTagSFUtil.modifyBTagsWithSF(isBTagged,    jetBtagSF,     expEff);
          }
          if(isBTagged) varBJets.push_back(jp4);
          else          varLightJets.push_back(jp4);
       } 
      }
         if (Thad.M()!=0) plots2d["ttbar_mass_shapes_exp"]->Fill((Tlep+Thad).M(),2*ivar+isign,newWgt);
        if (metp4.Pt()>0)plots2d["metptshapes_exp"]->Fill(metp4.Pt(),2*ivar+isign,newWgt);
        plots2d["nbjetshapes_exp"]->Fill(n_bj,2*ivar+isign,newWgt);
        if (varBJets.size() > 0.)plots2d["bjetptshapes_exp"]->Fill(varBJets[0].Pt(),2*ivar+isign,newWgt);
        if (varBJets.size() > 0.)h_ptbjets_nom->Fill(varBJets[0].Pt(),w);
        if (jp44.size()>0)plots2d["jetptshapes_exp"]->Fill(jp44[0].Pt(),2*ivar+isign,newWgt);
      }
    }
    */
 // TTree weight_tree;
//  weight_tree = new TTree("tree","tree");
//  weight_tree->Fill();
        }//entries loop
        theFile->Write();
        theFile->Close();
        }//function loop
