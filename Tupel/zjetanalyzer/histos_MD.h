 /*
  TH1D*   = new TH1D("","",);
  ->Sumw2();
  */
  ///////////////////////////////////////
/* int nJetPt_Zinc1jet(11);
  double jetPt_Zinc1jet[12] = {30, 40, 52, 68, 88, 113, 144, 184, 234, 297, 377, 480};
  int nJetPt_Zinc2jet(10);
  double jetPt_Zinc2jet[11] = {30, 40, 52, 68, 88, 113, 144, 184, 234, 297, 377};
*/
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH1D *h_ratio_ptsl_zpt   = new TH1D("ratio_ptsl_zpt","ratio_ptsl_zpt",100,0,5);
  TH2D *h_2dgenjet_recocorrjet =new TH2D("2dgenjet_recocorrjet","2dgenjet_recocorrjet",50,0.,1000,50,0.,1000);
  TH2D *h_2dgenjet_recorawjet =new TH2D("2dgenjet_recorawjet","2dgenjet_recorawjet",50,0.,1000,50,0.,1000);
  TH1D *h_deltaphi_jet_Z= new TH1D("deltaphi_jet_Z","deltaphi_jet_Z",40,0,3.14);
  TH1D *h_match_gen_deltaphi_jet_Z= new TH1D("match_gen_deltaphi_jet_Z","match_gen_deltaphi_jet_Z",40,0,3.14);
  TH1D *h_nomatch_gen_deltaphi_jet_Z= new TH1D("nomatch_gen_deltaphi_jet_Z","nomatch_gen_deltaphi_jet_Z",40,0,3.14);

  TH1D *h_deltaphi_jet_Z_ptj_20_30   = new TH1D("deltaphi_jet_Z_ptj_20_30",   "deltaphi_jet_Z_ptj_20_30"   ,40,0,3.14);
  TH1D *h_deltaphi_jet_Z_ptj_30_50   = new TH1D("deltaphi_jet_Z_ptj_30_50",   "deltaphi_jet_Z_ptj_30_50"   ,40,0,3.14);
  TH1D *h_deltaphi_jet_Z_ptj_50_100  = new TH1D("deltaphi_jet_Z_ptj_50_100",  "deltaphi_jet_Z_ptj_50_100"  ,40,0,3.14);
  TH1D *h_deltaphi_jet_Z_ptj_100_1000= new TH1D("deltaphi_jet_Z_ptj_100_1000","deltaphi_jet_Z_ptj_100_1000",40,0,3.14);

  TH1D *h_deltaphi_jet_Z_ptj_20_24   = new TH1D("deltaphi_jet_Z_ptj_20_24",   "deltaphi_jet_Z_ptj_20_24"   ,40,0,3.14);
  TH1D *h_deltaphi_jet_Z_ptj_24_30   = new TH1D("deltaphi_jet_Z_ptj_24_30",   "deltaphi_jet_Z_ptj_24_30"   ,40,0,3.14);
  TH1D *h_deltaphi_jet_Z_ptj_30_40   = new TH1D("deltaphi_jet_Z_ptj_30_40",   "deltaphi_jet_Z_ptj_30_40"   ,40,0,3.14);
  TH1D *h_deltaphi_jet_Z_ptj_40_50   = new TH1D("deltaphi_jet_Z_ptj_40_50",   "deltaphi_jet_Z_ptj_40_50"   ,40,0,3.14);
  TH1D *h_deltaphi_jet_Z_ptj_50_68   = new TH1D("deltaphi_jet_Z_ptj_50_68",   "deltaphi_jet_Z_ptj_50_68"   ,40,0,3.14);


  TH1D *h_match_gen_deltaphi_jet_Z_ptj_20_30   = new TH1D("match_gen_deltaphi_jet_Z_ptj_20_30",   "match_gen_deltaphi_jet_Z_ptj_20_30"   ,40,0,3.14);
  TH1D *h_match_gen_deltaphi_jet_Z_ptj_30_50   = new TH1D("match_gen_deltaphi_jet_Z_ptj_30_50",   "match_gen_deltaphi_jet_Z_ptj_30_50"   ,40,0,3.14);
  TH1D *h_match_gen_deltaphi_jet_Z_ptj_50_100  = new TH1D("match_gen_deltaphi_jet_Z_ptj_50_100",  "match_gen_deltaphi_jet_Z_ptj_50_100"  ,40,0,3.14);
  TH1D *h_match_gen_deltaphi_jet_Z_ptj_100_1000= new TH1D("match_gen_deltaphi_jet_Z_ptj_100_1000","match_gen_deltaphi_jet_Z_ptj_100_1000",40,0,3.14);
  TH1D *h_nomatch_gen_deltaphi_jet_Z_ptj_20_30   = new TH1D("nomatch_gen_deltaphi_jet_Z_ptj_20_30",   "nomatch_gen_deltaphi_jet_Z_ptj_20_30"   ,40,0,3.14);
  TH1D *h_nomatch_gen_deltaphi_jet_Z_ptj_30_50   = new TH1D("nomatch_gen_deltaphi_jet_Z_ptj_30_50",   "nomatch_gen_deltaphi_jet_Z_ptj_30_50"   ,40,0,3.14);
  TH1D *h_nomatch_gen_deltaphi_jet_Z_ptj_50_100  = new TH1D("nomatch_gen_deltaphi_jet_Z_ptj_50_100",  "nomatch_gen_deltaphi_jet_Z_ptj_50_100"  ,40,0,3.14);
  TH1D *h_nomatch_gen_deltaphi_jet_Z_ptj_100_1000= new TH1D("nomatch_gen_deltaphi_jet_Z_ptj_100_1000","nomatch_gen_deltaphi_jet_Z_ptj_100_1000",40,0,3.14);


  TH1D *h_inc_deltaphi_jet_Z= new TH1D("inc_deltaphi_jet_Z","deltaphi_jet_Z",40,0,3.14);
  TH1D *h_inc_match_gen_deltaphi_jet_Z= new TH1D("inc_match_gen_deltaphi_jet_Z","inc_match_gen_deltaphi_jet_Z",40,0,3.14);
  TH1D *h_inc_nomatch_gen_deltaphi_jet_Z= new TH1D("inc_nomatch_gen_deltaphi_jet_Z","inc_nomatch_gen_deltaphi_jet_Z",40,0,3.14);
  TH1D *h_inc_deltaphi_jet_Z_ptj_20_30   = new TH1D("inc_deltaphi_jet_Z_ptj_20_30",   "inc_deltaphi_jet_Z_ptj_20_30"   ,40,0,3.14);
  TH1D *h_inc_deltaphi_jet_Z_ptj_30_50   = new TH1D("inc_deltaphi_jet_Z_ptj_30_50",   "inc_deltaphi_jet_Z_ptj_30_50"   ,40,0,3.14);
  TH1D *h_inc_deltaphi_jet_Z_ptj_50_100  = new TH1D("inc_deltaphi_jet_Z_ptj_50_100",  "inc_deltaphi_jet_Z_ptj_50_100"  ,40,0,3.14);
  TH1D *h_inc_deltaphi_jet_Z_ptj_100_1000= new TH1D("inc_deltaphi_jet_Z_ptj_100_1000","inc_deltaphi_jet_Z_ptj_100_1000",40,0,3.14);

  TH1D *h_inc_1j_deltaphi_jettot_Z= new TH1D("inc_1j_deltaphi_jettot_Z","deltaphi_jet_Z",40,0,3.14);
  TH1D *h_inc_1j_match_gen_deltaphi_jettot_Z= new TH1D("inc_1j_match_gen_deltaphi_jettot_Z","inc_match_gen_deltaphi_jet_Z",40,0,3.14);
  TH1D *h_inc_1j_nomatch_gen_deltaphi_jettot_Z= new TH1D("inc_1j_nomatch_gen_deltaphi_jettot_Z","inc_nomatch_gen_deltaphi_jet_Z",40,0,3.14);
  TH1D *h_inc_1j_deltaphi_jettot_Z_ptj_20_30   = new TH1D("inc_1j_deltaphi_jettot_Z_ptj_20_30",   "inc_deltaphi_jet_Z_ptj_20_30"   ,40,0,3.14);
  TH1D *h_inc_1j_deltaphi_jettot_Z_ptj_30_50   = new TH1D("inc_1j_deltaphi_jettot_Z_ptj_30_50",   "inc_deltaphi_jet_Z_ptj_30_50"   ,40,0,3.14);
  TH1D *h_inc_1j_deltaphi_jettot_Z_ptj_50_100  = new TH1D("inc_1j_deltaphi_jettot_Z_ptj_50_100",  "inc_deltaphi_jet_Z_ptj_50_100"  ,40,0,3.14);
  TH1D *h_inc_1j_deltaphi_jettot_Z_ptj_100_1000= new TH1D("inc_1j_deltaphi_jettot_Z_ptj_100_1000","inc_deltaphi_jet_Z_ptj_100_1000",40,0,3.14);




  TH1D *h_inc_match_gen_deltaphi_jet_Z_ptj_20_30   = new TH1D("inc_match_gen_deltaphi_jet_Z_ptj_20_30",   "inc_match_gen_deltaphi_jet_Z_ptj_20_30"   ,40,0,3.14);
  TH1D *h_inc_match_gen_deltaphi_jet_Z_ptj_30_50   = new TH1D("inc_match_gen_deltaphi_jet_Z_ptj_30_50",   "inc_match_gen_deltaphi_jet_Z_ptj_30_50"   ,40,0,3.14);
  TH1D *h_inc_match_gen_deltaphi_jet_Z_ptj_50_100  = new TH1D("inc_match_gen_deltaphi_jet_Z_ptj_50_100",  "inc_match_gen_deltaphi_jet_Z_ptj_50_100"  ,40,0,3.14);
  TH1D *h_inc_match_gen_deltaphi_jet_Z_ptj_100_1000= new TH1D("inc_match_gen_deltaphi_jet_Z_ptj_100_1000","inc_match_gen_deltaphi_jet_Z_ptj_100_1000",40,0,3.14);
  TH1D *h_inc_nomatch_gen_deltaphi_jet_Z_ptj_20_30   = new TH1D("inc_nomatch_gen_deltaphi_jet_Z_ptj_20_30",   "inc_nomatch_gen_deltaphi_jet_Z_ptj_20_30"   ,40,0,3.14);
  TH1D *h_inc_nomatch_gen_deltaphi_jet_Z_ptj_30_50   = new TH1D("inc_nomatch_gen_deltaphi_jet_Z_ptj_30_50",   "inc_nomatch_gen_deltaphi_jet_Z_ptj_30_50"   ,40,0,3.14);
  TH1D *h_inc_nomatch_gen_deltaphi_jet_Z_ptj_50_100  = new TH1D("inc_nomatch_gen_deltaphi_jet_Z_ptj_50_100",  "inc_nomatch_gen_deltaphi_jet_Z_ptj_50_100"  ,40,0,3.14);
  TH1D *h_inc_nomatch_gen_deltaphi_jet_Z_ptj_100_1000= new TH1D("inc_nomatch_gen_deltaphi_jet_Z_ptj_100_1000","inc_nomatch_gen_deltaphi_jet_Z_ptj_100_1000",40,0,3.14);

  TH1D *h_inc_1j_match_gen_deltaphi_jettot_Z_ptj_20_30   = new TH1D("inc_1j_match_gen_deltaphi_jettot_Z_ptj_20_30",   "inc_1j_match_gen_deltaphi_jettot_Z_ptj_20_30"   ,40,0,3.14);
  TH1D *h_inc_1j_match_gen_deltaphi_jettot_Z_ptj_30_50   = new TH1D("inc_1j_match_gen_deltaphi_jettot_Z_ptj_30_50",   "inc_1j_match_gen_deltaphi_jettot_Z_ptj_30_50"   ,40,0,3.14);
  TH1D *h_inc_1j_match_gen_deltaphi_jettot_Z_ptj_50_100  = new TH1D("inc_1j_match_gen_deltaphi_jettot_Z_ptj_50_100",  "inc_1j_match_gen_deltaphi_jettot_Z_ptj_50_100"  ,40,0,3.14);
  TH1D *h_inc_1j_match_gen_deltaphi_jettot_Z_ptj_100_1000= new TH1D("inc_1j_match_gen_deltaphi_jettot_Z_ptj_100_1000","inc_1j_match_gen_deltaphi_jettot_Z_ptj_100_1000",40,0,3.14);
  TH1D *h_inc_1j_nomatch_gen_deltaphi_jettot_Z_ptj_20_30   = new TH1D("inc_1j_nomatch_gen_deltaphi_jettot_Z_ptj_20_30",   "inc_1j_nomatch_gen_deltaphi_jettot_Z_ptj_20_30"   ,40,0,3.14);
  TH1D *h_inc_1j_nomatch_gen_deltaphi_jettot_Z_ptj_30_50   = new TH1D("inc_1j_nomatch_gen_deltaphi_jettot_Z_ptj_30_50",   "inc_1j_nomatch_gen_deltaphi_jettot_Z_ptj_30_50"   ,40,0,3.14);
  TH1D *h_inc_1j_nomatch_gen_deltaphi_jettot_Z_ptj_50_100  = new TH1D("inc_1j_nomatch_gen_deltaphi_jettot_Z_ptj_50_100",  "inc_1j_nomatch_gen_deltaphi_jettot_Z_ptj_50_100"  ,40,0,3.14);
  TH1D *h_inc_1j_nomatch_gen_deltaphi_jettot_Z_ptj_100_1000= new TH1D("inc_1j_nomatch_gen_deltaphi_jettot_Z_ptj_100_1000","inc_1j_nomatch_gen_deltaphi_jettot_Z_ptj_100_1000",40,0,3.14);

TH1D * tackman_tau_cm_sum[100];
TH1D * tackman_tau_sum[100];
TH1D * tackman_tauc_cm_sum[100];
TH1D * tackman_tauc_sum[100];

TH1D * tackman_tau_cm_max[100];
TH1D * tackman_tau_max[100];
TH1D * tackman_tauc_cm_max[100];
TH1D * tackman_tauc_max[100];

TProfile* h_ptj_ratio_before  = new TProfile("ptj_ratio_before","Reco/gen vs eta",50,-5.1,5.1,0.,3.0);
TProfile* h_ptj_ratio_after  = new TProfile("ptj_ratio_after","Reco/gen vs eta",50,-5.1,5.1,0.,3.0);
TProfile* h_ptj_ratio_before_50  = new TProfile("ptj_ratio_before_50","Reco/gen vs eta",50,-5.1,5.1,0.,3.0);
TProfile* h_ptj_ratio_before_50_70  = new TProfile("ptj_ratio_before_50_70","Reco/gen vs eta",50,-5.1,5.1,0.,3.0);
TProfile* h_ptj_ratio_after_50  = new TProfile("ptj_ratio_after_50","Reco/gen vs eta",50,-5.1,5.1,0.,3.0);
TProfile* h_ptj_ratio_before_100  = new TProfile("ptj_ratio_before_100","Reco/gen vs eta",50,-5.1,5.1,0.,3.0);
TProfile* h_ptj_ratio_after_100  = new TProfile("ptj_ratio_after_100","Reco/gen vs eta",50,-5.1,5.1,0.,3.0);

TProfile* h_ptj_ratio2_before  = new TProfile("ptj_ratio2_before","Reco/gen vs pt",50,0.,1000.,0.,3.0);
TProfile* h_ptj_ratio2Raw_before  = new TProfile("ptj_ratio2Raw_before","Reco/gen vs pt",50,0.,1000.,0.,3.0);
TProfile* h_ptj_ratio2_after  = new TProfile("ptj_ratio2_after","Reco/gen vs pt",50,0.,1000.,0.,3.0);
TProfile* h_ptj_ratio2_before_50  = new TProfile("ptj_ratio2_before_50","Reco/gen vs pt",50,0.,1000.,0.,3.0);
TProfile* h_ptj_ratio2_after_50  = new TProfile("ptj_ratio2_after_50","Reco/gen vs pt",50,0.,1000.,0.,3.0);
TProfile* h_ptj_ratio2_before_100  = new TProfile("ptj_ratio2_before_100","Reco/gen vs pt",50,0.,1000.,0.,3.0);
TProfile* h_ptj_ratio2_after_100  = new TProfile("ptj_ratio2_after_100","Reco/gen vs pt",50,0.,1000.,0.,3.0);

TProfile* h_ptj_ratioRaw_before  = new TProfile("ptj_ratioRaw_before","Reco/gen vs pt",50,-5.1,5.1,0.,3.0);

TProfile* h_ptj_ratioRaw_after  = new TProfile("ptj_ratioRaw_after","Reco/gen vs pt",50,-5.1,5.1,0.,3.0);
TProfile* h_ptj_ratioRaw_before_50  = new TProfile("ptj_ratioRaw_before_50","Reco/gen vs pt",50,-5.1,5.1,0.,3.0);
TProfile* h_ptj_ratioRaw_before_50_70  = new TProfile("ptj_ratioRaw_before_50_70","Reco/gen vs eta",50,-5.1,5.1,0.,3.0);
TProfile* h_ptj_ratioRaw_after_50  = new TProfile("ptj_ratioRaw_after_50","Reco/gen vs pt",50,-5.1,5.1,0.,3.0);
TProfile* h_ptj_ratioRaw_before_100  = new TProfile("ptj_ratioRaw_before_100","Reco/gen vs pt",50,-5.1,5.1,0.,3.0);
TProfile* h_ptj_ratioRaw_after_100  = new TProfile("ptj_ratioRaw_after_100","Reco/gen vs pt",50,-5.1,5.1,0.,3.0);

  TH1D* h_weight_plus_sherpa= new TH1D("weight_plus_sherpa","weight_plus_sherpa",150,0.,15.);
  TH1D* h_weight_minus_sherpa= new TH1D("weight_minus_sherpa","weight_minus_sherpa",150,0.,15.);
  TH1D* h_ljet_rapidity= new TH1D("ljet_rapidity","ljet_rapidity",50,-5.1,5.1);
  TH1D* h_ljet_rapidity_50= new TH1D("ljet_rapidity_50","ljet_rapidity_50",50,-5.1,5.1);

  TH1D* h_ljet_eta= new TH1D("ljet_eta","ljet_eta",50,-5.1,5.1);
  TH1D* h_sljet_rapidity= new TH1D("sljet_rapidity","sljet_rapidity",50,-5.1,5.1);
  TH1D* h_sljet_eta= new TH1D("sljet_eta","sljet_eta",50,-5.1,5.1);

  TH1D* h_ljet_rapidity_0= new TH1D("ljet_rapidity_0","ljet_rapidity_0",50,-5.1,5.1);
  TH1D* h_ljet_eta_0= new TH1D("ljet_eta_0","ljet_eta_0",50,-5.1,5.1);
  TH1D* h_ljet_rapidity_1= new TH1D("ljet_rapidity_1","ljet_rapidity_1",50,-5.1,5.1);
  TH1D* h_ljet_eta_1= new TH1D("ljet_eta_1","ljet_eta_1",50,-5.1,5.1);
  TH1D* h_ljet_rapidity_2= new TH1D("ljet_rapidity_2","ljet_rapidity_2",50,-5.1,5.1);
  TH1D* h_ljet_eta_2= new TH1D("ljet_eta_2","ljet_eta_2",50,-5.1,5.1);
  TH1D* h_ljet_rapidity_3= new TH1D("ljet_rapidity_3","ljet_rapidity_3",50,-5.1,5.1);
  TH1D* h_ljet_eta_3= new TH1D("ljet_eta_3","ljet_eta_3",50,-5.1,5.1);
  TH1D* h_ljet_rapidity_4= new TH1D("ljet_rapidity_4","ljet_rapidity_4",50,-5.1,5.1);
  TH1D* h_ljet_eta_4= new TH1D("ljet_eta_4","ljet_eta_4",50,-5.1,5.1);
  TH1D* h_ljet_rapidity_5= new TH1D("ljet_rapidity_5","ljet_rapidity_5",50,-5.1,5.1);
  TH1D* h_ljet_eta_5= new TH1D("ljet_eta_5","ljet_eta_5",50,-5.1,5.1);
  TH1D* h_ljet_rapidity_6= new TH1D("ljet_rapidity_6","ljet_rapidity_6",50,-5.1,5.1);
  TH1D* h_ljet_eta_6= new TH1D("ljet_eta_6","ljet_eta_6",50,-5.1,5.1);
  TH1D* h_ljet_rapidity_7= new TH1D("ljet_rapidity_7","ljet_rapidity_7",50,-5.1,5.1);
  TH1D* h_ljet_eta_7= new TH1D("ljet_eta_7","ljet_eta_7",50,-5.1,5.1);
  TH1D* h_ljet_rapidity_8= new TH1D("ljet_rapidity_8","ljet_rapidity_8",50,-5.1,5.1);
  TH1D* h_ljet_eta_8= new TH1D("ljet_eta_8","ljet_eta_8",50,-5.1,5.1);
  TH1D* h_ljet_rapidity_9= new TH1D("ljet_rapidity_9","ljet_rapidity_9",50,-5.1,5.1);
  TH1D* h_ljet_eta_9= new TH1D("ljet_eta_9","ljet_eta_9",50,-5.1,5.1);
  TH1D* h_ljet_rapidity_10= new TH1D("ljet_rapidity_10","ljet_rapidity_10",50,-5.1,5.1);
  TH1D* h_ljet_eta_10= new TH1D("ljet_eta_10","ljet_eta_10",50,-5.1,5.1);
  TH1D* h_ljet_rapidity_11= new TH1D("ljet_rapidity_11","ljet_rapidity_11",50,-5.1,5.1);
  TH1D* h_ljet_eta_11= new TH1D("ljet_eta_11","ljet_eta_11",50,-5.1,5.1);
  TH1D* h_ljet_rapidity_12= new TH1D("ljet_rapidity_12","ljet_rapidity_12",50,-5.1,5.1);
  TH1D* h_ljet_eta_12= new TH1D("ljet_eta_12","ljet_eta_12",50,-5.1,5.1);



  TH1D* h_deltar = new TH1D("h_deltar","h_deltar",50,0.0,2.0);
  int nJetPt_Zinc1jet(22);
  double jetPt_Zinc1jet[23] = {20, 24, 30, 39, 49, 60, 72, 85, 100, 117, 136, 157, 187, 220, 258, 300, 350, 400, 450, 500, 590, 700, 1000};

 int nJetPt_Zinc2jet(21);
  double jetPt_Zinc2jet[22] = {20, 24, 30, 39, 49, 60, 72, 85, 100, 117, 136, 157, 187, 220, 258, 300, 350, 400, 450, 500, 590, 800};


  TH1D* h_nup = new TH1D("nup","nup",6,-0.5,5.5);
  TH1D* h_nup_elec = new TH1D("nup_elec","nup",6,-0.5,5.5);
  TH1D* h_nup_mu = new TH1D("nup_mu","nup",6,-0.5,5.5);
 TH1D* h_nup_elec_with_cuts_noweigt = new TH1D("nup_elec_with_cuts_noweigt","nup",6,-0.5,5.5);
  TH1D* h_nup_mu_with_cuts_noweigt = new TH1D("nup_mu_with_cuts_noweigt","nup",6,-0.5,5.5);

 TH1D* h_nup_elec_with_cuts = new TH1D("nup_elec_with_cuts","nup",6,-0.5,5.5);
  TH1D* h_nup_mu_with_cuts = new TH1D("nup_mu_with_cuts","nup",6,-0.5,5.5);
  TH1D* h_dE_gen_raw = new TH1D("dE_gen_raw","dE_gen_raw",50,-30,30);
  TH1D* h_dE_gen_corr = new TH1D(" dE_gen_corr","dE_gen_corr",50,-30,30);
  TH1D* h_dE_gen_corr_rescorr= new TH1D(" dE_gen_corr_rescorr","dE_gen_corr_rescorr",50,-30,30);
  TH1D* h_mjjjj = new TH1D("mjjjj","M(jjjj)",100, 10.,200.);
  TH1D* h_mjjjj_BEF = new TH1D("mjjjj_BEF","M(jjjj) BEF",100, 10.,200.);
 // TH1D* h_hist_efe  = new TH1D("h_hist_efe","#Delta R (Gj1,J1)",50,0.,1.);
 // TH1D* h_hist_efee  = new TH1D("h_hist_efee","#Delta R (Gj,J)",50,0.,1.);
 TH1D* h_deltar_reco_gen_jet  = new TH1D("deltar_reco_gen_jet","#Delta R (Gj,J)",50,0.,1.);
 TH1D* h_deltar_reco_gen_jet_matched  = new TH1D("deltar_reco_gen_jet_matched","#Delta R (Gj,J)",50,0.,1.);

  TH1D* h_nvtx = new TH1D("nvtx","nvtx",40,0,40);

  TH1D* h_nvtx_noWeight = new TH1D("nvtx_noWeight","nvtx_noWeight",40,0,40);
  
  TH1D* h_lpt = new TH1D("lpt","lpt",200,0,1000);
  TH1D*h_leta   = new TH1D("leta","leta",50,-5,5);
  TH1D*h_lphi   = new TH1D("lphi","lphi",50,-5,5);

  TH1D* h_slpt = new TH1D("slpt","slpt",200,0,1000);
  TH1D*h_sleta   = new TH1D("sleta","sleta",50,-5,5);
  TH1D*h_slphi   = new TH1D("slphi","slphi",50,-5,5);


  TH1D* h_MZ = new TH1D("MZ","M(Z)#rightarrow #mu#mu",40, 71,111.);
  TH1D* h_Zpt = new TH1D("Zpt","Zpt",200,0,1000);
  TH1D* h_Zrapidity = new TH1D("Zrapidity","Zrapidity",50,-2.5,2.5);
  TH1D* h_Zphi = new TH1D("Zphi","Zphi",12,-3.2,3.2);

  TH1D* h_MZ_0 = new TH1D("MZ_0","M(Z)#rightarrow #mu#mu",40, 71,111.);
  TH1D* h_Zpt_0 = new TH1D("Zpt_0","Zpt",200,0,1000);
  TH1D* h_Zrapidity_0 = new TH1D("Zrapidity_0","Zrapidity",50,-2.5,2.5);
  TH1D* h_Zphi_0 = new TH1D("Zphi_0","Zphi",12,-3.2,3.2);

  TH1D* h_MZ_1 = new TH1D("MZ_1","M(Z)#rightarrow #mu#mu",40, 71,111.);
  TH1D* h_Zpt_1 = new TH1D("Zpt_1","Zpt",200,0,1000);
  TH1D* h_Zrapidity_1 = new TH1D("Zrapidity_1","Zrapidity",50,-2.5,2.5);
  TH1D* h_Zphi_1 = new TH1D("Zphi_1","Zphi",12,-3.2,3.2);

  TH1D* h_MZ_2 = new TH1D("MZ_2","M(Z)#rightarrow #mu#mu",40, 71,111.);
  TH1D* h_Zpt_2 = new TH1D("Zpt_2","Zpt",200,0,1000);
  TH1D* h_Zrapidity_2 = new TH1D("Zrapidity_2","Zrapidity",50,-2.5,2.5);
  TH1D* h_Zphi_2 = new TH1D("Zphi_2","Zphi",12,-3.2,3.2);

  TH1D* h_MZ_3 = new TH1D("MZ_3","M(Z)#rightarrow #mu#mu",40, 71,111.);
  TH1D* h_Zpt_3 = new TH1D("Zpt_3","Zpt",200,0,1000);
  TH1D* h_Zrapidity_3 = new TH1D("Zrapidity_3","Zrapidity",50,-2.5,2.5);
  TH1D* h_Zphi_3 = new TH1D("Zphi_3","Zphi",12,-3.2,3.2);

  TH1D* h_MZ_4 = new TH1D("MZ_4","M(Z)#rightarrow #mu#mu",40, 71,111.);
  TH1D* h_Zpt_4 = new TH1D("Zpt_4","Zpt",200,0,1000);
  TH1D* h_Zrapidity_4 = new TH1D("Zrapidity_4","Zrapidity",50,-2.5,2.5);
  TH1D* h_Zphi_4 = new TH1D("Zphi_4","Zphi",12,-3.2,3.2);

 TH1D* h_Zpt_pc = new TH1D("Zpt_pc","Zpt_pc",200,0,1000);
  TH1D* h_Zphi_pc = new TH1D("Zphi_pc","Zphi_pc",12,-3.2,3.2);
  TH1D* h_Zpt_p0 = new TH1D("Zpt_p0","Zpt_p0",200,0,1000);
  TH1D* h_Zphi_p0 = new TH1D("Zphi_p0","Zphi_p0",12,-3.2,3.2);
  TH1D* h_Zpt_p1 = new TH1D("Zpt_p1","Zpt_p1",200,0,1000);
  TH1D* h_Zphi_p1 = new TH1D("Zphi_p1","Zphi_p1",12,-3.2,3.2);
  TH1D* h_Zpt_p2 = new TH1D("Zpt_p2","Zpt_p2",200,0,1000);
  TH1D* h_Zphi_p2 = new TH1D("Zphi_p2","Zphi_p2",12,-3.2,3.2);


 TH1D* h_Zpt_mc = new TH1D("Zpt_mc","Zpt_mc",200,0,1000);
  TH1D* h_Zphi_mc = new TH1D("Zphi_mc","Zphi_mc",12,-3.2,3.2);
 TH1D* h_Zpt_m0 = new TH1D("Zpt_m0","Zpt_m0",200,0,1000);
  TH1D* h_Zphi_m0 = new TH1D("Zphi_m0","Zphi_m0",12,-3.2,3.2);
  TH1D* h_Zpt_m1 = new TH1D("Zpt_m1","Zpt_m1",200,0,1000);
  TH1D* h_Zphi_m1 = new TH1D("Zphi_m1","Zphi_m1",12,-3.2,3.2);
  TH1D* h_Zpt_m2 = new TH1D("Zpt_m2","Zpt_m2",200,0,1000);
  TH1D* h_Zphi_m2 = new TH1D("Zphi_m2","Zphi_m2",12,-3.2,3.2);


  TH1D*h_felec_pt_p0    = new TH1D("felec_pt_p0 ","felec_pt_p0",200,0,1000);
  TH1D*h_felec_eta_p0   = new TH1D("felec_eta_p0","felec_eta_p0",50,-5,5);
  TH1D*h_felec_phi_p0   = new TH1D("felec_phi_p0","felec_phi_p0",50,-5,5);
  TH1D*h_celec_pt_p0    = new TH1D("celec_pt_p0 ","celec_pt_p0",200,0,1000);
  TH1D*h_celec_eta_p0   = new TH1D("celec_eta_p0","celec_eta_p0",50,-5,5);
  TH1D*h_celec_phi_p0   = new TH1D("celec_phi_p0","celec_phi_p0",50,-5,5);

  TH1D*h_lelec_pt_p0    = new TH1D("lelec_pt_p0  ","lelec_pt_p0  ",200,0,1000);
  TH1D*h_lelec_eta_p0   = new TH1D("lelec_eta_p0 ","lelec_eta_p0 ",50,-5,5);
  TH1D*h_lelec_phi_p0   = new TH1D("lelec_phi_p0 ","lelec_phi_p0 ",50,-5,5);
  TH1D*h_slelec_pt_p0   = new TH1D("slelec_pt_p0 ","slelec_pt_p0 ",200,0,1000);
  TH1D*h_slelec_eta_p0  = new TH1D("slelec_eta_p0","slelec_eta_p0",50,-5,5);
  TH1D*h_slelec_phi_p0  = new TH1D("slelec_phi_p0","slelec_phi_p0",50,-5,5);

  TH1D*h_felec_pt_p1    = new TH1D("felec_pt_p1 ","felec_pt_p1",200,0,1000);
  TH1D*h_felec_eta_p1   = new TH1D("felec_eta_p1","felec_eta_p1",50,-5,5);
  TH1D*h_felec_phi_p1   = new TH1D("felec_phi_p1","felec_phi_p1",50,-5,5);
  TH1D*h_celec_pt_p1    = new TH1D("celec_pt_p1 ","celec_pt_p1",200,0,1000);
  TH1D*h_celec_eta_p1   = new TH1D("celec_eta_p1","celec_eta_p1",50,-5,5);
  TH1D*h_celec_phi_p1   = new TH1D("celec_phi_p1","celec_phi_p1",50,-5,5);

  TH1D*h_lelec_pt_p1    = new TH1D("lelec_pt_p1  ","lelec_pt_p1  ",200,0,1000);
  TH1D*h_lelec_eta_p1   = new TH1D("lelec_eta_p1 ","lelec_eta_p1 ",50,-5,5);
  TH1D*h_lelec_phi_p1   = new TH1D("lelec_phi_p1 ","lelec_phi_p1 ",50,-5,5);
  TH1D*h_slelec_pt_p1   = new TH1D("slelec_pt_p1 ","slelec_pt_p1 ",200,0,1000);
  TH1D*h_slelec_eta_p1  = new TH1D("slelec_eta_p1","slelec_eta_p1",50,-5,5);
  TH1D*h_slelec_phi_p1  = new TH1D("slelec_phi_p1","slelec_phi_p1",50,-5,5);

  TH1D*h_felec_pt_p2    = new TH1D("felec_pt_p2 ","felec_pt_p2",200,0,1000);
  TH1D*h_felec_eta_p2   = new TH1D("felec_eta_p2","felec_eta_p2",50,-5,5);
  TH1D*h_felec_phi_p2   = new TH1D("felec_phi_p2","felec_phi_p2",50,-5,5);
  TH1D*h_celec_pt_p2    = new TH1D("celec_pt_p2 ","celec_pt_p2",200,0,1000);
  TH1D*h_celec_eta_p2   = new TH1D("celec_eta_p2","celec_eta_p2",50,-5,5);
  TH1D*h_celec_phi_p2   = new TH1D("celec_phi_p2","celec_phi_p2",50,-5,5);

  TH1D*h_lelec_pt_p2    = new TH1D("lelec_pt_p2  ","lelec_pt_p2  ",200,0,1000);
  TH1D*h_lelec_eta_p2   = new TH1D("lelec_eta_p2 ","lelec_eta_p2 ",50,-5,5);
  TH1D*h_lelec_phi_p2   = new TH1D("lelec_phi_p2 ","lelec_phi_p2 ",50,-5,5);
  TH1D*h_slelec_pt_p2   = new TH1D("slelec_pt_p2 ","slelec_pt_p2 ",200,0,1000);
  TH1D*h_slelec_eta_p2  = new TH1D("slelec_eta_p2","slelec_eta_p2",50,-5,5);
  TH1D*h_slelec_phi_p2  = new TH1D("slelec_phi_p2","slelec_phi_p2",50,-5,5);

  TH1D*h_felec_pt_pc    = new TH1D("felec_pt_pc ","felec_pt_pc",200,0,1000);
  TH1D*h_felec_eta_pc   = new TH1D("felec_eta_pc","felec_eta_pc",50,-5,5);
  TH1D*h_felec_phi_pc   = new TH1D("felec_phi_pc","felec_phi_pc",50,-5,5);
  TH1D*h_celec_pt_pc    = new TH1D("celec_pt_pc ","celec_pt_pc",200,0,1000);
  TH1D*h_celec_eta_pc   = new TH1D("celec_eta_pc","celec_eta_pc",50,-5,5);
  TH1D*h_celec_phi_pc   = new TH1D("celec_phi_pc","celec_phi_pc",50,-5,5);

  TH1D*h_lelec_pt_pc    = new TH1D("lelec_pt_pc  ","lelec_pt_pc  ",200,0,1000);
  TH1D*h_lelec_eta_pc   = new TH1D("lelec_eta_pc ","lelec_eta_pc ",50,-5,5);
  TH1D*h_lelec_phi_pc   = new TH1D("lelec_phi_pc ","lelec_phi_pc ",50,-5,5);
  TH1D*h_slelec_pt_pc   = new TH1D("slelec_pt_pc ","slelec_pt_pc ",200,0,1000);
  TH1D*h_slelec_eta_pc  = new TH1D("slelec_eta_pc","slelec_eta_pc",50,-5,5);
  TH1D*h_slelec_phi_pc  = new TH1D("slelec_phi_pc","slelec_phi_pc",50,-5,5);

  TH1D*h_felec_pt_mc    = new TH1D("felec_pt_mc ","felec_pt_mc",200,0,1000);
  TH1D*h_felec_eta_mc   = new TH1D("felec_eta_mc","felec_eta_mc",50,-5,5);
  TH1D*h_felec_phi_mc   = new TH1D("felec_phi_mc","felec_phi_mc",50,-5,5);
  TH1D*h_celec_pt_mc    = new TH1D("celec_pt_mc ","celec_pt_mc",200,0,1000);
  TH1D*h_celec_eta_mc   = new TH1D("celec_eta_mc","celec_eta_mc",50,-5,5);
  TH1D*h_celec_phi_mc   = new TH1D("celec_phi_mc","celec_phi_mc",50,-5,5);

  TH1D*h_lelec_pt_mc    = new TH1D("lelec_pt_mc  ","lelec_pt_mc  ",200,0,1000);
  TH1D*h_lelec_eta_mc   = new TH1D("lelec_eta_mc ","lelec_eta_mc ",50,-5,5);
  TH1D*h_lelec_phi_mc   = new TH1D("lelec_phi_mc ","lelec_phi_mc ",50,-5,5);
  TH1D*h_slelec_pt_mc   = new TH1D("slelec_pt_mc ","slelec_pt_mc ",200,0,1000);
  TH1D*h_slelec_eta_mc  = new TH1D("slelec_eta_mc","slelec_eta_mc",50,-5,5);
  TH1D*h_slelec_phi_mc  = new TH1D("slelec_phi_mc","slelec_phi_mc",50,-5,5);  

  TH1D*h_felec_pt_m0    = new TH1D("felec_pt_m0 ","felec_pt_m0",200,0,1000);
  TH1D*h_felec_eta_m0   = new TH1D("felec_eta_m0","felec_eta_m0",50,-5,5);
  TH1D*h_felec_phi_m0   = new TH1D("felec_phi_m0","felec_phi_m0",50,-5,5);
  TH1D*h_celec_pt_m0    = new TH1D("celec_pt_m0 ","celec_pt_m0",200,0,1000);
  TH1D*h_celec_eta_m0   = new TH1D("celec_eta_m0","celec_eta_m0",50,-5,5);
  TH1D*h_celec_phi_m0   = new TH1D("celec_phi_m0","celec_phi_m0",50,-5,5);

  TH1D*h_lelec_pt_m0    = new TH1D("lelec_pt_m0  ","lelec_pt_m0  ",200,0,1000);
  TH1D*h_lelec_eta_m0   = new TH1D("lelec_eta_m0 ","lelec_eta_m0 ",50,-5,5);
  TH1D*h_lelec_phi_m0   = new TH1D("lelec_phi_m0 ","lelec_phi_m0 ",50,-5,5);
  TH1D*h_slelec_pt_m0   = new TH1D("slelec_pt_m0 ","slelec_pt_m0 ",200,0,1000);
  TH1D*h_slelec_eta_m0  = new TH1D("slelec_eta_m0","slelec_eta_m0",50,-5,5);
  TH1D*h_slelec_phi_m0  = new TH1D("slelec_phi_m0","slelec_phi_m0",50,-5,5);

  TH1D*h_felec_pt_m1    = new TH1D("felec_pt_m1 ","felec_pt_m1",200,0,1000);
  TH1D*h_felec_eta_m1   = new TH1D("felec_eta_m1","felec_eta_m1",50,-5,5);
  TH1D*h_felec_phi_m1   = new TH1D("felec_phi_m1","felec_phi_m1",50,-5,5);
  TH1D*h_celec_pt_m1    = new TH1D("celec_pt_m1 ","celec_pt_m1",200,0,1000);
  TH1D*h_celec_eta_m1   = new TH1D("celec_eta_m1","celec_eta_m1",50,-5,5);
  TH1D*h_celec_phi_m1   = new TH1D("celec_phi_m1","celec_phi_m1",50,-5,5);

  TH1D*h_lelec_pt_m1    = new TH1D("lelec_pt_m1  ","lelec_pt_m1  ",200,0,1000);
  TH1D*h_lelec_eta_m1   = new TH1D("lelec_eta_m1 ","lelec_eta_m1 ",50,-5,5);
  TH1D*h_lelec_phi_m1   = new TH1D("lelec_phi_m1 ","lelec_phi_m1 ",50,-5,5);
  TH1D*h_slelec_pt_m1   = new TH1D("slelec_pt_m1 ","slelec_pt_m1 ",200,0,1000);
  TH1D*h_slelec_eta_m1  = new TH1D("slelec_eta_m1","slelec_eta_m1",50,-5,5);
  TH1D*h_slelec_phi_m1  = new TH1D("slelec_phi_m1","slelec_phi_m1",50,-5,5);

  TH1D*h_felec_pt_m2    = new TH1D("felec_pt_m2 ","felec_pt_m2",200,0,1000);
  TH1D*h_felec_eta_m2   = new TH1D("felec_eta_m2","felec_eta_m2",50,-5,5);
  TH1D*h_felec_phi_m2   = new TH1D("felec_phi_m2","felec_phi_m2",50,-5,5);
  TH1D*h_celec_pt_m2    = new TH1D("celec_pt_m2 ","celec_pt_m2",200,0,1000);
  TH1D*h_celec_eta_m2   = new TH1D("celec_eta_m2","celec_eta_m2",50,-5,5);
  TH1D*h_celec_phi_m2   = new TH1D("celec_phi_m2","celec_phi_m2",50,-5,5);

  TH1D*h_lelec_pt_m2    = new TH1D("lelec_pt_m2  ","lelec_pt_m2  ",200,0,1000);
  TH1D*h_lelec_eta_m2   = new TH1D("lelec_eta_m2 ","lelec_eta_m2 ",50,-5,5);
  TH1D*h_lelec_phi_m2   = new TH1D("lelec_phi_m2 ","lelec_phi_m2 ",50,-5,5);
  TH1D*h_slelec_pt_m2   = new TH1D("slelec_pt_m2 ","slelec_pt_m2 ",200,0,1000);
  TH1D*h_slelec_eta_m2  = new TH1D("slelec_eta_m2","slelec_eta_m2",50,-5,5);
  TH1D*h_slelec_phi_m2  = new TH1D("slelec_phi_m2","slelec_phi_m2",50,-5,5);



  TH1D* h_GMZElec = new TH1D("GMZElec","M(Z)#rightarrow #mu#mu",40, 71,111.);
  TH1D* h_GdiElecpt = new TH1D("GdiElecpt","diElecpt",200,0,1000);
  TH1D* h_GdiElecrapidity = new TH1D("GdiElecrapidity","diElecrapidity",50,-2.5,2.5);
  TH1D* h_GdiElecphi = new TH1D("GdiElecphi","diElecphi",12,-3.2,3.2);
 TH1D* h_GMZElec_mu = new TH1D("GMZElec_mu","M(Z)#rightarrow #mu#mu",40, 71,111.);
  TH1D* h_GdiElecpt_mu = new TH1D("GdiElecpt_mu","diElecpt",200,0,1000);
  TH1D* h_GdiElecrapidity_mu = new TH1D("GdiElecrapidity_mu","diElecrapidity",50,-2.5,2.5);
  TH1D* h_GdiElecphi_mu = new TH1D("GdiElecphi_mu","diElecphi",12,-3.2,3.2);
 TH1D* h_GMZElec_st3mu = new TH1D("GMZElec_st3mu","M(Z)#rightarrow #mu#mu",40, 71,111.);
  TH1D* h_GdiElecpt_st3mu = new TH1D("GdiElecpt_st3mu","diElecpt",200,0,1000);
  TH1D* h_GdiElecrapidity_st3mu = new TH1D("GdiElecrapidity_st3mu","diElecrapidity",50,-2.5,2.5);
  TH1D* h_GdiElecphi_st3mu = new TH1D("GdiElecphi_st3mu","diElecphi",12,-3.2,3.2);
 TH1D* h_GMZElec_st3e = new TH1D("GMZElec_st3e","M(Z)#rightarrow #mu#mu",40, 71,111.);
  TH1D* h_GdiElecpt_st3e = new TH1D("GdiElecpt_st3e","diElecpt",200,0,1000);
  TH1D* h_GdiElecrapidity_st3e = new TH1D("GdiElecrapidity_st3e","diElecrapidity",50,-2.5,2.5);
  TH1D* h_GdiElecphi_st3e = new TH1D("GdiElecphi_st3e","diElecphi",12,-3.2,3.2);


  TH1D *h_ljet_pt_y[100];
  TH1D *h_ljet_y_pt[100];
  TH1D *h_ljet_pt_y_z0[100];
  TH1D *h_ljet_pt_y_z1[100];
  TH1D *h_ljet_y_pt_z0[100];
  TH1D *h_ljet_y_pt_z1[100];

  TH1D *h_ljet_pt_y_Gen[100];
  TH1D *h_ljet_y_pt_Gen[100];
  TH1D *h_ljet_pt_y_z0_Gen[100];
  TH1D *h_ljet_pt_y_z1_Gen[100];
  TH1D *h_ljet_y_pt_z0_Gen[100];
  TH1D *h_ljet_y_pt_z1_Gen[100];


  TH2* h_ljet_pt_y_Gen_VS_reco[100];
  TH2* h_ljet_y_pt_Gen_VS_reco[100];
  TH2* h_ljet_pt_y_z0_Gen_VS_reco[100];
  TH2* h_ljet_pt_y_z1_Gen_VS_reco[100];
  TH2* h_ljet_y_pt_z0_Gen_VS_reco[100];
  TH2* h_ljet_y_pt_z1_Gen_VS_reco[100];

  TH1D *h_ratio_jetpt30_50_zpt[100];
  TH1D *h_ratio_jetpt50_inf_zpt[100];

  TH1D *h_deltaphi_jet_Z_ptj_20_30_jrap[100];
  TH1D *h_deltaphi_jet_Z_ptj_30_50_jrap[100];
  TH1D *h_deltaphi_jet_Z_ptj_50_100_jrap[100];
  TH1D *h_deltaphi_jet_Z_ptj_100_1000_jrap[100];

  TH1D *h_nomatch_gen_deltaphi_jet_Z_ptj_20_30_jrap[100];
  TH1D *h_nomatch_gen_deltaphi_jet_Z_ptj_30_50_jrap[100];
  TH1D *h_nomatch_gen_deltaphi_jet_Z_ptj_50_100_jrap[100];
  TH1D *h_nomatch_gen_deltaphi_jet_Z_ptj_100_1000_jrap[100];

  TH1D *h_match_gen_deltaphi_jet_Z_ptj_20_30_jrap[100];
  TH1D *h_match_gen_deltaphi_jet_Z_ptj_30_50_jrap[100];
  TH1D *h_match_gen_deltaphi_jet_Z_ptj_50_100_jrap[100];
  TH1D *h_match_gen_deltaphi_jet_Z_ptj_100_1000_jrap[100];

  TH1D *h_deltaphi_jet_Z_ptj_20_24_jrap[100];
  TH1D *h_deltaphi_jet_Z_ptj_24_30_jrap[100];
  TH1D *h_deltaphi_jet_Z_ptj_30_40_jrap[100];
  TH1D *h_deltaphi_jet_Z_ptj_40_50_jrap[100];
  TH1D *h_deltaphi_jet_Z_ptj_50_68_jrap[100];

  TH1D *h_nomatch_gen_deltaphi_jet_Z_ptj_20_24_jrap[100];
  TH1D *h_nomatch_gen_deltaphi_jet_Z_ptj_24_30_jrap[100];
  TH1D *h_nomatch_gen_deltaphi_jet_Z_ptj_30_40_jrap[100];
  TH1D *h_nomatch_gen_deltaphi_jet_Z_ptj_40_50_jrap[100];
  TH1D *h_nomatch_gen_deltaphi_jet_Z_ptj_50_68_jrap[100];

  TH1D *h_match_gen_deltaphi_jet_Z_ptj_20_24_jrap[100];
  TH1D *h_match_gen_deltaphi_jet_Z_ptj_24_30_jrap[100];
  TH1D *h_match_gen_deltaphi_jet_Z_ptj_30_40_jrap[100];
  TH1D *h_match_gen_deltaphi_jet_Z_ptj_40_50_jrap[100];
  TH1D *h_match_gen_deltaphi_jet_Z_ptj_50_68_jrap[100];


  TH1D *h_inc_deltaphi_jet_Z_ptj_20_30_jrap[100];
  TH1D *h_inc_deltaphi_jet_Z_ptj_30_50_jrap[100];
  TH1D *h_inc_deltaphi_jet_Z_ptj_50_100_jrap[100];
  TH1D *h_inc_deltaphi_jet_Z_ptj_100_1000_jrap[100];

  TH1D *h_inc_nomatch_gen_deltaphi_jet_Z_ptj_20_30_jrap[100];
  TH1D *h_inc_nomatch_gen_deltaphi_jet_Z_ptj_30_50_jrap[100];
  TH1D *h_inc_nomatch_gen_deltaphi_jet_Z_ptj_50_100_jrap[100];
  TH1D *h_inc_nomatch_gen_deltaphi_jet_Z_ptj_100_1000_jrap[100];

  TH1D *h_inc_match_gen_deltaphi_jet_Z_ptj_20_30_jrap[100];
  TH1D *h_inc_match_gen_deltaphi_jet_Z_ptj_30_50_jrap[100];
  TH1D *h_inc_match_gen_deltaphi_jet_Z_ptj_50_100_jrap[100];
  TH1D *h_inc_match_gen_deltaphi_jet_Z_ptj_100_1000_jrap[100];


  TH1D *h_inc_1j_deltaphi_jettot_Z_ptj_20_30_jrap[100];
  TH1D *h_inc_1j_deltaphi_jettot_Z_ptj_30_50_jrap[100];
  TH1D *h_inc_1j_deltaphi_jettot_Z_ptj_50_100_jrap[100];
  TH1D *h_inc_1j_deltaphi_jettot_Z_ptj_100_1000_jrap[100];

  TH1D *h_inc_1j_nomatch_gen_deltaphi_jettot_Z_ptj_20_30_jrap[100];
  TH1D *h_inc_1j_nomatch_gen_deltaphi_jettot_Z_ptj_30_50_jrap[100];
  TH1D *h_inc_1j_nomatch_gen_deltaphi_jettot_Z_ptj_50_100_jrap[100];
  TH1D *h_inc_1j_nomatch_gen_deltaphi_jettot_Z_ptj_100_1000_jrap[100];

  TH1D *h_inc_1j_match_gen_deltaphi_jettot_Z_ptj_20_30_jrap[100];
  TH1D *h_inc_1j_match_gen_deltaphi_jettot_Z_ptj_30_50_jrap[100];
  TH1D *h_inc_1j_match_gen_deltaphi_jettot_Z_ptj_50_100_jrap[100];
  TH1D *h_inc_1j_match_gen_deltaphi_jettot_Z_ptj_100_1000_jrap[100];


  TH1D *h_ratio_Zpt0_30_zpt_dphi_28[100];
  TH1D *h_ratio_Zpt30_50_zpt_dphi_28[100];
  TH1D *h_ratio_Zpt50_inf_zpt_dphi_28[100];
  TH1D *h_ratio_Zpt0_30_zpt_dphi_28_up[100];
  TH1D *h_ratio_Zpt30_50_zpt_dphi_28_up[100];
  TH1D *h_ratio_Zpt50_inf_zpt_dphi_28_up[100];
  TH1D *h_ratio_Zpt0_30_zpt_dphi_28_dn[100];
  TH1D *h_ratio_Zpt30_50_zpt_dphi_28_dn[100];
  TH1D *h_ratio_Zpt50_inf_zpt_dphi_28_dn[100];

  TH1D *h_ratio_jetpt30_50_zpt_dphi_28[100];
  TH1D *h_ratio_jetpt50_inf_zpt_dphi_28[100];
  TH1D *h_ratio_jetpt30_50_zpt_dphi_28_up[100];
  TH1D *h_ratio_jetpt50_inf_zpt_dphi_28_up[100];
  TH1D *h_ratio_jetpt30_50_zpt_dphi_28_dn[100];
  TH1D *h_ratio_jetpt50_inf_zpt_dphi_28_dn[100];
  TH1D *h_inc_ljet_pt_y[100];
  TH1D *h_inc_ljet_Z_y[100];
  TH1D *h_inc_ljet_pt_y_fb[100];
  TH1D *h_inc_ljet_pt_y_fb_nomatch[100];
  TH1D *h_inc_ljet_pt_y_fb_match[100];
  TH1D *h_inc_ljet_pt_y_fb_dphi_lt_1pt5[100];
  TH1D *h_inc_ljet_pt_y_fb_dphi_lt_1pt0[100];
  TH1D *h_inc_ljet_pt_y_fixbin[100];
 
double range_for_incl[]= {20,24,30,40,50,68,85,104,126,148,175,200,230,400};
  TH1D *h_inc_ljet_pt_y_fb_alljetrap_z0_pos = new TH1D ("inc_ljet_pt_y_fb_alljetrap_pos_z0","inc_ljet_pt_y_fb_alljetrap_pos_z0", 13,  range_for_incl);
  TH1D *h_inc_ljet_pt_y_fb_alljetrap_z1_pos = new TH1D ("inc_ljet_pt_y_fb_alljetrap_pos_z1","inc_ljet_pt_y_fb_alljetrap_pos_z1", 13,  range_for_incl);
  TH1D *h_inc_ljet_pt_y_fb_alljetrap_z0_neg = new TH1D ("inc_ljet_pt_y_fb_alljetrap_neg_z0","inc_ljet_pt_y_fb_alljetrap_neg_z0", 13,  range_for_incl);
  TH1D *h_inc_ljet_pt_y_fb_alljetrap_z1_neg = new TH1D ("inc_ljet_pt_y_fb_alljetrap_neg_z1","inc_ljet_pt_y_fb_alljetrap_neg_z1", 13,  range_for_incl);

  TH1D *h_inc_ljet_pt_y_fb_alljetrap_z0_pos_Gen = new TH1D ("inc_ljet_pt_y_fb_alljetrap_pos_Gen_z0","inc_ljet_pt_y_fb_alljetrap_pos_Gen_z0", 13,  range_for_incl);
  TH1D *h_inc_ljet_pt_y_fb_alljetrap_z1_pos_Gen = new TH1D ("inc_ljet_pt_y_fb_alljetrap_pos_Gen_z1","inc_ljet_pt_y_fb_alljetrap_pos_Gen_z1", 13,  range_for_incl);
  TH1D *h_inc_ljet_pt_y_fb_alljetrap_z0_neg_Gen = new TH1D ("inc_ljet_pt_y_fb_alljetrap_neg_Gen_z0","inc_ljet_pt_y_fb_alljetrap_neg_Gen_z0", 13,  range_for_incl);
  TH1D *h_inc_ljet_pt_y_fb_alljetrap_z1_neg_Gen = new TH1D ("inc_ljet_pt_y_fb_alljetrap_neg_Gen_z1","inc_ljet_pt_y_fb_alljetrap_neg_Gen_z1", 13,  range_for_incl);

  TH2D *h_inc_ljet_pt_y_fb_alljetrap_z0_pos_Gen_VS_reco = new TH2D ("inc_ljet_pt_y_fb_alljetrap_pos_Gen_VS_reco_z0","inc_ljet_pt_y_fb_alljetrap_pos_Gen_VS_reco_z0", 13,  range_for_incl,13,  range_for_incl);
  TH2D *h_inc_ljet_pt_y_fb_alljetrap_z1_pos_Gen_VS_reco = new TH2D ("inc_ljet_pt_y_fb_alljetrap_pos_Gen_VS_reco_z1","inc_ljet_pt_y_fb_alljetrap_pos_Gen_VS_reco_z1", 13,  range_for_incl,13,  range_for_incl);
  TH2D *h_inc_ljet_pt_y_fb_alljetrap_z0_neg_Gen_VS_reco = new TH2D ("inc_ljet_pt_y_fb_alljetrap_neg_Gen_VS_reco_z0","inc_ljet_pt_y_fb_alljetrap_neg_Gen_VS_reco_z0", 13,  range_for_incl,13,  range_for_incl);
  TH2D *h_inc_ljet_pt_y_fb_alljetrap_z1_neg_Gen_VS_reco = new TH2D ("inc_ljet_pt_y_fb_alljetrap_neg_Gen_VS_reco_z1","inc_ljet_pt_y_fb_alljetrap_neg_Gen_VS_reco_z1", 13,  range_for_incl,13,  range_for_incl);


  TH1D *h_inc_ljet_pt_y_fb_z0[100];
  TH1D *h_inc_ljet_pt_y_fb_z1[100];
  TH1D *h_inc_ljet_pt_y_fb_z0_pos[100];
  TH1D *h_inc_ljet_pt_y_fb_z1_pos[100];
  TH1D *h_inc_ljet_pt_y_fb_z0_neg[100];
  TH1D *h_inc_ljet_pt_y_fb_z1_neg[100];
  TH1D *h_inc_ljet_pt_sljetpt[100];
  TH1D *h_inc_ljet_y_pt[100];
  TH1D *h_inc_ljet_pt_y_z0[100];
  TH1D *h_inc_ljet_pt_y_z1[100];

  TH1D *h_inc_ljet_pt_sljetpt_z0[100];
  TH1D *h_inc_ljet_pt_sljetpt_z1[100];
  TH1D *h_inc_ljet_y_pt_z0[100];
  TH1D *h_inc_ljet_y_pt_z1[100];

  TH1D *h_inc_ljet_pt_y_Gen[100];
  TH1D *h_inc_ljet_Z_y_Gen[100];
  TH1D *h_inc_ljet_pt_sljetpt_Gen[100];
  TH1D *h_inc_ljet_y_pt_Gen[100];
  TH1D *h_inc_ljet_pt_y_z0_Gen[100];
  TH1D *h_inc_ljet_pt_y_z1_Gen[100];
  TH1D *h_inc_ljet_pt_sljetpt_z0_Gen[100];
  TH1D *h_inc_ljet_pt_sljetpt_z1_Gen[100];
  TH1D *h_inc_ljet_y_pt_z0_Gen[100];
  TH1D *h_inc_ljet_y_pt_z1_Gen[100];
  TH1D *h_inc_ljet_pt_y_fb_Gen[100];
  TH1D *h_inc_ljet_pt_y_fb_z0_Gen[100];
  TH1D *h_inc_ljet_pt_y_fb_z1_Gen[100];
  TH1D *h_inc_ljet_pt_y_fb_z0_pos_Gen[100];
  TH1D *h_inc_ljet_pt_y_fb_z1_pos_Gen[100];
  TH1D *h_inc_ljet_pt_y_fb_z0_neg_Gen[100];
  TH1D *h_inc_ljet_pt_y_fb_z1_neg_Gen[100];

  TH2* h_inc_ljet_pt_y_Gen_VS_reco[100];
  TH2 *h_inc_ljet_Z_y_Gen_VS_reco[100];
  TH2 *h_inc_ljet_pt_sljetpt_Gen_VS_reco[100];
  TH2* h_inc_ljet_y_pt_Gen_VS_reco[100];
  TH2* h_inc_ljet_pt_y_z0_Gen_VS_reco[100];
  TH2* h_inc_ljet_pt_y_z1_Gen_VS_reco[100];
  TH2 *h_inc_ljet_pt_sljetpt_z0_Gen_VS_reco[100];
  TH2 *h_inc_ljet_pt_sljetpt_z1_Gen_VS_reco[100];
  TH2* h_inc_ljet_y_pt_z0_Gen_VS_reco[100];
  TH2* h_inc_ljet_y_pt_z1_Gen_VS_reco[100];
  TH2 *h_inc_ljet_pt_y_fb_Gen_VS_reco[100];
  TH2 *h_inc_ljet_pt_y_fb_z0_Gen_VS_reco[100];
  TH2 *h_inc_ljet_pt_y_fb_z1_Gen_VS_reco[100];
  TH2 *h_inc_ljet_pt_y_fb_z0_pos_Gen_VS_reco[100];
  TH2 *h_inc_ljet_pt_y_fb_z1_pos_Gen_VS_reco[100];
  TH2 *h_inc_ljet_pt_y_fb_z0_neg_Gen_VS_reco[100];
  TH2 *h_inc_ljet_pt_y_fb_z1_neg_Gen_VS_reco[100];

  TH1D *h_inc_sljet_pt_y[100];
  TH1D *h_inc_sljet_y_pt[100];
  TH1D *h_inc_sljet_pt_y_z0[100];
  TH1D *h_inc_sljet_pt_y_z1[100];
  TH1D *h_inc_sljet_y_pt_z0[100];
  TH1D *h_inc_sljet_y_pt_z1[100];

  TH1D *h_inc_sljet_pt_y_Gen[100];
  TH1D *h_inc_sljet_y_pt_Gen[100];
  TH1D *h_inc_sljet_pt_y_z0_Gen[100];
  TH1D *h_inc_sljet_pt_y_z1_Gen[100];
  TH1D *h_inc_sljet_y_pt_z0_Gen[100];
  TH1D *h_inc_sljet_y_pt_z1_Gen[100];


  TH2* h_inc_sljet_pt_y_Gen_VS_reco[100];
  TH2* h_inc_sljet_y_pt_Gen_VS_reco[100];
  TH2* h_inc_sljet_pt_y_z0_Gen_VS_reco[100];
  TH2* h_inc_sljet_pt_y_z1_Gen_VS_reco[100];
  TH2* h_inc_sljet_y_pt_z0_Gen_VS_reco[100];
  TH2* h_inc_sljet_y_pt_z1_Gen_VS_reco[100];

  TH1D *h_z_pt_y[100][100];
  TH1D *h_z_y_pt[100][100];
  TH1D *h_z_pt_y_Gen[100][100];
  TH1D *h_z_y_pt_Gen[100][100];
  TH2* h_z_pt_y_Gen_VS_reco[100][100];
  TH2* h_z_y_pt_Gen_VS_reco[100][100];

  char name[100];


  double j_pT_range[9]={30.,50.,70.,100.,130.,170.,220.,350.,1000.};
  double j_pT_range2[21] = { 30, 39, 49, 60, 72, 85, 100, 117, 136, 157, 187, 220, 258, 300, 350, 400, 450, 500, 590, 700, 1000};
/*double j_pT_range3[7][23] = {
			     {30,35,40,45,50,55,60,70,80,90,100,110,120,130,145,160,180,200,220,240,270,310,450},  
  			     {30,35,40,45,50,55,60,70,80,90,100,110,120,130,145,160,180,200,220,240,270,310,420},  
                             {30,35,40,45,50,55,60,70,80,90,100,110,120,130,145,160,180,200,220,240,270,390},
                             {30,35,40,45,50,55,60,70,80,90,100,110,120,130,145,160,180,200,220,360},
                             {30,35,40,45,50,55,60,70,80,90,100,110,120,130,145,160,180,275},
                             {30,35,40,45,50,55,60,70,80,90,100,110,120,130,145,220},
                             {30,35,40,45,50,55,60,70,80,145}};
  double range_2d[7]={22,22,21,19,17,15,9};*/
  
/*double j_pT_range3[7][16] = {
{30,40,50,60,75,90,105,120,135,155,180,210,240,270,310,450},  
{30,40,50,60,75,90,105,120,135,155,180,210,240,270,310,420}, 
{30,40,50,60,75,90,105,120,135,155,180,210,240,270,390}, 
{30,40,50,60,75,90,105,120,135,155,180,210,360}, 
{30,40,50,60,75,90,105,120,135,155,180,280}, 
{30,40,50,60,75,90,105,120,135,220}, 
{30,40,50,60,75,90,145}};  
  double range_2d[7]={15,15,14,12,11,9,6};*/

double j_pT_range3[7][18] = {
{20,24,30,40,50,68,85,104,126,148,175,200,230,270,320,550},  
{20,24,30,40,50,68,85,104,126,148,175,200,230,270,320,500}, 
{20,24,30,40,50,68,85,104,126,148,175,200,230,270,450}, 
{20,24,30,40,50,68,85,104,126,148,175,200,230,400}, 
{20,24,30,40,50,68,85,104,126,148,175,200,320}, 
{20,24,30,40,50,68,85,104,126,148,240}, 
{20,24,30,40,50,68,85,104,180}};  

  double range_2d[7]={15,15,14,13,12,10,8};

  double j_Y_range[8]={0.,0.5,1.0,1.5,2.0,2.5,3.2,4.7};

  double Z_pT_range[6]={20.,30.,50.,100.,200.,1000.};
  double Z_Y_range[6]={0.,0.5,1.0,1.5,2.0,2.5};


double j_pT_range33[3][18] = {

{20,24,30,40,50,68,85,104,126,148,175,200,230,270,320,500}, 
{20,24,30,40,50,68,85,104,126,148,175,200,230,400}, 
{20,24,30,40,50,68,85,104,126,148,240}};  

  double range_3d[7]={15,13,10};

 //   Z_y_range[5]={0,0.9,1.5,2.1,2.5};//
 //   threed_j_y_range[4]={0,1.5,2.4,4.7};//
 // pt_range[6]={20.,30.,50.,100.,200.,1000.};//
for(int i = 0; i<5;i++){//Z rapidity range
  for(int j = 0; j<3;j++){//jet rapidity range
//Z PT vrt jet 
    sprintf(name,"Z_pt_y_%i_%i",i,j);
    h_z_pt_y[i][j] = new TH1D (name,name, 5, Z_pT_range);
    sprintf(name,"Z_y_pt_%i_%i",i,j);
    h_z_y_pt[i][j] = new TH1D (name,name,5, Z_Y_range);
  
    sprintf(name,"Z_pt_y_Gen_%i_%i",i,j);
    h_z_pt_y_Gen[i][j] = new TH1D (name,name, 5, Z_pT_range);
    sprintf(name,"Z_y_pt_Gen_%i_%i",i,j);
    h_z_y_pt_Gen[i][j] = new TH1D (name,name, 5, Z_Y_range);

    sprintf(name,"Z_pt_y_Gen_VS_reco_%i_%i",i,j);
    h_z_pt_y_Gen_VS_reco[i][j] = new TH2D (name,name, 5, Z_pT_range, 5, Z_pT_range);
    sprintf(name,"Z_y_pt_Gen_VS_reco_%i_%i",i,j);
    h_z_y_pt_Gen_VS_reco[i][j] = new TH2D (name,name, 5, Z_Y_range, 5, Z_Y_range);
     
  }
}
for(int i = 0; i<9;i++){//j rapidity range


sprintf(name,"tackman_tau_cm_sum_%i",i);
tackman_tau_cm_sum[i] = new TH1D (name,name, 100,0, 100);
sprintf(name,"tackman_tau_sum_%i",i);
tackman_tau_sum[i] = new TH1D (name,name, 100,0, 100);
sprintf(name,"tackman_tauc_cm_sum_%i",i);
tackman_tauc_cm_sum[i] = new TH1D (name,name, 100,0, 100);
sprintf(name,"tackman_tauc_sum_%i",i);
tackman_tauc_sum[i] = new TH1D (name,name, 100,0, 100);

sprintf(name,"tackman_tau_cm_max_%i",i);
tackman_tau_cm_max[i] = new TH1D (name,name, 100,0, 100);
sprintf(name,"tackman_tau_max_%i",i);
tackman_tau_max[i] = new TH1D (name,name, 100,0, 100);
sprintf(name,"tackman_tauc_cm_max_%i",i);
tackman_tauc_cm_max[i] = new TH1D (name,name, 100,0, 100);
sprintf(name,"tackman_tauc_max_%i",i);
tackman_tauc_max[i] = new TH1D (name,name, 100,0, 100);

  sprintf(name,"ljet_pt_y_%i",i);
  h_ljet_pt_y[i] = new TH1D (name,name, 8, j_pT_range);
  sprintf(name,"ljet_y_pt_%i",i);
  h_ljet_y_pt[i] = new TH1D (name,name, 7,j_Y_range);
  sprintf(name,"ljet_pt_y_z0_%i",i);
  h_ljet_pt_y_z0[i] = new TH1D (name,name, 8, j_pT_range);
  sprintf(name,"ljet_pt_y_z1_%i",i);
  h_ljet_pt_y_z1[i] = new TH1D (name,name, 8, j_pT_range);
  sprintf(name,"ljet_y_pt_z0_%i",i);
  h_ljet_y_pt_z0[i] = new TH1D (name,name, 7,j_Y_range);
  sprintf(name,"ljet_y_pt_z1_%i",i);
  h_ljet_y_pt_z1[i] = new TH1D (name,name, 7,j_Y_range);


  sprintf(name,"ljet_pt_y_Gen_%i",i);
  h_ljet_pt_y_Gen[i] = new TH1D(name,name,8,j_pT_range);
  sprintf(name,"ljet_y_pt_Gen_%i",i);
  h_ljet_y_pt_Gen[i] = new TH1D(name,name,7,j_Y_range);
  sprintf(name,"ljet_pt_y_Gen_z0_%i",i);
  h_ljet_pt_y_z0_Gen[i] = new TH1D(name,name,8,j_pT_range);
  sprintf(name,"ljet_pt_y_Gen_z1_%i",i);
  h_ljet_pt_y_z1_Gen[i] = new TH1D(name,name,8,j_pT_range);
  sprintf(name,"ljet_y_pt_Gen_z0_%i",i);
  h_ljet_y_pt_z0_Gen[i] = new TH1D(name,name,7,j_Y_range);
  sprintf(name,"ljet_y_pt_Gen_z1_%i",i);
  h_ljet_y_pt_z1_Gen[i] = new TH1D(name,name,7,j_Y_range);

  sprintf(name,"ljet_pt_y_Gen_VS_reco_%i",i);
  h_ljet_pt_y_Gen_VS_reco[i] = new TH2D(name,name,8,j_pT_range,8,j_pT_range);
  sprintf(name,"ljet_y_pt_Gen_VS_reco_%i",i);
  h_ljet_y_pt_Gen_VS_reco[i] = new TH2D(name,name,7,j_Y_range,7,j_Y_range);
  sprintf(name,"ljet_pt_y_Gen_VS_reco_z0_%i",i);
  h_ljet_pt_y_z0_Gen_VS_reco[i] = new TH2D(name,name,8,j_pT_range,8,j_pT_range);
  sprintf(name,"ljet_pt_y_Gen_VS_reco_z1_%i",i);
  h_ljet_pt_y_z1_Gen_VS_reco[i] = new TH2D(name,name,8,j_pT_range,8,j_pT_range);
  sprintf(name,"ljet_y_pt_Gen_VS_reco_z0_%i",i);
  h_ljet_y_pt_z0_Gen_VS_reco[i] = new TH2D(name,name,7,j_Y_range,7,j_Y_range);
  sprintf(name,"ljet_y_pt_Gen_VS_reco_z1_%i",i);
  h_ljet_y_pt_z1_Gen_VS_reco[i] = new TH2D(name,name,7,j_Y_range,7,j_Y_range);

  sprintf(name,"inc_ljet_pt_y_%i",i);
  h_inc_ljet_pt_y[i] = new TH1D (name,name, 8, j_pT_range);

  sprintf(name,"inc_ljet_Z_y_%i",i);
  h_inc_ljet_Z_y[i] = new TH1D (name,name, 20,-2.5,2.5);



  sprintf(name,"ratio_jetpt30_50_zpt_%i",i);
  h_ratio_jetpt30_50_zpt[i]= new TH1D (name,name, 30,0.,3.);
  sprintf(name,"ratio_jetpt50_inf_zpt_%i",i);
  h_ratio_jetpt50_inf_zpt[i]= new TH1D (name,name, 30,0.,3.);

  sprintf(name,"deltaphi_jet_Z_ptj_20_30_jrap_%i",i);
  h_deltaphi_jet_Z_ptj_20_30_jrap[i]= new TH1D (name,name, 40,0.,3.14);
  sprintf(name,"deltaphi_jet_Z_ptj_30_50_jrap_%i",i);
  h_deltaphi_jet_Z_ptj_30_50_jrap[i]= new TH1D (name,name, 40,0.,3.14);
  sprintf(name,"deltaphi_jet_Z_ptj_50_100_jrap_%i",i);
  h_deltaphi_jet_Z_ptj_50_100_jrap[i]= new TH1D (name,name, 40,0.,3.14);
  sprintf(name,"deltaphi_jet_Z_ptj_100_1000_jrap_%i",i);
  h_deltaphi_jet_Z_ptj_100_1000_jrap[i]= new TH1D (name,name, 40,0.,3.14);

  sprintf(name,"match_gen_deltaphi_jet_Z_ptj_20_30_jrap_%i",i);
  h_match_gen_deltaphi_jet_Z_ptj_20_30_jrap[i]= new TH1D (name,name, 40,0.,3.14);
  sprintf(name,"match_gen_deltaphi_jet_Z_ptj_30_50_jrap_%i",i);
  h_match_gen_deltaphi_jet_Z_ptj_30_50_jrap[i]= new TH1D (name,name, 40,0.,3.14);
  sprintf(name,"match_gen_deltaphi_jet_Z_ptj_50_100_jrap_%i",i);
  h_match_gen_deltaphi_jet_Z_ptj_50_100_jrap[i]= new TH1D (name,name, 40,0.,3.14);
  sprintf(name,"match_gen_deltaphi_jet_Z_ptj_100_1000_jrap_%i",i);
  h_match_gen_deltaphi_jet_Z_ptj_100_1000_jrap[i]= new TH1D (name,name, 40,0.,3.14);

  sprintf(name,"nomatch_gen_deltaphi_jet_Z_ptj_20_30_jrap_%i",i);
  h_nomatch_gen_deltaphi_jet_Z_ptj_20_30_jrap[i]= new TH1D (name,name, 40,0.,3.14);
  sprintf(name,"nomatch_gen_deltaphi_jet_Z_ptj_30_50_jrap_%i",i);
  h_nomatch_gen_deltaphi_jet_Z_ptj_30_50_jrap[i]= new TH1D (name,name, 40,0.,3.14);
  sprintf(name,"nomatch_gen_deltaphi_jet_Z_ptj_50_100_jrap_%i",i);
  h_nomatch_gen_deltaphi_jet_Z_ptj_50_100_jrap[i]= new TH1D (name,name, 40,0.,3.14);
  sprintf(name,"nomatch_gen_deltaphi_jet_Z_ptj_100_1000_jrap_%i",i);
  h_nomatch_gen_deltaphi_jet_Z_ptj_100_1000_jrap[i]= new TH1D (name,name, 40,0.,3.14);

  sprintf(name,"deltaphi_jet_Z_ptj_20_24_jrap_%i",i);
  h_deltaphi_jet_Z_ptj_20_24_jrap[i]= new TH1D (name,name, 40,0.,3.14);
  sprintf(name,"deltaphi_jet_Z_ptj_24_30_jrap_%i",i);
  h_deltaphi_jet_Z_ptj_24_30_jrap[i]= new TH1D (name,name, 40,0.,3.14);
  sprintf(name,"deltaphi_jet_Z_ptj_30_40_jrap_%i",i);
  h_deltaphi_jet_Z_ptj_30_40_jrap[i]= new TH1D (name,name, 40,0.,3.14);
  sprintf(name,"deltaphi_jet_Z_ptj_40_50_jrap_%i",i);
  h_deltaphi_jet_Z_ptj_40_50_jrap[i]= new TH1D (name,name, 40,0.,3.14);
  sprintf(name,"deltaphi_jet_Z_ptj_50_68_jrap_%i",i);
  h_deltaphi_jet_Z_ptj_50_68_jrap[i]= new TH1D (name,name, 40,0.,3.14);

  sprintf(name,"match_gen_deltaphi_jet_Z_ptj_20_24_jrap_%i",i);
  h_match_gen_deltaphi_jet_Z_ptj_20_24_jrap[i]= new TH1D (name,name, 40,0.,3.14);
  sprintf(name,"match_gen_deltaphi_jet_Z_ptj_24_30_jrap_%i",i);
  h_match_gen_deltaphi_jet_Z_ptj_24_30_jrap[i]= new TH1D (name,name, 40,0.,3.14);
  sprintf(name,"match_gen_deltaphi_jet_Z_ptj_30_40_jrap_%i",i);
  h_match_gen_deltaphi_jet_Z_ptj_30_40_jrap[i]= new TH1D (name,name, 40,0.,3.14);
  sprintf(name,"match_gen_deltaphi_jet_Z_ptj_40_50_jrap_%i",i);
  h_match_gen_deltaphi_jet_Z_ptj_40_50_jrap[i]= new TH1D (name,name, 40,0.,3.14);
  sprintf(name,"match_gen_deltaphi_jet_Z_ptj_50_68_jrap_%i",i);
  h_match_gen_deltaphi_jet_Z_ptj_50_68_jrap[i]= new TH1D (name,name, 40,0.,3.14);

  sprintf(name,"nomatch_gen_deltaphi_jet_Z_ptj_20_24_jrap_%i",i);
  h_nomatch_gen_deltaphi_jet_Z_ptj_20_24_jrap[i]= new TH1D (name,name, 40,0.,3.14);
  sprintf(name,"nomatch_gen_deltaphi_jet_Z_ptj_24_30_jrap_%i",i);
  h_nomatch_gen_deltaphi_jet_Z_ptj_24_30_jrap[i]= new TH1D (name,name, 40,0.,3.14);
  sprintf(name,"nomatch_gen_deltaphi_jet_Z_ptj_30_40_jrap_%i",i);
  h_nomatch_gen_deltaphi_jet_Z_ptj_30_40_jrap[i]= new TH1D (name,name, 40,0.,3.14);
  sprintf(name,"nomatch_gen_deltaphi_jet_Z_ptj_40_50_jrap_%i",i);
  h_nomatch_gen_deltaphi_jet_Z_ptj_40_50_jrap[i]= new TH1D (name,name, 40,0.,3.14);
  sprintf(name,"nomatch_gen_deltaphi_jet_Z_ptj_50_68_jrap_%i",i);
  h_nomatch_gen_deltaphi_jet_Z_ptj_50_68_jrap[i]= new TH1D (name,name, 40,0.,3.14);




  sprintf(name,"inc_deltaphi_jet_Z_ptj_20_30_jrap_%i",i);
  h_inc_deltaphi_jet_Z_ptj_20_30_jrap[i]= new TH1D (name,name, 40,0.,3.14);
  sprintf(name,"inc_deltaphi_jet_Z_ptj_30_50_jrap_%i",i);
  h_inc_deltaphi_jet_Z_ptj_30_50_jrap[i]= new TH1D (name,name, 40,0.,3.14);
  sprintf(name,"inc_deltaphi_jet_Z_ptj_50_100_jrap_%i",i);
  h_inc_deltaphi_jet_Z_ptj_50_100_jrap[i]= new TH1D (name,name, 40,0.,3.14);
  sprintf(name,"inc_deltaphi_jet_Z_ptj_100_1000_jrap_%i",i);
  h_inc_deltaphi_jet_Z_ptj_100_1000_jrap[i]= new TH1D (name,name, 40,0.,3.14);

  sprintf(name,"inc_match_gen_deltaphi_jet_Z_ptj_20_30_jrap_%i",i);
  h_inc_match_gen_deltaphi_jet_Z_ptj_20_30_jrap[i]= new TH1D (name,name, 40,0.,3.14);
  sprintf(name,"inc_match_gen_deltaphi_jet_Z_ptj_30_50_jrap_%i",i);
  h_inc_match_gen_deltaphi_jet_Z_ptj_30_50_jrap[i]= new TH1D (name,name, 40,0.,3.14);
  sprintf(name,"inc_match_gen_deltaphi_jet_Z_ptj_50_100_jrap_%i",i);
  h_inc_match_gen_deltaphi_jet_Z_ptj_50_100_jrap[i]= new TH1D (name,name, 40,0.,3.14);
  sprintf(name,"inc_match_gen_deltaphi_jet_Z_ptj_100_1000_jrap_%i",i);
  h_inc_match_gen_deltaphi_jet_Z_ptj_100_1000_jrap[i]= new TH1D (name,name, 40,0.,3.14);

  sprintf(name,"inc_nomatch_gen_deltaphi_jet_Z_ptj_20_30_jrap_%i",i);
  h_inc_nomatch_gen_deltaphi_jet_Z_ptj_20_30_jrap[i]= new TH1D (name,name, 40,0.,3.14);
  sprintf(name,"inc_nomatch_gen_deltaphi_jet_Z_ptj_30_50_jrap_%i",i);
  h_inc_nomatch_gen_deltaphi_jet_Z_ptj_30_50_jrap[i]= new TH1D (name,name, 40,0.,3.14);
  sprintf(name,"inc_nomatch_gen_deltaphi_jet_Z_ptj_50_100_jrap_%i",i);
  h_inc_nomatch_gen_deltaphi_jet_Z_ptj_50_100_jrap[i]= new TH1D (name,name, 40,0.,3.14);
  sprintf(name,"inc_nomatch_gen_deltaphi_jet_Z_ptj_100_1000_jrap_%i",i);
  h_inc_nomatch_gen_deltaphi_jet_Z_ptj_100_1000_jrap[i]= new TH1D (name,name, 40,0.,3.14);

  sprintf(name,"inc_1j_deltaphi_jettot_Z_ptj_20_30_jrap_%i",i);
  h_inc_1j_deltaphi_jettot_Z_ptj_20_30_jrap[i]= new TH1D (name,name, 40,0.,3.14);
  sprintf(name,"inc_1j_deltaphi_jettot_Z_ptj_30_50_jrap_%i",i);
  h_inc_1j_deltaphi_jettot_Z_ptj_30_50_jrap[i]= new TH1D (name,name, 40,0.,3.14);
  sprintf(name,"inc_1j_deltaphi_jettot_Z_ptj_50_100_jrap_%i",i);
  h_inc_1j_deltaphi_jettot_Z_ptj_50_100_jrap[i]= new TH1D (name,name, 40,0.,3.14);
  sprintf(name,"inc_1j_deltaphi_jettot_Z_ptj_100_1000_jrap_%i",i);
  h_inc_1j_deltaphi_jettot_Z_ptj_100_1000_jrap[i]= new TH1D (name,name, 40,0.,3.14);

  sprintf(name,"inc_1j_match_gen_deltaphi_jettot_Z_ptj_20_30_jrap_%i",i);
  h_inc_1j_match_gen_deltaphi_jettot_Z_ptj_20_30_jrap[i]= new TH1D (name,name, 40,0.,3.14);
  sprintf(name,"inc_1j_match_gen_deltaphi_jettot_Z_ptj_30_50_jrap_%i",i);
  h_inc_1j_match_gen_deltaphi_jettot_Z_ptj_30_50_jrap[i]= new TH1D (name,name, 40,0.,3.14);
  sprintf(name,"inc_1j_match_gen_deltaphi_jettot_Z_ptj_50_100_jrap_%i",i);
  h_inc_1j_match_gen_deltaphi_jettot_Z_ptj_50_100_jrap[i]= new TH1D (name,name, 40,0.,3.14);
  sprintf(name,"inc_1j_match_gen_deltaphi_jettot_Z_ptj_100_1000_jrap_%i",i);
  h_inc_1j_match_gen_deltaphi_jettot_Z_ptj_100_1000_jrap[i]= new TH1D (name,name, 40,0.,3.14);

  sprintf(name,"inc_1j_nomatch_gen_deltaphi_jettot_Z_ptj_20_30_jrap_%i",i);
  h_inc_1j_nomatch_gen_deltaphi_jettot_Z_ptj_20_30_jrap[i]= new TH1D (name,name, 40,0.,3.14);
  sprintf(name,"inc_1j_nomatch_gen_deltaphi_jettot_Z_ptj_30_50_jrap_%i",i);
  h_inc_1j_nomatch_gen_deltaphi_jettot_Z_ptj_30_50_jrap[i]= new TH1D (name,name, 40,0.,3.14);
  sprintf(name,"inc_1j_nomatch_gen_deltaphi_jettot_Z_ptj_50_100_jrap_%i",i);
  h_inc_1j_nomatch_gen_deltaphi_jettot_Z_ptj_50_100_jrap[i]= new TH1D (name,name, 40,0.,3.14);
  sprintf(name,"inc_1j_nomatch_gen_deltaphi_jettot_Z_ptj_100_1000_jrap_%i",i);
  h_inc_1j_nomatch_gen_deltaphi_jettot_Z_ptj_100_1000_jrap[i]= new TH1D (name,name, 40,0.,3.14);

  sprintf(name,"ratio_jetpt30_50_zpt_dphi_28_%i",i);
  h_ratio_jetpt30_50_zpt_dphi_28[i]= new TH1D (name,name, 30,0.,3.);
  sprintf(name,"ratio_jetpt50_inf_zpt_dphi_28_%i",i);
  h_ratio_jetpt50_inf_zpt_dphi_28[i]= new TH1D (name,name, 30,0.,3.);

  sprintf(name,"ratio_jetpt30_50_zpt_dphi_28_up_%i",i);
  h_ratio_jetpt30_50_zpt_dphi_28_up[i]= new TH1D (name,name, 30,0.,3.);
  sprintf(name,"ratio_jetpt50_inf_zpt_dphi_28_up_%i",i);
  h_ratio_jetpt50_inf_zpt_dphi_28_up[i]= new TH1D (name,name, 30,0.,3.);

  sprintf(name,"ratio_jetpt30_50_zpt_dphi_28_dn_%i",i);
  h_ratio_jetpt30_50_zpt_dphi_28_dn[i]= new TH1D (name,name, 30,0.,3.);
  sprintf(name,"ratio_jetpt50_inf_zpt_dphi_28_dn_%i",i);
  h_ratio_jetpt50_inf_zpt_dphi_28_dn[i]= new TH1D (name,name, 30,0.,3.);



 sprintf(name,"ratio_Zpt0_30_zpt_dphi_28_%i",i);
  h_ratio_Zpt0_30_zpt_dphi_28[i]= new TH1D (name,name, 30,0.,3.);
  sprintf(name,"ratio_Zpt30_50_zpt_dphi_28_%i",i);
  h_ratio_Zpt30_50_zpt_dphi_28[i]= new TH1D (name,name, 30,0.,3.);
  sprintf(name,"ratio_Zpt50_inf_zpt_dphi_28_%i",i);
  h_ratio_Zpt50_inf_zpt_dphi_28[i]= new TH1D (name,name, 30,0.,3.);

  sprintf(name,"ratio_Zpt0_30_zpt_dphi_28_up_%i",i);
  h_ratio_Zpt0_30_zpt_dphi_28_up[i]= new TH1D (name,name, 30,0.,3.);
  sprintf(name,"ratio_Zpt30_50_zpt_dphi_28_up_%i",i);
  h_ratio_Zpt30_50_zpt_dphi_28_up[i]= new TH1D (name,name, 30,0.,3.);
  sprintf(name,"ratio_Zpt50_inf_zpt_dphi_28_up_%i",i);
  h_ratio_Zpt50_inf_zpt_dphi_28_up[i]= new TH1D (name,name, 30,0.,3.);

  sprintf(name,"ratio_Zpt0_30_zpt_dphi_28_dn_%i",i);
  h_ratio_Zpt0_30_zpt_dphi_28_dn[i]= new TH1D (name,name, 30,0.,3.);
  sprintf(name,"ratio_Zpt30_50_zpt_dphi_28_dn_%i",i);
  h_ratio_Zpt30_50_zpt_dphi_28_dn[i]= new TH1D (name,name, 30,0.,3.);
  sprintf(name,"ratio_Zpt50_inf_zpt_dphi_28_dn_%i",i);
  h_ratio_Zpt50_inf_zpt_dphi_28_dn[i]= new TH1D (name,name, 30,0.,3.);


  sprintf(name,"inc_ljet_pt_y_fb_%i",i);
  h_inc_ljet_pt_y_fb[i] = new TH1D (name,name, range_2d[i], j_pT_range3[i]);

  sprintf(name,"inc_ljet_pt_y_fb_dphi_lt_1p5_%i",i);
  h_inc_ljet_pt_y_fb_dphi_lt_1pt5[i] = new TH1D (name,name, range_2d[i], j_pT_range3[i]);

  sprintf(name,"inc_ljet_pt_y_fb_dphi_lt_1pt0_%i",i);
  h_inc_ljet_pt_y_fb_dphi_lt_1pt0[i] = new TH1D (name,name, range_2d[i], j_pT_range3[i]);

 sprintf(name,"inc_ljet_pt_y_fb_match_%i",i);
  h_inc_ljet_pt_y_fb_match[i] = new TH1D (name,name, range_2d[i], j_pT_range3[i]);

 sprintf(name,"inc_ljet_pt_y_fb_nomatch_%i",i);
  h_inc_ljet_pt_y_fb_nomatch[i] = new TH1D (name,name, range_2d[i], j_pT_range3[i]);


  sprintf(name,"inc_ljet_pt_y_fixbin_%i",i);
  h_inc_ljet_pt_y_fixbin[i] = new TH1D (name,name, 200,0.,1000. );


  sprintf(name,"inc_ljet_y_pt_%i",i);
  h_inc_ljet_y_pt[i] = new TH1D (name,name, 7,j_Y_range);
  sprintf(name,"inc_ljet_pt_y_z0_%i",i);
  h_inc_ljet_pt_y_z0[i] = new TH1D (name,name, 8, j_pT_range);
  sprintf(name,"inc_ljet_pt_y_z1_%i",i);
  h_inc_ljet_pt_y_z1[i] = new TH1D (name,name, 8, j_pT_range);

  sprintf(name,"inc_ljet_pt_y_fb_z0_%i",i);
  h_inc_ljet_pt_y_fb_z0[i] = new TH1D (name,name, range_2d[i], j_pT_range3[i]);
  sprintf(name,"inc_ljet_pt_y_fb_z1_%i",i);
  h_inc_ljet_pt_y_fb_z1[i] = new TH1D (name,name, range_2d[i], j_pT_range3[i]);

  sprintf(name,"inc_ljet_pt_y_fb_pos_z0_%i",i);
  h_inc_ljet_pt_y_fb_z0_pos[i] = new TH1D (name,name, range_3d[i], j_pT_range33[i]);
  sprintf(name,"inc_ljet_pt_y_fb_pos_z1_%i",i);
  h_inc_ljet_pt_y_fb_z1_pos[i] = new TH1D (name,name, range_3d[i], j_pT_range33[i]);

  sprintf(name,"inc_ljet_pt_y_fb_neg_z0_%i",i);
  h_inc_ljet_pt_y_fb_z0_neg[i] = new TH1D (name,name, range_3d[i], j_pT_range33[i]);
  sprintf(name,"inc_ljet_pt_y_fb_neg_z1_%i",i);
  h_inc_ljet_pt_y_fb_z1_neg[i] = new TH1D (name,name, range_3d[i], j_pT_range33[i]);

  sprintf(name,"inc_ljet_pt_sljetpt_%i",i);
  h_inc_ljet_pt_sljetpt[i] = new TH1D (name,name, 8, j_pT_range);
  sprintf(name,"inc_ljet_pt_sljetpt_z0_%i",i);
  h_inc_ljet_pt_sljetpt_z0[i] = new TH1D (name,name, 8, j_pT_range);
  sprintf(name,"inc_ljet_pt_sljetpt_z1_%i",i);
  h_inc_ljet_pt_sljetpt_z1[i] = new TH1D (name,name, 8, j_pT_range);
  sprintf(name,"inc_ljet_y_pt_z0_%i",i);
  h_inc_ljet_y_pt_z0[i] = new TH1D (name,name, 7,j_Y_range);
  sprintf(name,"inc_ljet_y_pt_z1_%i",i);
  h_inc_ljet_y_pt_z1[i] = new TH1D (name,name, 7,j_Y_range);


  sprintf(name,"inc_ljet_pt_y_Gen_%i",i);
  h_inc_ljet_pt_y_Gen[i] = new TH1D(name,name,8,j_pT_range);
  sprintf(name,"inc_ljet_Z_y_Gen_%i",i);
  h_inc_ljet_Z_y_Gen[i] = new TH1D(name,name, 20,-2.5,2.5);

  sprintf(name,"inc_ljet_pt_y_fb_Gen_%i",i);
  h_inc_ljet_pt_y_fb_Gen[i] = new TH1D(name,name,range_2d[i], j_pT_range3[i]);

  sprintf(name,"inc_ljet_y_pt_Gen_%i",i);
  h_inc_ljet_y_pt_Gen[i] = new TH1D(name,name,7,j_Y_range);
  sprintf(name,"inc_ljet_pt_y_Gen_z0_%i",i);
  h_inc_ljet_pt_y_z0_Gen[i] = new TH1D(name,name,8,j_pT_range);
  sprintf(name,"inc_ljet_pt_y_Gen_z1_%i",i);
  h_inc_ljet_pt_y_z1_Gen[i] = new TH1D(name,name,8,j_pT_range);

  sprintf(name,"inc_ljet_pt_y_fb_Gen_z0_%i",i);
  h_inc_ljet_pt_y_fb_z0_Gen[i] = new TH1D(name,name,range_2d[i], j_pT_range3[i]);
  sprintf(name,"inc_ljet_pt_y_fb_Gen_z1_%i",i);
  h_inc_ljet_pt_y_fb_z1_Gen[i] = new TH1D(name,name,range_2d[i], j_pT_range3[i]);

  sprintf(name,"inc_ljet_pt_y_fb_pos_Gen_z0_%i",i);
  h_inc_ljet_pt_y_fb_z0_pos_Gen[i] = new TH1D(name,name,range_3d[i], j_pT_range33[i]);
  sprintf(name,"inc_ljet_pt_y_fb_pos_Gen_z1_%i",i);
  h_inc_ljet_pt_y_fb_z1_pos_Gen[i] = new TH1D(name,name,range_3d[i], j_pT_range33[i]);

  sprintf(name,"inc_ljet_pt_y_fb_neg_Gen_z0_%i",i);
  h_inc_ljet_pt_y_fb_z0_neg_Gen[i] = new TH1D(name,name,range_3d[i], j_pT_range33[i]);
  sprintf(name,"inc_ljet_pt_y_fb_neg_Gen_z1_%i",i);
  h_inc_ljet_pt_y_fb_z1_neg_Gen[i] = new TH1D(name,name,range_3d[i], j_pT_range33[i]);


  sprintf(name,"inc_ljet_pt_sljetpt_Gen_%i",i);
  h_inc_ljet_pt_sljetpt_Gen[i] = new TH1D (name,name, 8, j_pT_range);
  sprintf(name,"inc_ljet_pt_sljetpt_Gen_z0_%i",i);
  h_inc_ljet_pt_sljetpt_z0_Gen[i] = new TH1D (name,name, 8, j_pT_range);
  sprintf(name,"inc_ljet_pt_sljetpt_Gen_z1_%i",i);
  h_inc_ljet_pt_sljetpt_z1_Gen[i] = new TH1D (name,name, 8, j_pT_range);

  sprintf(name,"inc_ljet_y_pt_Gen_z0_%i",i);
  h_inc_ljet_y_pt_z0_Gen[i] = new TH1D(name,name,7,j_Y_range);
  sprintf(name,"inc_ljet_y_pt_Gen_z1_%i",i);
  h_inc_ljet_y_pt_z1_Gen[i] = new TH1D(name,name,7,j_Y_range);


  sprintf(name,"inc_ljet_y_pt_Gen_VS_reco_%i",i);
  h_inc_ljet_y_pt_Gen_VS_reco[i] = new TH2D(name,name,7,j_Y_range,7,j_Y_range);
  sprintf(name,"inc_ljet_pt_y_Gen_VS_reco_%i",i);
  h_inc_ljet_pt_y_Gen_VS_reco[i] = new TH2D(name,name,8,j_pT_range,8,j_pT_range);

  sprintf(name,"inc_ljet_Z_y_Gen_VS_reco_%i",i);
  h_inc_ljet_Z_y_Gen_VS_reco[i] = new TH2D(name,name,20,-2.5,2.5, 20,-2.5,2.5);

  sprintf(name,"inc_ljet_pt_y_Gen_VS_reco_z0_%i",i);
  h_inc_ljet_pt_y_z0_Gen_VS_reco[i] = new TH2D(name,name,8,j_pT_range,8,j_pT_range);
  sprintf(name,"inc_ljet_pt_y_Gen_VS_reco_z1_%i",i);
  h_inc_ljet_pt_y_z1_Gen_VS_reco[i] = new TH2D(name,name,8,j_pT_range,8,j_pT_range);

  sprintf(name,"inc_ljet_pt_y_fb_Gen_VS_reco_%i",i);
  h_inc_ljet_pt_y_fb_Gen_VS_reco[i] = new TH2D(name,name,range_2d[i], j_pT_range3[i],range_2d[i], j_pT_range3[i]);
  sprintf(name,"inc_ljet_pt_y_fb_Gen_VS_reco_z0_%i",i);
  h_inc_ljet_pt_y_fb_z0_Gen_VS_reco[i] = new TH2D(name,name,range_2d[i], j_pT_range3[i],range_2d[i], j_pT_range3[i]);
  sprintf(name,"inc_ljet_pt_y_fb_Gen_VS_reco_z1_%i",i);
  h_inc_ljet_pt_y_fb_z1_Gen_VS_reco[i] = new TH2D(name,name,range_2d[i], j_pT_range3[i],range_2d[i], j_pT_range3[i]);

  sprintf(name,"inc_ljet_pt_y_fb_pos_Gen_VS_reco_z0_%i",i);
  h_inc_ljet_pt_y_fb_z0_pos_Gen_VS_reco[i] = new TH2D(name,name,range_3d[i], j_pT_range33[i],range_3d[i], j_pT_range33[i]);
  sprintf(name,"inc_ljet_pt_y_fb_pos_Gen_VS_reco_z1_%i",i);
  h_inc_ljet_pt_y_fb_z1_pos_Gen_VS_reco[i] = new TH2D(name,name,range_3d[i], j_pT_range33[i],range_3d[i], j_pT_range33[i]);

  sprintf(name,"inc_ljet_pt_y_fb_neg_Gen_VS_reco_z0_%i",i);
  h_inc_ljet_pt_y_fb_z0_neg_Gen_VS_reco[i] = new TH2D(name,name,range_3d[i], j_pT_range33[i],range_3d[i], j_pT_range33[i]);
  sprintf(name,"inc_ljet_pt_y_fb_neg_Gen_VS_reco_z1_%i",i);
  h_inc_ljet_pt_y_fb_z1_neg_Gen_VS_reco[i] = new TH2D(name,name,range_3d[i], j_pT_range33[i],range_3d[i], j_pT_range33[i]);

  sprintf(name,"inc_ljet_pt_sljetpt_Gen_VS_reco_%i",i);
  h_inc_ljet_pt_sljetpt_Gen_VS_reco[i] = new TH2D(name,name,8,j_pT_range,8,j_pT_range);
  sprintf(name,"inc_ljet_pt_sljetpt_Gen_VS_reco_z0_%i",i);
  h_inc_ljet_pt_sljetpt_z0_Gen_VS_reco[i] = new TH2D(name,name,8,j_pT_range,8,j_pT_range);
  sprintf(name,"inc_ljet_pt_sljetpt_Gen_VS_reco_z1_%i",i);
  h_inc_ljet_pt_sljetpt_z1_Gen_VS_reco[i] = new TH2D(name,name,8,j_pT_range,8,j_pT_range);
  sprintf(name,"inc_ljet_y_pt_Gen_VS_reco_z0_%i",i);

  sprintf(name,"inc_ljet_y_pt_Gen_VS_reco_z0_%i",i);
  h_inc_ljet_y_pt_z0_Gen_VS_reco[i] = new TH2D(name,name,7,j_Y_range,7,j_Y_range);
  sprintf(name,"inc_ljet_y_pt_Gen_VS_reco_z1_%i",i);
  h_inc_ljet_y_pt_z1_Gen_VS_reco[i] = new TH2D(name,name,7,j_Y_range,7,j_Y_range);

sprintf(name,"inc_sljet_pt_y_%i",i);
  h_inc_sljet_pt_y[i] = new TH1D (name,name, 8, j_pT_range);
  sprintf(name,"inc_sljet_y_pt_%i",i);
  h_inc_sljet_y_pt[i] = new TH1D (name,name, 7,j_Y_range);
  sprintf(name,"inc_sljet_pt_y_z0_%i",i);
  h_inc_sljet_pt_y_z0[i] = new TH1D (name,name, 8, j_pT_range);
  sprintf(name,"inc_sljet_pt_y_z1_%i",i);
  h_inc_sljet_pt_y_z1[i] = new TH1D (name,name, 8, j_pT_range);
  sprintf(name,"inc_sljet_y_pt_z0_%i",i);
  h_inc_sljet_y_pt_z0[i] = new TH1D (name,name, 7,j_Y_range);
  sprintf(name,"inc_sljet_y_pt_z1_%i",i);
  h_inc_sljet_y_pt_z1[i] = new TH1D (name,name, 7,j_Y_range);


  sprintf(name,"inc_sljet_pt_y_Gen_%i",i);
  h_inc_sljet_pt_y_Gen[i] = new TH1D(name,name,8,j_pT_range);
  sprintf(name,"inc_sljet_y_pt_Gen_%i",i);
  h_inc_sljet_y_pt_Gen[i] = new TH1D(name,name,7,j_Y_range);
  sprintf(name,"inc_sljet_pt_y_Gen_z0_%i",i);
  h_inc_sljet_pt_y_z0_Gen[i] = new TH1D(name,name,8,j_pT_range);
  sprintf(name,"inc_sljet_pt_y_Gen_z1_%i",i);
  h_inc_sljet_pt_y_z1_Gen[i] = new TH1D(name,name,8,j_pT_range);
  sprintf(name,"inc_sljet_y_pt_Gen_z0_%i",i);
  h_inc_sljet_y_pt_z0_Gen[i] = new TH1D(name,name,7,j_Y_range);
  sprintf(name,"inc_sljet_y_pt_Gen_z1_%i",i);
  h_inc_sljet_y_pt_z1_Gen[i] = new TH1D(name,name,7,j_Y_range);

  sprintf(name,"inc_sljet_pt_y_Gen_VS_reco_%i",i);
  h_inc_sljet_pt_y_Gen_VS_reco[i] = new TH2D(name,name,8,j_pT_range,8,j_pT_range);
  sprintf(name,"inc_sljet_y_pt_Gen_VS_reco_%i",i);
  h_inc_sljet_y_pt_Gen_VS_reco[i] = new TH2D(name,name,7,j_Y_range,7,j_Y_range);
  sprintf(name,"inc_sljet_pt_y_Gen_VS_reco_z0_%i",i);
  h_inc_sljet_pt_y_z0_Gen_VS_reco[i] = new TH2D(name,name,8,j_pT_range,8,j_pT_range);
  sprintf(name,"inc_sljet_pt_y_Gen_VS_reco_z1_%i",i);
  h_inc_sljet_pt_y_z1_Gen_VS_reco[i] = new TH2D(name,name,8,j_pT_range,8,j_pT_range);
  sprintf(name,"inc_sljet_y_pt_Gen_VS_reco_z0_%i",i);
  h_inc_sljet_y_pt_z0_Gen_VS_reco[i] = new TH2D(name,name,7,j_Y_range,7,j_Y_range);
  sprintf(name,"inc_sljet_y_pt_Gen_VS_reco_z1_%i",i);
  h_inc_sljet_y_pt_z1_Gen_VS_reco[i] = new TH2D(name,name,7,j_Y_range,7,j_Y_range);

}

//X-Check histograms
  double n_jet_range[9]={-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5};

  TH1D* h_n_jet_30_unfold  = new TH1D("n_jet_30_unfold","n_jet_30_unfold",8,n_jet_range);
  TH1D* h_n_jet_30_gen  = new TH1D("n_jet_30_gen","n_jet_30_gen",8,n_jet_range);
  TH2* h_n_jet_30_Gen_VS_reco = new TH2D("n_jet_30_Gen_VS_reco","n_jet_30_Gen_VS_reco",8,n_jet_range,8,n_jet_range);

  TH1D* h_n_jet_30_unselected = new TH1D("n_jet_30_unselected","n_jet_30_unselected",9,n_jet_range);
//  TH2* h_jetpt_Gen_VS_reco= new TH2D("jetpt_Gen_VS_reco","jetpt_Gen_VS_reco",17,j_pT_range_unf,17,j_pT_range_unf);
//  TH2* h_jetY_Gen_VS_reco= new TH2D("h_jetY_Gen_VS_reco","h_jetY_Gen_VS_reco",25,-2.5,2.5,25,-2.5,2.5);
 TH1D* h_1jetpt  = new TH1D("jetpt_1","1jetpt",nJetPt_Zinc1jet,jetPt_Zinc1jet);
  h_1jetpt->Sumw2(); 
  TH1D* h_2jetpt  = new TH1D("jetpt_2","2jetpt",nJetPt_Zinc2jet,jetPt_Zinc2jet);
  h_2jetpt->Sumw2();
  //TH1D* h_jetY = new TH1D("jetY","jetY",24,-2.4,2.4);
  //h_jetY->Sumw2();
  TH1D* h_1jetY  = new TH1D("jetY_1","1jetY",40,0,4.7);
  h_1jetY->Sumw2();
  TH1D* h_2jetY  = new TH1D("jetY_2","2jetY",40,0,4.7);
  h_2jetY->Sumw2();

  //TH1D* h_Gjetpt  = new TH1D("Gjetpt","Gjetpt",nJetPt_Zinc1jet,jetPt_Zinc1jet);
  //h_Gjetpt->Sumw2();
  TH1D* h_1Gjetpt  = new TH1D("Gjetpt_1","G1jetpt",nJetPt_Zinc1jet,jetPt_Zinc1jet);
  h_1Gjetpt->Sumw2(); 
  TH1D* h_2Gjetpt  = new TH1D("Gjetpt_2","G2jetpt",nJetPt_Zinc2jet,jetPt_Zinc2jet);
  h_2Gjetpt->Sumw2();
  //TH1D* h_GjetY = new TH1D("GjetY","GjetY",24,-2.4,2.4);
  //h_GjetY->Sumw2();
  TH1D* h_1GjetY  = new TH1D("GjetY_1","G1jetY",40,0,4.7);
  h_1GjetY->Sumw2();
  TH1D* h_2GjetY  = new TH1D("GjetY_2","G2jetY",40,0,4.7);
  h_2GjetY->Sumw2();

 TH1D* h_3Gjetpt  = new TH1D("Gjetpt_3","G3jetpt",nJetPt_Zinc2jet,jetPt_Zinc2jet);
  h_3Gjetpt->Sumw2();
  TH1D* h_4Gjetpt  = new TH1D("Gjetpt_4","G4jetpt",nJetPt_Zinc2jet,jetPt_Zinc2jet);
  h_4Gjetpt->Sumw2();
  TH1D* h_3GjetY  = new TH1D("GjetY_3","G3jetY",24,-2.4,2.4);
  h_3GjetY->Sumw2();
  TH1D* h_4GjetY  = new TH1D("GjetY_4","G4jetY",24,-2.4,2.4);
  h_4GjetY->Sumw2();

  TH1D* h_5Gjetpt  = new TH1D("Gjetpt_5","G5jetpt",nJetPt_Zinc2jet,jetPt_Zinc2jet);
  h_5Gjetpt->Sumw2();
  TH1D* h_6Gjetpt  = new TH1D("Gjetpt_6","G6jetpt",nJetPt_Zinc2jet,jetPt_Zinc2jet);
  h_6Gjetpt->Sumw2();
  TH1D* h_5GjetY  = new TH1D("GjetY_5","G5jetY",24,-2.4,2.4);
  h_5GjetY->Sumw2();
  TH1D* h_6GjetY  = new TH1D("GjetY_6","G6jetY",24,-2.4,2.4);
  h_6GjetY->Sumw2();



  //TH2* h_jetpt_Gen_VS_reco= new TH2D("jetpt_Gen_VS_reco","jetpt_Gen_VS_reco",nJetPt_Zinc1jet,jetPt_Zinc1jet,nJetPt_Zinc1jet,jetPt_Zinc1jet);
  //TH2* h_jetY_Gen_VS_reco= new TH2D("h_jetY_Gen_VS_reco","h_jetY_Gen_VS_reco",24,-2.4,2.4,24,-2.4,2.4);
  TH2* h_1jetpt_Gen_VS_reco= new TH2D("jet1pt_Gen_VS_reco","jet1pt_Gen_VS_reco",nJetPt_Zinc1jet,jetPt_Zinc1jet,nJetPt_Zinc1jet,jetPt_Zinc1jet);
  TH2* h_1jetY_Gen_VS_reco = new TH2D("jet1Y_Gen_VS_reco","jet1Y_Gen_VS_reco",40,0,4.7,40,0,4.7);
  TH2* h_2jetpt_Gen_VS_reco= new TH2D("jet2pt_Gen_VS_reco","2jet2pt_Gen_VS_reco",nJetPt_Zinc2jet,jetPt_Zinc2jet,nJetPt_Zinc2jet,jetPt_Zinc2jet);
  TH2* h_2jetY_Gen_VS_reco= new TH2D("jet2Y_Gen_VS_reco","jet2Y_Gen_VS_reco",40,0,4.7,40,0,4.7);


