#!/bin/sh
#cd /storage_mnt/storage/user/mgul/Higgs_tottbar/anlyzer765/CMSSW_7_6_5/src/backgourd/ttbar
#/user/mgul/Higgs_tottbar/analyzer8024/CMSSW_8_0_24/src/Tupel/Tupel/Permutaion_check_ntuple.root
#/pnfs/iihe/cms/store/user/mgul/final_Tuple_reminiAOD/muon/MC13TeV_TTJets/MergednTuple_MC13TeV_TTJets_76.root
#/user/mgul/Higgs_tottbar/analyzer8024/CMSSW_8_0_24/src/Tupel/Tupel/ntuple.root
#source $VO_CMS_SW_DIR/cmsset_default.sh 
#/user/mgul/Higgs_tottbar/analyzer8024/CMSSW_8_0_24/src/Tupel/Tupel/ntuple.root
#/pnfs/iihe/cms/store/user/mgul/Tuple_reminiAOD/muon/data_mc/fed7288/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_MC13TeV_TTJets/170306_002531/0000/ntuple_990.root
eval `scram runtime -sh` 
root -l -b <<EOF 
.x /user/mgul/Higgs_tottbar/analyzer8024/CMSSW_8_0_24/src/Tupel/Tupel/qcd_estim/el_qcdEstimator.cc+ ("/pnfs/iihe/cms/store/user/mgul/Tuple_reminiAOD/muon/data_mc/fed7288/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_MC13TeV_TTJets/170306_002531/0000/ntuple_986.root","testing/qcdntuple") 
.q; 
EOF
