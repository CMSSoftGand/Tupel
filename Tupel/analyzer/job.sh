#!/bin/sh
#cd /storage_mnt/storage/user/mgul/Higgs_tottbar/anlyzer765/CMSSW_7_6_5/src/backgourd/ttbar
#/pnfs/iihe/cms/store/user/mgul/Tupel_morion17_v3/data_mc/fbe1e1d/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/crab_MC13TeV_ST_sch/170124_210100/0000/ntuple_50.root
#/pnfs/iihe/cms/store/user/mgul/Tupel_morion17_v3/data_mc/fbe1e1d/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_MC13TeV_TTJets/170124_211257/0000/ntuple_992.root
#/pnfs/iihe/cms/store/user/mgul/final_Tupel_morion17_v2/muon/MC13TeV_TTJets/MergednTuple_MC13TeV_TTJets_90.root
#/user/mgul/Higgs_tottbar/analyzer8024/CMSSW_8_0_24/src/Tupel/Tupel/Permutaion_check_ntuple.root
#source $VO_CMS_SW_DIR/cmsset_default.sh 
#.x simpleReader.cc+ ("dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/mgul/tuple_09_10_16/8d7def9/TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_MC13TeV_TTJets_powheg/161009_154416/0000/ntuple_989.root","testing/testoutput_0000_5") 
eval `scram runtime -sh` 
root -l -b <<EOF 
.x /afs/cern.ch/work/m/mgul/public/Hto_ttbar/cmssw8026_patch1/CMSSW_8_0_26_patch1/src/Tupel/Tupel/analyzer/simpleReader.cc+ ("/afs/cern.ch/work/m/mgul/public/Hto_ttbar/cmssw8026_patch1/CMSSW_8_0_26_patch1/src/Tupel/Tupel/Permutaion_check_ntuple_v2.root","testing/ntuple") 
.q; 
EOF
