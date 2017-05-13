#!/bin/sh
#cd /storage_mnt/storage/user/mgul/Higgs_tottbar/anlyzer765/CMSSW_7_6_5/src/backgourd/ttbar
#/pnfs/iihe/cms/store/user/mgul/final_Tuple_reminiAOD/muon/MC13TeV_TTJets/MergednTuple_MC13TeV_TTJets_86.root
#/pnfs/iihe/cms/store/user/mgul/Tuple_reminiAOD/muon/data_mc/new_topStatus/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_MC13TeV_TTJets_backup_new_topStatus/170415_122815/0000/ntuple_977.root
#source $VO_CMS_SW_DIR/cmsset_default.sh 
#.x simpleReader.cc+ ("dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/mgul/tuple_09_10_16/8d7def9/TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_MC13TeV_TTJets_powheg/161009_154416/0000/ntuple_989.root","testing/testoutput_0000_5") 
eval `scram runtime -sh` 
root -l -b <<EOF 
.x /user/mgul/Higgs_tottbar/analyzer8024/CMSSW_8_0_24/src/Tupel/Tupel/template/template_08_02_17.cc+ ("/pnfs/iihe/cms/store/user/mgul/final_Tuple_reminiAOD/muon/MC13TeV_TTJets/MergednTuple_MC13TeV_TTJets_86.root","testing/testoutput_0000_5") 
.q; 
EOF
