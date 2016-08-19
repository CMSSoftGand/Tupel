#!/bin/sh
#cd /storage_mnt/storage/user/mgul/Higgs_tottbar/anlyzer765/CMSSW_7_6_5/src/backgourd/ttbar
#/pnfs/iihe/cms/store/user/mgul/Tuples765/bkg/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/160806_130240/0000/simulation_ntuple_999.root
#/pnfs/iihe/cms/store/user/mgul/Tuples765/bkg/TT_TuneCUETP8M1_13TeV-powheg-pythia8/TT_TuneCUETP8M1_13TeV-powheg-pythia8/160806_130135/0000/simulation_ntuple_1.root
source $VO_CMS_SW_DIR/cmsset_default.sh 
eval `scram runtime -sh` 
root -l -b <<EOF 
.x simpleReader.C+ ("/pnfs/iihe/cms/store/user/mgul/Tupel_8_0_11/1c09f31/MC13TeV_SingleTbar_t/MergednTuple_0.root","testoutput_0000_1" ) 
.q; 
EOF
