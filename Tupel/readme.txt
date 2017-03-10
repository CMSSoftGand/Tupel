Instructions for ntuple Producer
ssh -Y <uname>@lxplus.cern.ch
#execute the cmssw.sh script in scripts dir
https://github.com/UGent/Tupel/tree/Tuple_MiniAOD_TTbar_8x/Tupel/scripts/cmssw.sh
./cmssw.sh or copy paste the lines
//to fetch the most recent pseudotop producer.
git clone -b Tuple_MiniAOD_TTbar_8x git@github.com:UGent/Tupel
scram b -j8

// eos mount (active for ~24 hours)
mkdir ~/eos
eosmount ~/eos

cd Tupel/Tupel/
// Make symbolic link if need correction from db
ln -s data/moriond17/Summer16_23Sep2016AllV2_DATA.db
ln -s data/moriond17/Summer16_23Sep2016V2_MC.db
// To run the analyzer locally
cmsRun scripts/simple_run_80X_cfg.py
// Make a test for crab using test_samples.json;
// and then submit the whole samples via crab run the following command
python scripts/submitToGrid.py -j data/era2016/test_samples.json -c ${CMSSW_BASE}/src/Tupel/Tupel/scripts/simple_run_80X_cfg.py --lfn /store/user/mgul/test -s

// See the crab status first using the following command and then resubmit with the second command
   tree -d -L 1 grid/ | awk '{printf("crab status -d grid/%s\n",$NF);}' >crab_status.sh && chmod 777 crab_status.sh && ./crab_status.sh
   tree -d -L 1 grid/ | awk '{printf("crab resubmit -d grid/%s\n",$NF);}' >crab_resubmit.sh && chmod 777 crab_resubmit.sh && ./crab_resubmit.sh

//  get crab report with the following script which will print the resultant files used for merging.
./scripts/crab_report.sh
//  Merge json files
mergeJSON.py grid/crab_Data13TeV_SingleMuon_2016B/results/processedLumis.json grid/crab_Data13TeV_SingleMuon_2016C/results/processedLumis.json grid/crab_Data13TeV_SingleMuon_2016D/results/processedLumis.json --output data/era2016/Data13TeV_SingleMuon_lumis.json

// You can then run the brilcalc tool to get the integrated luminosity in total and per run (see https://twiki.cern.ch/twiki/bin/view/CMS/2015LumiNormtag for more details).
export PATH=$HOME/.local/bin:/afs/cern.ch/cms/lumi/brilconda-1.0.3/bin:$PATH
brilcalc lumi --normtag /afs/cern.ch/user/l/lumipro/public/normtag_file/normtag_DATACERT.json -i data/era2016/Data13TeV_SingleMuon_lumis.json

// To add output files and reduce the numbers of files use the following [1],[2] command for iihe batch.Use checkProductionIntegrity.py for /pnfs/... It creates a tmp_combine directory which consists of all the .sh files for datasets. For scratch area a seperate file is available name for_scratchArea_checkProductionIntegrity.py.
  First create proxy: voms-proxy-init --voms cms --valid=168:00
[1]: python scripts/submitCheckProductionIntegrity.py -i /store/user/mgul/tuple_8011_new1/1c09f31/ -o /store/user/mgul/Htottbar/1c09f31 
if you want to use /scratch -->[2]: python scripts/submitCheckProductionIntegrity.py -i /store/user/mgul/tuple_8011_new1/1c09f31/ -o /scratch/mgul/Htottbar/1c09f31 

// Pileup weighting. To update the pileup distributions run the script below. It will store the data pileup distributions for different min.bias cross section in data/pileupWgts.root
 python scripts/runPileupEstimation.py --json data/era2016/Data13TeV_SingleMuon_lumisBCD.json --out data/era2016/pileupWgts.root --mbXsec 69200

// B-tagging. To apply corrections to the simulation one needs the expected efficiencies stored somwewhere. The script below will project the jet pT spectrum from the TTbar sample before and after applying b-tagging, to compute the expecte efficiencies (this is done only for one file but could to done by many files giving the directory path). The result will be stored in data/era2016/expTageff.root
>> python scripts/saveExpectedBtagEff.py -i /store/user/mgul/Htottbar/b2c6591/MC13TeV_WJets -o data/era2016/expTageff.root 

// to run the analyzer locally [1] and submit the for the whole sample is [2]:
[1] ./job.sh 
[2] python scripts/bash_ttbar.py -i /store/user/mgul/Tupel_8_0_11/1c09f31

// The output will be saved in root_files directory. To merge the output use the command:
 find root_files/ -type d -maxdepth 1 | awk '{printf("python scripts/mergeOutputs.py %s\n",$NF);}' >merge_output.sh && chmod 777 merge_output.sh && sh merge_output.sh

// copy the output files to a single directory ~/work/root_files
 find root_files/ -type f -maxdepth 2 | awk '{printf(" cp %s ~/work/root_files/\n",$NF);}' >copy_output.sh && printf '0a\nif [ ! -d ~/work/root_files ] \n   then mkdir ~/work/root_files \n   else echo "File exists" \n   fi \n.\nw\n' | ed copy_output.sh && chmod 777 copy_output.sh && sh copy_output.sh

// SF values can be calculated using script calc_SF.cpp

extra info:

  channel    = cms.string( 'smu' ), 

is for skimming the trees at production level. Options:
dimu, smu, dielec,selec, 
for no skimming, use noselection

  keepparticlecoll    = cms.bool(False),
to keep the particle flow objects and gen level particles. Needed for UE studies. To keep, switch to True
