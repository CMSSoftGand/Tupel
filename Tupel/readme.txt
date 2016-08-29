Instructions for ntuple Producer
ssh -Y <uname>@lxplus.cern.ch
cmsrel CMSSW_8_0_11
cd CMSSW_8_0_11/src
cmsenv
git-cms-merge-topic 13960
//to fetch the most recent pseudotop producer.
git clone -b Tuple_MiniAOD_TTbar_8x git@github.com:UGent/Tupel
scram b -j8
// met corrections
git cms-merge-topic cms-met:metTool80X
scram b -j8

// eos mount (active for ~24 hours)
mkdir ~/eos
eosmount ~/eos

cd Tupel/Tupel/
// Make symbolic link
ln -s data/era2016/Spring16_25nsV6_DATA.db
ln -s data/era2016/Spring16_25nsV6_MC.db
// To run the analyzer locally
cmsRun scripts/simple_run_80X_cfg.py
// Make a test for crab using test_samples.json;
// and then submit the whole samples via crab run the following command
python scripts/submitToGrid.py -j data/era2016/test_samples.json -c ${CMSSW_BASE}/src/Tupel/Tupel/scripts/simple_run_80X_cfg.py --lfn /store/user/mgul/test -s

//  get crab report
crab report grid/crab_Data13TeV_SingleMuon_2016B/
//  Merge json files
mergeJSON.py grid/crab_Data13TeV_SingleMuon_2016B/results/processedLumis.json grid/crab_Data13TeV_SingleMuon_2016C/results/processedLumis.json grid/crab_Data13TeV_SingleMuon_2016D/results/processedLumis.json --output data/era2016/Data13TeV_SingleMuon_lumis.json

// You can then run the brilcalc tool to get the integrated luminosity in total and per run (see https://twiki.cern.ch/twiki/bin/view/CMS/2015LumiNormtag for more details).
export PATH=$HOME/.local/bin:/afs/cern.ch/cms/lumi/brilconda-1.0.3/bin:$PATH
brilcalc lumi --normtag /afs/cern.ch/user/l/lumipro/public/normtag_file/normtag_DATACERT.json -i data/era2016/Data13TeV_SingleMuon_lumis.json

// To add output files and reduce the numbers of files use the following [1] command for iihe batch but there is a proxy problem during jot submission which should to be fixed. It creates a tmp_combine directory which consists of all the .sh files for datasets. (You can submit one-by-one. ) 
[1]: python scripts/submitCheckProductionIntegrity.py -i /store/user/mgul/tuple_8011_new1/1c09f31/ -o /store/user/mgul/Htottbar/1c09f31 

// Pileup weighting. To update the pileup distributions run the script below. It will store the data pileup distributions for different min.bias cross section in data/pileupWgts.root
python scripts/runPileupEstimation.py --json data/era2016/Data13TeV_SingleMuon_lumis.json --out data/era2016/pileupWgts.root

// B-tagging. To apply corrections to the simulation one needs the expected efficiencies stored somwewhere. The script below will project the jet pT spectrum from the TTbar sample before and after applying b-tagging, to compute the expecte efficiencies (this is done only for one file but could to done by many files giving the directory path). The result will be stored in data/era2016/expTageff.root
>> python scripts/saveExpectedBtagEff.py -i /store/user/mgul/ntuple.root -o data/era2016/expTageff.root 

// to run the analyzer locally [1] and submit the for the whole sample is [2]:
[1] ./job.sh 
[2] python scripts/bash_ttbar.py -i /store/user/mgul/Tupel_8_0_11/1c09f31

// SF values can be calculated using script calc_SF.cpp

extra info:

  channel    = cms.string( 'smu' ), 

is for skimming the trees at production level. Options:
dimu, smu, dielec,selec, 
for no skimming, use noselection

  keepparticlecoll    = cms.bool(False),
to keep the particle flow objects and gen level particles. Needed for UE studies. To keep, switch to True
