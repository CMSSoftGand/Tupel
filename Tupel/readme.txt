Instructions for ntuple Producer
ssh -Y <uname>@lxplus.cern.ch
cmsrel CMSSW_8_0_11
cd CMSSW_8_0_11/src
cmsenv
git-cms-merge-topic 13960
//to fetch the most recent pseudotop producer.
git clone -b Tupel_80X git@github.com:UGent/Tupel
scram b -j8
// met corrections
git cms-merge-topic cms-met:metTool80X
scram b -j8

// eos mount (active for ~24 hours)
mkdir ~/eos
eosmount ~/eos

cd /Tupel/Tupel/
// Make symbolic link
ln -s data/era2016/Spring16_25nsV6_DATA.db
ln -s data/era2016/Spring16_25nsV6_MC.db
// To run the analyzer locally
cmsRun scripts/simple_run_80X_cfg.py
// To submit the whole samples via crab run the following command
python scripts/submitToGrid.py -j data/era2016/samples.json -c ${CMSSW_BASE}/src/Tupel/Tupel/scripts/simple_run_80X_cfg.py --lfn /store/user/mgul/test -s


extra info:

  channel    = cms.string( 'smu' ), 

is for skimming the trees at production level. Options:
dimu, smu, dielec,selec, 
for no skimming, use noselection

  keepparticlecoll    = cms.bool(False),
to keep the particle flow objects and gen level particles. Needed for UE studies. To keep, switch to True
