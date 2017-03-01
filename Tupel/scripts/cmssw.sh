#!/bin/bash

# Recent 80X release
release="CMSSW_8_0_26_patch1"
scramv1 project CMSSW $release  # cmsrel alias expanded

cd $release/src
eval `scramv1 runtime -sh`  # this is cmsenv alias expanded


git cms-init

# Updates for MET [1]
# [1] https://twiki.cern.ch/twiki/bin/view/CMS/MissingETUncertaintyPrescription?rev=57#Instructions_for_8_0_X_X_20_for
git cms-merge-topic -u cms-met:METRecipe_8020
git cms-merge-topic -u cms-met:METRecipe_80X_part2

#electron ID summer-16
#https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2?rev=39#Recipe_for_regular_users_for_8_0
git cms-merge-topic ikrav:egm_id_80X_v2

#deepCSV
git cms-merge-topic -u mverzett:DeepFlavour-from-CMSSW_8_0_21
mkdir RecoBTag/DeepFlavour/data/
cd RecoBTag/DeepFlavour/data/
wget http://home.fnal.gov/~verzetti//DeepFlavour/training/DeepFlavourNoSL.json
cd -
scram b -j8
