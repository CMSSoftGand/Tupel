

#include "Tupel/Tupel/interface/correction.h"



std::map<BTagEntry::JetFlavor,BTagCalibrationReader *> getBTVcalibrationReaders(TString era,BTagEntry::OperatingPoint btagOP)
{
    //start the btag calibration
      TString btagUncUrl(era+"/btagSFactors.csv");
       gSystem->ExpandPathName(btagUncUrl);
          BTagCalibration btvcalib("csvv2", btagUncUrl.Data());
    
    //start calibration readers for b,c and udsg separately including the up/down variations
std::map<BTagEntry::JetFlavor,BTagCalibrationReader *> btvCalibReaders;
  btvCalibReaders[BTagEntry::FLAV_B]=new BTagCalibrationReader(btagOP, "central", {"up", "down"});
    btvCalibReaders[BTagEntry::FLAV_B]->load(btvcalib,BTagEntry::FLAV_B,"mujets");
      btvCalibReaders[BTagEntry::FLAV_C]=new BTagCalibrationReader(btagOP, "central", {"up", "down"});
        btvCalibReaders[BTagEntry::FLAV_C]->load(btvcalib,BTagEntry::FLAV_C,"mujets");
          btvCalibReaders[BTagEntry::FLAV_UDSG]=new BTagCalibrationReader(btagOP, "central", {"up", "down"});
            btvCalibReaders[BTagEntry::FLAV_UDSG]->load(btvcalib,BTagEntry::FLAV_UDSG,"incl");

              //all done
                return btvCalibReaders;
                }
                
