#include "../interface/Permutation.h"
#include "../interface/TTBarSolver.h"
#include "../interface/NeutrinoSolver.h"
#include <iostream>
#include <limits>


TTBarSolver::TTBarSolver(bool active){  
	if(!active) return;

}

TTBarSolver::~TTBarSolver()
{}

void TTBarSolver::Solve(Permutation &hyp, bool lazy)
{
  if(!lazy && !hyp.IsComplete()) {                          
    Logger::log().fatal() << "The permutation you are trying to solve is not complete!" << std::endl;
    throw 42;
  }
  bool USEBTAG_(false),useptratios_(false), usewjetqgtag_(false);
  bool USEMASS_(true) ,USENS_(true);
  TH1 *N_right_;
  TH2 *WTmass_right_;
  TFile *probfile=TFile::Open("/afs/cern.ch/work/m/mgul/public/Hto_ttbar/cmssw8026_patch1/CMSSW_8_0_26_patch1/src/Tupel/Tupel/data/moriond17/htt_permutations.root");
  TDirectory *dir = (TDirectory*) probfile->Get("nosys");
  gDirectory->Get("nosys");
  N_right_     =(TH1 *)gDirectory->Get("nosys/nusolver_chi2_right");
  N_right_->Scale(1./N_right_->Integral(),"width");
//  N_right_->Scale(1./N_right_->Integral("width"));
  WTmass_right_=(TH2 *)gDirectory->Get("nosys/mWhad_vs_mtophad_right");
  WTmass_right_->Scale(1./WTmass_right_->Integral("width"));
	double nschi    = numeric_limits<double>::max();
	double res      = numeric_limits<double>::max();
	double nstest   = numeric_limits<double>::max();
	double masstest = numeric_limits<double>::max();
	double btagtest = numeric_limits<double>::max();
	double ptratios = numeric_limits<double>::max();
	double qgtest   = numeric_limits<double>::max();
	NeutrinoSolver NS(hyp.L(), hyp.BLep(), mw_, mtop_);
	hyp.Nu(NS.GetBest(hyp.MET()->Px(), hyp.MET()->Py(), 1, 1, 0., nschi)); //ignore MET covariance matrix, take bare distance
	double mwhad = hyp.WHad().M();
	double mthad = hyp.THad().M();
//  auto min_max = [] (const double& a, const double& b) -> pair<const double,const double> {
//        return (b<a) ? std::make_pair(b, a) : std::make_pair(a, b); };
//----------------------New Method-------------
    if(nschi > 0 && N_right_->GetXaxis()->GetXmin()>=0. && N_right_->GetXaxis()->GetXmax()<=200.){
        int nubin = N_right_->FindFixBin(Sqrt(nschi));
        if(nubin <= N_right_->GetNbinsX())
            nstest = -1.*std::log(N_right_->GetBinContent(nubin));
//cout<<"this is nstest:  "<<"  nschi:  "<<nschi<<",  nstest:  "<<nstest<<endl;
        }
   
    if( WTmass_right_->GetXaxis()->GetXmax()<=500 && WTmass_right_->GetYaxis()->GetXmax()<=500 ) {
        int massbin = WTmass_right_->FindFixBin(mwhad, mthad);
        double massdisval_ = 0.;
//cout<<"mass bin: "<<massbin<<",  WTmass_right_->GetNbinsX():  "<<WTmass_right_->GetNbinsX()<<",  WTmass_right_->GetNbinsY():  "<<WTmass_right_->GetNbinsY()<<endl;
//        if(massbin <= WTmass_right_->GetNbinsX() && massbin <= WTmass_right_->GetNbinsY()) 
         if (!WTmass_right_->IsBinOverflow(massbin)) massdisval_ = WTmass_right_->GetBinContent(massbin);

          if(massdisval_ > 1.0E-10) masstest = -1.*Log(massdisval_);
        }

//---------------------------------------
  probfile->Close();

	res = 0.;
	if(USEMASS_) {res += masstest;}
	if(USENS_  ) {res += nstest;}
	hyp.Prob(res);
	hyp.NuChisq(nschi);
	hyp.NuDiscr(nstest);
  hyp.MassDiscr(masstest);

//cout<<"nschi:  "<<nschi<<",  nsdisc:  "<<nstest<<",  mass disc: "<<masstest<<",  full disc:  "<<res<<", mwhad: "<<mwhad<<", mthad: "<<mthad<<endl;
}

