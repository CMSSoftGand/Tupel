#ifndef PERMUTATION_H
#define PERMUTATION_H
#include <vector>
#include <limits>
#include <TLorentzVector.h>
//#include "Analyses/URTTbar/interface/URStreamer.h"
//#include "Analyses/URTTbar/interface/IDJet.h"
#include <iostream>
#include <stdexcept>
#include <sstream>
using namespace std;

//class IDMet;

std::ostream & operator<<(std::ostream &os, const TLorentzVector& p);

class Permutation
{
	private:
		double prob_ = numeric_limits<double>::max();
		double nu_chisq_          = numeric_limits<double>::max();
		double nu_discriminant_   = numeric_limits<double>::max();
		double btag_discriminant_ = numeric_limits<double>::max();
		double mass_discriminant_ = numeric_limits<double>::max();
		double qgtag_discriminant_ = numeric_limits<double>::max();
		double jratio_discriminant_ = numeric_limits<double>::max();
		TLorentzVector* wja_ = 0;
		TLorentzVector* wjb_ = 0;
		TLorentzVector* bjh_ = 0;
		TLorentzVector* bjl_ = 0;

		TLorentzVector* lep_ = 0;
    int lepcharge_=0;
		TLorentzVector* met_ = 0;
		TLorentzVector nu_;
		bool kinfit_ = false;
    size_t perm_jets_ = 0; //number of jets used for permutations
	public:
		Permutation() {}
		Permutation(TLorentzVector* wja, TLorentzVector* wjb, TLorentzVector* bjh, TLorentzVector* bjl, TLorentzVector* lep, TLorentzVector* met, int lcharge=0);
    int LepCharge() {return lepcharge_;}
    void LepCharge(int c) {lepcharge_ = c;}
		void Reset();
		bool IsWHadComplete() const {return(wja_ != 0 && wjb_ != 0);}
		bool IsTHadComplete() const {return(IsWHadComplete() && bjh_ != 0);}
		bool IsTLepComplete() const {return(bjl_ != 0 && lep_ != 0 && met_ != 0);}
		bool IsComplete() const {return (IsTHadComplete() && IsTLepComplete());}
		int NumBJets() const {return((bjl_ != 0 ? 1 : 0) + (bjh_ != 0 ? 1 : 0));}
		int NumWJets() const {return((wja_ != 0 ? 1 : 0) + (wjb_ != 0 ? 1 : 0));}
		int NumTTBarJets() const {return(NumBJets() + NumWJets());}
		TLorentzVector* WJa() const {return(wja_);}
		TLorentzVector* WJb() const {return(wjb_);}
		TLorentzVector* BHad() const {return(bjh_);}
		TLorentzVector* BLep() const {return(bjl_);}
    int permutating_jets() const {return perm_jets_;}
    void permutating_jets(int njets) { perm_jets_=njets;}
		TLorentzVector* L() const {return(lep_);}
		TLorentzVector* MET() const {return(met_);}
		void SetMET(TLorentzVector* met) {met_ = met;}
		void WJa(TLorentzVector* wja){wja_=wja;}
		void WJb(TLorentzVector* wjb){wjb_=wjb;}
		void BHad(TLorentzVector* bjh){bjh_=bjh;}
		void BLep(TLorentzVector* bjl){bjl_=bjl;}
		void L(TLorentzVector* lep){lep_=lep;}
		void MET(TLorentzVector* met){met_=met;}
  friend std::ostream & operator<<(std::ostream &os, const Permutation& p);

		TLorentzVector Nu() const {return(nu_);}
		void Nu(TLorentzVector v) {nu_ = v;}
		const TLorentzVector* NuPtr() const {return(&nu_);}
		TLorentzVector WHad() const {
			if(!WJa()) {
				stringstream err;
				err << "WJa is not defined in the permutation: " << *this;
				throw std::runtime_error(err.str());
			}
			return (WJb()) ? (*WJa() + *WJb()) : *WJa();
		}
		TLorentzVector bhadt() const {return(*BHad());}
		TLorentzVector blept() const {return(*BLep());}
		TLorentzVector WLep() const {return((*L() + Nu()));}
		TLorentzVector THad() const {return((WHad() + *BHad()));}
		TLorentzVector TLep() const {return((WLep() + *BLep()));}
    TLorentzVector LVect() const {return THad()+TLep();}

		double Prob()      const {return prob_             ;}
		double NuChisq() 	 const {return nu_chisq_         ;}
		double NuDiscr() 	 const {return nu_discriminant_  ;}
		double BDiscr()  	 const {return btag_discriminant_;}
		double MassDiscr() const {return mass_discriminant_;}
		double QGDiscr()   const {return qgtag_discriminant_;}
		double JRatioDiscr() const {return jratio_discriminant_;}
		void Prob(     double val) {prob_              = val;}
		void NuChisq(  double val) {nu_chisq_          = val;}
		void NuDiscr(  double val) {nu_discriminant_   = val;}
		void BDiscr(   double val) {btag_discriminant_ = val;}
		void MassDiscr(double val) {mass_discriminant_ = val;}
		void QGDiscr(  double val) {qgtag_discriminant_= val;} 
		void JRatioDiscr(double val) {jratio_discriminant_= val;}
		bool IsValid() const
		{
			if(WJa() != 0 && WJa() == WJb()) {return(false);}
			if(BHad() != 0 && (BHad() == WJa() || BHad() == WJb())) {return(false);}
			if(BLep() != 0 && (BLep() == BHad() || BLep() == WJa() || BLep() == WJb())) {return(false);}
			return(true);
		}

		bool AreBsCorrect(const Permutation& other) const //bjets are selected correct, but not necessarily at the right position!!!!!!!!!!
		{
			return((BLep() == other.BLep() && BHad() == other.BHad()) || (BHad() == other.BLep() && BLep() == other.BHad()));
		}
		bool AreJetsCorrect(const Permutation& other) const
		{
			if(BLep() != other.BLep() && BLep() != other.BHad() && BLep() != other.WJa() && BLep() != other.WJb()){return(false);}
			if(BHad() != other.BLep() && BHad() != other.BHad() && BHad() != other.WJa() && BHad() != other.WJb()){return(false);}
			if(WJa() != other.BLep() && WJa() != other.BHad() && WJa() != other.WJa() && WJa() != other.WJb()){return(false);}
			if(WJb() != other.BLep() && WJb() != other.BHad() && WJb() != other.WJa() && WJb() != other.WJb()){return(false);}
			return(true);
		}

		bool AreHadJetsCorrect(const Permutation& other) const
		{
			if(BHad() != other.BHad() && BHad() != other.WJa() && BHad() != other.WJb()){return(false);}
			if(WJa() != other.BHad() && WJa() != other.WJa() && WJa() != other.WJb()){return(false);}
			if(WJb() != other.BHad() && WJb() != other.WJa() && WJb() != other.WJb()){return(false);}
			return(true);
		}

		bool IsBLepCorrect(const Permutation& other) const
		{
			return(BLep() == other.BLep());
		}
		bool IsBHadCorrect(const Permutation& other) const
		{
			return(BHad() == other.BHad());
		}
		bool IsWHadCorrect(const Permutation& other) const
		{
			return((WJa() == other.WJa() && WJb() == other.WJb()) || (WJa() == other.WJb() && WJb() == other.WJa()));
		}
		bool IsTHadCorrect(const Permutation& other) const
		{
			return(IsBHadCorrect(other) && IsWHadCorrect(other));
		}
		bool IsCorrect(const Permutation& other) const
		{
			return(IsTLepCorrect(other) && IsTHadCorrect(other));
		}
		bool IsJetIn(TLorentzVector* jet)
		{
			return(jet == WJa() || jet == WJb() || jet == BHad() || jet == BLep());
		}

		bool IsLooselyCorrect(const Permutation& other) const
		{
			return(IsTLepCorrect(other) && AreHadJetsCorrect(other));
		}
		bool IsTLepCorrect(const Permutation& other) const
		{
			return(IsBLepCorrect(other) && (L() == other.L()) );
		}		
};

bool operator<(const Permutation& A, const Permutation& B);
bool operator>(const Permutation& A, const Permutation& B);

#endif