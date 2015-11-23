#ifndef TGMCMC_H_
#define TGMCMC_H_

#include "nrmm_sampler.h"
#include "bhc.h"

namespace npbayes
{
	class TGMCMC : public NRMMSampler
	{
	public:
		TGMCMC(NRM *mu, const std::vector<SuffStats*> &ss,
			int num_sm = 20, int depth = 2);
		~TGMCMC(void) { }
		void Init(bool noisy = false);
		
		Node * SampleSub(Node *nd, bool draw_leaf, 
			double &log_trans) const;		
		double SampleSub(Node *nd, Node *anc) const;
		void SelectiveGibbs(int depth);

		void StocInsert(Node::Set &ndset, Node::Queue &ndque, 
			double &log_trans);

		void SplitMerge(void); 
		void Sweep(void);
		void PrintAdditionalInfo(std::ostream &os)
		{
			for (int i = 0; i < log_r_.size(); ++i)
				os << log_r_[i] << std::endl;
		}

		double LogJoint(void) const { return bhc_.LogJoint(ndset_); }
		int NumClusters(void) const { return ndset_.size(); }
		void GetLabels(std::vector<int> &labels) const;

	private:
		BHC bhc_;
		Node::Set ndset_;
		std::vector<double> log_r_;
		int num_sm_, depth_;
	};
}


#endif