#ifndef NRMM_SPLIT_MERGE_SAMPLER_H_
#define NRMM_SPLIT_MERGE_SAMPLER_H_

#include "nrmm_gibbs_sampler.h"

namespace npbayes
{
	class NRMMSplitMergeSampler : public NRMMGibbsSampler
	{
	public:
		NRMMSplitMergeSampler(NRM *mu, const std::vector<SuffStats*> &ss, 
			int subset_size = 0) : NRMMGibbsSampler(mu, ss, subset_size) { }
		static const int num_rgibbs_iter;
		static const int num_gibbs_iter;
		double RestrictedGibbs(std::map<int, int> &cind, 
			std::vector<Cluster*> &cl) const;
		double RestrictedGibbs(const std::map<int, int> &cind_target,
			std::map<int, int> &cind, std::vector<Cluster*> &cl) const;
		double log_h(Cluster *cl) const 
		{ 
			return mu_->LogKappa(cl->n) + mu_->LogMarginal(cl->ss);
		}
		double log_h(Cluster *cl0, Cluster *cl1) const 
		{
			return mu_->LogKappa(cl0->n + cl1->n) 
				+ mu_->LogMarginal(cl0->ss, cl1->ss);
		}
		void Sweep(void);
		void PrintAdditionalInfo(std::ostream &os)
		{
			for (int i = 0; i < log_r_.size(); ++i)
				os << log_r_[i] << std::endl;
		}
	private:
		std::vector<double> log_r_;
	};
}

#endif