#ifndef NRMM_GIBBS_SAMPLER_H_
#define NRMM_GIBBS_SAMPLER_H_

#include "nrmm_sampler.h"

namespace npbayes
{
	class NRMMGibbsSampler : public NRMMSampler
	{
	public:
		NRMMGibbsSampler(NRM *mu, const std::vector<SuffStats*> &ss, 
			int subset_size = 0);
		~NRMMGibbsSampler(void) { }
		void SampleCl(int ind);
		void Init(void);
		void Init(const std::vector<int> &labels);
		void Sweep(void);
		double LogJoint(void) const { return mu_->LogJoint(clset_); }
		int NumClusters(void) const { return clset_.size(); }
		void GetLabels(std::vector<int> &labels) const;		
	protected:
		std::vector<int> indices_;
		ClusterSet clset_;
		std::vector<Cluster*> ind_to_cl_;
		int subset_size_;
	};
}

#endif