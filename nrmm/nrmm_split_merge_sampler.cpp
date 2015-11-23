#include "nrmm_split_merge_sampler.h"

namespace npbayes
{
	const int NRMMSplitMergeSampler::num_rgibbs_iter = 10;
	const int NRMMSplitMergeSampler::num_gibbs_iter = 1;

	double NRMMSplitMergeSampler::RestrictedGibbs
		(std::map<int, int> &cind, std::vector<Cluster*> &cl) const
	{
		double log_trans = 0;
		std::vector<int> indices;
		double log_p0, log_p1;
		foreach(it, cind) indices.push_back(it.first);
		std::random_shuffle(indices.begin(), indices.end());
		foreach(ind, indices) {
			cl[cind[ind]]->Subtract(ss_[ind]);
			log_p0 = mu_->LogKappaJoin(cl[0]->n) + 
				mu_->LogPred(ss_[ind], cl[0]->ss);
			log_p1 = mu_->LogKappaJoin(cl[1]->n) +
				mu_->LogPred(ss_[ind], cl[1]->ss);
			double nc = LogSumExp(log_p0, log_p1);
			if (Randu() < exp(log_p0 - nc)) {
				log_trans += log_p0 - nc;
				cl[0]->Add(ss_[ind]);
				cind[ind] = 0;
			}
			else {
				log_trans += log_p1 - nc;
				cl[1]->Add(ss_[ind]);
				cind[ind] = 1;
			}
		}
		return log_trans;
	}

	double NRMMSplitMergeSampler::RestrictedGibbs
		(const std::map<int, int> &cind_target,
		std::map<int, int> &cind, std::vector<Cluster*> &cl) const
	{
		double log_trans = 0;
		std::vector<int> indices;
		double log_p0, log_p1;
		foreach(it, cind) indices.push_back(it.first);
		foreach(ind, indices) {
			cl[cind[ind]]->Subtract(ss_[ind]);
			log_p0 = mu_->LogKappaJoin(cl[0]->n) +
				mu_->LogPred(ss_[ind], cl[0]->ss);
			log_p1 = mu_->LogKappaJoin(cl[1]->n) +
				mu_->LogPred(ss_[ind], cl[1]->ss);
			double nc = LogSumExp(log_p0, log_p1);
			if (cind_target.at(ind) == 0) {
				log_trans += log_p0 - nc;
				cl[0]->Add(ss_[ind]);
			}
			else {
				log_trans += log_p1 - nc;
				cl[1]->Add(ss_[ind]);
			}
		}
		return log_trans;
	}

	void NRMMSplitMergeSampler::Sweep(void)
	{
		std::random_shuffle(indices_.begin(), indices_.end());
		int ind0 = indices_[0], ind1 = indices_[1];
		Cluster *cl0 = ind_to_cl_[ind0];
		Cluster *cl1 = ind_to_cl_[ind1];

		// construct the launch state
		std::vector<Cluster*> cl_copy(2, 0);
		cl_copy[0] = new Cluster(ss_[ind0]);		
		cl_copy[1] = new Cluster(ss_[ind1]);		

		std::map<int, int> cind;
		for (int i = 2; i < indices_.size(); ++i) {
			if (ind_to_cl_[indices_[i]] == cl0 ||
				ind_to_cl_[indices_[i]] == cl1) {
				int j = Randu() < 0.5 ? 0 : 1;
				cind[indices_[i]] = j;
				cl_copy[j]->Add(ss_[indices_[i]]);
			}
		}

		for (int i = 0; i < num_rgibbs_iter; ++i)
			RestrictedGibbs(cind, cl_copy);

		if (cl0 == cl1) { // split
			double log_trans = RestrictedGibbs(cind, cl_copy);
			double log_r = log_h(cl_copy[0]) + log_h(cl_copy[1])
				- log_h(cl0) - log_trans;
			log_r_.push_back(log_r);
			if (Randu() < std::min(1.0, exp(log_r))) {
				clset_.erase(cl0);
				delete cl0;
				clset_.insert(cl_copy[0]);
				clset_.insert(cl_copy[1]);
				ind_to_cl_[ind0] = cl_copy[0];
				ind_to_cl_[ind1] = cl_copy[1];
				foreach(it, cind)
					ind_to_cl_[it.first] = cl_copy[it.second];				
			}			
		}
		else {	// merge
			std::map<int, int> cind_target;
			foreach(it, cind) {
				if (ind_to_cl_[it.first] == cl0) cind_target[it.first] = 0;
				else cind_target[it.first] = 1;
			}
			double log_trans = RestrictedGibbs(cind_target, cind, cl_copy);			
			double log_r = log_trans + log_h(cl0, cl1)
				- log_h(cl0) - log_h(cl1);
			log_r_.push_back(log_r);
			if (Randu() < std::min(1.0, exp(log_r))) {
				clset_.erase(cl0);
				clset_.erase(cl1);
				cl0->Add(cl1);				
				clset_.insert(cl0);
				delete cl1;
				ind_to_cl_[ind0] = ind_to_cl_[ind1] = cl0;
				foreach(it, cind) ind_to_cl_[it.first] = cl0;				
			}			
		}
		for (int i = 0; i < num_gibbs_iter; ++i)
			NRMMGibbsSampler::Sweep();
	}
}