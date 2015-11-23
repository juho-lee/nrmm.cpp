#include "nrmm_gibbs_sampler.h"

namespace npbayes
{	
	NRMMGibbsSampler::NRMMGibbsSampler(NRM *mu, const std::vector<SuffStats*> &ss,
		int subset_size) : NRMMSampler(mu, ss), subset_size_(subset_size)
	{
		indices_.assign(ss_.size(), 0);
		for (int i = 0; i < ss_.size(); ++i)
			indices_[i] = i;
		ind_to_cl_.assign(ss_.size(), 0);
	}

	void NRMMGibbsSampler::SampleCl(int ind)
	{
		std::vector<double> p(clset_.size() + 1, 0);
		double nc = -INF;
		int j = 0;
		foreach(it, clset_) {
			p[j] = mu_->LogKappaJoin(it->n) + mu_->LogPred(ss_[ind], it->ss);
			nc = LogSumExp(nc, p[j++]);
		}
		p[j] = mu_->LogKappaNew() + mu_->LogMarginal(ss_[ind]);
		nc = LogSumExp(nc, p[j]);
		for (j = 0; j < p.size(); ++j)
			p[j] = exp(p[j] - nc);
		j = RandMult(p);
		if (j < clset_.size()) {
			Cluster *cl = *std::next(clset_.begin(), j);
			cl->Add(ss_[ind]);
			ind_to_cl_[ind] = cl;
		}
		else {
			Cluster *new_cl = new Cluster(ss_[ind]);
			clset_.insert(new_cl);
			ind_to_cl_[ind] = new_cl;
		}
	}

	void NRMMGibbsSampler::Init(void)
	{
		std::random_shuffle(indices_.begin(), indices_.end());
		foreach(it, indices_) SampleCl(it);
	}

	void NRMMGibbsSampler::Init(const std::vector<int> &labels)
	{
		std::map<int, Cluster*> cind_to_cl;
		for (int i = 0; i < labels.size(); ++i) {
			auto it = cind_to_cl.find(labels[i]);
			Cluster *cl;
			if (it == cind_to_cl.end()) {
				cl = new Cluster();
				cind_to_cl[labels[i]] = cl;
			}
			else cl = it->second;
			clset_.insert(cl);
			cl->Add(ss_[i]);
			ind_to_cl_[i] = cl;
		}
	}

	void NRMMGibbsSampler::Sweep(void)
	{		
		std::random_shuffle(indices_.begin(), indices_.end());
		std::set<int> indices;
		if (subset_size_) {
			indices.insert(indices_.begin(),
				std::next(indices_.begin(), subset_size_));
		}
		else indices.insert(indices_.begin(), indices_.end());
		
		foreach(it, indices) {
			Cluster *old_cl = ind_to_cl_[it];
			old_cl->Subtract(ss_[it]);
			if (old_cl->n == 0) {
				clset_.erase(old_cl);
				delete old_cl;
			}
			SampleCl(it);
		}
		mu_->SampleU(clset_);
		mu_->SampleHyperparams(clset_);
	}

	void NRMMGibbsSampler::GetLabels(std::vector<int> &labels) const
	{
		labels.assign(ss_.size(), 0);
		int l = 0;
		std::map<Cluster*, int> lmap;
		foreach(it, clset_) lmap[it] = l++;
		for (int i = 0; i < ss_.size(); ++i)
			labels[i] = lmap[ind_to_cl_[i]];
	}
}