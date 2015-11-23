#include "mult_dir.h"

namespace npbayes
{
	MultDirSuffStats::MultDirSuffStats(const std::map<int, int> &hist)
		: sum_(0), log_norm_(0), hist_(hist)
	{		
		foreach(it, hist) {
			sum_ += it.second;
			log_norm_ -= LogGamma(it.second + 1.0);
		}
		log_norm_ += LogGamma(sum_ + 1.0);
	}

	void MultDirSuffStats::Add(SuffStats *ss)
	{
		MultDirSuffStats *mdss = (MultDirSuffStats *)(ss);
		foreach(it, mdss->hist_) 
			hist_[it.first] += it.second;
		sum_ += mdss->sum_;
		log_norm_ += mdss->log_norm_;
	}

	void MultDirSuffStats::Subtract(SuffStats *ss)
	{
		MultDirSuffStats *mdss = (MultDirSuffStats *)(ss);
		foreach(it, mdss->hist_) {
			hist_[it.first] -= it.second;
			if (hist_[it.first] == 0) hist_.erase(it.first);
		}
		sum_ -= mdss->sum_;
		log_norm_ -= mdss->log_norm_;
	}

	double MultDir::LogMarginal(SuffStats *ss) const
	{
		MultDirSuffStats *mdss = (MultDirSuffStats *)(ss);
		double y = lg_sum_alpha_ + mdss->log_norm_ 
			- LogGamma(sum_alpha_ + mdss->sum_);
		foreach(it, mdss->hist_)
			y += LogGamma(alpha_ + it.second) - lg_alpha_;
		return y;
	}

	double MultDir::LogMarginal(SuffStats *ss0, SuffStats *ss1) const
	{
		MultDirSuffStats *mdss0 = (MultDirSuffStats *)(ss0);
		MultDirSuffStats *mdss1 = (MultDirSuffStats *)(ss1);
		double y = lg_sum_alpha_ + mdss0->log_norm_ + mdss1->log_norm_
			- LogGamma(sum_alpha_ + mdss0->sum_ + mdss1->sum_);
		std::map<int, int> hist = mdss0->hist_;
		foreach(it, mdss1->hist_) hist[it.first] += it.second;
		foreach(it, hist)
			y += LogGamma(alpha_ + it.second) - lg_alpha_;
		return y;
	}
}