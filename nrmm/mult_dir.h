#ifndef MULT_DIR_H_
#define MULT_DIR_H_

#include <map>
#include "Math.h"
#include "base_measure.h"

namespace npbayes
{
	class MultDirSuffStats : public SuffStats
	{
	public:
		MultDirSuffStats(void) : sum_(0), log_norm_(0) { }
		MultDirSuffStats(const std::map<int, int> &hist);
		~MultDirSuffStats(void) { }
		SuffStats * Copy(void) const
		{
			return new MultDirSuffStats(*this);
		}
		void Add(SuffStats *ss);
		void Subtract(SuffStats *ss);
		friend class MultDir;
	private:
		int sum_;
		double log_norm_;
		std::map<int, int> hist_;
	};

	class MultDir : public BaseMeasure
	{
	public:
		MultDir(int dim, double alpha)
			: dim_(dim), alpha_(alpha), lg_alpha_(LogGamma(alpha_)),
			sum_alpha_(dim*alpha), lg_sum_alpha_(LogGamma(sum_alpha_)) { }
		double LogMarginal(SuffStats *ss) const;
		double LogMarginal(SuffStats *ss0, SuffStats *ss1) const;
	private:
		int dim_;
		double alpha_, lg_alpha_;
		double sum_alpha_, lg_sum_alpha_;		
	};
}

#endif