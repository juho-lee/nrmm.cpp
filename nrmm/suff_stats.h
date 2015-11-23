#ifndef SUFF_STATS_H_
#define SUFF_STATS_H_

namespace npbayes
{
	class SuffStats
	{
	public:
		SuffStats(void) { }
		virtual ~SuffStats(void) { }
		virtual SuffStats * Copy(void) const = 0;
		virtual void Add(SuffStats *ss) = 0;
		virtual void Subtract(SuffStats *ss) = 0;
	};
}

#endif