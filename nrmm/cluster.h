#ifndef CLUSTER_H_
#define CLUSTER_H_

#include <set>
#include "suff_stats.h"

namespace npbayes
{
	class Cluster
	{
	public:
		Cluster(void) : n(0), ss(0) { }
		Cluster(SuffStats *ss) : n(1), ss(ss->Copy()) { }
		~Cluster(void) { if (ss) delete ss; }		
		void Add(SuffStats *ss, int n = 1);
		void Add(Cluster *cl);
		void Subtract(SuffStats *ss, int n = 1);
		void Subtract(Cluster *cl);
		int n;
		SuffStats *ss;		
	};	

	typedef std::set<Cluster*> ClusterSet;
}

#endif