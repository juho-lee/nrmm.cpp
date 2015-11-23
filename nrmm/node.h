#ifndef NODE_H_
#define NODE_H_

#include <queue>
#include "cluster.h"
#include "phash.h"
#include "math.h"

namespace npbayes
{
	class Node : public Cluster
	{
	public:
		typedef struct Stats
		{
			double log_k;
			double log_m;
			double log_d;
			double log_h(void) { return log_k + log_m; }
			double log_t(void) {
				return log_k + log_m + LogSumExp(0, log_d);
			}
		} stats;

		typedef struct DistLessThan
		{
			bool operator()(const Node *nd0, const Node *nd1) const
			{
				if (nd0->st.log_d == nd1->st.log_d) {
					if (nd0->hv == nd1->hv) return nd0 < nd1;
					else return nd0->hv < nd1->hv;
				}
				else return nd0->st.log_d < nd1->st.log_d;
			}
		} DistLessThan;
		typedef std::set<Node*, DistLessThan> SortedSet;

		typedef struct LessThan
		{
			bool operator()(const Node *nd0, const Node *nd1) const
			{
				if (nd0->hv == nd1->hv) return nd0 < nd1;
				else return nd0->hv < nd1->hv;
			}
		} LessThan;
		typedef std::set<Node*, LessThan> Set;

		typedef std::queue<Node*> Queue;
		static void Shuffle(const Set &ndset, Queue &queue);

		Node(void) : Cluster(), ind(-1), par(0), ch0(0), ch1(0) { }
		Node(SuffStats *ss) : Cluster(ss), ind(-1), par(0), ch0(0), ch1(0) { }
		~Node(void) { }

		bool IsLeaf(void) const { return !ch0 && !ch1 ; }
		void GetLeaves(Set &leaves);
		void GetLeaves(Queue &leaves);
		Node * sib(void) const 
		{ 
			if (par) {
				if (par->ch0 == this) return par->ch1;
				else return par->ch0;
			}
			return 0;
		}
		Node * GetRoot(void);
		
		int ind;
		size_t hv;
		Stats st;
		Node *par;
		Node *ch0;
		Node *ch1;
	};

	void Delete(Node *nd);
	Node * Copy(Node *src);
}


#endif