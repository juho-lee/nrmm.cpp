#ifndef BHC_H_
#define BHC_H_

#include <stack>
#include "nrm.h"
#include "node.h"

namespace npbayes
{
	class BHC
	{
	public:
		BHC(NRM *mu) : mu_(mu) { }
		~BHC(void) { }

		Node * MakeLeaf(int ind, SuffStats *ss) const;

		void ComputeStats(Node *ch0, Node *ch1, Node::Stats &st) const;
		void ComputeStats(Node *nd, bool compute_log_h = true) const;

		Node * MakePair(Node *ch0, Node *ch1, const Node::Stats &st) const;
		void MergePair(Node *pair) const;
		void GreedyConstruction(Node::Set &ndset) const;

		void SeqMergePair(Node *pair) const;
		void Destroy(Node *nd, Node::Queue &ndque) const;
		void Update(Node *&nd, Node::Queue &ndque) const;
		void Update(Node *&nd) const;
		void Insert(Node::Set &ndset, Node *nd0, Node *nd1,
			const Node::Stats &st, Node::Queue &ndque) const;
		void Insert
			(Node::Set &ndset, Node *nd0, Node *nd1, Node::Queue &ndque) const
		{
			Node::Stats st;
			ComputeStats(nd0, nd1, st);
			Insert(ndset, nd0, nd1, st, ndque);
		}
		void Insert(Node::Set &ndset, Node *nd0, Node *nd1,
			const Node::Stats &st) const;
		void Insert(Node::Set &ndset, Node *nd0, Node *nd1) const
		{
			Node::Stats st;
			ComputeStats(nd0, nd1, st);
			Insert(ndset, nd0, nd1, st);
		}
		void SequentialConstruction(Node::Set &ndset, bool noisy = false) const;

		void Detach(Node::Set &ndset, Node *nd) const;

		void Recompute(Node *nd) const;

		double LogJoint(const Node::Set &ndset) const
		{
			ClusterSet clset(ndset.begin(), ndset.end());
			return mu_->LogJoint(clset);
		}

		bool FindError(Node *nd) const;
		bool FindError(const Node::Set &ndset, int n) const;
	private:
		NRM *mu_;
	};
}

#endif