#include "bhc.h"

namespace npbayes
{
	Node * BHC::MakeLeaf(int ind, SuffStats *ss) const
	{
		static const double log_k = mu_->LogKappa(1);
		Node *nd = new Node(ss);
		nd->ind = ind;
		nd->hv = phash<size_t>(ind, ind);
		nd->st.log_k = log_k;
		nd->st.log_m = mu_->LogMarginal(ss);
		nd->st.log_d = -INF;
		return nd;
	}

	void BHC::ComputeStats(Node *ch0, Node *ch1, Node::Stats &st) const
	{
		st.log_k = mu_->LogKappa(ch0->n + ch1->n);
		st.log_m = mu_->LogMarginal(ch0->ss, ch1->ss);
		st.log_d = ch0->st.log_t() + ch1->st.log_t() - st.log_h();
	}

	void BHC::ComputeStats(Node *nd, bool compute_log_h) const
	{
		if (compute_log_h) {
			nd->st.log_k = mu_->LogKappa(nd->n);
			nd->st.log_m = mu_->LogMarginal(nd->ss);
		}
		if (!nd->IsLeaf()) {
			nd->st.log_d = nd->ch0->st.log_t() + nd->ch1->st.log_t()
				- nd->st.log_h();
		}
	}

	Node * BHC::MakePair(Node *ch0, Node *ch1, const Node::Stats &st) const
	{
		Node *pair = new Node();
		pair->hv = phash<size_t>(ch0->ind, ch1->ind);
		pair->st = st;
		pair->ch0 = ch0;
		pair->ch1 = ch1;
		return pair;
	}

	void BHC::MergePair(Node *pair) const
	{
		if (pair->ch0 && pair->ch1) {			
			pair->Add(pair->ch0);
			pair->Add(pair->ch1);
			pair->ch0->par = pair->ch1->par = pair;
		}
	}

	void BHC::GreedyConstruction(Node::Set & ndset) const
	{
		Node::SortedSet ndsset;
		Node::Stats st;
		for (auto it = ndset.begin(); it != ndset.end(); ++it) {
			auto jt = it;
			for (++jt; jt != ndset.end(); ++jt) {
				ComputeStats(*it, *jt, st);
				if (st.log_d < 0)
					ndsset.insert(MakePair(*it, *jt, st));
			}
		}
		while (!ndsset.empty())
		{
			Node *pair = *ndsset.begin();
			ndsset.erase(ndsset.begin());
			Node *ch0 = pair->ch0, *ch1 = pair->ch1;
			ndset.erase(ch0);
			ndset.erase(ch1);
			for (auto it = ndsset.begin(); it != ndsset.end();) {
				if ((*it)->ch0 == ch0 || (*it)->ch1 == ch0 ||
					(*it)->ch0 == ch1 || (*it)->ch1 == ch1) {
					delete *it;
					it = ndsset.erase(it);
				}
				else ++it;
			}
			MergePair(pair);
			foreach(it, ndset) {
				ComputeStats(it, pair, st);
				if (st.log_d < 0)
					ndsset.insert(MakePair(it, pair, st));
			}
			ndset.insert(pair);
		}
	}

	void BHC::SeqMergePair(Node *pair) const
	{
		Node *nd = pair->ch0, *leaf = pair->ch1;
		if (nd->IsLeaf()) std::swap(nd, leaf);
		Node::Stats st0, st1;
		while (!nd->IsLeaf()) {			
			ComputeStats(nd->ch0, leaf, st0);
			ComputeStats(nd->ch1, leaf, st1);
			if (nd->st.log_d < st0.log_d && nd->st.log_d < st1.log_d)
				break;
			else {
				nd->Add(leaf);
				nd->st = pair->st;
				if (st0.log_d < st1.log_d) {
					pair->st = st0;
					nd = nd->ch0;
				}
				else {
					pair->st = st1;
					nd = nd->ch1;
				}
			}
		}		
		if (nd->par) {
			pair->par = nd->par;
			Node *sib = nd->sib();
			nd->par->ch0 = pair;
			nd->par->ch1 = sib;			
		}		
		pair->ch0 = nd;
		pair->ch1 = leaf;
		MergePair(pair);
	}

	void BHC::Destroy(Node *nd, Node::Queue &ndque) const
	{
		Node *par = nd->par, *sib = nd->sib(), *temp;
		nd->par = 0;
		while (par) {
			temp = par;
			sib->par = 0;
			ndque.push(sib);
			sib = par->sib();
			par = par->par;
			delete temp;
		}
	}

	void BHC::Update(Node *&nd, Node::Queue &ndque) const
	{
		Node *par = nd->par;
		while (par) {
			ComputeStats(par, false);
			if (par->st.log_d > 0) {
				Destroy(nd, ndque);
				break;
			}
			nd = par;
			par = nd->par;
		}
	}

	void BHC::Update(Node *&nd) const
	{
		Node *par = nd->par;
		while (par) {
			ComputeStats(par, false);
			nd = par;
			par = nd->par;
		}
	}

	void BHC::Insert(Node::Set &ndset, Node *nd0, Node *nd1,
		const Node::Stats &st, Node::Queue &ndque) const
	{
		ndset.erase(nd0);
		Node *pair = MakePair(nd0, nd1, st);
		SeqMergePair(pair);
		Update(pair, ndque);
		ndset.insert(pair);
	}

	void BHC::Insert(Node::Set &ndset, Node *nd0, Node *nd1,
		const Node::Stats &st) const
	{
		ndset.erase(nd0);
		Node *pair = MakePair(nd0, nd1, st);
		SeqMergePair(pair);
		Update(pair);
		ndset.insert(pair);
	}

	void BHC::SequentialConstruction(Node::Set &ndset, bool noisy) const
	{
		Node::Queue ndque;
		Node::Shuffle(ndset, ndque);
		ndset.clear();
		while (!ndque.empty()) {
			Node *nd = ndque.front();
			ndque.pop();
			if (ndset.empty()) ndset.insert(nd);
			else {
				Node *nnd = 0;
				Node::Stats st, temp;
				foreach(it, ndset) {
					ComputeStats(nd, it, temp);
					if (nnd == 0 || temp.log_d < st.log_d) {
						nnd = it;
						st = temp;
					}
				}
				if (st.log_d > 0) ndset.insert(nd);
				else {
					if (noisy) {
						ndset.erase(nnd);
						nnd = MakePair(nd, nnd, st);
						MergePair(nnd);
						ndset.insert(nnd);
					}
					else Insert(ndset, nnd, nd, st, ndque);
				}
			}
		}
	}

	void BHC::Detach(Node::Set &ndset, Node *nd) const
	{
		Node *par = nd->par;
		Node *sib = nd->sib();
		nd->par = 0;
		if (par) {
			if (par->par) {
				sib->par = par->par;
				if (par->par->ch0 == par) par->par->ch0 = sib;
				else par->par->ch1 = sib;										
				delete par;
				Node *temp = sib;
				par = temp->par;
				sib = temp->sib();
				while (par) {
					ComputeStats(temp, sib, par->st);
					par->Subtract(nd);
					temp = par;
					sib = temp ? temp->sib() : 0;
					par = temp ? temp->par : 0;
				}
			}
			else {
				ndset.erase(par);
				sib->par = 0;
				delete par;
				ndset.insert(sib);
			}
		}
		else ndset.erase(nd);
	}


	void BHC::Recompute(Node *nd) const
	{
		std::stack<Node*> ndstack;
		do {
			while (nd) {
				if (nd->ch1) ndstack.push(nd->ch1);
				ndstack.push(nd);
				nd = nd->ch0;
			}
			nd = ndstack.top();
			ndstack.pop();
			if (!ndstack.empty() && (nd->ch1 && ndstack.top() == nd->ch1)) {
				ndstack.pop();
				ndstack.push(nd);
				nd = nd->ch1;
			}
			else {
				ComputeStats(nd);
				nd = 0;
			}
		} while (!ndstack.empty());
	}

	bool BHC::FindError(Node *nd) const
	{
		Node::Queue ndque;
		ndque.push(nd);
		while (!ndque.empty()) {
			nd = ndque.front();
			ndque.pop();
			if (!nd->IsLeaf()) {
				if (nd->ch0->par != nd || nd->ch1->par != nd) 
					return true;
				double k1 = mu_->LogKappa(nd->n);
				double k2 = nd->st.log_k;
				double m1 = mu_->LogMarginal(nd->ss), m2 = nd->st.log_m;
				double d1 = nd->ch0->st.log_t() + nd->ch1->st.log_t() -k1 - m1;				
				double d2 = nd->st.log_d;
				double eps = 1.0e-7;
				if (abs(k1 - k2) > eps || abs(m1 - m2) > eps || abs(d1 - d2) > eps)
					return true;
				ndque.push(nd->ch0);
				ndque.push(nd->ch1);
			}
		}
		return false;
	}

	bool BHC::FindError(const Node::Set &ndset, int n) const
	{
		Node *err = 0;
		int sum = 0;
		foreach(it, ndset) {
			if (FindError(it)) return true;
			sum += it->n;
		}
		if (sum != n)
			return true;
		else return false;
	}
}