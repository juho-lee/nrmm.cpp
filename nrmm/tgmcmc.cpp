#include "tgmcmc.h"

namespace npbayes
{
	TGMCMC::TGMCMC(NRM *mu, const std::vector<SuffStats*> &ss,
		int num_sm, int depth) : NRMMSampler(mu, ss), bhc_(BHC(mu)), 
		num_sm_(num_sm), depth_(depth)
	{
		for (int i = 0; i < ss.size(); ++i)
			ndset_.insert(bhc_.MakeLeaf(i, ss[i]));
	}

	void TGMCMC::Init(bool noisy)
	{
		bhc_.SequentialConstruction(ndset_, noisy);
	}

	Node * TGMCMC::SampleSub(Node *nd, bool draw_leaf, 
		double &log_trans) const
	{
		if (nd->IsLeaf()) return nd;
		std::vector<Node*> ndvec;
		std::vector<double> p;
		double max_log_d = -INF;
		Node::Queue ndque;
		ndque.push(nd);
		while (!ndque.empty()) {
			nd = ndque.front();
			ndque.pop();
			if (draw_leaf || !nd->IsLeaf()) {
				ndvec.push_back(nd);
				p.push_back(nd->st.log_d);
				max_log_d = std::max(max_log_d, nd->st.log_d);
				if (!nd->IsLeaf()) {
					ndque.push(nd->ch0);
					ndque.push(nd->ch1);
				}
			}
		}
		double nc = -INF;
		for (int j = 0; j < p.size(); ++j) {
			p[j] = LogSumExp(p[j], max_log_d);
			nc = LogSumExp(p[j], nc);
		}
		for (int j = 0; j < p.size(); ++j)
			p[j] = exp(p[j] - nc);
		int j = RandMult(p);
		log_trans += log(p[j]);
		return ndvec[j];
	}

	double TGMCMC::SampleSub(Node *nd, Node *target) const
	{
		if (nd->IsLeaf()) return 0;
		int pos = 0;
		std::vector<double> p;
		double max_log_d = -INF;
		Node::Queue ndque;
		ndque.push(nd);
		while (!ndque.empty()) {
			nd = ndque.front();
			ndque.pop();
			if (!nd->IsLeaf()) {
				p.push_back(nd->st.log_d);
				if (nd == target) pos = p.size() - 1;
				max_log_d = std::max(max_log_d, nd->st.log_d);
				if (!nd->IsLeaf()) {
					ndque.push(nd->ch0);
					ndque.push(nd->ch1);
				}
			}
		}
		double nc = -INF;
		for (int j = 0; j < p.size(); ++j) {
			p[j] = LogSumExp(p[j], max_log_d);
			nc = LogSumExp(p[j], nc);
		}
		return p[pos] - nc;
	}

	void TGMCMC::SelectiveGibbs(int depth)
	{
		Node::Queue ndque;
		double log_trans = 0;
		foreach(it, ndset_) {
			Node *nd = it;
			for (int i = 0; i < depth; ++i) {
				nd = SampleSub(nd, true, log_trans);
				if (nd->IsLeaf()) break;
			}
			nd->GetLeaves(ndque);
		}				
		while (!ndque.empty()) {
			Node *leaf = ndque.front();
			ndque.pop();
			Node *rt = leaf->GetRoot();
			if (rt == leaf) ndset_.erase(rt);						
			std::vector<double> p(ndset_.size() + 1, 0);
			int j = 0;
			double nc = -INF;
			foreach(it, ndset_) {				
				if (it == rt) {
					SuffStats *ss = rt->ss->Copy();
					ss->Subtract(leaf->ss);
					p[j] = mu_->LogKappaJoin(rt->n - 1) + rt->st.log_m
						- mu_->LogMarginal(ss);
					delete ss;								
				}
				else {
					p[j] = mu_->LogKappaJoin(it->n) - it->st.log_m
						+ mu_->LogMarginal(leaf->ss, it->ss);
				}
				nc = LogSumExp(p[j++], nc);
			}
			p[j] = mu_->LogKappaNew() + leaf->st.log_m;
			nc = LogSumExp(p[j], nc);
			for (j = 0; j < p.size(); ++j)
				p[j] = exp(p[j] - nc);
			j = RandMult(p);
			if (j == ndset_.size()) {
				bhc_.Detach(ndset_, leaf);
				ndset_.insert(leaf);
			}
			else {
				Node *nd = *std::next(ndset_.begin(), j);
				if (nd != rt) {
					bhc_.Detach(ndset_, leaf);
					bhc_.Insert(ndset_, nd, leaf);
				}
			}
		}
	}
	
	void TGMCMC::StocInsert(Node::Set &ndset, Node::Queue &ndque, 
		double &log_trans)
	{
		while (!ndque.empty()) {
			Node *nd = ndque.front();
			ndque.pop();
			if (ndset.empty()) ndset.insert(nd);
			else {
				std::vector<double> p(ndset.size(), 0);
				std::vector<Node::Stats> st(ndset.size(), Node::Stats());
				p.push_back(0);
				double nc = 0;
				int j = 0;
				foreach(it, ndset) {
					bhc_.ComputeStats(nd, it, st[j]);
					p[j] = -st[j].log_d;
					nc = LogSumExp(p[j++], nc);
				}
				for (j = 0; j < p.size(); ++j)
					p[j] = exp(p[j] - nc);
				j = RandMult(p);
				log_trans += log(p[j]);
				if (j == ndset.size()) ndset.insert(nd);
				else {
					bhc_.Insert(ndset, *std::next(ndset.begin(), j),
						nd, st[j]);
				}
			}
		}
	}

	void TGMCMC::SplitMerge(void)
	{
		double log_trans = 0;
		double log_rtrans = 0;

		// step 1: pick a node
		Node *nd = *std::next(ndset_.begin(), rand() % ndset_.size());
		log_trans -= log(ndset_.size());
		ndset_.erase(nd);

		// collect the set M
		Node::Set M;
		for (auto it = ndset_.begin(); it != ndset_.end(); ) {
			Node::Stats st;
			bhc_.ComputeStats(nd, *it, st);
			double log_p = -LogSumExp(0, st.log_d);
			if (log(Randu()) < log_p) {
				M.insert(*it);
				log_trans += log_p;
				it = ndset_.erase(it);
			}	
			else {
				log_trans += st.log_d + log_p;
				++it;
			}
		}

		if (M.empty()) { // split
			if (nd->IsLeaf()) {
				ndset_.insert(nd);
				return;
			}
			Node *nd_copy = Copy(nd);
			Node *anc = SampleSub(nd, false, log_trans);
			Node::Queue Q;
			bhc_.Destroy(anc, Q);
			anc->ch0->par = 0;
			anc->ch1->par = 0;
			Node::Set S;
			S.insert(anc->ch0);
			S.insert(anc->ch1);
			delete anc;
			StocInsert(S, Q, log_trans);
			
			// compute reverse transition prob
			foreach(it, S) Q.push(it);
	
			int k = ndset_.size() + S.size();
			double log_rtrans_sub = -INF;
			while (!Q.empty()) {
				double log_rtrans_sub_sub = -log(k);
				Node *temp = Q.front();
				Q.pop();
				Node::Stats st;
				foreach(it, ndset_) {
					bhc_.ComputeStats(temp, it, st);
					log_rtrans_sub_sub += st.log_d - LogSumExp(0, st.log_d);
				}
				foreach(it, S) {
					if (it != temp) {
						bhc_.ComputeStats(temp, it, st);
						log_rtrans_sub_sub += -LogSumExp(0, st.log_d);
					}
				}
				log_rtrans_sub = LogSumExp(log_rtrans_sub, log_rtrans_sub_sub);
			}
			log_rtrans += log_rtrans_sub;

			double log_r = log_rtrans -log_trans - nd_copy->st.log_h();
			foreach(it, S) log_r += it->st.log_h();
			log_r_.push_back(log_r);
			if (Randu() < std::min(1.0, exp(log_r))) {
				Delete(nd_copy);
				ndset_.insert(S.begin(), S.end());
			}
			else {
				foreach(it, S) Delete(it);
				ndset_.insert(nd_copy);
			}
		}
		else { // merge
			Node *anc = 0;
			Node *merged = nd;
			Node::Set S;

			foreach(it, M) {
				Node::Stats st;
				bhc_.ComputeStats(merged, it, st);
				merged = bhc_.MakePair(merged, it, st);
				bhc_.MergePair(merged);
				if (!anc) {
					anc = merged;
					S.insert(anc->ch0);
					S.insert(anc->ch1);
				}
				else {
					double nc = 0;
					foreach(jt, S) {
						bhc_.ComputeStats(jt, it, st);
						nc = LogSumExp(nc, -st.log_d);
					}
					log_rtrans -= nc;
					S.insert(it);
				}
			}
			
			log_rtrans += -log(ndset_.size() + 1);
			foreach(it, ndset_) {
				Node::stats st;
				bhc_.ComputeStats(merged, it, st);
				log_rtrans += st.log_d - LogSumExp(0, st.log_d);
			}
			log_rtrans += SampleSub(merged, anc);

			double log_r = log_rtrans - log_trans + merged->st.log_h()
				- nd->st.log_h();
			foreach(it, M) log_r -= it->st.log_h();
			log_r_.push_back(log_r);
			
			if (Randu() < std::min(1.0, exp(log_r))) 
				ndset_.insert(merged);
			else {
				Node::Queue Q;
				bhc_.Destroy(anc, Q);
				anc->ch0->par = 0;
				anc->ch1->par = 0;
				ndset_.insert(anc->ch0);
				ndset_.insert(anc->ch1);
				delete(anc);
				while (!Q.empty()) {
					ndset_.insert(Q.front());
					Q.pop();
				}
			}
		}
	}

	void TGMCMC::Sweep(void)
	{
		// int num_sm = 20;
		for (int i = 0; i < num_sm_; ++i) SplitMerge();
		
		// int depth = 2;
		SelectiveGibbs(depth_);
	
		ClusterSet clset(ndset_.begin(), ndset_.end());
		mu_->SampleU(clset);
		mu_->SampleHyperparams(clset);		
		foreach(it, ndset_) bhc_.Recompute(it);

		//if (bhc_.FindError(ndset_, ss_.size())) {
		//	std::cout << "something is wrong" << std::endl;
		//}
	}

	void TGMCMC::GetLabels(std::vector<int> &labels) const
	{
		labels.assign(ss_.size(), 0);
		int l = 0;
		foreach(it, ndset_) {
			Node::Set leaves;
			it->GetLeaves(leaves);
			foreach(jt, leaves)
				labels[jt->ind] = l;
			++l;
		}
	}
}

