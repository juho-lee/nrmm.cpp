#include "cluster.h"

namespace npbayes
{
	void Cluster::Add(SuffStats *ss, int n)
	{
		this->n += n;
		if (this->ss) this->ss->Add(ss);
		else this->ss = ss->Copy();
	}

	void Cluster::Add(Cluster *cl)
	{
		n += cl->n;
		if (ss) ss->Add(cl->ss);
		else ss = cl->ss->Copy();
	}

	void Cluster::Subtract(SuffStats *ss, int n)
	{
		if (this->ss) {
			this->n -= n;
			if (this->n) this->ss->Subtract(ss);
			else {
				delete this->ss;
				this->ss = 0;
			}			
		}		
	}

	void Cluster::Subtract(Cluster *cl)
	{
		if (ss) {
			n -= cl->n;
			if (n) ss->Subtract(cl->ss);
			else {
				delete ss;
				ss = 0;
			}
		}
	}
}