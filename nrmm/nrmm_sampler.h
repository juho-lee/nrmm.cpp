#ifndef NRMM_SAMPLER_H_
#define NRMM_SAMPLER_H_

#include <fstream>
#include <iostream>
#include <algorithm>
#include "random.h"
#include "math.h"
#include "nrm.h"

namespace npbayes
{
	class NRMMSampler
	{
	public:
		NRMMSampler(NRM *mu, const std::vector<SuffStats*> &ss)
			: mu_(mu), ss_(ss) { }
		~NRMMSampler(void) { }
		virtual void Init(void) { }
		virtual void Init(const std::vector<int> &labels) { Init(); }
		virtual void Sweep(void) { }
		virtual double LogJoint(void) const { return 0; }
		virtual int NumClusters(void) const { return 0; }
		virtual void PrintAdditionalInfo(std::ostream &os) { }
		virtual void GetLabels(std::vector<int> &labels) const { }

		NRM * mu(void) const { return mu_; }
		const std::vector<SuffStats*> & ss(void) const { return ss_; }

	protected:
		NRM * mu_;
		const std::vector<SuffStats*> &ss_;
	};
}

#endif