#ifndef NormalWishart_H_
#define NormalWishart_H_

#include "Vec.h"
#include "Math.h"
#include "base_measure.h"

namespace npbayes
{
	class NormalWishartSuffStats : public SuffStats
	{
	public:
		NormalWishartSuffStats(void) { }
		NormalWishartSuffStats(const Vec &x)
			: n_(1), x_(x), X_(x*x.transpose()) { }
		~NormalWishartSuffStats(void) { }
		SuffStats * Copy(void) const 
		{ 
			return new NormalWishartSuffStats(*this); 
		}
		void Add(SuffStats *ss);
		void Subtract(SuffStats *ss);
		friend class NormalWishart;
	private:
		int n_;
		Vec x_;
		Mat X_;
	};

	class NormalWishart : public BaseMeasure
	{
	public:
		NormalWishart(void) { }
		NormalWishart(const Mat &X, double c = 0.1);
		~NormalWishart(void) { }
		static double LogMvnGammaRatio(int nu, int n, int d);
		double LogMarginal(SuffStats *ss) const;
		double LogMarginal(SuffStats *ss0, SuffStats *ss1) const;
	private:
		int d_;
		double r_;
		int nu_;
		Vec xi_;
		Mat chol_Psi_inv_;
		Mat Psi_rxixit_;
		double nu_log_det_Psi_;
	};		
}

#endif